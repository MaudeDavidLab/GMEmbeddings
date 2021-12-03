

# This function takes in the best_hits datatable from the blast output, and creates an asv by embedding database asv
# matrix where each element is 1 over the number of hits for that ASV in the entire database. Instead of breaking ties arbitrarily
# we allow each ASV to hit multiple embedding ASV
getHitTransformation <- function(best_hits){
  asv2embed_split <- split(best_hits, best_hits$qseqid)
  asv2embed_list <- lapply(asv2embed_split, function(x) {
    num <- data.frame(matrix(rep(1, nrow(x)), nrow = 1))
    colnames(num) <- x$sseqid
    return(num)
  })
  
  asv2embed_short <- rbind.fill(asv2embed_list)
  asv2embed_short[is.na(asv2embed_short)] <- 0
  asv2embed_short <- t(apply(asv2embed_short, 1, function(x) return(x / sum(x))))
  rownames(asv2embed_short) <- names(asv2embed_list)
  return(asv2embed_short)
}

#This function runs blast to find the closest hits to each query sequence in the embedding database
#If blast has been run externally for increased speed, the function may take in the blast output table
#If this function is taking more than a few minutes to run, consider running BLAST externally
getBestHits <- function(blast_hits, id_thresh = 99){
  #save best hits per query ASV
  print("Filtering best hits from blast output")
  best_hits_list <- lapply(unique(blast_hits$qseqid), function(asv_name){
    tmp <- blast_hits[blast_hits$qseqid == asv_name, ]
    tmp <- tmp[tmp$evalue == min(tmp$evalue), ]
    tmp <- tmp[tmp$pident >= id_thresh, ]
    return(tmp)
  })
  best_hits <- do.call(rbind, best_hits_list)
  return(best_hits)
}

#For each query sequence, assigns the name of the closest embedding sequence. If there are multiple closest hits, splits the abundance
#of the query sequence evenly among all closest hits
transformSeqtab <- function(seqtab, best_hits){
  
  #sample by ASV
  print('converting ids')
  
  num_seqs_aligned <- sum(colnames(seqtab) %in% as.character(best_hits$qseqid))
  percent_seqs_aligned <- num_seqs_aligned / ncol(seqtab)
  cat("Number of sequences from query dataset that aligned: ", num_seqs_aligned, "\n")
  cat("Percent of sequences from query dataset that aligned: ", percent_seqs_aligned, "\n")
  
  if(num_seqs_aligned == 0){
    stop("No reads aligned to database. Check that column names of sequence table are ids, and match the queryid column of blast_hits")
  }
  seqtab <- seqtab[, colnames(seqtab) %in% as.character(best_hits$qseqid)]
  
  # Split counts among all best hits
  print('renaming ASVs by embedding database')
  asv2embed_short <- getHitTransformation(best_hits)
  seqtab <- seqtab[ , colnames(seqtab) %in% rownames(asv2embed_short)]
  asv2embed_short <- asv2embed_short[colnames(seqtab), ]
  colnames(seqtab) == rownames(asv2embed_short)
  
  seqtab_transformed <- as.matrix(seqtab) %*% as.matrix(asv2embed_short)
  return(list(seqtab_transformed, num_seqs_aligned, percent_seqs_aligned))
}



EmbedAsvTable <- function(seqtab, blast_hits, embedding_matrix, id_thresh = 99){
  best_hits = getBestHits(blast_hits = blast_hits, id_thresh = id_thresh)
  
  # There are just a hand ful of sequences that got thrown away during the embedding process, so are present in the best_hits hits column, but not the embedding matrix itself
  best_hits <- best_hits[best_hits$sseqid %in% rownames(embedding_matrix), ]
  tmp <- transformSeqtab(seqtab = seqtab, best_hits = best_hits)
  seqtab_transformed <- tmp[[1]]
  num_seqs_aligned <- tmp[[2]]
  percent_seqs_aligned <- tmp[[3]]
  seqtab_transformed <- seqtab_transformed[ , colnames(seqtab_transformed) %in% rownames(embedding_matrix)]
  embedding_matrix <- embedding_matrix[colnames(seqtab_transformed), ]
  embedded <- as.matrix(seqtab_transformed) %*% as.matrix(embedding_matrix)
  return(list("embedded_matrix" = embedded, "num_seqs_aligned" = num_seqs_aligned, "percent_seqs_aligned" = percent_seqs_aligned))
}
