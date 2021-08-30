
#' @export
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
  #run if column names are in full sequence form. Don't run if the column names are ASV
  seqtab <- seqtab[, colnames(seqtab) %in% as.character(best_hits$qseqid)]
  
  # Split counts among all best hits
  print('renaming ASVs by embedding database')
  asv2embed_short <- getHitTransformation(best_hits)
  seqtab <- seqtab[ , colnames(seqtab) %in% rownames(asv2embed_short)]
  asv2embed_short <- asv2embed_short[colnames(seqtab), ]
  colnames(seqtab) == rownames(asv2embed_short)
  
  seqtab_transformed <- as.matrix(seqtab) %*% as.matrix(asv2embed_short)
  return(seqtab_transformed)
}

###################

embedSeqtab <- function(seqtab, best_hits, embedding_file_name){
  seqtab_transformed <- seqtab
    # transformSeqtab(seqtab = seqtab, best_hits = best_hits)
  embedding_matrix <- getEmbeddingMatrix(embedding_file_name)
  embedding_matrix <- embedding_matrix[colnames(seqtab_transformed), ]
  embedded <- as.matrix(seqtab_transformed) %*% as.matrix(embedding_matrix)
  return(embedded)
}

EmbedAsvTable <- function(seqtab, blast_hits, embedding_file_name){
  best_hits = getBestHits(blast_hits = blast_hits, id_thresh = 99)
  seqtab <- embedSeqtab(seqtab, best_hits = best_hits, embedding_file_name)
  seqtab
}