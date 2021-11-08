#' #' @export
#' getExampleSeqtab <- function(){
#'   seqtab <- utils::read.csv(seqtab)
#'   rownames(seqtab) <- seqtab$X
#'   rownames(seqtab) <- gsub("X", "", rownames(seqtab))
#'   seqtab <- seqtab[, 2:ncol(seqtab)]
#'   return(seqtab)
#' }
#' 
#' getExampleSeqtab_AG <- function(){
#'   seqtab <- data.table::fread("data/seqtab_final_filter.07.txt", header = F, sep = '\t')
#'   seqs <- scan("data/sequences_.07.txt", what = character())
#'   seqtab <- as.data.frame(seqtab)
#'   rownames(seqtab) <- seqtab$V1
#'   rownames(seqtab) <- gsub("X", "", rownames(seqtab))
#'   seqtab <- seqtab[, 2:ncol(seqtab)]
#'   colnames(seqtab) <- seqs
#' }
#' 
#' getExampleSeqtab_Halfvarson <- function(){
#'   seqtab <- data.table::fread("data/halfvarson/seqtab.txt", header = F, sep = '\t')
#'   fasta <- seqinr::read.fasta("data/halfvarson/sequences.fasta")
#'   seqs <- toupper(unlist(seqinr::getSequence(fasta, as.string = T)))
#'   seqtab <- as.data.frame(seqtab)
#'   rownames(seqtab) <- seqtab$V1
#'   rownames(seqtab) <- gsub("X", "", rownames(seqtab))
#'   seqtab <- seqtab[, 2:ncol(seqtab)]
#'   colnames(seqtab) <- seqs
#'   return(seqtab)
#' }


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

# getEmbeddingHits <- function(fasta_file, blast_hits_file = NA, id_thresh = 99, out_dir = "data/blast_output/"){
#   if(is.na(blast_hits_file)){
#     runBlast(fasta_file = fasta_file, out_dir = out_dir)
#     blast_hits_file = file.path(out_dir, "blast_hits.tsv")
#     blast_hits <- read.delim(blast_hits_file, header= TRUE, sep = "\t")
#   }
#   else{
#     blast_hits <- read.delim(blast_hits_file, header = TRUE, sep = "\t")
#   }
#   return(getBestHits(blast_hits))
# }


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
  return(seqtab_transformed)
}

###################
# convertIDs <- function(ids, from_id, to_id, fasta_file){
#   fasta_list <- getFastaList(fasta_file)
#   if(from_id == "full_length_seq" && to_id == "ASV_label"){
#     fasta_df <- data.frame(fasta_list$embed_ids, row.names = fasta_list$seqs)
#   }
#   if(from_id == "ASV_label" && to_id == "full_length_seq"){
#     fasta_df <- data.frame(fasta_list$seqs, row.names = fasta_list$embed_ids)
#   }
#   return(as.character(fasta_df[ids, 1]))
# }
# 
# getFastaList <- function(fasta_file){
#   fasta <- seqinr::read.fasta(fasta_file)
#   embed_ids <- names(fasta)
#   seqs <- toupper(unlist(getSequence(fasta, as.string = T)))
#   return(list(embed_ids = embed_ids, seqs=seqs))
# }


EmbedAsvTable <- function(seqtab, blast_hits, embedding_matrix, id_thresh = 99){
  best_hits = getBestHits(blast_hits = blast_hits, id_thresh = id_thresh)
  
  # There are just a hand ful of sequences that got thrown away during the embedding process, so are present in the best_hits hits column, but not the embedding matrix itself
  best_hits <- best_hits[best_hits$sseqid %in% rownames(embedding_matrix), ]
  seqtab_transformed <- transformSeqtab(seqtab = seqtab, best_hits = best_hits)
  embedding_matrix <- embedding_matrix[colnames(seqtab_transformed), ]
  embedded <- as.matrix(seqtab_transformed) %*% as.matrix(embedding_matrix)
  seqtab
}

#Read in the dataframe blast_hits from that data folder
#blast_hits <- read.delim("data/blast_output/blast_hits.tsv")

#Set the variable embedding_file_name to be the file path (string) of the embedding transformation matrix
#embedding_file_name <- "data/embed_matrices/embed_.07_100dim.txt"

#Call you function
#EmbedAsvTable(blast_hits, embedding_file_name)