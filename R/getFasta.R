### Ahh! Why this no work?
getExampleFasta <- function(fasta_file_name){
  fasta_file <- fasta_file_name
  return(fasta_file)
}

getFastaDF <- function(fasta_file){
  fasta <- read.fasta(fasta_file)
  asv_ids <- names(fasta)
  asv_seqs <- toupper(unlist(seqinr::getSequence(fasta, as.string = T)))
  fasta_df <- data.frame(asv_ids, row.names = asv_seqs)
  return(fasta_df)
}

getFasta <- function(fasta_file_name){
  fasta_file <- getExampleFasta(fasta_file_name)
  fasta_df <- getFastaDF(fasta_file)
  return(fasta_df)
}
