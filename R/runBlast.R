#' @export
runBlast <- function(fasta_file, out_dir, database =  "inst/extdata/blastdb/embedding_db_.07"){

  ###########################
  ### Read in query seqs ####
  ###########################
  print("Reading in query fasta")
  query_seqs <- Biostrings::readDNAStringSet(fasta_file)

  #######################
  ## Read in blast db ###
  #######################
  print("Reading blast database")
  bl <- rBLAST::blast(db= database)

  ########################
  ### Run Blast ##########
  ########################
  print("Running Blast")
  blast_hits <- predict(bl, query_seqs, custom_format = "qseqid sseqid qseq sseq evalue bitscore length pident" )
  write.table(blast_hits, paste(out_dir, "/blast_hits.tsv", sep=""),
              sep = "\t", quote = FALSE, row.names = F)

}

#' @export
filterBlastHits <- function(out_dir){
  ### Put blast hits in the following order:
    #1. Alphabetically by query sequence id
    #2. Lowest e-value hit
    #3. Highest length
    #4. Highest percent identity
  blast_hits_file = paste(out_dir, "/blast_hits.tsv", sep="")
  best_hits_file <- paste(out_dir, "/best_hits.tsv", sep =)

  blast_hits_ord <- blast_hits[
    #id, evalue, length, perc_id
    order(as.character(blast_hits[ , 1]), blast_hits[ , 5], rev(blast_hits[ , 7]), rev(blast_hits[ , 8])),
  ]

}
#ec5f31d6-8daa-4b53-b38d-3603fcaa2f14ec5f31d6-8daa-4b53-b38d-3603fcaa2f14
