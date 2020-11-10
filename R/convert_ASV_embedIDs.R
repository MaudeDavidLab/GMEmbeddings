getFastaList <- function(fasta_file){
  fasta <- seqinr::read.fasta(fasta_file)
  embed_ids <- names(fasta)
  seqs <- toupper(unlist(getSequence(fasta, as.string = T)))
  return(list(embed_ids = embed_ids, seqs=seqs))

}


#'Convert ids from full ASV labels to embed_ids as assigned in label_qualvec_transform_mat which outputs embed/data/seqs_.07_embed linking the two id types
#' @export
#' @param
convertIDs <- function(ids, from_id, to_id, fasta_file= "data/embed_matrices/seqs_.07.fasta"){
  fasta_list <- getFastaList(fasta_file)

  if(from_id == "ASV" && to_id == "embedIDs"){
    fasta_df <- data.frame(fasta_list$embed_ids, row.names = fasta_list$seqs)
  }
  if(from_id == "embedIDs" && to_id == "ASV"){
    fasta_df <- data.frame(fasta_list$seqs, row.names = fasta_list$embed_ids)
  }
  return(as.character(fasta_df[ids, 1]))
}

