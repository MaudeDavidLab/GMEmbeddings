getEmbeddingMatrix <- function(embedding_file_name){
  qual_vecs <- read.csv(embedding_file_name, row.names=1, sep="")
  return(qual_vecs)
}
