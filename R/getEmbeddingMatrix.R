getEmbeddingMatrix <- function(embedding_file_name){
  embedding_matrix <- read.csv(embedding_file_name, row.names=1, sep="")
  return(embedding_matrix)
}
