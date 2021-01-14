library(GMEmbeddings)
library(plyr)

blast_hits <- read.delim("C:/Users/ctata/Documents/Lab/gut_microbiome_embeddings_package//data/blast_output/blast_hits.tsv")
embedding_file_name <- "C:/Users/ctata/Documents/Lab/gut_microbiome_embeddings_package//data/embed_matrices/embed_.07_100dim.txt"
seqtab <- getExampleSeqtab()
fasta_file <- getExampleFasta()

tmp <- EmbedAsvTable(seqtab, fasta_file, blast_hits, embedding_file_name)
tmp2 <- EmbedAsvTable(seqtab, fasta_file, blast_hits, embedding_file_name)
