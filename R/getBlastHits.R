getBlastHits <- function(blast_hits_file){
  blast_hits <- read.delim(blast_hits_file)
  return(blast_hits)
}
