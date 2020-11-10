#source("GME/R/embed_asv_table.R")
#runBlast(fasta_file = "GME/data/test.fasta", out_dir = "GME/test/")
#filterBlastHits(out_dir = "GME/test/")

#'A test function
#' @export
#' @param string The string you want output. Defaults to Hi
test <- function(string = 'Hi'){
  print(string)
}

#embedSeqtab(getExampleSeqtab(), getExampleFasta())
