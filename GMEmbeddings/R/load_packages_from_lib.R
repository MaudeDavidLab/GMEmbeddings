# library(Biostrings)
# library(rBLAST)
# library(seqinr)
# library(plyr)
# library(data.table)
#

# install_or_load_pack <- function(pack){
#   create.pkg <- pack[!(pack %in% installed.packages()[, "Package"])]
#   if (length(create.pkg))
#     install.packages(create.pkg, dependencies = TRUE)
#   sapply(pack, require, character.only = TRUE)
# }
#
#
# packages <- c("Biostrings", "rBlast",  "seqinr", "plyr", "data.table")
# install_or_load_pack(packages)

load_packages_from_lib <- function(){
  library(Biostrings)
  library(rBLAST)
  library(seqinr)
  library(plyr)
  library(data.table)
  library(pbapply)
}

load_packages_from_lib()
