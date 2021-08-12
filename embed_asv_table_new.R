
#' @export
getExampleSeqtab <- function(){
  #seqtab <- utils::read.csv("/Library/Frameworks/R.framework/Versions/4.0/Resources/library/GMEmbeddings/extdata/seqtab_test.csv")#######################
  seqtab <- utils::read.csv(seqtab)
  rownames(seqtab) <- seqtab$X
  rownames(seqtab) <- gsub("X", "", rownames(seqtab))
  seqtab <- seqtab[, 2:ncol(seqtab)]
  return(seqtab)
}

getExampleSeqtab_AG <- function(){
  seqtab <- data.table::fread("data/seqtab_final_filter.07.txt", header = F, sep = '\t')
  seqs <- scan("data/sequences_.07.txt", what = character())
  seqtab <- as.data.frame(seqtab)
  rownames(seqtab) <- seqtab$V1
  rownames(seqtab) <- gsub("X", "", rownames(seqtab))
  seqtab <- seqtab[, 2:ncol(seqtab)]
  colnames(seqtab) <- seqs
}

getExampleSeqtab_Halfvarson <- function(){
  seqtab <- data.table::fread("data/halfvarson/seqtab.txt", header = F, sep = '\t')
  fasta <- seqinr::read.fasta("data/halfvarson/sequences.fasta")
  seqs <- toupper(unlist(seqinr::getSequence(fasta, as.string = T)))
  seqtab <- as.data.frame(seqtab)
  rownames(seqtab) <- seqtab$V1
  rownames(seqtab) <- gsub("X", "", rownames(seqtab))
  seqtab <- seqtab[, 2:ncol(seqtab)]
  colnames(seqtab) <- seqs
  return(seqtab)
}


#' #' @export
# getExampleFasta <- function(){
#   fasta_file <- fasta_file_name
#   return(fasta_file)
# }

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
#"$blast_software_dir/blastn" -db "$blast_db/embedding_db_.07" -query "$data_dir/repseqs.fasta" -out "$out_dir/blast_hits.tsv"  -outfmt "6 qseqid sseqid qseq sseq evalue bitscore length pident"
#cat "$out_dir/blast_hits.tsv" | sort -k1,1 -k5,5g -k6,6nr | sort -u -k1,1 --merge > "$out_dir/best_hits.tsv"
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

getEmbeddingHits <- function(fasta_file, blast_hits_file = NA, id_thresh = 99, out_dir = "data/blast_output/"){
  if(is.na(blast_hits_file)){
    runBlast(fasta_file = fasta_file, out_dir = out_dir)
    blast_hits_file = file.path(out_dir, "blast_hits.tsv")
    blast_hits <- read.delim(blast_hits_file, header= TRUE, sep = "\t")
  }
  else{
    blast_hits <- read.delim(blast_hits_file, header = TRUE, sep = "\t")
  }
  return(getBestHits(blast_hits))
}


# getFastaDF <- function(fasta_file){
#   fasta <- read.fasta(fasta_file)
#   asv_ids <- names(fasta)
#   asv_seqs <- toupper(unlist(seqinr::getSequence(fasta, as.string = T)))
#   fasta_df <- data.frame(asv_ids, row.names = asv_seqs)
#   return(fasta_df)
# }

#For each query sequence, assigns the name of the closest embedding sequence. If there are multiple closest hits, splits the abundance
#of the query sequence evenly among all closest hits
transformSeqtab <- function(seqtab, fasta_file, best_hits){

  #sample by ASV
  print('converting ids')
  #run if column names are in full sequence form. Don't run if the column names are ASV
  #colnames(seqtab) <- convertIDs(ids = colnames(seqtab), from_id = "full_length_seq", to_id = "ASV_label", fasta_file) #fasta should contain all sequences in seqtab
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

embedSeqtab <- function(seqtab, fasta_file, best_hits, embedding_file_name){
  seqtab_transformed <- transformSeqtab(seqtab = seqtab, fasta_file = fasta_file, best_hits = best_hits)
  embedding_matrix <- getEmbeddingMatrix(embedding_file_name)
  embedding_matrix <- embedding_matrix[colnames(seqtab_transformed), ]
  embedded <- as.matrix(seqtab_transformed) %*% as.matrix(embedding_matrix)
  return(embedded)
}


convertIDs <- function(ids, from_id, to_id, fasta_file){
  fasta_list <- getFastaList(fasta_file)
  if(from_id == "full_length_seq" && to_id == "ASV_label"){
    fasta_df <- data.frame(fasta_list$embed_ids, row.names = fasta_list$seqs)
  }
  if(from_id == "ASV_label" && to_id == "full_length_seq"){
    fasta_df <- data.frame(fasta_list$seqs, row.names = fasta_list$embed_ids)
  }
  return(as.character(fasta_df[ids, 1]))
}

getFastaList <- function(fasta_file){
  fasta <- seqinr::read.fasta(fasta_file)
  embed_ids <- names(fasta)
  seqs <- toupper(unlist(getSequence(fasta, as.string = T)))
  return(list(embed_ids = embed_ids, seqs=seqs))
}


EmbedAsvTable <- function(seqtab, fasta_file_name, blast_hits, embedding_file_name){

  best_hits = getBestHits(blast_hits = blast_hits, id_thresh = 99)

  seqtab <- embedSeqtab(seqtab, fasta_file = fasta_file_name, best_hits = best_hits, embedding_file_name)

  seqtab
}

#Read in the dataframe blast_hits from that data folder
#blast_hits <- read.delim("data/blast_output/blast_hits.tsv")

#Set the variable embedding_file_name to be the file path (string) of the embedding transformation matrix
#embedding_file_name <- "data/embed_matrices/embed_.07_100dim.txt"

#Call you function
#EmbedAsvTable(blast_hits, embedding_file_name)
