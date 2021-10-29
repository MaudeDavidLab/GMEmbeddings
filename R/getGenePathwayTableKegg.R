library(KEGGREST)
library(pheatmap)
library(plotly)
library(htmlwidgets)

# Function imports the output of picrust to make an ASV x pathway table
# showing which pathway is present in which taxa (ASV)
# Does not need to be run, find output object asv_by_pathway.rds in inst/extdata/pathways
prepareASVbyPathwayTable <- function(){
  # Step 1. BLAST sequences in embedding database to KEGG database to assign a genome code to each
  # Run picrust on server at christine/gut_microbiome_embeddings/picrust/ to get an ASV by KO gene table
  
  asv_by_gene <- read.delim("inst/extdata/pathways/KO_predicted.tsv", sep = "\t", header = T, row.names = 1)
  
  
  # Step 2. Build pathway by KO gene table directly from kegg
  possible_pathways <- readRDS("inst/extdata/pathways/prokaryote_pathway_ids.rds")
  
  
  pathway_by_genes_list <- lapply(possible_pathways, function(pathway){
    
    possible_genes <- rep(0, ncol(asv_by_gene))
    names(possible_genes) <- colnames(asv_by_gene)
    
    tmp <- keggGet(pathway)
    genes <- names(tmp[[1]]$ORTHOLOGY)
    
    genes <- genes[genes %in% names(possible_genes)]
    possible_genes[genes] <- 1
    
    return(possible_genes)
  })
  
  pathway_by_genes <- do.call(rbind, pathway_by_genes_list)
  rownames(pathway_by_genes) <- possible_pathways
  
  
  # Step 3. Build pathway by ASV table using matrix multiplication
  
  pathway_by_asv <- pathway_by_genes %*% t(asv_by_gene)
  
  # Step 4. Rename ASVs with their full length sequence for clarity
  
  fasta <- read.table("inst/extdata/pathways/glove_emb_fullseq.fasta", quote="\"", comment.char="")
  headers <- gsub(">", "", fasta[seq(1, nrow(fasta), 2), ])
  seqs <- fasta[seq(2, nrow(fasta), 2), ]
  df <- data.frame(seqs, row.names = headers)
  colnames(pathway_by_asv) <- df[colnames(pathway_by_asv), ]
  
  # Step 5. Transform to asv by pathway
  asv_by_pathway <- t(pathway_by_asv)
  rownames(asv_by_pathway) <- colnames(pathway_by_asv)
  
  # Step 6. Save
  saveRDS(asv_by_pathway, "inst/extdata/pathways/asv_by_pathway.rds")
}



# Internal function. Take correlation between every column of the embedded table (sample x property) and every column of the sample x pathway table
getCorMat <- function(embedding_table, pathway_table){
  cor_list <- list()
  for(i in seq(1, ncol(embedding_table))){
    cor = apply(pathway_table, 2, function(pathway_vec) return(cor(pathway_vec, embedding_table[ ,i])))
    cor_list[[i]] <- cor
  }
  cor_mat <- data.frame(matrix(unlist(cor_list), byrow = T, nrow = length(cor_list)))
  colnames(cor_mat) <- colnames(pathway_table)
  rownames(cor_mat) <- colnames(embedding_table)
  return(cor_mat)
}

# Internal function. Return the pathway name that is most correlation with every dimension in the provided embed_table
getMaxCorr <- function(embed_table, pathway_table){
  cor_mat <- getCorMat(embed_table, pathway_table)
  max_corr_inx <- apply(cor_mat, 1, function(cor_vec){
    return(which(cor_vec == max(cor_vec))[1])
  })
  return(colnames(pathway_table)[max_corr_inx])
}


#' @export
getPathwayCorrelationMatrix <- function(embed){
  asv_by_pathway <- readRDS("inst/extdata/pathways/asv_by_pathway.rds")
  
  # Align ASVs between embedding transformation matrix and asv x pathway table
  asvs <- intersect(rownames(embed), rownames(asv_by_pathway))
  embed <- embed[asvs, ]
  asv_by_pathway <- asv_by_pathway[asvs, ]
  
  corrmat <- getCorMat(embed, asv_by_pathway)
  corrmat <- corrmat[, !is.na(colSums(corrmat))]
  
  return(corrmat)
}

#' @export
getHeatmap <- function(corrmat){
  full_names <- lapply(colnames(corrmat), function(x){ return(keggGet(x)[[1]]$NAME)})
  plot_names <- paste(colnames(corrmat), full_names, sep = ": ")
  
  fig <- plot_ly(z = as.matrix(corrmat), type = "heatmap", y = rownames(corrmat), x = plot_names)
  return(fig)
}

makeHeatmapFigures <- function(){
  # 50 dimensions
  embed <- read.table("inst/extdata/embedding_transformation_matrices/glove_emb_fullseq_50.txt", quote="\"", comment.char="", row.names = 1)
  corrmat <- getPathwayCorrelationMatrix(embed)
  fig <- getHeatmap(corrmat)
  saveWidget(fig, "results/pathways/correlation_50dimensions_pathways.html")
  
  # 100 dimensions
  embed <- read.table("inst/extdata/embedding_transformation_matrices/glove_emb_fullseq_100.txt", quote="\"", comment.char="", row.names = 1)
  corrmat <- getPathwayCorrelationMatrix(embed)
  fig <- getHeatmap(corrmat)
  saveWidget(fig, "results/pathways/correlation_100dimensions_pathways.html")
  
  # 250 dimensions
  embed <- read.table("inst/extdata/embedding_transformation_matrices/glove_emb_fullseq_250.txt", quote="\"", comment.char="", row.names = 1)
  corrmat <- getPathwayCorrelationMatrix(embed)
  fig <- getHeatmap(corrmat)
  saveWidget(fig, "results/pathways/correlation_250dimensions_pathways.html")
  
}
