
library(reshape2)
library(ggplot2)
library(pheatmap)
library(tibble)
library(dplyr)


# define starting variables:
project <- "hgsoc_repeats"
expName <- "exp9"
STypes <- c("arPT", "erPT", "mrPT", "msST",  "prPT", "rcAF","rfPT")
#HGSOCtypes <- c("acquired_resistant", "extreme_responder", "multiple_responder", primary_resistant", "recurrent", "primary_refractory")
annot <- "custom3"

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
#homeDir <- "/Users/jamestorpy/Documents/Garvan/phd"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
                     expName, "/HGSOC_vs_HGSOC_bowtell_cats/")
plotDir <- paste0(resultsDir, "/R/", expName,
                  "/plots/DEplots/HGSOC_vs_HGSOC_bowtell_cats/")


system(paste0("mkdir -p ", plotDir))


### 1. Load in inputs and arrange into dfs ###

allGene <- readRDS(file=paste0(RobjectDir, "/", annot, "_DEepiGenes.rds"))

allNames <- unlist(lapply(allGene, function(x) {
  return(paste0(strsplit(names(x)[1], "_")[[1]][1], "_as_ctl"))
}))

allGene <- lapply(allGene, function(x) {
  return(do.call("rbind", x))
})

names(allGene) <- allNames


### 2. Create logFC and FDR heatmaps:

# fetch sample number to add to colnames:
sample_nos <- readRDS(paste0(RobjectDir, "/sample_no_per_cat.rds"))

library(grid)

draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

for (i in 1:length(allGene)) {
  
  fcDF <- as.data.frame(allGene[[i]]$logFC, row.names = rownames(allGene[[i]]))
  colnames(fcDF) <- "logFC"
  fcDF$sample <- gsub("\\..*$", "", rownames(fcDF))
  fcDF$gene_id <- gsub("^.*\\.", "", rownames(fcDF))
  fcDF <- fcDF %>%
    dcast(gene_id ~ sample, value.var = "logFC") %>%
    column_to_rownames("gene_id")
  
  fdrDF <- as.data.frame(allGene[[i]]$FDR, row.names = rownames(allGene[[i]]))
  colnames(fdrDF) <- "FDR"
  fdrDF$sample <- gsub("\\..*$", "", rownames(fdrDF))
  fdrDF$gene_id <- gsub("^.*\\.", "", rownames(fdrDF))
  fdrDF <- fdrDF %>%
    dcast(gene_id ~ sample, value.var = "FDR") %>%
    column_to_rownames("gene_id")
  
  
  # resize margins for plots:
  par(mar=c(4,4,4,4))
  
  pheatmap(fcDF, color = colorRampPalette(c("#08519C", "white", "firebrick3"))(50),
           display_numbers = as.matrix(ifelse(fdrDF < 0.1, "*", "")), fontsize = 8, scale="column",
           cluster_cols=F)
  
  # FC pheatmap with column clustering including 'both' DE:
  log2FC <- pheatmap(fcDF, color = colorRampPalette(c("#08519C", "white", "firebrick3"))(50),
                     display_numbers = as.matrix(ifelse(fdrDF < 0.1, "*", "")), fontsize = 8,
                     scale="column")
  
  pdf(file=paste0(plotDir, "/log2FC_HGSOC_heatmap_epi_", names(allGene)[i], ".pdf"), width=10, height=10)
  pheatmap(fcDF, color = colorRampPalette(c("#08519C", "white", "firebrick3"))(50),
           display_numbers = as.matrix(ifelse(fdrDF < 0.1, "*", "")), fontsize = 8,
           scale="column")
  dev.off()
  
  
  # FDR pheatmap with column clustering including 'both' DE:
  logfdrDF <- log10(fdrDF)
  
  pdf(file=paste0(plotDir, "/log10_FDR_HGSOC_heatmap_epi_", names(allGene)[i], ".pdf"), width=10, height=10)
  pheatmap(logfdrDF, color = colorRampPalette(c("#08519C", "white"))(50), 
           display_numbers = as.matrix(ifelse(fdrDF < 0.1, "*", "")), fontsize = 8)
  dev.off()
  
  
  # split fc DF into more manageble parts:
  #spl_no <- nrow(fcDF)/8
  
  #for (i in 1:8) {
  #  print(i)
  #  vec <- (rep(i, spl_no))
  #  if (i==1) {
  #    spl_vec <- c(vec)
  #  } else {
  #    spl_vec <- append(spl_vec, vec)
  #  }
  #}
  
  #spl_fc <- split(fcDF, spl_vec)
  #spl_fdr <- split(fdrDF, spl_vec)
  
  #for (j in 1:length(spl_fc)) {
  #  print(j)
  #  pdf(file=paste0(plotDir, "/log2FC_HGSOC_heatmap_", names(allGene)[j], ".pdf"), width=10, height=10)
  #  pheatmap(spl_fc[[j]], color = colorRampPalette(c("#08519C", "white", "firebrick3"))(50),
  #           display_numbers = as.matrix(ifelse(spl_fdr[[j]] < 0.1, "*", "")), fontsize = 8,
  #           scale="column")
  #  dev.off()
  #}
  
  interesting_fc <- fcDF %>%
    rownames_to_column('gene_id') %>%
    dplyr::filter(rowSums(fdrDF < 0.1) >= 1) %>%
    column_to_rownames('gene_id')
  
  interesting_fdr <-  fdrDF[rownames(interesting_fc),]
  
  if (nrow(interesting_fc) > 1) {
    pdf(file=paste0(plotDir, "/log2FC_HGSOC_heatmap_interesting_epi_", names(allGene)[i], ".pdf"), width=10, height=10)
    pheatmap(interesting_fc, color = colorRampPalette(c("#08519C", "white", "firebrick3"))(50),
             display_numbers = as.matrix(ifelse(interesting_fdr < 0.1, "*", "")), fontsize = 8,
             scale="column")
    dev.off()
  }

}

