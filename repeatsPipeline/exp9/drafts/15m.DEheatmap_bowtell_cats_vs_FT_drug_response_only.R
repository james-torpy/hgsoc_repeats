
library(reshape2)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tibble)


# define starting variables:
project <- "hgsoc_repeats"
expName <- "exp9"
compTypes <- c("extreme_responders_vs_FT", "multiple_responders_vs_FT", 
               "acquired_resistant_vs_FT", "primary_resistant_vs_FT", 
               "refractory_vs_FT", "metastatic_vs_FT")
annot <- "custom3"

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
#homeDir <- "/Users/jamestorpy/Documents/Garvan/phd"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
#RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/", expName, "/test/", groupage, "/")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
                     expName, "/drug_response_only/")
plotDir <- paste0(resultsDir, "/R/", expName,
                  "/plots/DEplots/drug_response_only/")

system(paste0("mkdir -p ", plotDir))


### 1. Load in inputs and arrange into df ###

allGene <- readRDS(file=paste0(RobjectDir, "/", annot, "_DEreps.rds"))

### 2. Fetch vector of all, top and str genes: ###

allName <- rownames(allGene[[1]])


### 3. Create logFC and FDR heatmaps:

# fetch sample number to add to colnames:
sample_nos <- readRDS(paste0(RobjectDir, "/sample_no_per_cat.rds"))

  
allFC <- do.call("rbind", allGene)
allFC$sample <- gsub("\\..*$", "", rownames(allFC))
allFC$repeat_id <- gsub("^.*\\.*\\.", "", rownames(allFC))
allFC <- subset(allFC, select=c(logFC, sample, repeat_id))
allFC[,4] <- allFC$logFC
allFC <- subset(allFC, select=-logFC)
names(allFC)[3] <- "logFC"

fcDF <- dcast(allFC, repeat_id ~ sample)
rownames(fcDF) <- fcDF$repeat_id
fcDF <- subset(fcDF, select=-repeat_id)
fcDF <- fcDF[,compTypes]

allFDR <- do.call("rbind", allGene)
allFDR$sample <- gsub("\\..*$", "", rownames(allFDR))
allFDR$repeat_id <- gsub("^.*\\.*\\.", "", rownames(allFDR))
allFDR <- subset(allFDR, select=c(FDR, sample, repeat_id))
allFDR[,4] <- allFDR$FDR
allFDR <- subset(allFDR, select=-FDR)
names(allFDR)[3] <- "FDR"

fdrDF <- dcast(allFDR, repeat_id ~ sample)
rownames(fdrDF) <- fdrDF$repeat_id
fdrDF <- subset(fdrDF, select=-repeat_id)
fdrDF <- fdrDF[,compTypes]

# resize margins for plots:
par(mar=c(4,4,4,4))

pheatmap(fcDF, color = colorRampPalette(c("#08519C", "white", "firebrick3"))(50),
         display_numbers = as.matrix(ifelse(fdrDF < 0.05, "*", "")), fontsize = 8, scale="column",
         cluster_cols=F)

library(grid)

draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

# FC pheatmap with column clustering including 'both' DE:
log2FC <- pheatmap(fcDF, color = colorRampPalette(c("#08519C", "white", "firebrick3"))(50),
                   display_numbers = as.matrix(ifelse(fdrDF < 0.05, "*", "")), fontsize = 8,
                   scale="column", cluster_cols = F)

pdf(file=paste0(plotDir, "/log2FC_HGSOC_heatmap_with_both.pdf"), width=10, height=10)
pheatmap(fcDF, color = colorRampPalette(c("#08519C", "white", "firebrick3"))(50),
         display_numbers = as.matrix(ifelse(fdrDF < 0.05, "*", "")), fontsize = 8,
         scale="column", cluster_cols = F)
dev.off()


# FDR pheatmap with column clustering including 'both' DE:
logfdrDF <- log10(fdrDF)

pdf(file=paste0(plotDir, "/log10_FDR_HGSOC_heatmap_with_both.pdf"), width=10, height=10)
pheatmap(logfdrDF, color = colorRampPalette(c("#08519C", "white"))(50), 
         display_numbers = as.matrix(ifelse(fdrDF < 0.1, "*", "")), fontsize = 8, cluster_cols = F)
dev.off()
 


interesting_fc <- fcDF %>%
  rownames_to_column('gene_id') %>%
  dplyr::filter(rowSums(fcDF > 2.3 | fcDF < -2.3) >= 1) %>%
  column_to_rownames('gene_id')

# filter out rows of interesting_fc with no samples having fdr < 0.05
sig_fdr <- rownames(fdrDF %>%
  rownames_to_column('gene_id') %>%
  dplyr::filter(rowSums(fdrDF < 0.05) >= 1) %>%
  column_to_rownames('gene_id'))
interesting_fc <- interesting_fc[rownames(interesting_fc) %in% sig_fdr,]

interesting_fdr <-  fdrDF[rownames(interesting_fc),]

pdf(file=paste0(plotDir, "/log2FC_HGSOC_heatmap_interesting.pdf"), width=10, height=10)
pheatmap(interesting_fc, color = colorRampPalette(c("#08519C", "white", "firebrick3"))(50),
         display_numbers = as.matrix(ifelse(interesting_fdr < 0.05, "*", "")), fontsize = 8,
         scale="column", cluster_cols = F)
dev.off()
