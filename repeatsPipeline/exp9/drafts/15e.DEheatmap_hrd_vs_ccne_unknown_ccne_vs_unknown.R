library(reshape2)
library(ggplot2)
library(pheatmap)
library(tibble)
library(dplyr)


# define starting variables:
project <- "hgsoc_repeats"
expName <- "exp9"
sTypes <- c("bothDrivers", "CCNEAmp", "HRD", "unknownDrivers")
annot <- "custom3"

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
#homeDir <- "/Users/jamestorpy/Documents/Garvan/phd"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
                     expName, "/hrd_vs_ccne_unknown_and_ccne_vs_hrd_unknown/")
plotDir <- paste0(resultsDir, "/R/", expName, "/plots/DEplots/hrd_vs_ccne_unknown_ccne_vs_unknown/")

system(paste0("mkdir -p ", plotDir))


### 1. Load in inputs and arrange into df ###

allGene <- readRDS(file=paste0(RobjectDir, "/", annot, "_DEepiGenes.rds"))


### 2. Fetch vector of all, top and str genes: ###

allName <- rownames(allGene[[1]])


### 3. Create logFC and FDR heatmaps:

# fetch sample number to add to colnames:
sample_nos <- readRDS(paste0(RobjectDir, "/sample_no_per_cat.rds"))

Names <- names(allGene)
for (i in 1:length(Names)) {
  names <- c(strsplit(Names[i], "_")[[1]][1], strsplit(Names[i], "_")[[1]][3])
  nos <- c(sample_nos[grep(names[1], names(sample_nos))], sample_nos[grep(names[2], names(sample_nos))])
  if (i==1) {
    newNames <- c(paste0(names[1], "_vs_", names[2], " (n=", nos[1], ",", nos[2], ")"))
  } else {
    newNames[i] <- paste0(names[1], "_vs_", names[2], " (n=", nos[1], ",", nos[2], ")")
  }
}

names(allGene) <- newNames
  
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

# resize margins for plots:
par(mar=c(4,4,4,4))

library(grid)

draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

pheatmap(fcDF, color = colorRampPalette(c("#08519C", "white", "firebrick3"))(50),
         display_numbers = as.matrix(ifelse(fdrDF < 0.1, "*", "")), fontsize = 8, scale="column",
         cluster_cols=F)

# FC pheatmap with column clustering including 'both' DE:
log2FC <- pheatmap(fcDF, color = colorRampPalette(c("#08519C", "white", "firebrick3"))(50),
                   display_numbers = as.matrix(ifelse(fdrDF < 0.1, "*", "")), fontsize = 8,
                   scale="column")

pdf(file=paste0(plotDir, "/log2FC_HGSOC_heatmap_epi_.pdf"), width=10, height=10)
pheatmap(fcDF, color = colorRampPalette(c("#08519C", "white", "firebrick3"))(50),
         display_numbers = as.matrix(ifelse(fdrDF < 0.1, "*", "")), fontsize = 8,
         scale="column")
dev.off()


# FDR pheatmap with column clustering including 'both' DE:
logfdrDF <- log10(fdrDF)

pdf(file=paste0(plotDir, "/log10_FDR_HGSOC_heatmap_epi_.pdf"), width=10, height=10)
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
  pdf(file=paste0(plotDir, "/log2FC_HGSOC_heatmap_interesting_epi_.pdf"), width=10, height=10)
  pheatmap(interesting_fc, color = colorRampPalette(c("#08519C", "white", "firebrick3"))(50),
           display_numbers = as.matrix(ifelse(interesting_fdr < 0.1, "*", "")), fontsize = 8,
           scale="column", cluster_cols = F)
  dev.off()
}




