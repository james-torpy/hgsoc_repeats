library(reshape2)
library(ggplot2)
library(pheatmap)
library(tibble)
library(dplyr)


# define starting variables:
project <- "hgsoc_repeats"
expName <- "exp9"
Type <- "custom3"
descrip <- "htseq_HGSOC_unknown_vs_HRD_CCNEamp"

# specify what combination of repeat genes (repeats), epigenetic modulators (epiMods),
# RNAi genes (RNAi) and protein-coding genes (pCoding) should contribute to the results:
resultTypes <- c("repeats")

# specify what FDR and log2 fold change thresholds were used:
FDRthresh <- 0.05
FCthresh <- 0

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
#homeDir <- "/Users/jamestorpy/Documents/Garvan/phd"
projectDir <- paste0(homeDir, "/projects/", project, "/RNA-seq/")
resultsDir <- paste0(projectDir, "/results")
RobjectDir <- paste0(projectDir, "/Robjects/",
                     expName, "/", descrip, "/")
plotDir <- paste0(resultsDir, "/R/", expName,
                  "/plots/DEplots/", descrip, "/")

system(paste0("mkdir -p ", RobjectDir))
system(paste0("mkdir -p ", plotDir))


### 1. Load in inputs ###

# load each comparison:
comparisons <- list(readRDS(file=paste0(projectDir, "/Robjects/exp9/htseq_HGSOC_HRD_unknown_vs_CCNEamp/", Type, "_DEsigReps.rds")))
comparisons[[2]] <- readRDS(file=paste0(projectDir, "/Robjects/exp9/htseq_HGSOC_CCNEamp_unknown_vs_HRD/", Type, "_DEsigReps.rds"))

# fetch unknown_driver DFs and add to list:
for (n in 1:length(comparisons)) {
  ind <- grep("unknown_driver", names(comparisons[[n]]))
  if (n==1) {
    allGene <- list(comparisons[[n]][[ind]])
    names(allGene)[n] <- names(comparisons[[n]])[ind]
  } else {
    allGene[[n]] <- comparisons[[n]][[ind]]
    names(allGene)[n] <- names(comparisons[[n]])[ind]
  }
}

# fetch sample number to add to colnames:
sample_nos <- readRDS(paste0(projectDir, "/Robjects/exp9/htseq_HGSOC_HRD_unknown_vs_CCNEamp//sample_no_per_cat.rds"))

# add sample numbers per group to names:
Names <- names(allGene)
for (i in 1:length(Names)) {
  names <- c(strsplit(Names[i], "_vs_")[[1]][1], strsplit(Names[i], "_vs_")[[1]][2])
  nos <- c(sample_nos[grep(names[1], names(sample_nos))], sample_nos[grep(names[2], names(sample_nos))])
  if (i==1) {
    newNames <- c(paste0(names[1], "_vs_", names[2], " (n=", nos[1], ",", nos[2], ")"))
  } else {
    newNames[i] <- paste0(names[1], "_vs_", names[2], " (n=", nos[1], ",", nos[2], ")")
  }
}
names(allGene) <- newNames


### 2. Split into logFC and FDR dfs:

allFC <- as.data.frame(do.call("rbind", allGene))
allFC$sample <- gsub("\\..*$", "", rownames(allFC))
allFC$repeat_id <- gsub("^.*\\.*\\.", "", rownames(allFC))
allFC <- subset(allFC, select=c(logFC, sample, repeat_id))
allFC[,4] <- allFC$logFC
allFC <- subset(allFC, select=-logFC)
names(allFC)[3] <- "logFC"

fcDF <- dcast(allFC, repeat_id ~ sample)
rownames(fcDF) <- fcDF$repeat_id
fcDF <- subset(fcDF, select=-repeat_id)
fcDF <- na.omit(fcDF)

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
fdrDF <- na.omit(fdrDF)


### 3. Create logFC and FDR heatmaps:

##### fix below this point #####

# resize margins for plots:
par(mar=c(4,4,4,4))

library(grid)

draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust =  0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

# create FC pheatmap without column clustering:
#pdf(file=paste0(plotDir, "/log2FC_HGSOC_heatmap_no_colclustering_", names( allGene)[j], ".pdf"), width=10, height=10)
#pheatmap(fcDF, color = colorRampPalette(c("#08519C", "white",  "firebrick3"))(50),
#         display_numbers = as.matrix(ifelse(fdrDF < 0.1, "*", "")),  fontsize = 8,
#         cluster_cols=F)
#dev.off()

# create FC pheatmap with column clustering:
#pdf(file=paste0(plotDir, "/log2FC_HGSOC_heatmap_colclustering_", names(  allGene)[j], ".pdf"), width=10, height=10)
#pheatmap(fcDF, color = colorRampPalette(c("#08519C", "white",  "firebrick3"))(50),
#         display_numbers = as.matrix(ifelse(fdrDF < 0.1, "*", "")),  fontsize = 8,
#         scale="column")
#dev.off()

# create FDR pheatmap without column clustering:
logfdrDF <- log10(fdrDF)

#pdf(file=paste0(plotDir, "/log10_FDR_HGSOC_heatmap_no_colclustering_", names (allGene)[j], ".pdf"), width=10, height=10)
#pheatmap(logfdrDF, color = colorRampPalette(c("#08519C", "white"))(50), 
#         display_numbers = as.matrix(ifelse(fdrDF < 0.1, "*", "")),  fontsize = 8, 
#         scale="column")
#dev.off()

# create FDR pheatmap with column clustering:
#pdf(file=paste0(plotDir, "/log10_FDR_HGSOC_heatmap_colclustering_", names( allGene)[j], ".pdf"), width=10, height=10)
#pheatmap(logfdrDF, color = colorRampPalette(c("#08519C", "white"))(50), 
#         display_numbers = as.matrix(ifelse(fdrDF < 0.1, "*", "")),  fontsize = 8)
#dev.off()

# isolate most significant genes and create FC heatmap:
interesting_fc <- fcDF %>%
  rownames_to_column('gene_id') %>%
  dplyr::filter(rowSums(fdrDF < FDRthresh & fcDF > FCthresh|fcDF <  -FCthresh) >= 1) %>%
  column_to_rownames('gene_id')

# create fdr for interesting repeats:
interesting_fdr <-  fdrDF[rownames(interesting_fc),]

# save interesting repeats:
saveRDS(rownames(interesting_fc), file=paste0(RobjectDir, "/interesting_fc.rds"))


if (nrow(interesting_fc) > 1) {
  # FC without clustering:
  #pdf(file=paste0(plotDir, "/log2FC_HGSOC_heatmap_FDR", FDRthresh,   "_colclustering.pdf"), width=10, height=10)
  #pheatmap(interesting_fc, color = colorRampPalette(c("#08519C", "white",  "firebrick3"))(50),
  #         display_numbers = as.matrix(ifelse(interesting_fdr < 0.1, "*",  "")), fontsize = 8,
  #         cluster_cols = F)
  #dev.off()
  
  # FC with clustering:
  #pdf(file=paste0(plotDir, "/log2FC_HGSOC_heatmap_FDR", FDRthresh,   "_no_colclustering.pdf"), width=10, height=10)
  #pheatmap(interesting_fc, color = colorRampPalette(c("#08519C", "white",  "firebrick3"))(50),
  #         display_numbers = as.matrix(ifelse(interesting_fdr < 0.1, "*",  "")), fontsize = 8,
  #         cluster_cols = T)
  #dev.off()
  
  # FC with row clustering, manual breaks:
  interesting_fc <- interesting_fc[, -1][c(2, 1, 3)]
  interesting_fdr <- interesting_fdr[, -1][c(2, 1, 3)]
  
  paletteLength <- 50
  myBreaks <- c(seq(min(fcDF), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(fcDF)/paletteLength, max(fcDF), length.out=floor( paletteLength/2)))
  
  pdf(file=paste0(plotDir, "/log2FC_HGSOC_heatmap_FDR", FDRthresh, "_FC",   FCthresh, "_colclustering_", 
                  ".pdf"), width=10, height=20)
  pheatmap(interesting_fc, color = colorRampPalette(c("#08519C", "white"))(50),
           display_numbers = as.matrix(ifelse(interesting_fdr < 0.1, "*",   "")), fontsize = 8,
           cluster_cols = F)
  dev.off()
  
  
  
  

#  pdf(file=paste0(plotDir, "/log2FC_HGSOC_sats_heatmap_FDR", FDRthresh, "_FC",   FCthresh, "_colclustering_", 
#                  ".pdf"), width=10, height=20)
#  pheatmap(int_sats, color = colorRampPalette(c("#08519C", "white"))(50),
#           display_numbers = as.matrix(ifelse(int_sat_fdr < 0.1, "*",   "")), fontsize = 8,
#           cluster_cols = T)
#  dev.off()
  
  
  
  
  # FDR without clustering:
  #log_interesting_fdr <- log10(interesting_fdr)
  
  #pdf(file=paste0(plotDir, "/log10_FDR_HGSOC_heatmapFDR", FDRthresh,   "_colclustering_", names(allGene)[j], ".pdf"), width=10, height=10)
  #pheatmap(logfdrDF, color = colorRampPalette(c("#08519C", "white"))(50), 
  #         display_numbers = as.matrix(ifelse(fdrDF < 0.1, "*", "")),  fontsize = 8,
  #         cluster_cols = F)
  #dev.off()
  
  # FDR with clustering:
  #pdf(file=paste0(plotDir, "/log10_FDR_HGSOC_heatmapFDR", FDRthresh, "_FC",   FCthresh, "colclustering_", names(allGene)[j], ".pdf"), width=10,   height=10)
  #pheatmap(logfdrDF, color = colorRampPalette(c("#08519C", "white"))(50), 
  #         display_numbers = as.matrix(ifelse(fdrDF < 0.1, "*", "")),   fontsize = 8)
  #dev.off()
}

