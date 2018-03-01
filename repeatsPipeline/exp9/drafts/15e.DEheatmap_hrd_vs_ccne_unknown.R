
library(reshape2)
library(ggplot2)
library(pheatmap)


# define starting variables:
project <- "hgsoc_repeats"
expName <- "exp9"
sTypes <- c("FT", "unknownDrivers", "HRD", "CCNEAmp", "bothDrivers")
annot <- "custom3"

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
#homeDir <- "/Users/jamestorpy/Documents/Garvan/phd"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/", expName, "/hrd_vs_ccne_unknown/")
plotDir <- paste0(resultsDir, "/R/", expName, "/plots/DEplots/hrd_vs_ccne_unknown/")

system(paste0("mkdir -p ", plotDir))


### 1. Load in inputs ###

allGene <- readRDS(file=paste0(RobjectDir, "/", annot, "_DEsigReps.rds"))

### 2. Fetch vector of all, top and str genes: ###

allName <- rownames(allGene[[1]])


### 3. Create logFC and FDR heatmaps:
  
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


par(mar=c(4,4,4,4))

pheatmap(fcDF, color = colorRampPalette(c("#08519C", "white", "firebrick3"))(50),
         display_numbers = as.matrix(ifelse(fdrDF < 0.1, "*", "")), fontsize = 8, scale="column",
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

# FC pheatmap with column clustering:
log2FC <- pheatmap(fcDF, color = colorRampPalette(c("#08519C", "white", "firebrick3"))(50),
                   display_numbers = as.matrix(ifelse(fdrDF < 0.1, "*", "")), fontsize = 8,
                   scale="column")

pdf(file=paste0(plotDir, "/log2FC_HGSOCvsFT_heatmap.pdf"), width=10, height=10)
pheatmap(fcDF, color = colorRampPalette(c("#08519C", "white", "firebrick3"))(50),
         display_numbers = as.matrix(ifelse(fdrDF < 0.1, "*", "")), fontsize = 8,
         scale="column")
dev.off()


# FC pheatmap with column clustering:
logfdrDF <- log10(fdrDF)
colnames(logfdrDF) <- Names
#myBreaks <- seq(min(logfdrDF), 1.26, length.out=ceiling(50/2) + 1)
#myBreaks <- unique(c(brk, sort((brk*-1))))
pdf(file=paste0(plotDir, "/log10_FDR_HGSOCvsFT_heatmap.pdf"), width=10, height=10)
pheatmap(logfdrDF, color = colorRampPalette(c("#08519C", "white"))(50), display_numbers = as.matrix(ifelse(fdrDF < 0.1, "*", "")), fontsize = 8)
dev.off()
 


