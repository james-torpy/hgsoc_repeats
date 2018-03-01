
library(reshape2)
library(ggplot2)
library(pheatmap)


# define starting variables:
project <- "hgsoc_repeats"
expName <- "exp9"
STypes <- c("FT", "HGSOC")
annot <- "custom3"

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
#homeDir <- "/Users/jamestorpy/Documents/Garvan/phd"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/", expName, "/allHGSOC/")
plotDir <- paste0(resultsDir, "/R/", expName, "/plots/DEplots/test/allHGSOC/")

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
pdf(file=paste0(plotDir, "/log2FC_heatmap.pdf"), width=10, height=10)
pheatmap(fcDF, fontsize = 7, cluster_rows = F, cluster_cols = F)
dev.off()

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
pdf(file=paste0(plotDir, "/log2FDR_heatmap.pdf"), width=10, height=10)
pheatmap(fdrDF, fontsize = 7)
dev.off()


#colnames(fcDF) <- c("acquired_resistance", "extreme_response", "multiple_response", "metastatic", "primary_ascites", "primary_resistant", "refractory_ascites", "refractory")
Names <- gsub(
  "FT_vs_", "", gsub(
    "prPT", "primary_resistant", gsub(
      "rfPT", "primary_refractory", gsub(
        "typeF", "", colnames(fcDF)
       )
     )
  )
)

colnames(fcDF) <- Names
fcDF <- cbind(fcDF[,grep("resistant", colnames(fcDF))], fcDF[,grep("refractory", colnames(fcDF))])

colnames(fdrDF) <- Names
fdrDF <- cbind(fdrDF[,grep("resistant", colnames(fdrDF))], fdrDF[,grep("refractory", colnames(fdrDF))])


######

fake_cat <- c("AOCS_055_primary_resistant", "AOCS_079_primary_resistant", "AOCS_081_primary_resistant")
ind <- colnames(fcDF) %in% fake_cat

# arrange dfs to group annotation together:
fcDF <- cbind(fcDF[,which(ind)], fcDF[,which(!ind)])
fdrDF <- cbind(fdrDF[,which(ind)], fdrDF[,which(!ind)])

# create annotation around fake cat:
Annotation <- data.frame(fakeCat = factor(ind))
rownames(Annotation) <- colnames(fcDF)
Annotation$fakeCat <- as.factor(
  gsub(
    "TRUE", "mutated", gsub(
      "FALSE", "wt", Annotation$fakeCat
    )
  )
)

# change the colors of annotation:
fakeCat <- c("navy", "darkgreen")
names(fakeCat) <- c("wt", "mutated")
anno_cols <- list(fakeCat = fakeCat)

# create heatmap with annotation bar:
pheatmap(fcDF, color = colorRampPalette(c("#08519C", "white", "firebrick3"))(50),
         display_numbers = as.matrix(ifelse(fdrDF < 0.1, "*", "")), fontsize = 8, scale="column",
         cluster_cols = F, annotation = Annotation, annotation_colors = anno_cols)

######

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

# FC pheatmap without column clustering:


log2FC <- pheatmap(fcDF, color = colorRampPalette(c("#08519C", "white", "firebrick3"))(50),
                   display_numbers = as.matrix(ifelse(fdrDF < 0.1, "*", "")), fontsize = 8,
                   scale="column", cluster_cols = F)

pdf(file=paste0(plotDir, "/log2FC_HGSOCvsFT_heatmap_nocolclust.pdf"), width=10, height=10)
pheatmap(fcDF, color = colorRampPalette(c("#08519C", "white", "firebrick3"))(50),
         display_numbers = as.matrix(ifelse(fdrDF < 0.1, "*", "")), fontsize = 8,
         scale="column", cluster_cols = F)
dev.off()

#create the breaks
logfdrDF <- log10(fdrDF)
colnames(logfdrDF) <- Names
#myBreaks <- seq(min(logfdrDF), 1.26, length.out=ceiling(50/2) + 1)
#myBreaks <- unique(c(brk, sort((brk*-1))))
pdf(file=paste0(plotDir, "/log10_FDR_HGSOCvsFT_heatmap.pdf"), width=10, height=10)
pheatmap(logfdrDF, color = colorRampPalette(c("#08519C", "white"))(50), display_numbers = as.matrix(ifelse(fdrDF < 0.1, "*", "")), fontsize = 8)
dev.off()
 


