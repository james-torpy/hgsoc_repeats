
library(reshape2)
library(ggplot2)
library(pheatmap)


# define starting variables:
project <- "hgsoc_repeats"
methodName <- ""
refName <- "human-89.repeats.tab"
expName <- "exp1"
STypes <- c("FT", "erPT", "mrPT", "arPT", "prPT", "rfPT", "msST", "pAF", "rcAF")
CTypes <- c("arPT", "erPT", "FT", "mrPT", "msST", "pAF", "prPT", "rcAF", "rfPT")
annot <- "custom3"

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
#homeDir <- "/Users/jamestorpy/Documents/Garvan/phd"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/", expName, "/")
plotDir <- paste0(resultsDir, "/R/", expName, "/plots/DEplots/")

system(paste0("mkdir -p ", plotDir))


### 1. Load in inputs ###

# topRep is top 15 with FDR<0.1, strRep is stricter - top 15 with FDR<0.1 and -1 < logFC < 1:
allGene <- readRDS(file=paste0(RobjectDir, "/", annot, "_DEreps.rds"))
topGene <- readRDS(file=paste0(RobjectDir, "/", annot, "_topRep.rds"))
strGene <- readRDS(file=paste0(RobjectDir, "/", annot, "_strRep.rds"))


### 2. Fetch vector of all, top and str genes: ###

allName <- rownames(allGene[[1]])

for (i in 1:length(topGene)) {
  if (i==1) {
    topName <- list(rownames(topGene[[i]]))
  } else {
    topName[[i]] <- rownames(topGene[[i]])
  }
}
topName <- gsub("^.*\\.", "", unique(unlist(topName)))

for (i in 1:length(strGene)) {
  if (i==1) {
    strName <- list(rownames(strGene[[i]]))
  } else {
    strName[[i]] <- rownames(strGene[[i]])
  }
}
strName <- gsub("^.*\\.", "", unique(unlist(strName)))


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
#pheatmap(fcDF, fontsize = 7)

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
#pheatmap(fdrDF, fontsize = 7)

strFDR <- allFDR[allFDR$repeat_id %in% strName,]
strFDRdf <- dcast(strFDR, repeat_id ~ sample)
rownames(strFDRdf) <- strFDRdf$repeat_id
strFDRdf <- subset(strFDRdf, select=-repeat_id)

strFC <- allFC[allFC$repeat_id %in% strName,]
strFCdf <- dcast(strFC, repeat_id ~ sample)
rownames(strFCdf) <- strFCdf$repeat_id
strFCdf <- subset(strFCdf, select=-repeat_id)
colnames(strFCdf) <- c("acquired_resistance", "extreme_response", "multiple_response", "metastatic", "primary_ascites", "primary_resistant", "refractory_ascites", "refractory")
pheatmap(strFCdf, color = colorRampPalette(c("#08519C", "white", "firebrick3"))(50), display_numbers = as.matrix(ifelse(strFDRdf < 0.1, "*", "")), fontsize = 8, scale="column")

library(grid)

draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

log2FC <- pheatmap(strFCdf, color = colorRampPalette(c("#08519C", "white", "firebrick3"))(50),
                   display_numbers = as.matrix(ifelse(strFDRdf < 0.1, "*", "")), fontsize = 8,
                   scale="column")

pdf(file=paste0(plotDir, "/log2FC_HGSOCvsFT_heatmap.pdf"))
pheatmap(strFCdf, color = colorRampPalette(c("#08519C", "white", "firebrick3"))(50),
         display_numbers = as.matrix(ifelse(strFDRdf < 0.1, "*", "")), fontsize = 8,
         scale="column")
dev.off()

#create the breaks
logstrFDRdf <- log10(strFDRdf)
colnames(logstrFDRdf) <- c("acquired_resistance", "extreme_response", "multiple_response", "metastatic", "primary_ascites", "primary_resistant", "refractory_ascites", "refractory")
myBreaks <- seq(min(logstrFDRdf), 1.26, length.out=ceiling(50/2) + 1)
#myBreaks <- unique(c(brk, sort((brk*-1))))
pdf(file=paste0(plotDir, "/FDR_HGSOCvsFT_heatmap.pdf"))
pheatmap(logstrFDRdf, breaks=myBreaks, color = colorRampPalette(c("#08519C", "white", "firebrick3"))(50), display_numbers = as.matrix(ifelse(strFDRdf < 0.1, "*", "")), fontsize = 8)
dev.off()
 


