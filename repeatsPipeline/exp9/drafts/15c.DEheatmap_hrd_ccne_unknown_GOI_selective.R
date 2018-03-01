
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
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/", expName, "/hrd_ccne_unknown/")
plotDir <- paste0(resultsDir, "/R/", expName, "/plots/DEplots/hrd_ccne_unknown/")

system(paste0("mkdir -p ", plotDir))


### 1. Load in inputs ###

allGene <- readRDS(file=paste0(RobjectDir, "/", annot, "_DEsigReps.rds"))
allGene <- allGene[2:4]
names(allGene) <- c("CCNEAmp_vs_FT", "HRD_vs_FT", "unknownDrivers_vs_FT")

# 2. Isolate custom genes of interest:
cenGenes <- lapply(allGene, function(x) {
  return(x[c("ACRO1", "(CATTC)n", "GSATX", "HSATII", "CER", 
             "BSR_Beta", "GSAT"),])
})

rtGenes <- lapply(allGene, function(x) {
  return(x[c("AluYi6", "CER", "L1M3a", "L1M3c", "L1MC", "L1P4", "L1PA5", 
             "FAM", "L1MD3", "L1PA2", "REP522"),])
})

groupList <- list(cenGenes, rtGenes)
names(groupList) <- c("cenGenes", "rtGenes")


### 3. Create logFC and FDR heatmaps:

# fetch sample number to add to colnames:
sample_nos <- readRDS(paste0(RobjectDir, "/sample_no_per_cat.rds"))
temp_names <- paste0(rep("FT_vs_", length(sample_nos)), names(sample_nos))
sample_nos <- paste0(rep( sample_nos[grep("FT", names(sample_nos))], length(sample_nos) ), ",", sample_nos)
names(sample_nos) <- temp_names
sample_nos <- c(sample_nos[1:grep("_FT", names(sample_nos))-1], sample_nos[(grep("_FT", names(sample_nos))+1):length(sample_nos)])
sample_nos <- sample_nos[2:4]
names(sample_nos) <- c("CCNEAmp_vs_FT", "HRD_vs_FT", "unknownDrivers_vs_FT")

for (i in 1:length(groupList)) {
  allFC <- do.call("rbind", groupList[[i]])
  allFC$sample <- gsub("\\..*$", "", rownames(allFC))
  allFC$repeat_id <- gsub("^.*\\.*\\.", "", rownames(allFC))
  allFC <- subset(allFC, select=c(logFC, sample, repeat_id))
  allFC[,4] <- allFC$logFC
  allFC <- subset(allFC, select=-logFC)
  names(allFC)[3] <- "logFC"
  
  fcDF <- dcast(allFC, repeat_id ~ sample)
  rownames(fcDF) <- fcDF$repeat_id
  fcDF <- subset(fcDF, select=-repeat_id)
  
  colnames(fcDF) <- paste0(colnames(fcDF), " (n=", sample_nos, ")")
  
  # change order of columns:
  fcDF <- fcDF[,c(3, 1, 2)]
  
  allFDR <- do.call("rbind", groupList[[i]])
  allFDR$sample <- gsub("\\..*$", "", rownames(allFDR))
  allFDR$repeat_id <- gsub("^.*\\.*\\.", "", rownames(allFDR))
  allFDR <- subset(allFDR, select=c(FDR, sample, repeat_id))
  allFDR[,4] <- allFDR$FDR
  allFDR <- subset(allFDR, select=-FDR)
  names(allFDR)[3] <- "FDR"
  
  fdrDF <- dcast(allFDR, repeat_id ~ sample)
  rownames(fdrDF) <- fdrDF$repeat_id
  fdrDF <- subset(fdrDF, select=-repeat_id)

  # add sample numbers to column names:
  colnames(fdrDF) <- paste0(colnames(fdrDF), " (n=", sample_nos, ")")
  
  # change order of columns:
  fdrDF <- fdrDF[,c(3, 1, 2)]
  
  if (i==1) {
    fc <- list(fcDF)
    fdr <- list(fdrDF)
  } else {
    fc[[i]] <- fcDF
    fdr[[i]] <- fdrDF
  }
  
  par(mar=c(4,4,4,4))
  
  library(grid)
  
  draw_colnames_45 <- function (coln, gaps, ...) {
    coord = pheatmap:::find_coordinates(length(coln), gaps)
    x = coord$coord - 0.5 * coord$size
    res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust =   0.5, hjust = 1, rot = 45, gp = gpar(...))
    return(res)
  }
  
  # 'Overwrite' default draw_colnames with my own version:
  assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                    ns=asNamespace("pheatmap"))
  
  # set up colour and breaks for heatmaps:
  paletteLength <- 50
  
  #if (i==1) {
  #  myColor <- colorRampPalette(c("firebrick3", "white", "#08519C"))(paletteLength)
  #} else {
    myColor <- colorRampPalette(c("#08519C", "white", "firebrick3"))(paletteLength)
  #}
  
  myBreaks <- c(seq(min(fcDF), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(fcDF)/paletteLength, max(fcDF), length.out=floor(paletteLength/2)))
  
  # save heatmap of selected repeats:
  pdf(file=paste0(plotDir, "/log2FC_HGSOCvsFT_heatmap_no_both_", names(groupList)[i], ".pdf"), width=10,   height=10)
  pheatmap(fcDF, display_numbers = as.matrix(ifelse(fdrDF < 0.1, "*", "")), 
           color = myColor, fontsize = 8, scale="none", 
           breaks = myBreaks, cluster_cols=F)
  dev.off()
  
  #pheatmap(fcDF, display_numbers = as.matrix(ifelse(fdrDF < 0.1, "*", "")), 
  #         color = colorRampPalette(c(hmCol))(50), fontsize = 8, scale="none", 
  #         breaks = myBreaks, cluster_cols=F)
  
  # FC pheatmap without column clustering:
  logfdrDF <- log10(fdrDF)
  
  pdf(file=paste0(plotDir, "/log10_FDR_HGSOCvsFT_heatmap_no_both_", names(groupList)[i], ".pdf"), width=10,    height=10)
  pheatmap(logfdrDF, display_numbers = as.matrix(ifelse(fdrDF < 0.1, "*", "")), fontsize = 8,
           color = colorRampPalette(c("white", "#08519C"))(50),
           scale="none", cluster_cols=F)
  dev.off()
}






