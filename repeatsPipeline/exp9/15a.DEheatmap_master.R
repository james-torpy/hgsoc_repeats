library(reshape2)
library(ggplot2)
library(pheatmap)
library(tibble)
library(dplyr)


# define starting variables:
project <- "hgsoc_repeats"
expName <- "exp9"
Type <- "custom3"
descrip <- "htseq_hgsoc_split_more_by_drug_response_vs_FT"

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
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
                     expName, "/", descrip, "/")
plotDir <- paste0(resultsDir, "/R/", expName,
                  "/plots/DEplots/", descrip, "/")


### 1. Load in inputs ###

for (i in 1:length(resultTypes)) {
  if (resultTypes[i] == "repeats") {
    if (i==1) {
      allGene <- list(readRDS(file=paste0(RobjectDir, "/", Type, "_DEsigReps.rds")))
    } else {
      allGene[[i]] <- readRDS(file=paste0(RobjectDir, "/", Type, "_DEsigReps.rds"))
    }
    names(allGene)[i] <- "sigRep"
  } else if (resultTypes[i] == "epiMods") {
    if (i==1) {
      allGene <- list(readRDS(file=paste0(RobjectDir, "/", Type, "_DEsigReps.rds")))
    } else {
      allGene[[i]] <- readRDS(file=paste0(RobjectDir, "/", Type, "_DEepiGenes.rds"))
    }
    names(allGene)[i] <- "allEpi"
  }
}

# fetch sample number to add to colnames:
sample_nos <- readRDS(paste0(RobjectDir, "/sample_no_per_cat.rds"))

allGene <- lapply(allGene, function(x) {
  Names <- names(x)
  for (i in 1:length(Names)) {
    names <- c(strsplit(Names[i], "_vs_")[[1]][1], strsplit(Names[i], "_vs_")[[1]][2])
    nos <- c(sample_nos[grep(names[1], names(sample_nos))], sample_nos[grep(names[2], names(sample_nos))])
    if (i==1) {
      newNames <- c(paste0(names[1], "_vs_", names[2], " (n=", nos[1], ",", nos[2], ")"))
    } else {
      newNames[i] <- paste0(names[1], "_vs_", names[2], " (n=", nos[1], ",", nos[2], ")")
    }
  }
  names(x) <- newNames
  return(x)
})


### 2. Split into logFC and FDR dfs:

for (j in 1:length(allGene)) {
  allFC <- do.call("rbind", allGene[[j]])
  allFC$sample <- gsub("\\..*$", "", rownames(allFC))
  allFC$repeat_id <- gsub("^.*\\.*\\.", "", rownames(allFC))
  allFC <- subset(allFC, select=c(logFC, sample, repeat_id))
  allFC[,4] <- allFC$logFC
  allFC <- subset(allFC, select=-logFC)
  names(allFC)[3] <- "logFC"
  
  fcDF <- dcast(allFC, repeat_id ~ sample)
  rownames(fcDF) <- fcDF$repeat_id
  fcDF <- subset(fcDF, select=-repeat_id)
  
  allFDR <- do.call("rbind", allGene[[j]])
  allFDR$sample <- gsub("\\..*$", "", rownames(allFDR))
  allFDR$repeat_id <- gsub("^.*\\.*\\.", "", rownames(allFDR))
  allFDR <- subset(allFDR, select=c(FDR, sample, repeat_id))
  allFDR[,4] <- allFDR$FDR
  allFDR <- subset(allFDR, select=-FDR)
  names(allFDR)[3] <- "FDR"
  
  fdrDF <- dcast(allFDR, repeat_id ~ sample)
  rownames(fdrDF) <- fdrDF$repeat_id
  fdrDF <- subset(fdrDF, select=-repeat_id)
  
  
  ### 3. Create logFC and FDR heatmaps:
  
  ##### fix below this point #####
  
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
  
  # create FC pheatmap without column clustering:
  #pdf(file=paste0(plotDir, "/log2FC_HGSOC_heatmap_no_colclustering_", names(allGene)[j], ".pdf"), width=10, height=10)
  #pheatmap(fcDF, color = colorRampPalette(c("#08519C", "white", "firebrick3"))(50),
  #         display_numbers = as.matrix(ifelse(fdrDF < 0.1, "*", "")), fontsize = 8,
  #         cluster_cols=F)
  #dev.off()
  
  # create FC pheatmap with column clustering:
  #pdf(file=paste0(plotDir, "/log2FC_HGSOC_heatmap_colclustering_", names(allGene)[j], ".pdf"), width=10, height=10)
  #pheatmap(fcDF, color = colorRampPalette(c("#08519C", "white", "firebrick3"))(50),
  #         display_numbers = as.matrix(ifelse(fdrDF < 0.1, "*", "")), fontsize = 8,
  #         scale="column")
  #dev.off()
  
  # create FDR pheatmap without column clustering:
  logfdrDF <- log10(fdrDF)
  
  #pdf(file=paste0(plotDir, "/log10_FDR_HGSOC_heatmap_no_colclustering_", names(allGene)[j], ".pdf"), width=10, height=10)
  #pheatmap(logfdrDF, color = colorRampPalette(c("#08519C", "white"))(50), 
  #         display_numbers = as.matrix(ifelse(fdrDF < 0.1, "*", "")), fontsize = 8, 
  #         scale="column")
  #dev.off()
  
  # create FDR pheatmap with column clustering:
  #pdf(file=paste0(plotDir, "/log10_FDR_HGSOC_heatmap_colclustering_", names(allGene)[j], ".pdf"), width=10, height=10)
  #pheatmap(logfdrDF, color = colorRampPalette(c("#08519C", "white"))(50), 
  #         display_numbers = as.matrix(ifelse(fdrDF < 0.1, "*", "")), fontsize = 8)
  #dev.off()
  
  # isolate most significant genes and create FC heatmap:
  interesting_fc <- fcDF %>%
    rownames_to_column('gene_id') %>%
    dplyr::filter(rowSums(fdrDF < FDRthresh & fcDF > FCthresh|fcDF < -FCthresh) >= 1) %>%
    column_to_rownames('gene_id')
  
  # create fdr for interesting repeats:
  interesting_fdr <-  fdrDF[rownames(interesting_fc),]

  # save interesting repeats:
  if (j==1) {
    saveRDS(rownames(interesting_fc), file=paste0(RobjectDir, "/interesting_fc.rds"))
  }
  
  if (nrow(interesting_fc) > 1) {
    # FC without clustering:
    #pdf(file=paste0(plotDir, "/log2FC_HGSOC_heatmap_FDR", FDRthresh, "_colclustering.pdf"), width=10, height=10)
    #pheatmap(interesting_fc, color = colorRampPalette(c("#08519C", "white", "firebrick3"))(50),
    #         display_numbers = as.matrix(ifelse(interesting_fdr < 0.1, "*", "")), fontsize = 8,
    #         cluster_cols = F)
    #dev.off()
    
    # FC with clustering:
    #pdf(file=paste0(plotDir, "/log2FC_HGSOC_heatmap_FDR", FDRthresh, "_no_colclustering.pdf"), width=10, height=10)
    #pheatmap(interesting_fc, color = colorRampPalette(c("#08519C", "white", "firebrick3"))(50),
    #         display_numbers = as.matrix(ifelse(interesting_fdr < 0.1, "*", "")), fontsize = 8,
    #         cluster_cols = T)
    #dev.off()
    
    # FC with clustering, manual breaks:
    paletteLength <- 50
    myBreaks <- c(seq(min(fcDF), 0, length.out=ceiling(paletteLength/2) + 1), 
                  seq(max(fcDF)/paletteLength, max(fcDF), length.out=floor(paletteLength/2)))
    
    pdf(file=paste0(plotDir, "/log2FC_HGSOC_heatmap_FDR", FDRthresh, "_FC", FCthresh, "_colclustering_", 
                    names(allGene)[j], ".pdf"), width=10, height=20)
    pheatmap(interesting_fc, color = colorRampPalette(c("#08519C", "white", "firebrick3"))(50),
             display_numbers = as.matrix(ifelse(interesting_fdr < 0.1, "*", "")), fontsize = 8,
             cluster_cols = T,  breaks=myBreaks)
    dev.off()
    
    # FDR without clustering:
    #log_interesting_fdr <- log10(interesting_fdr)
    
    #pdf(file=paste0(plotDir, "/log10_FDR_HGSOC_heatmapFDR", FDRthresh, "_colclustering_", names(allGene)[j], ".pdf"), width=10, height=10)
    #pheatmap(logfdrDF, color = colorRampPalette(c("#08519C", "white"))(50), 
    #         display_numbers = as.matrix(ifelse(fdrDF < 0.1, "*", "")), fontsize = 8,
    #         cluster_cols = F)
    #dev.off()
    
    # FDR with clustering:
    pdf(file=paste0(plotDir, "/log10_FDR_HGSOC_heatmapFDR", FDRthresh, "_FC", FCthresh, "colclustering_", names(allGene)[j], ".pdf"), width=10, height=10)
    pheatmap(logfdrDF, color = colorRampPalette(c("#08519C", "white"))(50), 
             display_numbers = as.matrix(ifelse(fdrDF < 0.1, "*", "")), fontsize = 8)
    dev.off()
  }
}
