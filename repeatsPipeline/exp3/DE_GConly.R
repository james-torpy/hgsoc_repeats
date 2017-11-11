### normaliseCounts.R ###

# This script takes a list of dfs of different classes and types of repeat counts and
# normalises using RUVseq, EdgeR and DEseq, then compares methods

### 0. Define variables/paths ###

# load packages needed:
library(tibble)
library(dplyr)
library(RColorBrewer)
library(RUVSeq)
library(ggplot2)
library(ggrepel)
library(preprocessCore)

# define starting variables:
project <- "hgsoc_repeats"
methodName <- ""
refName <- "human-89.repeats.tab"
expName <- "exp3"
STypes <- c("FT", "prPT")
#annot <- "custom3"

# define directories:
#homeDir <- "/Users/jamestorpy/clusterHome/"
homeDir <- "/share/ScratchGeneral/jamtor/"
#homeDir <- "/Users/jamestorpy/Documents/Garvan/phd"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/", expName, "/")
plotDir <- paste0(resultsDir, "/R/", expName, "/plots/DEplots/")

system(paste0("mkdir -p ", plotDir))


### 1. Load in counts ###

#if (file.exists(paste0(RobjectDir, "/",  "_DEreps.rds"))) {
#  repGeneL <- readRDS(file=paste0(RobjectDir, "/",  "_DEreps.rds"))
#} else {
#  
#  if (file.exists(paste0(RobjectDir, "/",  "_counts.RData"))) {
#    Counts <- load(file=paste0(RobjectDir, "/",  "_counts.RData"))
#  } else {
#    Counts <- readRDS(file=paste0(RobjectDir,  "_RepeatCounts/all_",  "RepeatCountDFs.rds"))
#    if (annot=="custom3") {
#      Counts2 <- readRDS(file=paste0(RobjectDir, "c_RepeatCounts/all_cRepeatCountDFs.rds"))
#    }
#    
#    if (annot=="custom3") {
#      Counts <- c(Counts[1], Counts[3:length(Counts)])
#    }
#    
#    countsDF <- do.call("rbind", Counts)
#    
#    if (annot=="custom3") {
#      countsDF <- rbind(countsDF, Counts2$other[2:3,][,-grep("prPT8", colnames(Counts2$other))])
#    }
#    # simplify the row names:
#    #rownames(countsDF) <- gsub("^.*\\.", "", rownames(countsDF))
#    
#    # remove duplicate samples:
#    countsDF <- countsDF[,-which(duplicated(colnames(countsDF)))]
    
    
    ### 2. Load in GCcounts ###
    
    gcFiles <- list.files(RobjectDir, pattern="GCfeatureCounts", full.names=T)
    gcFiles <- grep("subset", gcFiles, value=T, invert=T)
    gcNames <- gsub("_.*$", "", basename(gcFiles))
    
    gcL <- list()
    j=1
    for (f in gcFiles) {
      temp <- readRDS(file = f)
      gcL[[j]] <- as.data.frame(temp$counts)
      GCrownames <- rownames(temp$counts)
      j=j+1
    }
    
    GCcountsDF <- do.call("cbind", gcL)
    colnames(GCcountsDF) <- gcNames
    GCcountsDF$gene_id <- GCrownames
    GCcountsDF$gene_id <- gsub("\\..*$", "", GCcountsDF$gene_id)
    GCcountsDF <- aggregate(.~gene_id, GCcountsDF, mean)
    rownames(GCcountsDF) <- GCcountsDF$gene_id
    GCcountsDF <- subset(GCcountsDF, select=-gene_id)
    
    
    # remove duplicate samples:
    #GCcountsDF <- GCcountsDF[-which(duplicated(colnames(GCcountsDF)))]
    
    # adjust so GCcountsDF only has columns of samples present in countsDF:
    #GCcountsDF <- GCcountsDF[,colnames(GCcountsDF) %in% colnames(countsDF)]
    
#    # append GCcountsDF to each GRanges object of custom1Counts:
#    Counts <- rbind(countsDF, GCcountsDF)
#    
#    # remove lrcT sample:
#    Counts <- as.data.frame(Counts[,-(grep("lrcT", colnames(Counts)))])

    Counts <- GCcountsDF

    saveRDS(Counts, file=paste0(RobjectDir, "/", "GCcounts.RData"))
  
  
  
  ### 3. Perform pre-normalisation PCA and RLE plots ###
  
  # eliminate lowly expressed genes (rows where there are less than 3 counts where df > 4):
  print(paste0("No. rows before filtering is: ", nrow(Counts)))
  Counts <- Counts %>%
    rownames_to_column('repeat_id') %>%
    dplyr::filter(rowSums(Counts > 5) >= 3) %>%
    column_to_rownames('repeat_id')
  print(paste0("No. rows after  filtering: ", nrow(Counts)))
  
  # create pre-normalised PCA plot from counts and plot:
  Counts <- apply(Counts, 2, unlist)
  if (ncol(Counts) > nrow(Counts)) {
    pca <- prcomp(Counts)
  } else {
    pca <- princomp(Counts)	  
  }
  
  if (file.exists(paste0(plotDir,  "_pcaCompsPrenormGC.pdf"))) {
    print(paste0(plotDir,  "_pcaCompsPrenormGC.pdf already exists, no need to create"))
  } else {
    print(paste0("Creating ", plotDir,  "_pcaCompsPrenormGC.pdf"))
    pdf(file = paste0(plotDir,  "_pcaCompsPrenormGC.pdf"))
    plot(pca)
    dev.off()
  }
  
  splt <- unlist(lapply(split(colnames(Counts), gsub("[0-9]", "", colnames(Counts))), length))
  for (i in 1:length(splt)) {
    if (i==1) {
      typeF <- c(rep(names(splt)[i], splt[i]))
    } else {
      typeF <- c(typeF, rep(names(splt)[i], splt[i]))
    }
  }
  levels(typeF) <- STypes
  
  # convert matrix into SeqExpressionSet:
  set <- newSeqExpressionSet(Counts, phenoData = data.frame(typeF, row.names=colnames(Counts)))
  
  # create pre-norm RLE plot:
  if (file.exists(paste0(plotDir,  "_RLEPrenormGC.pdf"))) {
    print(paste0(plotDir,  "_RLEPrenormGC.pdf already exists, no need to create"))
  } else {
    print(paste0("Creating ", plotDir,  "_RLEPrenormGC.pdf"))
    par(mar=c(1,1,1,1))
    pdf(file = paste0(plotDir,  "_RLEPrenormGC.pdf"))
    plotRLE(set)
    dev.off()
  }
  
  # create RUVseq pre-norm PCA:
  if (file.exists(paste0(plotDir,  "_pcaPrenormGC.pdf"))) {
    print(paste0(plotDir,  "_pcaPrenormGC.pdf already exists, no need to create"))
  } else {
    print(paste0("Creating ", plotDir,  "_pcaPrenormGC.pdf"))
    pdf(file = paste0(plotDir,  "_pcaPrenormGC.pdf"), height = 15, width = 20)
    plotPCA(set, cex=0.7)
    dev.off()
  }
  
  
  ### 4. perform normalisation on counts using RUVseq:
  
  # perform between lane full normalisation:
  nSet <- betweenLaneNormalization(set, which="full")
  pdf(file = paste0(plotDir,  "_RLElaneNormGC.pdf"))
  plotRLE(nSet, outline=FALSE, ylim=c(-4, 4))
  dev.off()
  
  #grps <- split(colnames(Counts), gsub("\\d.*$", "", colnames(Counts)))
  #for (i in 1:length(grps)) {
  #  if (i==1) {
  #    grpVec <- c(rep(i, length(grps[[i]])))
  #  } else {
  #    grpVec <- c(grpVec, rep(i, length(grps[[i]])))
  #  }
  #}
  
  pdf(file = paste0(plotDir,  "_pcalaneNormGC.pdf"), height = 15, width = 20)
  plotPCA(nSet, cex=0.7)
  dev.off()
  
  
  ### 5. Perform differential expression comparing normalised FT controls to cancer samples ###
  
  genes <- rownames(Counts)
  # normalise using upper-quartile normalisation (http://inaykmittal.blogspot.com.au/2013/10/	pkmrpkm-normalization-caveat-and-upper.html)
  # design matrix specifying all samples as the thing to be compared to:
  design <- model.matrix(~0+typeF, data=pData(nSet))
  # convert set into a DGElist format which specifies library size and normalisation factors:
  y <- DGEList(counts=counts(nSet), group=typeF)
  
  # create MDS plot:
  pdf(file = paste0(plotDir,  "_mdslaneNormGC.pdf"), height = 15, width = 20)
  plotMDS(y)
  dev.off()
  
  
  # calculate the deviance residuals from a first-pass GLM regression of the counts on the co-variates 	of interest (p8 RUVseq manual):
  # estimate dispersion:
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  # adjust values using dispersion:
  fit <- glmFit(y, design)
  
  saveRDS(fit, file=paste0(RobjectDir, "/",  "DEfitGC.rds"))
  
  con <- c(-1, 1)


  comp <- paste0("FT_vs_prPT")
  
  # perform likelihood ratio test:
  lrt <- glmLRT(fit, contrast = con)
  
  # determine the top DE genes:
  topTags(lrt)
  
  # save lrt as RDS:
  saveRDS(lrt, file=paste0(RobjectDir, "/", expName, "/FT_vs_prPT_lrtGC.rds"))

  
  ### 5. Calculate differential expression values ###
  
  # PCA and RLE plots of RUVseq RUVr-normalised data looked best with etween lane normalisation = 	'full', will go with this #
  
  # fetch summary of differentially expressed genes (those with FDR =< 0.05:
  DEs <- summary(result <- decideTestsDGE((lrt)))
  
  # fetch all gene DE info, 
  allGenes <- as.data.frame(topTags(lrt, n=Inf))
  allGenes$threshold <- as.factor(allGenes$FDR < 0.1)
  sigGenes <- filter(allGenes, FDR<0.1)
  print(sigGenes)
  
  # plot on volcano plot:
  # include the control genes for labelling:
  # add positive and negative control genes CD47, CCNE1, GAPDH, b-actin:
  posGenes <- rbind(allGenes["ENSG00000091831",], allGenes["ENSG00000196776",], allGenes["ENSG00000105173",])
  negGenes <- rbind(allGenes["ENSG00000111640",], allGenes["ENSG00000075624",], allGenes["ENSG00000204574",],
                    allGenes["ENSG00000104904",], allGenes["ENSG00000089157",])
  rownames(posGenes) <- c("ESR1", "CD47", "CCNE1")
  rownames(negGenes) <- c("GAPDH", "beta-actin", "ABCF1", "OAZ1", "RPLP0")
 
   posGenes$threshold <- "POSITIVE"
  if (nrow(posGenes[posGenes$FDR<0.1,])>0) {
    posGenes[posGenes$FDR<0.1,]$threshold <- "POSSIG"
  }
  
  negGenes$threshold = "NEGATIVE"
  if (nrow(negGenes[negGenes$FDR<0.1,])>0) {
    negGenes[negGenes$FDR<0.1,]$threshold <- "NEGSIG"
  }
 
  lab <- rbind(posGenes, negGenes)
  allGenes <- rbind(rbind(allGenes, posGenes), negGenes)
  lab$genes <- rownames(lab)
  
  p <- ggplot(data=allGenes, aes(x=logFC, y=-log10(FDR), color=threshold))
  p <- p + geom_point(data=allGenes)
  p <- p + geom_text_repel(data=lab, aes(label=genes))
  p <- p + theme(legend.position = "none")
  p <- p + labs(x="log2 fold change vs FT control", y="-log10 FDR")
  p <- p +  xlim(c(-18, 18))
  if (file.exists(paste0(plotDir, "/",  "_volcano_FDR0.1_", comp, "GC.pdf"))) {
    print(paste0(plotDir, "/",  "_volcano_FDR0.1_", comp, ".GCpdf"))
    p
  } else {
    print(paste0("Creating ",plotDir, "/",  "_volcano_FDR0.1_", comp, "GC.pdf"))
    pdf(file = paste0(plotDir, "/",  "_volcano_FDR0.1_", comp, "GC.pdf"))
    print(p)
    dev.off()
  }
      
  # create plot with only DE expresed genes with FDR<0.05
  # and logFC > 10:
  allGenes <- as.data.frame(topTags(lrt, n=Inf))
  allGenes$threshold <- as.factor(allGenes$FDR < 0.1)
  allGenes$genes <- rownames(allGenes)
  sigGenes <- filter(allGenes, FDR<0.1)
  sigGenes <- sigGenes[sigGenes$logFC > 10,]
  print(sigGenes)
  allGenes <- subset(allGenes, select=-genes)
      
  p <- ggplot(data=allGenes, aes(x=logFC, y=-log10(FDR), color=threshold))
  p <- p + geom_point(data=allGenes)
  p <- p + geom_text_repel(data=sigGenes, aes(label=genes))
  p <- p + theme(legend.position = "none")
  p <- p + labs(x="log2 fold change vs FT control", y="-log10 FDR")
  p <- p +  xlim(c(-18, 18))
  if (file.exists(paste0(plotDir, "/",  "volcano_FDR0.05_logFC10_", comp, "GC.pdf"))) {
    print(paste0(plotDir, "/",  "volcano_FDR0.05_logFC10_", comp, ".GCpdf"))
    p
  } else {
    print(paste0("Creating ",plotDir, "/",  "volcano_FDR0.05_logFC10_", comp, "GC.pdf"))
    pdf(file = paste0(plotDir, "/",  "volcano_FDR0.05_logFC10_", comp, "GC.pdf"))
    print(p)
    dev.off()
  }

saveRDS(allGenes, file=paste0(RobjectDir, "/",  "_allGenesGC.rds"))

if (!file.exists(paste0(RobjectDir, "DEImg_", expName, "GC.RData"))) {
  save.image(file = paste0(RobjectDir, "DEImg_", expName, "GC.RData"))
}



