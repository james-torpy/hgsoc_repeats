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
methodName <- "method1"
refName <- "human-89.repeats.tab"
sampleName <- ""

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
#homeDir <- "/Users/jamestorpy/Documents/Garvan/phd"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/")
outDir <- paste0(resultsDir, "/R/", methodName, "/")


### 1. Load in counts ###

# fetch inFiles:
inFiles <- grep("subset", list.files(RobjectDir, pattern = "RepeatsCounts", full.names = T), value=T, invert=T)

# define the sampleNames:
sampleNames <- strsplit(inFiles, "/")
sampleNames <- lapply(sampleNames, function(x) {
  return(x[12])
})
sampleNames <- unlist(sampleNames)

# Load in and cbind the sc counts for each sample:
scFiles <- grep("sc_", inFiles, value=T)
scCountsL <- list()
i=1
for (f in scFiles) {
  scCountsL[[i]] <- readRDS(f)
  i=i+1
}
scCounts <- do.call("cbind", scCountsL)
colnames(scCounts) <- sampleNames

# Load in and cbind the c counts for each sample:
cFiles <- grep("c_", inFiles, value=T)
cFiles <- grep("sc_", cFiles, value=T, invert=T)
cCountsL <- list()
i=1
for (f in cFiles) {
  cCountsL[[i]] <- readRDS(f)
  i=i+1
}
cCounts <- do.call("cbind", cCountsL)
colnames(cCounts) <- sampleNames

# Load in and cbind the t counts for each sample:
tFiles <- grep("t_", inFiles, value=T)
tCountsL <- list()
i=1
for (f in tFiles) {
  tCountsL[[i]] <- readRDS(f)
  i=i+1
}
tCounts <- do.call("cbind", tCountsL)
colnames(tCounts) <- sampleNames


### 2. Perform pre-normalisation PCA and RLE plots ###

# load in GCcountsDF:
GCcountsDF <- readRDS(file = paste0(RobjectDir, "GCcountsDF.rds"))
rownames(GCcountsDF) <- GCcountsDF$gene_id
GCcountsDF <- subset(GCcountsDF, select=-gene_id)

# append GCcounts to each GRanges object of sc, cCounts:
scGCcounts <- rbind(scCounts, GCcountsDF)
cGCcounts <- rbind(cCounts, GCcountsDF)
tGCcounts <- rbind(tCounts, GCcountsDF)

types <- c("sc", "c", "t")
i=1
for (x in list(scGCcounts, cGCcounts, tGCcounts)) {
	Type <- types[i]
	# eliminate lowly expressed genes (rows where there are less than 3 counts where df > 4):
	print(paste0("No. rows before filtering is: ", nrow(x)))
  Counts <- x %>%
  		rownames_to_column('repeat_id') %>%
  		dplyr::filter(rowSums(x > 4) >= 3) %>%
  		column_to_rownames('repeat_id')
  print(paste0("No. rows after  filtering: ", nrow(Counts)))

	# define plot colour scheme:
	cols <- brewer.pal(6, "Paired")

	# create pre-normalised PCA plot from counts and plot:
	Counts <- apply(Counts, 2, unlist)
	if (ncol(Counts) > nrow(Counts)) {
	  pca <- prcomp(Counts)
	} else {
	  pca <- princomp(Counts)	  
	}
  
	if (file.exists(paste0(outDir, Type, "_pcaCompsSubsetPrenormGC.pdf"))) {
	  print(paste0(outDir, Type, "_pcaCompsSubsetPrenormGC.pdf already exists, no need to create"))
	} else {
	  print(paste0("Creating ", outDir, Type, "_pcaCompsSubsetPrenormGC.pdf"))
	  pdf(file = paste0(outDir, Type, "_pcaCompsSubsetPrenormGC.pdf"))
	  plot(pca)
	  dev.off()
	}

	x <- as.factor(rep(STypes, each=3))
	
	######
	# try preprocessCore normalize.quantiles
	#nCount <- normalizeQuantiles(cCounts)
	######
		
	# convert matrix into SeqExpressionSet:
	set <- newSeqExpressionSet(Counts, phenoData = data.frame(x, row.names = colnames(Counts)))
	
	# create pre-norm RLE plot:
	if (file.exists(paste0(outDir, Type, "_RLESubsetPrenormGC.pdf"))) {
	  print(paste0(outDir, Type, "_RLESubsetPrenormGC.pdf already exists, no need to create"))
	} else {
	  print(paste0("Creating ", outDir, Type, "_RLESubsetPrenormGC.pdf"))
	  par(mar=c(1,1,1,1))
	  pdf(file = paste0(outDir, Type, "_RLESubsetPrenormGC.pdf"))
	  plotRLE(set)
	  dev.off()
	}
	
	# create RUVseq pre-norm PCA:
	if (file.exists(paste0(outDir, Type, "_pcaSubsetPrenormGC.pdf"))) {
	  print(paste0(outDir, Type, "_pcaSubsetPrenormGC.pdf already exists, no need to create"))
	} else {
	  print(paste0("Creating ", outDir, Type, "_pcaSubsetPrenormGC.pdf"))
	  pdf(file = paste0(outDir, Type, "_pcaSubsetPrenormGC.pdf"))
	  plotPCA(set, cex=0.7)
	  dev.off()
	}
	
	
	### 3. perform normalisation on counts using RUVseq:
	
	######
	
	#set1 <- betweenLaneNormalization(set, which="upper")
	#plotRLE(set1, outline=FALSE, ylim=c(-4, 4))
	#plotPCA(set1, cex=1.2)
	
	#set2 <- betweenLaneNormalization(set, which="full")
	#plotRLE(set2, outline=FALSE, ylim=c(-4, 4))
	#plotPCA(set2, cex=1.2)
	
	#set3 <- betweenLaneNormalization(set, which="median")
	#plotRLE(set3, outline=FALSE, ylim=c(-4, 4))
	#plotPCA(set3, cex=1.2)
	
	######
	
	# perform between lane full normalisation:
	nSet <- betweenLaneNormalization(set, which="full")
	pdf(file = paste0(outDir, Type, "_RLESubsetlaneNormGC.pdf"))
	plotRLE(nSet, outline=FALSE, ylim=c(-4, 4))
	dev.off()
	
	pdf(file = paste0(outDir, Type, "_pcaSubsetlaneNormGC.pdf"))
	plotPCA(nSet, cex=0.7)
	dev.off()
	
	
	### 4. Perform differential expression comparing normalised FT controls to bowtell PRs ###
	
	genes <- rownames(Counts)
	# normalise using upper-quartile normalisation (http://vinaykmittal.blogspot.com.au/2013/10/fpkmrpkm-normalization-caveat-and-upper.html)
	# design matrix specifying all samples as the thing to be compared to:
	design <- model.matrix(~x, data=pData(nSet))
	# convert set into a DGElist format which specifies library size and normalisation factors:
	y <- DGEList(counts=counts(nSet), group=x)
	# calculate the deviance residuals from a first-pass GLM regression of the counts on the co-variates of interest (p8 RUVseq manual):
	# estimate dispersion:
	y <- estimateGLMCommonDisp(y, design)
	y <- estimateGLMTagwiseDisp(y, design)
	# adjust values using dispersion:
	fit <- glmFit(y, design)
	# perform likelihood ratio test:
	lrt <- glmLRT(fit, coef=2)
	
	# determine the top DE genes:
	topTags(lrt)
	
	create_rds <- function(x, rds) {
	  if (file.exists(rds)) {
	    print(paste0(rds, " exists, no need to create"))
	  } else {
	    print(paste0("Creating ", rds))
	    saveRDS(x, file = rds)
	  }
	}
	
	create_rds(lrt, paste0(RobjectDir, "bowtellFTvsPRSubsetGC_RUVlrt.rds"))
	#lrt <- readRDS(file = paste0(RobjectDir, "bowtellFTvsPRSubset_RUVlrt.rds"))
	
	### 5. Calculate differential expression values ###
	
	# PCA and RLE plots of RUVseq RUVr-normalised data looked best with between lane normalisation = 'full', will go with this #
	
	# fetch summary of differentially expressed genes (those with FDR =< 0.05:
	DEs <- summary(result <- decideTestsDGE((lrt)))
	
	# fetch all gene DE info, 
	allGenes <- as.data.frame(topTags(lrt, n=Inf))
	repGenes <- allGenes[grep("ENS", rownames(allGenes), invert = T),]
	sigGenes <- filter(repGenes, FDR<0.05)
	
	# add negative log p-value column to allGenes:
	#allGenes$negLog10PValue <- -log10(allGenes$PValue)
	
	# plot on volcano plot:
	if (Type == "sc" | Type == "c") {
	  repGenes$threshold <- as.factor(repGenes$FDR < 0.3)
	  sig <- subset(repGenes, threshold == T)
	  sig$genes <- rownames(sig)
	  p <- ggplot(data=repGenes, aes(x=logFC, y=-log10(FDR), color=threshold))
	  p <- p + geom_point(data=repGenes)
	  p <- p + geom_text_repel(data=sig, aes(label=genes))
	  p <- p + theme(legend.position = "none")
	  p <- p + labs(x="log2 fold change vs FT control", y="-log10 FDR")
	  p <- p +  xlim(c(-5, 5))
	  pdf(file = paste0(outDir, "plots/", Type, "_volcanoFDR0.3.pdf"))
	  print(p)
	  dev.off()
	} else if (Type == "t") {
	  repGenes$threshold <- as.factor(repGenes$FDR < 0.1)
	  sig <- subset(repGenes, threshold == T)
	  sig$genes <- rownames(sig)
	  p <- ggplot(data=repGenes, aes(x=logFC, y=-log10(FDR), color=threshold))
	  p <- p + geom_point(data=repGenes)
	  p <- p + geom_text_repel(data=sig, aes(label=genes))
	  p <- p + theme(legend.position = "none")
	  p <- p + labs(x="log2 fold change vs FT control", y="-log10 FDR")
	  p <- p +  xlim(c(-5, 5))
	  pdf(file = paste0(outDir, "plots/", Type, "_volcanoFDR0.1.pdf"))
	  print(p)
	  dev.off()
	}
	
	i=i+1
}

if (!file.exists(paste0(RobjectDir, "DEsubsetsImg.RData"))) {
  save.image(file = paste0(RobjectDir, "DEsubsetsImg.RData"))
}



