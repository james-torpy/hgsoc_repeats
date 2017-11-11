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
STypes <- c("bowtell_FT", "bowtell_PR", "grant_endo", "grant_endosis", "gtx_ft", "kur_v")
annot <- "custom3"
comp <- "bowtell_FT_vs_PR_vs_gtx_ft_vs_PR"

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
#homeDir <- "/Users/jamestorpy/Documents/Garvan/phd"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/")
outDir <- paste0(resultsDir, "/R/", methodName, "/")


### 1. Load in counts ###

Counts <- readRDS(file=paste0(RobjectDir, annot, "_RepeatCounts/all_", annot, "RepeatCountDFs.rds"))
#  remove duplicate list element:
Counts <- c(Counts[1], Counts[3:length(Counts)])
# bind into df:
countsDF <- do.call("rbind", Counts)
# simplify the row names:
rownames(countsDF) <- gsub("^.*\\.", "", rownames(countsDF))


### 2. Load in GCcounts ###

gcFiles <- list.files(paste0(RobjectDir, annot, "_RepeatCounts"), pattern="GCcountsDF", full.names=T)
gcFiles <- grep("subset", gcFiles, value=T, invert=T)

gcL <- list()
j=1
for (f in gcFiles) {
	gcL[[j]] <- as.data.frame(readRDS(file = f)[,2])
	GCrownames <- readRDS(file = f)[,1]
	j=j+1
}

GCcountsDF <- do.call("cbind", gcL)
rownames(GCcountsDF) <- GCrownames
colnames(GCcountsDF) <- colnames(countsDF)

# append GCcountsDF to each GRanges object of custom1Counts:
Counts <- rbind(countsDF, GCcountsDF)


### 3. Perform pre-normalisation PCA and RLE plots ###

# eliminate lowly expressed genes (rows where there are less than 3 counts where df > 4):
print(paste0("No. rows before filtering is: ", nrow(Counts)))
  Counts <- Counts %>%
    rownames_to_column('repeat_id') %>%
    dplyr::filter(rowSums(Counts > 6) >= 9) %>%
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

if (file.exists(paste0(outDir, annot, "_pcaCompsPrenormGC.pdf"))) {
  print(paste0(outDir, annot, "_pcaCompsPrenormGC.pdf already exists, no eed to create"))
} else {
  print(paste0("Creating ", outDir, annot, "_pcaCompsPrenormGC.pdf"))
  pdf(file = paste0(outDir, annot, "_pcaCompsPrenormGC.pdf"))
  plot(pca)
  dev.off()
}
typeF <- as.factor(rep(STypes, each=3))

# convert matrix into SeqExpressionSet:
set <- newSeqExpressionSet(Counts, phenoData = data.frame(typeF, row.names=colnames(Counts)))

# create pre-norm RLE plot:
if (file.exists(paste0(outDir, annot, "_RLEPrenormGC.pdf"))) {
  print(paste0(outDir, annot, "_RLEPrenormGC.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", outDir, annot, "_RLEPrenormGC.pdf"))
  par(mar=c(1,1,1,1))
  pdf(file = paste0(outDir, annot, "_RLEPrenormGC.pdf"))
  plotRLE(set)
  dev.off()
}

# create RUVseq pre-norm PCA:
if (file.exists(paste0(outDir, annot, "_pcaPrenormGC.pdf"))) {
  print(paste0(outDir, annot, "_pcaPrenormGC.pdf already exists, no need to reate"))
} else {
  print(paste0("Creating ", outDir, annot, "_pcaPrenormGC.pdf"))
  pdf(file = paste0(outDir, annot, "_pcaPrenormGC.pdf"))
  plotPCA(set, cex=0.7)
  dev.off()
}


### 4. perform normalisation on counts using RUVseq:

# perform between lane full normalisation:
nSet <- betweenLaneNormalization(set, which="full")
pdf(file = paste0(outDir, annot, "_RLElaneNormGC.pdf"))
plotRLE(nSet, outline=FALSE, ylim=c(-4, 4))
dev.off()

pdf(file = paste0(outDir, annot, "_pcalaneNormGC.pdf"))
plotPCA(nSet, cex=0.7)
dev.off()


### 5. Perform differential expression comparing normalised FT controls to owtell PRs ###

genes <- rownames(Counts)
# design matrix specifying all samples as the thing to be compared to:
design <- model.matrix(~0+typeF, data=pData(nSet))
# convert set into a DGElist format which specifies library size and normalisation factors:
y <- DGEList(counts=counts(nSet), group=typeF)

# plot y in MDS plot:
pdf(file = paste0(outDir, annot, "_mdslaneNormGC.pdf"))
plotMDS(y)
dev.off()

temp <- y

# calculate the deviance residuals from a first-pass GLM regression of the counts on the co-variates of interest (p8 RUVseq manual):
# estimate dispersion:
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

# create BCV plot:
pdf(file = paste0(outDir, annot, "_BCVlaneNormGC.pdf"))
plotBCV(y)
dev.off()

# adjust values using dispersion:
fit <- glmFit(y, design)
# perform likelihood ratio test:
lrt_bowtell <- glmLRT(fit, contrast = c(-1, 1, rep(0, 4)))
lrt_gtx <- glmLRT(fit, contrast = c(0, 1, 0, 0, -1, 0))


### 5. Calculate differential expression values ###

# fetch all gene DE info, 
allGenes_bowtell <- as.data.frame(topTags(lrt_bowtell, n=Inf))
allGenes_gtx <- as.data.frame(topTags(lrt_gtx, n=Inf))
repGenes_bowtell <- allGenes_bowtell[grep("ENS", rownames(allGenes_bowtell), invert = T),]
repGenes_gtx <- allGenes_gtx[grep("ENS", rownames(allGenes_gtx), invert = T),]
sigGenes_bowtell <- filter(repGenes_bowtell, FDR<0.05)
sigGenes_gtx <- filter(repGenes_gtx, FDR<0.05)

# add negative log p-value column to allGenes:
#allGenes$negLog10PValue <- -log10(allGenes$PValue)

# plot on volcano plot:

  repGenes$threshold <- as.factor(repGenes$FDR < 0.1)
  sig <- subset(repGenes, threshold == T)
  sig$genes <- rownames(sig)
  p <- ggplot(data=repGenes, aes(x=logFC, y=-log10(FDR), color=threshold))
  p <- p + geom_point(data=repGenes)
  p <- p + geom_text_repel(data=sig, aes(label=genes))
  p <- p + theme(legend.position = "none")
  p <- p + labs(x="log2 fold change vs FT control", y="-log10 FDR")
  p <- p +  xlim(c(-6, 6))
  if (file.exists(paste0(outDir, "plots/", annot, "_volcanoFDR0.1", comp, ".pdf"))) {
    print(paste0(outDir, "plots/", annot, "_volcanoFDR0.1", comp, ".pdf already exists"))
    p
  } else {
    print(paste0("Creating ", outDir, "plots/", annot, "_volcanoFDR0.1", comp, ".pdf"))
    pdf(file = paste0(outDir, "plots/", annot, "_volcanoFDR0.1", comp, ".pdf"))
    print(p)
    dev.off()
  }



if (!file.exists(paste0(RobjectDir, "DEImg.RData"))) {
  save.image(file = paste0(RobjectDir, "DEImg.RData"))
}



