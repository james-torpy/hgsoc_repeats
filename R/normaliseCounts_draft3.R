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

# define starting variables:
project <- "hgsoc_repeats"
methodName <- "method1"
refName <- "human-89.repeats.tab"
sampleName <- "subset"
classType <- "c"
types <- c("bowtell_FT", "bowtell_PR", "grant_endo", "grant_endosis", "gtx", "kur")

# define directories:
#homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
homeDir <- "/Users/jamestorpy/Documents/Garvan/phd"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/")
outDir <- paste0(resultsDir, "/R/", methodName, "/")


### 1. Load in counts ###

cCounts <- readRDS(file = paste0(RobjectDir, classType, "_subsetRepeatCountsDF.rds"))
# define sampleNames:
sampleNames = colnames(cCounts[[1]])

### 2. Perform pre-normalisation PCA and RLE plots ###

# eliminate lowly expressed genes (rows where there are less than 3 counts where df > 4):
i=1
cCounts <- lapply(cCounts, function(x) {
  print(paste0("No. rows before filtering in ", names(cCounts)[i], " :", nrow(x)))
  x <- x %>%
  rownames_to_column('repeat_id') %>%
  dplyr::filter(rowSums(x > 4) >= 3) %>%
  column_to_rownames('repeat_id')
  print(paste0("No. rows after  filtering: ", nrow(x)))
  i <<- i+1
  return(x)
})

# convert cCounts list into df containing all data points for PCA analysis - this creates data.frame rows, so convert to matrix:
cCountsAll <- as.matrix(do.call("rbind", cCounts))

# define plot colour scheme:
cols <- brewer.pal(6, "Paired")

# create pre-normalised PCA plot from counts and plot:
pca <- princomp(cCountsAll)
if (file.exists(paste0(outDir, "pca1CompsSubsetPrenorm.pdf"))) {
  print(paste0(outDir, "pca1CompsSubsetPrenorm.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", outDir, "pca1CompsSubsetPrenorm.pdf"))
  pdf(file = paste0(outDir, "pca1CompsSubsetPrenorm.pdf"))
  plot(pca)
  dev.off()
}

# plot pca for components 1 and 2:
if (file.exists(paste0(outDir, "pca1SubsetPrenorm.pdf"))) {
  print(paste0(outDir, "pca1SubsetPrenorm.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", outDir, "pca1SubsetPrenorm.pdf"))
  pdf(file = paste0(outDir, "pca1SubsetPrenorm.pdf"))
  plot(pca$loading, pch=19, cex=1, col=cols)
  text(pca$loading, colnames(cCountsAll), pos=1, cex=0.7)
  dev.off()
}

# plot pca for components 1 and 3:
if (file.exists(paste0(outDir, "pca1SubsetComp1and3Prenorm.pdf"))) {
  print(paste0(outDir, "pca1SubsetComp1and3Prenorm.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", outDir, "pcaSubsetComp1and3Prenorm.pdf"))
  pdf(file = paste0(outDir, "pca1SubsetComp1and3Prenorm.pdf"))
  plot(x = pca$loading[,1], y = pca$loading[,3], pch=19, cex=1, col=cols)
  text(x = pca$loading[,1], y =  pca$loading[,3], colnames(cCountsAll), pos=1, cex=0.7)
  dev.off()
}

# make factor vector of sample types:
x <- as.factor(rep(types, each=3))
# convert matrix into SeqExpressionSet:
set <- newSeqExpressionSet(cCountsAll, phenoData = data.frame(x, row.names = colnames(cCountsAll)))

# create pre-norm RLE plot:
if (file.exists(paste0(outDir, "RLESubsetPrenorm.pdf"))) {
  print(paste0(outDir, "RLESubsetPrenorm.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", outDir, "RLESubsetPrenorm.pdf"))
  par(mar=c(1,1,1,1))
  pdf(file = paste0(outDir, "RLESubsetPrenorm.pdf"))
  plotRLE(set)
  dev.off()
}

# create RUVseq pre-norm PCA:
if (file.exists(paste0(outDir, "pca2SubsetPrenorm.pdf"))) {
  print(paste0(outDir, "pca2SubsetPrenorm.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", outDir, "pca2SubsetPrenorm.pdf"))
  pdf(file = paste0(outDir, "pca2SubsetPrenorm.pdf"))
  plotPCA(set, cex=0.7)
  dev.off()
}


### 3. perform normalisation on counts using RUVseq:


genes <- rownames(cCountsAll)
# normalise using upper-quartile normalisation (http://vinaykmittal.blogspot.com.au/2013/10/fpkmrpkm-normalization-caveat-and-upper.html)
# design matrix specifying all samples as the thing to be compared to:
design <- model.matrix(~x, data=pData(set))
# convert set into a DGElist format which specifies library size and normalisation factors:
y <- DGEList(counts=counts(set), group=x)
# calculate the normalisation factors for each sample:
y <- calcNormFactors(y, method="upperquartile")
# calculate the deviance residuals from a first-pass GLM regression of the counts on the co-variates of interest (p8 RUVseq manual):
# estimate dispersion:
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
# adjust values using dispersion:
fit <- glmFit(y, design)
# calculate residuals:
res <- residuals(fit, type="deviance")
# use all the genes to estimate the factors of unwanted variation:
cCountsNorm <- RUVr(set, genes, k=1, res)

# create pre-normalised PCA plot from counts and plot:
pcaNorm <- princomp(fit$fitted.values)

# plot normalised pca:
if (file.exists(paste0(outDir, "pca2SubsetRUVSeq.pdf"))) {
  print(paste0(outDir, "pca2SubsetRUVseq.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", outDir, "pca2SubsetRUVseq.pdf"))
  pdf(file = paste0(outDir, "pca2SubsetRUVseq.pdf"))
  plotPCA(cCountsNorm, cex=0.7)
  dev.off()
}

# create normalised RLE plot:
if (file.exists(paste0(outDir, "RLESubsetRUVseq.pdf"))) {
  print(paste0(outDir, "RLESubsetRUVseq.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", outDir, "RLESubsetRUVseq.pdf"))
  par(mar=c(1,1,1,1))
  pdf(file = paste0(outDir, "RLESubsetRUVseq.pdf"))
  plotRLE(cCountsNorm, ces=0.7)
  dev.off()
}


### 4. Perform differential expression comparing normalised FT controls to bowtell PRs ###

# remodel the matrix with the factors of unwanted variation:
design <- model.matrix(~x + W_1, data=pData(cCountsNorm))
# convert set into a DGElist format which specifies library size and normalisation factors:
y <- DGEList(counts=counts(cCountsNorm), group=x)
# calculate the normalisation factors for each sample:
y <- calcNormFactors(y, method="upperquartile")
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

create_rds(lrt, paste0(outDir, "bowtellFTvsPRSubset_RUVlrt.rds"))

# plot PCA on fitted values of RUVseq normalisation:
pca <- princomp(lrt$fitted.values)
if (file.exists(paste0(outDir, "pca1SubsetRUVseq.pdf"))) {
  print(paste0(outDir, "pca1SubsetRUVseq.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", outDir, "pca1SubsetRUVseq.pdf"))
  pdf(file = paste0(outDir, "pca1SubsetRUVseq.pdf"))
  plot(pca$loading, pch=19, cex=1, col=cols)
  text(pca$loading, colnames(cCountsAll), pos=1, cex=0.7)
  dev.off()
}


### 5. Perform differential expression comparing non-RUVseq normalised bowtell FTs to PRs ###

design <- model.matrix(~x, data=pData(set))
# convert set into a DGElist format which specifies library size and normalisation factors:
y <- DGEList(counts=counts(set), group=x)
# calculate the normalisation factors for each sample:
y <- calcNormFactors(y, method="upperquartile")
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

# plot PCA on fitted values of RUVseq normalisation:
pca <- princomp(lrt$fitted.values)
if (file.exists(paste0(outDir, "pca1SubsetEdgeR.pdf"))) {
  print(paste0(outDir, "pca1SubsetEdgeR.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", outDir, "pca1SubsetEdgeR.pdf"))
  pdf(file = paste0(outDir, "pca1SubsetEdgeR.pdf"))
  plot(pca$loading, pch=19, cex=1, col=cols)
  text(pca$loading, colnames(cCountsAll), pos=1, cex=0.7)
  dev.off()
}