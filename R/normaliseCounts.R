### normaliseCounts.R ###

# This script takes a list of dfs of different classes and types of repeat counts and
# normalises using RUVseq, EdgeR and DEseq, then compares methods

### 0. Define variables/paths ###

# load packages needed:
library(tibble)
library(dplyr)
library(RColorBrewer)
library(ggplot2)

# define starting variables:
project <- "hgsoc_repeats"
methodName <- "method1"
refName <- "human-89.repeats.tab"
sampleName <- "subset"
classType <- "c"

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

### 2. Perform pre-normalisation PCA plots ###

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
pdf(file = paste0(outDir, "pcaCompsSubset.pdf"))
pca <- princomp(cCountsAll)
plot(pca)
dev.off()

# plot pca for components 1 and 2:
pdf(file = paste0(outDir, "pcaSubset.pdf"))
plot(pca$loading, pch=19, cex=1, col=cols)
text(pca$loading, colnames(cCountsAll), pos=1, cex=0.7)
dev.off()

# plot pca for components 1 and 3:
pdf(file = paste0(outDir, "pcaSubset1and3.pdf"))
plot(x = pca$loading[,1], y = pca$loading[,3], pch=19, cex=1, col=cols)
text(x = pca$loading[,1], y =  pca$loading[,3], colnames(cCountsAll), pos=1, cex=0.7)
dev.off()





# select one df for initial testing:
xCounts <- cCounts[[1]]