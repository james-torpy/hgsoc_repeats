### 13.DE_FT_vs_HGSOC.R ###

# This script takes a list of dfs of different classes and types of repeat counts for
# HGSOC and FT control RNA-seq data sets and performs DE analysis:

### 0. Define variables/paths ###

# load packages needed:
library(DESeq2)
library(tibble)
library(dplyr)
library(RColorBrewer)
library(RUVSeq)
library(ggplot2)
library(ggrepel)
library("BiocParallel")
register(MulticoreParam(workers=12))
#request 12 cores on cluster for parallel processing with 10 workers

# define starting variables:
project <- "hgsoc_repeats"
expName <- "exp9"

Type <- "custom3"
descrip <- "htseq_single_HGSOCs_vs_FT"


# define sample group to use as control:
ctl <- "FT"

# specify what combination of repeat genes (repeats), epigenetic modulators (epiMods),
# RNAi genes (RNAi) and protein-coding genes (pCoding) should contribute to the results:
#resultTypes <- c("repeats", "epiMods")
resultTypes <- c("repeats")

# specify what FDR and log2 fold change thresholds to use:
FDRthresh <- 0.05
FCthresh <- 0

# specify control genes to include:
posGeneIDs <- c("ENSG00000111640", "ENSG00000196776")
posGeneNames <- c("GAPDH", "CD47")
negGeneIDs <- c("ENSG00000075624", "ENSG00000169919")
negGeneNames <- c("beta-actin", "GUSB")


# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project)
rawDir <- paste0("/Users/jamestorpy/clusterHome2/projects/hgsoc_repeats/RNA-seq/raw_files")
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
                     expName, "/")
newRobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
                        expName, "/", descrip, "/")
plotDir <- paste0(resultsDir, "/R/", expName,
                  "/plots/DEplots/", descrip, "/")

system(paste0("mkdir -p ", plotDir))
system(paste0("mkdir -p ", newRobjectDir))


### 1. Load in all counts ###

custom3Counts <- readRDS(paste0(RobjectDir, "/", Type,
                                "_allcounts.htseq.rds"))
gcCounts <- readRDS(paste0(RobjectDir,
                           "/gc_allcounts.htseq.rds"))

# append gcCounts to custom3Counts:
Counts <- rbind(custom3Counts, gcCounts)

# make rownames gene_id, get rid of latter column and change
# storage mode from factor to integer:
rownames(Counts) <- Counts$gene_id
Counts <- subset(Counts, select=-gene_id)

# save Counts:
saveRDS(Counts, file=paste0(newRobjectDir, "/", Type, "_counts.RData"))


### 2. Perform pre-normalisation PCA and RLE plots ###

# eliminate lowly expressed genes (rows where there are less than 3 counts where df > 4):
print(paste0("No. rows before filtering is: ", nrow(Counts)))
Counts <- Counts %>%
  rownames_to_column('gene_id') %>%
  dplyr::filter(rowSums(Counts > 5) >= (ncol(Counts)/3)) %>%
  column_to_rownames('gene_id')
print(paste0("No. rows after  filtering: ", nrow(Counts)))

# create pre-normalised PCA plot from counts and plot:
if (ncol(Counts) > nrow(Counts)) {
  pca <- prcomp(Counts)
} else {
  pca <- princomp(Counts)   
}

if (file.exists(paste0(plotDir, "/", Type, "_pcaCompsPrenormGC.pdf"))) {
  print(paste0(plotDir, "/", Type, "_pcaCompsPrenormGC.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", plotDir, "/", Type, "_pcaCompsPrenormGC.pdf"))
  pdf(file = paste0(plotDir, "/", Type, "_pcaCompsPrenormGC.pdf"))
  plot(pca)
  dev.off()
}

# remove pAF and remove categorisations of HGSOC samples:
Counts <- Counts[,-(grep("pAF", colnames(Counts)))]
#colnames(Counts) <- gsub("_[a-z][a-z][A-Z][A-Z]$", "", colnames(Counts))
# append '.2' onto duplicate patient IDs:
#colnames(Counts)[duplicated(colnames(Counts))] <- paste0(colnames(Counts)[duplicated(colnames(Counts))], ".2")

# change the order of columns of Counts to 0.1betical order of subtypes:
Counts <- Counts[,order(
  gsub(
    "AOCS.*_[0-9][0-9][0-9]_", "", colnames(Counts)
  )
)]



# define sample groups:
splt <- unlist(
  lapply(
    split(
      colnames(Counts), gsub(
        "AOCS.*_[0-9][0-9][0-9]_FT", "FT", colnames(Counts)
      )
    ), length
  )
)

typeF <- factor(gsub("^.*FT", "FT", colnames(Counts)), levels = c("FT", colnames(Counts)[!(colnames(Counts) %in% "FT")]))

# save number of samples in each group:
saveRDS(splt, file = paste0(newRobjectDir, "/sample_no_per_cat.rds"))

# delist elements need to be delisted and change to integers:
Counts <- apply(Counts, 2, unlist)
storage.mode(Counts) <- "integer"
# add '.2' to IDs of samples with duplicate names:
colnames(Counts)[duplicated(colnames(Counts))] <- gsub("_HGSOC", ".2_HGSOC", colnames(Counts)[duplicated(colnames(Counts))])
# convert Counts into SeqExpressionSet
set <- newSeqExpressionSet(Counts, phenoData = data.frame(typeF, row.names=colnames(Counts)))

# create pre-norm RLE plot:
if (file.exists(paste0(plotDir, "/", Type, "_RLEPrenormGC.pdf"))) {
  print(paste0(plotDir, "/", Type, "_RLEPrenormGC.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", plotDir, "/", Type, "_RLEPrenormGC.pdf"))
  par(mar=c(1,1,1,1))
  pdf(file = paste0(plotDir, "/", Type, "_RLEPrenormGC.pdf"))
  plotRLE(set)
  dev.off()
}

# create RUVseq pre-norm PCA:
if (file.exists(paste0(plotDir, "/", Type, "_pcaPrenormGC.pdf"))) {
  print(paste0(plotDir, "/", Type, "_pcaPrenormGC.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", plotDir, "/", Type, "_pcaPrenormGC.pdf"))
  pdf(file = paste0(plotDir, "/", Type, "_pcaPrenormGC.pdf"), height = 10, width = 12)
  plotPCA(set, cex=0.7)
  dev.off()
}


### 3. perform normalisation on counts using RUVseq:

# perform within lane normalisation:
dataWithin <- withinLaneNormalization(set,"gc", which="full")

# perform between lane full normalisation:
nSet <- betweenLaneNormalization(set, which="full")

# create post-norm RLE plot:
if (file.exists(paste0(plotDir, "/", Type, "_RLElaneNormGC.pdf"))) {
  print(paste0(plotDir, "/", Type, "_RLElaneNormGC.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", plotDir, "/", Type, "_RLElaneNormGC.pdf"))
  pdf(file = paste0(plotDir, "/", Type, "_RLElaneNormGC.pdf"))
  plotRLE(nSet, outline=FALSE, ylim=c(-4, 4))
  dev.off()
}

# create RUVseq post-norm PCA:
if (file.exists(paste0(plotDir, "/", Type, "_pcalaneNormGC.pdf"))) {
  print(paste0(plotDir, "/", Type, "_pcalaneNormGC.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", plotDir, "/", Type, "_pcalaneNormGC.pdf"))
  pdf(file = paste0(plotDir, "/", Type, "_pcalaneNormGC.pdf"), height = 15, width = 20)
  plotPCA(nSet, cex=0.7)
  dev.off()
}


### 4. Perform DESeq2 on counts ###

# define condition groups:
#Counts <- assayData(nSet)$normalizedCounts
Counts <- counts(nSet)

condition <- gsub("^.*FT", "FT",colnames(Counts))

typ <- rep("paired-end", ncol(Counts))

coldata <- data.frame(condition, typ)
rownames(coldata) <- colnames(Counts)
colnames(coldata) <- c("condition", "type")
  
# check rownames of coldata are the same as colnames of Counts:
paste0("Are rownames of coldata are identical to colnames of Counts? ", identical(rownames(coldata), colnames(Counts)))


### 2.  ###
  
dds <- DESeqDataSetFromMatrix(countData = Counts,
                              colData = coldata,
                              design = ~ condition)


# set factor levels - whatever is first will act as the control group:
dds$condition <- factor(dds$condition, levels = c("FT", grep("FT", condition, value=T, invert=T)))

ddsDE <- DESeq(dds, parallel = TRUE)
pdf(paste0(plotDir, "/postDESeq2_RLE.pdf"))
boxplot(log10(assays(ddsDE)[["cooks"]]), range=0, las=2, cex.axis=0.3)
dev.off()

save.image(file=paste0(RobjectDir, "/post_DESeq.RData"))

register(MulticoreParam(workers=10))
for (i in 3:length(resultsNames(ddsDE))) {
  print(i)
  if (i==2) {
    resLFC <- list(lfcShrink(ddsDE, coef=i, parallel=TRUE))
  } else if (i>2) {
    resLFC[[i]] <- lfcShrink(ddsDE, coef=i, parallel=TRUE)
  }
}
#names(resLFC) <- resultsNames(ddsDE)[1:length(resultsNames(ddsDE))]





#res <- results(dds)

# filter out significant genes based on FDR adjusted p-values
filtered <- resLFC[(resLFC$padj < 0.1) & !is.infinite(resLFC$log2FoldChange) & !is.nan(resLFC$log2FoldChange) & !is.na(resLFC$padj),]
# order by p-value, and print out only the gene name, mean count, and log2 fold change
sorted <- filtered[order(filtered$padj),c(1,2,6)]


