### 17.DESeq2.R ###

# This script takes a list of htseq-count outputs of different classes and types of repeat counts for
# HGSOC and FT control RNA-seq data sets and performs DE analysis:


### 0. Define variables/paths ###

# load packages needed:
library(DESeq2)
library(tibble)
library(dplyr)
library(RUVSeq)
library(ggplot2)
library("BiocParallel")
register(MulticoreParam(workers=10))
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
#homeDir <- "/Users/jamestorpy/clusterHome/"
homeDir <- "/share/ScratchGeneral/jamtor/"
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


### 2. Perform normalisation of counts by RUVseq ###

rnames <- rownames(Counts)
Counts <- as.data.frame(sapply(Counts, function(x) {
  storage.mode(x) <- "integer"
  return(x)
}))
rownames(Counts) <- rnames

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


# define condition groups:
condition <- gsub(
  "AOCS.*$", "cancer", gsub(
    "^.*_FT", "normal", colnames(Counts)
  )
)

typ <- rep("paired-end", ncol(Counts))

coldata <- data.frame(condition, typ)
rownames(coldata) <- colnames(Counts)
colnames(coldata) <- c("condition", "type")
  
# check rownames of coldata are the same as colnames of Counts:
paste0("Are rownames of coldata are identical to colnames of Counts? ", identical(rownames(coldata), colnames(Counts)))




### 2.  ###
  
dds <- DESeqDataSetFromMatrix(countData = counts(nSet),
                              colData = coldata,
                              design = ~ condition)

# eliminate lowly expressed genes (rows where there are less than 3 counts where df > 4):
print(paste0("No. rows before filtering is: ", nrow(Counts)))
keep <- rowSums(Counts > 5) >= (ncol(Counts)/3)
dds <- dds[keep,]
print(paste0("No. rows after  filtering: ", length(rownames(dds))))

# set factor levels - whatever is first will act as the control group:
dds$condition <- factor(dds$condition, levels = c("normal", "cancer"))

ddsDE <- DESeq(dds, parallel = TRUE)
pdf(paste0(plotDir, "/prenorm_RLE.pdf"))
boxplot(log10(assays(ddsDE)[["cooks"]]), range=0, las=2)
dev.off()

save.image(file=paste0(newRobjectDir, "/postDE.RData"))
#load(file=paste0(newRobjectDir, "/postDE.RData"))

resAsh <- lfcShrink(ddsDE, coef=2, type="ashr")
pdf(paste0(plotDir, "/reAsh_norm_RLE.pdf"))
boxplot(log10(assays(resAsh)[["cooks"]]), range=0, las=2)
dev.off()







# create PCA plot from dds:
vsdPre <- vst(dds)

prePCA <- plotPCA(vsdPre, intgroup=c("condition", "type"), returnData=TRUE)

p <- ggplot(prePCA, aes(PC1, PC2, color=condition))
p <- p + geom_point(size=1)

if (file.exists(paste0(plotDir, "/", Type, "_pcaPrenormGC.pdf"))) {
  print(paste0(plotDir, "/", Type, "_pcaPrenormGC.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", plotDir, "/", Type, "_pcaPrenormGC.pdf"))
  pdf(file = paste0(plotDir, "/", Type, "_pcaPrenormGC.pdf"))
  p
  dev.off()
}

# perform DE with DESeq2:
#parallel = TRUE

#ddsDEsh <- DESeq(dds, betaPrior=T)
#parallel = FALSE

res <- results(object = ddsDE)
summary(res)





# create PCA plot from dds:
vsdPost <- vst(dds)

postPCA <- plotPCA(vsdPost, intgroup=c("condition", "type"), returnData=TRUE)

p <- ggplot(postPCA, aes(PC1, PC2, color=condition))
p <- p + geom_point(size=1)

if (file.exists(paste0(plotDir, "/", Type, "_pcaPrenormGC.pdf"))) {
  print(paste0(plotDir, "/", Type, "_pcaPrenormGC.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", plotDir, "/", Type, "_pcaPrenormGC.pdf"))
  pdf(file = paste0(plotDir, "/", Type, "_pcaPrenormGC.pdf"))
  p
  dev.off()
}
