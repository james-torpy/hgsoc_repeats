### 5.DE_FT_vs_prPT.R ###

# This script takes a list of dfs of different classes and types of repeat counts for
# primary resistant and FT control RNA-seq data sets and performs DE analysis:


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
expName <- "exp7"
STypes <- c("FT", "prPT")
Type <- "custom3"

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
#homeDir <- "/Users/jamestorpy/Documents/Garvan/phd"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/", expName, "/")
plotDir <- paste0(resultsDir, "/R/", expName, "/plots/DEplots/")

system(paste0("mkdir -p ", plotDir))


### 1. Load in all counts ###

custom3Counts <- readRDS(paste0(RobjectDir, "/", Type, "_allcounts.htseq.rds"))
gcCounts <- readRDS(paste0(RobjectDir, "/gc_allcounts.htseq.rds"))

# append GCcountsDF to each GRanges object of custom1Counts:
Counts <- rbind(custom3Counts, gcCounts)

# make rownames gene_id, get rid of latter column and change
# storage mode from factor to integer:
Rnames <- Counts$gene_id
Counts <- sapply(subset(Counts, select=-gene_id), unlist)
storage.mode(Counts) <- "integer"
Counts <- as.data.frame(Counts)
rownames(Counts) <- Rnames

# save Counts:
saveRDS(Counts, file=paste0(RobjectDir, "/", Type, "_counts.RData"))


### 3. Perform pre-normalisation PCA and RLE plots ###

# eliminate lowly expressed genes (rows where there are less than 3 counts where df > 4):
print(paste0("No. rows before filtering is: ", nrow(Counts)))
Counts <- Counts %>%
  rownames_to_column('gene_id') %>%
  dplyr::filter(rowSums(Counts > 5) >= (ncol(Counts)/3)) %>%
  column_to_rownames('gene_id')
print(paste0("No. rows after  filtering: ", nrow(Counts)))

# create pre-normalised PCA plot from counts and plot:
Counts <- apply(Counts, 2, unlist)
if (ncol(Counts) > nrow(Counts)) {
  pca <- prcomp(Counts)
} else {
  pca <- princomp(Counts)	  
}

if (file.exists(paste0(plotDir, "htseq_", Type, "_pcaCompsPrenormGC.pdf"))) {
  print(paste0(plotDir, "htseq_", Type, "_pcaCompsPrenormGC.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", plotDir, "htseq_", Type, "_pcaCompsPrenormGC.pdf"))
  pdf(file = paste0(plotDir, "htseq_", Type, "_pcaCompsPrenormGC.pdf"))
  plot(pca)
  dev.off()
}

splt <- unlist(lapply(split(colnames(Counts), gsub("^.*_*_", "", colnames(Counts))), length))
for (i in 1:length(splt)) {
  if (i==1) {
    typeF <- c(rep(names(splt)[i], splt[i]))
  } else {
    typeF <- c(typeF, rep(names(splt)[i], splt[i]))
  }
}
levels(typeF) <- c("FT", "prPT")

# convert matrix into SeqExpressionSet:
set <- newSeqExpressionSet(Counts, phenoData = data.frame(typeF, row.names=colnames(Counts)))

# create pre-norm RLE plot:
if (file.exists(paste0(plotDir, "htseq_", Type, "_RLEPrenormGC.pdf"))) {
  print(paste0(plotDir, "htseq_", Type, "_RLEPrenormGC.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", plotDir, "htseq_", Type, "_RLEPrenormGC.pdf"))
  par(mar=c(1,1,1,1))
  pdf(file = paste0(plotDir, "htseq_", Type, "_RLEPrenormGC.pdf"))
  plotRLE(set)
  dev.off()
}

# create RUVseq pre-norm PCA:
if (file.exists(paste0(plotDir, "htseq_", Type, "_pcaPrenormGC.pdf"))) {
  print(paste0(plotDir, "htseq_", Type, "_pcaPrenormGC.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", plotDir, "htseq_", Type, "_pcaPrenormGC.pdf"))
  pdf(file = paste0(plotDir, "htseq_", Type, "_pcaPrenormGC.pdf"), height = 10, width = 12)
  plotPCA(set, cex=0.7)
  dev.off()
}


### 4. perform normalisation on counts using RUVseq:

# perform between lane full normalisation:
nSet <- betweenLaneNormalization(set, which="full")
pdf(file = paste0(plotDir, "htseq_", Type, "_RLElaneNormGC.pdf"))
plotRLE(nSet, outline=FALSE, ylim=c(-4, 4))
dev.off()

######
dCounts <- Counts
storage.mode(dCounts) <- "double"
nCounts <- normalize.quantiles(dCounts)
colnames(nCounts) <- colnames(Counts)
rownames(nCounts) <- rownames(Counts)
pdf(file = paste0(plotDir, "htseq_", Type, "_RLElaneNormGC.pdf"))
plotRLE(nCounts, outline=FALSE, ylim=c(-4, 4))
dev.off()

storage.mode(nCounts) <- "integer"
nset <- newSeqExpressionSet(nCounts, phenoData = data.frame(typeF, row.names=colnames(nCounts)))
pdf(file = paste0(plotDir, "htseq_c_pcalaneNormGC.pdf"), height = 15, width = 20)
plotPCA(nCounts, cex=0.7)
dev.off()


# remove datasets that didn't normalise properly and re-normalise:
plotDir <- paste0(plotDir, "/strict_normalisation/")
system(paste0("mkdir -p ", plotDir))

Counts <- Counts[,!colnames(Counts) %in% c("AOCS_175_FT", "AOCS_083_prPT", "AOCS_085_prPT", "AOCS_107_prPT", "AOCS_111_prPT", "AOCS_113_prPT", 
                                           "AOCS_116_prPT", "AOCS_163_prPT")]
dCounts <- Counts
storage.mode(dCounts) <- "double"
nCounts <- normalize.quantiles.robust(dCounts)
colnames(nCounts) <- colnames(Counts)
rownames(nCounts) <- rownames(Counts)
pdf(file = paste0(plotDir, "htseq_", Type, "_RLElaneNormGC.pdf"))
plotRLE(nCounts, outline=FALSE, ylim=c(-4, 4))
dev.off()

storage.mode(nCounts) <- "integer"
nset <- newSeqExpressionSet(nCounts, phenoData = data.frame(typeF, row.names=colnames(nCounts)))
pdf(file = paste0(plotDir, "htseq_", Type, "_pcalaneNormGC.pdf"), height = 15, width = 20)
plotPCA(nCounts, cex=0.7)
dev.off()
######

### 5. Perform differential expression comparing normalised FT controls to cancer samples ###

genes <- rownames(Counts)
# normalise using upper-quartile normalisation (http://inaykmittal.blogspot.com.au/2013/10/	pkmrpkm-normalization-caveat-and-upper.html)
# design matrix specifying all samples as the thing to be compared to:
design <- model.matrix(~0+typeF, data=pData(nSet))
# convert set into a DGElist format which specifies library size and normalisation factors:
y <- DGEList(counts=counts(nSet), group=typeF)

# create MDS plot:
pdf(file = paste0(plotDir, Type, "_mdslaneNormGC.pdf"), height = 15, width = 20)
plotMDS(y)
dev.off()


# calculate the deviance residuals from a first-pass GLM regression of the counts on the co-variates 	of interest (p8 RUVseq manual):
# estimate dispersion:
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
# adjust values using dispersion:
fit <- glmFit(y, design)

saveRDS(fit, file=paste0(RobjectDir, "/", Type, "DEfit.rds"))

con <- c(-1, 1)

comp <- paste0("FT_vs_PT")

# perform likelihood ratio test:
con[i] <- 1
lrt <- glmLRT(fit, contrast = con)

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

create_rds(lrt, paste0(RobjectDir, "FT_vs_", CTypes[i],   "_lrt.rds"))


### 5. Calculate differential expression values with FDR threshold only ###

# PCA and RLE plots of RUVseq RUVr-normalised data looked   best with etween lane normalisation = 	'full', will go  with this #

# fetch summary of differentially expressed genes (those  with FDR =< 0.05:
DEs <- summary(result <- decideTestsDGE((lrt)))

# fetch all gene DE info, 
thresh <- 0.05
allGenes <- as.data.frame(topTags(lrt, n=Inf))
repGenes <- allGenes[grep("ENS", rownames(allGenes), invert   = T),]
sigGenes <- filter(repGenes, FDR<thresh)
print(repGenes)


# add negative log p-value column to allGenes:
#allGenes$negLog10PValue <- -log10(allGenes$PValue)

# plot on volcano plot:
repGenes$threshold <- as.factor(repGenes$FDR < thresh)
sig <- subset(repGenes, threshold == T)
# include the control genes for labelling:
# add positive and negative control genes CD47, CCNE1,  GAPDH, b-actin:

posGenes <- rbind(allGenes["ENSG00000091831",],   allGenes["ENSG00000196776",], allGenes["ENSG00000105173",])
negGenes <- rbind(allGenes["ENSG00000111640",],   allGenes["ENSG00000075624",], allGenes["ENSG00000204574",],
                  allGenes["ENSG00000104904",],   allGenes["ENSG00000089157",])
rownames(posGenes) <- c("ESR1", "CD47", "CCNE1")
rownames(negGenes) <- c("GAPDH", "beta-actin", "ABCF1",   "OAZ1", "RPLP0")

posGenes$threshold <- "POSITIVE"
if (nrow(posGenes[posGenes$FDR<thresh,])>0) {
  posGenes[posGenes$FDR<thresh,]$threshold <- "POSSIG"
}

negGenes$threshold = "NEGATIVE"
if (nrow(negGenes[negGenes$FDR<thresh,])>0) {
  negGenes[negGenes$FDR<thresh,]$threshold <- "NEGSIG"
}

lab <- rbind(rbind(sig, posGenes), negGenes)
repGenes <- rbind(rbind(repGenes, posGenes), negGenes)
lab$genes <- rownames(lab)


p <- ggplot(data=repGenes, aes(x=logFC, y=-log10(FDR),  color=threshold))
p <- p + geom_point(data=repGenes)
p <- p + geom_text_repel(data=lab, aes(label=genes))
p <- p + theme(legend.position = "none")
p <- p + labs(x="log2 fold change vs FT control", y="-log10   FDR")
p <- p +  xlim(c(-7, 7))
if (file.exists(paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, comp, ".pdf"))) {
  print(paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, comp, ".pdf"))
  p
} else {
  print(paste0("Creating ",plotDir, "/", Type,  "_volcano_FDR", thresh, comp, ".pdf"))
  pdf(file = paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, comp, ".pdf"))
  print(p)
  dev.off()
}


### 6. Calculate differential expression values with FDR and FC thresholds ###

# PCA and RLE plots of RUVseq RUVr-normalised data looked   best with etween lane normalisation = 	'full', will go  with this #

# fetch summary of differentially expressed genes (those  with FDR =< 0.05:
DEs <- summary(result <- decideTestsDGE((lrt)))

# fetch all gene DE info, 
thresh <- 0.05
FCthresh <- 2
allGenes <- as.data.frame(topTags(lrt, n=Inf))
repGenes <- allGenes[grep("ENS", rownames(allGenes), invert = T),]
sigGenes <- filter(repGenes, (FDR < thresh & logFC < -(FCthresh))|(FDR < thresh & logFC > FCthresh))
print(repGenes)


# add negative log p-value column to allGenes:
#allGenes$negLog10PValue <- -log10(allGenes$PValue)

# plot on volcano plot:
repGenes$threshold <- as.factor((repGenes$FDR < thresh & repGenes$logFC < -(FCthresh))|(repGenes$FDR < thresh & repGenes$logFC > FCthresh))
sig <- subset(repGenes, threshold == T)
# include the control genes for labelling:
# add positive and negative control genes CD47, CCNE1,  GAPDH, b-actin:

posGenes <- rbind(allGenes["ENSG00000091831",],   allGenes["ENSG00000196776",], allGenes["ENSG00000105173",])
negGenes <- rbind(allGenes["ENSG00000111640",],   allGenes["ENSG00000075624",], allGenes["ENSG00000204574",],
                  allGenes["ENSG00000104904",],   allGenes["ENSG00000089157",])
rownames(posGenes) <- c("ESR1", "CD47", "CCNE1")
rownames(negGenes) <- c("GAPDH", "beta-actin", "ABCF1",   "OAZ1", "RPLP0")

posGenes$threshold <- "POSITIVE"
if (nrow(posGenes[posGenes$FDR<thresh,])>0) {
  posGenes[posGenes$FDR<thresh,]$threshold <- "POSSIG"
}

negGenes$threshold = "NEGATIVE"
if (nrow(negGenes[negGenes$FDR<thresh,])>0) {
  negGenes[negGenes$FDR<thresh,]$threshold <- "NEGSIG"
}

lab <- rbind(rbind(sig, posGenes), negGenes)
repGenes <- rbind(rbind(repGenes, posGenes), negGenes)
lab$genes <- rownames(lab)

gcGenes <- allGenes[grep("ENS", rownames(allGenes)),]
gcGenes$threshold <- "GC"

# choose colours:
cols <- c("#77AADD", "deeppink1", "darkgoldenrod3", "forestgreen", "#114477", "darkgoldenrod3")

p <- ggplot(data=repGenes, aes(x=logFC, y=-log10(FDR),  color=threshold))
p <- p + geom_point(data=gcGenes, color="#C6C6C6", alpha=0.5)
p <- p + geom_point(data=repGenes)
p <- p + scale_color_manual(values=cols)
p <- p + geom_text_repel(data=lab, aes(label=genes))
#p <- p + theme(legend.position = "none")
p <- p + labs(x="log2 fold change vs FT control", y="-log10   FDR")
p <- p + xlim(c(-7, 7)) + ylim(c(0, 10))
if (file.exists(paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, "_logFC", FCthresh, "_", comp, ".pdf"))) {
  print(paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, "_logFC", FCthresh, "_", comp, ".pdf"))
  p
} else {
  print(paste0("Creating ",plotDir, "/", Type,  "_volcano_FDR", thresh, "_logFC", FCthresh, "_", comp, ".pdf"))
  pdf(file = paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, "_logFC", FCthresh, "_", comp, ".pdf"))
  print(p)
  dev.off()
}


saveRDS(repGenes, file=paste0(RobjectDir, "/", Type, "_DEreps", thresh, "_logFC", FCthresh, "_", comp, ".rds"))


#topRep <- lapply(repGeneL, function(x) {
#  result <- head(x[x$FDR<0.1, ], 15)
#  return(rbind(result, x[c("ESR1", "CD47", "CCNE1", "GAPDH", "beta-actin", "ABCF1", "OAZ1", "RPLP0"),]))
#})


#strRep <- lapply(repGeneL, function(x) {
#  result <- head(x[x$FDR<0.1 & (x$logFC < -1 | x$logFC > 1), ], 15)
#  return(rbind(result, x[c("ESR1", "CD47", "CCNE1", "GAPDH", "beta-actin", "ABCF1", "OAZ1", "RPLP0"),]))
#})

#saveRDS(topRep, file=paste0(RobjectDir, "/", Type, "_topRep.rds"))
#saveRDS(strRep, file=paste0(RobjectDir, "/", Type, "_strRep.rds"))

if (!file.exists(paste0(RobjectDir, "DEImg_", expName, ".RData"))) {
  save.image(file = paste0(RobjectDir, "DEImg_", expName, ".RData"))
}

load(file = paste0(RobjectDir, "DEImg_", expName, ".RData"))

