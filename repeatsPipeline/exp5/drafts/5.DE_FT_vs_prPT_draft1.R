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
expName <- "exp5"
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
if (file.exists(paste0(plotDir, "htseq_RLEPrenormGC.pdf"))) {
  print(paste0(plotDir, "htseq_RLEPrenormGC.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", plotDir, "htseq_RLEPrenormGC.pdf"))
  par(mar=c(1,1,1,1))
  pdf(file = paste0(plotDir, "htseq_RLEPrenormGC.pdf"))
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
pdf(file = paste0(plotDir, "htseq_c_RLElaneNormGC.pdf"))
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

pdf(file = paste0(plotDir, "htseq_c_pcalaneNormGC.pdf"), height = 15, width = 20)
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


### 5. Calculate differential expression values ###

# PCA and RLE plots of RUVseq RUVr-normalised data looked   best with etween lane normalisation = 	'full', will go  with this #

# fetch summary of differentially expressed genes (those  with FDR =< 0.05:
DEs <- summary(result <- decideTestsDGE((lrt)))

# fetch all gene DE info, 
allGenes <- as.data.frame(topTags(lrt, n=Inf))
repGenes <- allGenes[grep("ENS", rownames(allGenes), invert   = T),]
sigGenes <- filter(repGenes, FDR<0.05)
print(repGenes)


# add negative log p-value column to allGenes:
#allGenes$negLog10PValue <- -log10(allGenes$PValue)

# plot on volcano plot:
repGenes$threshold <- as.factor(repGenes$FDR < 0.1)
sig <- subset(repGenes, threshold == T)
# include the control genes for labelling:
# add positive and negative control genes CD47, CCNE1,  GAPDH, b-actin:
posGenes <- rbind(allGenes["ENSG00000091831",],   allGenes["ENSG00000196776",], allGenes["ENSG00000105173",])
negGenes <- rbind(allGenes["ENSG00000111640",],   allGenes["ENSG00000075624",], allGenes["ENSG00000204574",],
                  allGenes["ENSG00000104904",],   allGenes["ENSG00000089157",])
rownames(posGenes) <- c("ESR1", "CD47", "CCNE1")
rownames(negGenes) <- c("GAPDH", "beta-actin", "ABCF1",   "OAZ1", "RPLP0")

posGenes$threshold <- "POSITIVE"
if (nrow(posGenes[posGenes$FDR<0.1,])>0) {
  posGenes[posGenes$FDR<0.1,]$threshold <- "POSSIG"
}

negGenes$threshold = "NEGATIVE"
if (nrow(negGenes[negGenes$FDR<0.1,])>0) {
  negGenes[negGenes$FDR<0.1,]$threshold <- "NEGSIG"
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
if (file.exists(paste0(plotDir, "/", Type,  "_volcano_FDR0.1_", comp, ".pdf"))) {
  print(paste0(plotDir, "/", Type, "_volcano_FDR0.1_",  comp, ".pdf"))
  p
} else {
  print(paste0("Creating ",plotDir, "/", Type,  "_volcano_FDR0.1_", comp, ".pdf"))
  pdf(file = paste0(plotDir, "/", Type, "_volcano_FDR0.1_",   comp, ".pdf"))
  print(p)
  dev.off()
}

repGenes <- allGenes[grep("ENS", rownames(allGenes), invert   = T),]
sigGenes <- filter(repGenes, FDR<0.05)
repGenes$threshold <- as.factor(repGenes$FDR < 0.05)
sig <- subset(repGenes, threshold == T)
posGenes$threshold <- "POSITIVE"

if (nrow(posGenes[posGenes$FDR<0.05,])>0) {
  posGenes[posGenes$FDR<0.05,]$threshold <- "POSSIG"
}

negGenes$threshold = "NEGATIVE"
if (nrow(negGenes[negGenes$FDR<0.05,])>0) {
  negGenes[negGenes$FDR<0.05,]$threshold <- "NEGSIG"
}

lab <- rbind(rbind(sig, posGenes), negGenes)
repGenes <- rbind(rbind(repGenes, posGenes), negGenes)
lab$genes <- rownames(lab)

p <- ggplot(data=repGenes, aes(x=logFC, y=-log10(FDR),  color=threshold))
p <- p + geom_point(data=repGenes)
p <- p + geom_text_repel(data=lab, aes(label=genes))
p <- p + theme(legend.position = "none")
p <- p + labs(x="log2 fold change vs FT control", y="-log10   FDR")
p <- p +  xlim(c(-5, 5))
if (file.exists(paste0(plotDir, "/", Type,  "_volcano_FDR0.05_", comp, ".pdf"))) {
  print(paste0(plotDir, "/", Type, "_volcano_FDR0.05_",   comp, ".pdf"))
  p
} else {
  print(paste0("Creating ", plotDir, "/", Type,   "_volcano_FDR0.05_", comp, ".pdf"))
  pdf(file = paste0(plotDir, "/", Type,   "_volcano_FDR0.05_", comp, ".pdf"))
  print(p)
  dev.off()
}

repGenes <- allGenes[grep("ENS", rownames(allGenes), invert   = T),]
repGenes$threshold <- as.factor(repGenes$FDR < 0.05 & ( repGenes$logFC>=1 | repGenes$logFC<=-1))
sig <- subset(repGenes, threshold == T)

posGenes$threshold <- "POSITIVE"
if (any(posGenes$FDR<0.05 & (posGenes$logFC>=1 |  posGenes$logFC<=-1))) {
  posGenes[posGenes$FDR<0.05 & (posGenes$logFC>=1 |   posGenes$logFC<=-1),]$threshold <- "POSSIG"
}
negGenes$threshold = "NEGATIVE"
if (any(negGenes$FDR<0.05 & (negGenes$logFC>=1 |  negGenes$logFC<=-1))) {
  negGenes[negGenes$FDR<0.05 & (negGenes$logFC>=1 |   negGenes$logFC<=-1),]$threshold <- "NEGSIG"
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
if (file.exists(paste0(plotDir, "/", Type,  "_volcano_FDR0.05_FC2_", comp, ".pdf"))) {
  print(paste0(plotDir, "/", Type, "_volcano_FDR0.05_FC2_",   comp, ".pdf"))
  p
} else {
  print(paste0("Creating ", plotDir, "/", Type,   "_volcano_FDR0.05_FC2_", comp, ".pdf"))
  pdf(file = paste0(plotDir, "/", Type,   "_volcano_FDR0.05_FC2_", comp, ".pdf")) 
  print(p)
  dev.off()
}

}  
con <- c(0, 0, -1, rep(0, 6))
}

names(repGeneL) <- c(CTypes[1:2], CTypes[4:length(CTypes)])

for (i in 1:length(repGeneL)) {
  repGeneL[[i]]$threshold <- as.factor(repGeneL[[i]]$FDR < 0.05   & (repGeneL[[i]]$logFC>=1.5 | 	repGeneL[[i]]$logFC<=-1.5))
  sig <- subset(repGeneL[[i]], threshold == T)
  sig$genes <- rownames(sig)
  p <- ggplot(data=repGeneL[[i]], aes(x=logFC, y=-log10(FDR),   color=threshold))
  p <- p + geom_point(data=repGeneL[[i]])
  p <- p + geom_text_repel(data=lab, aes(label=genes))
  p <- p + theme(legend.position = "none")
  p <- p + labs(x="log2 fold change vs FT control", y="-log10   FDR")
  p <- p +  xlim(c(-7, 7))
  if (file.exists(paste0(plotDir, "/", Type,  "_volcano_FDR0.05_FC_log1.5_", names(repGeneL)[i], 	"_vs_FT.pdf"))) {
    print(paste0(plotDir, "/", Type,  "_volcano_FDR0.05_FC_log1.5_", names(repGeneL)[i], 	"_vs_FT.pdf"))
    p
  } else {
    print(paste0("Creating ", plotDir, "/", Type,   "_volcano_FDR0.05_FC_log1.5_", names(repGeneL)[i], 	  "_vs_FT.pdf"))
    pdf(file = paste0(plotDir, "/", Type,   "_volcano_FDR0.05_FC_log1.5_", names(repGeneL)[i], 	  "_vs_FT.pdf"))
    print(p)
    dev.off()
  }  

saveRDS(repGeneL, file=paste0(RobjectDir, "/", Type, "_DEreps.rds"))


topRep <- lapply(repGeneL, function(x) {
  result <- head(x[x$FDR<0.1, ], 15)
  return(rbind(result, x[c("ESR1", "CD47", "CCNE1", "GAPDH", "beta-actin", "ABCF1", "OAZ1", "RPLP0"),]))
})


strRep <- lapply(repGeneL, function(x) {
  result <- head(x[x$FDR<0.1 & (x$logFC < -1 | x$logFC > 1), ], 15)
  return(rbind(result, x[c("ESR1", "CD47", "CCNE1", "GAPDH", "beta-actin", "ABCF1", "OAZ1", "RPLP0"),]))
})

saveRDS(topRep, file=paste0(RobjectDir, "/", Type, "_topRep.rds"))
saveRDS(strRep, file=paste0(RobjectDir, "/", Type, "_strRep.rds"))

if (!file.exists(paste0(RobjectDir, "DEImg_", expName, ".RData"))) {
  save.image(file = paste0(RobjectDir, "DEImg_", expName, ".RData"))
}



