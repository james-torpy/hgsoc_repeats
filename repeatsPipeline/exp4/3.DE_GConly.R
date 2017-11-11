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
expName <- "exp4"
STypes <- c("FT", "prPT")

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
#homeDir <- "/Users/jamestorpy/Documents/Garvan/phd"
projectDir <- paste0(homeDir, "/projects/", project, 
  "/RNA-seq")
resultsDir <- paste0(projectDir, "/results")
inDir <- paste0(resultsDir, "/rsem/", expName)
RobjectDir <- paste0(projectDir, "/Robjects/", expName, 
  "/")
plotDir <- paste0(resultsDir, "/R/", expName, 
  "/plots/DEplots/")

system(paste0("mkdir -p ", plotDir))
system(paste0("mkdir -p ", RobjectDir))

  
### 2. Load in GCcounts ###

GCfiles <- list.files(inDir, full.names=T, recursive=T, 
  pattern="genes")
GCfiles <- grep("subset", GCfiles, value=T, invert=T)

gcL <- list()
j=1
for (f in GCfiles) {
  df <- read.table(file = f, header=T)
  gcL[[j]] <- data.frame(df[,1], as.numeric(df[,5]))
  j=j+1
}

GCcountsDF <- do.call("cbind", gcL)
rownames(GCcountsDF) <- GCcountsDF[,1]    

Counts <- GCcountsDF[,seq(2, ncol(GCcountsDF), 
  2)]
colnames(Counts) <- gsub("\\..*$", "", 
  basename(GCfiles))

rownames(Counts) <- gsub("\\..*$", "", rownames(Counts))
Counts$gene_id <- rownames(Counts)
Counts <- aggregate(.~gene_id, Counts, mean)
rownames(Counts) <- Counts$gene_id
Counts <- subset(Counts, select=-gene_id)


if (file.exists(paste0(RobjectDir, "/GCcounts.RData"))) {
  Counts <- readRDS(paste0(RobjectDir, "/GCcounts.RData"))
} else {
  saveRDS(Counts, file=paste0(RobjectDir, "/GCcounts.RData"))
}


### 3. Perform pre-normalisation PCA and RLE plots ###

# eliminate lowly expressed genes (rows where there are less than 3 counts where df > 4):
print(paste0("No. rows before filtering is: ", nrow(Counts)))
Counts <- Counts %>%
  rownames_to_column('gene_id') %>%
  dplyr::filter(rowSums(Counts > 5) >= 3) %>%
  column_to_rownames('gene_id')
print(paste0("No. rows after  filtering: ", nrow(Counts)))

# create pre-normalised PCA plot from counts and plot:
Counts <- apply(Counts, 2, unlist)
if (ncol(Counts) > nrow(Counts)) {
  pca <- prcomp(Counts)
} else {
  pca <- princomp(Counts)	  
}

if (file.exists(paste0(plotDir, "pcaCompsPrenormGC.pdf"))) {
  print(paste0(plotDir, "pcaCompsPrenormGC.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", plotDir, "pcaCompsPrenormGC.pdf"))
  pdf(file = paste0(plotDir, "pcaCompsPrenormGC.pdf"))
  plot(pca)
  dev.off()
}

splt <- unlist(lapply(split(colnames(Counts), gsub("[1-3]", "", colnames(Counts))), length))
for (i in 1:length(splt)) {
  if (i==1) {
    typeF <- c(rep(names(splt)[i], splt[i]))
  } else {
    typeF <- c(typeF, rep(names(splt)[i], splt[i]))
  }
}
levels(typeF) <- STypes
  
# save original colnames:
#origColNames <- colnames(Counts)
# remove the sample numbers from colnames:
#colnames(Counts) <- gsub("[0-9]", "", colnames(Counts))

storage.mode(Counts) <- "integer"

# convert matrix into SeqExpressionSet:
set <- newSeqExpressionSet(Counts, phenoData = data.frame(typeF, row.names=colnames(Counts)))

# create pre-norm RLE plot:
if (file.exists(paste0(plotDir,  "RLEPrenormGC.pdf"))) {
  print(paste0(plotDir,  "RLEPrenormGC.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", plotDir,  "RLEPrenormGC.pdf"))
  par(mar=c(1,1,1,1))
  pdf(file = paste0(plotDir,  "RLEPrenormGC.pdf"))
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
pdf(file = paste0(plotDir,  "_pcalaneNormGC.pdf"), height = 15, width = 20)
plotPCA(nSet, cex=0.7)
dev.off()


### 5. Perform differential expression comparing normalised FT controls to cancer samples ###

genes <- rownames(Counts)
# normalise using upper-quartile normalisation (http://inaykmittal.blogspot.com.au/2013/10/ pkmrpkm-normalization-caveat-and-upper.html)
# design matrix specifying all samples as the thing to be compared to:
design <- model.matrix(~0+typeF, data=pData(nSet))
# convert set into a DGElist format which specifies library size and normalisation factors:
y <- DGEList(counts=counts(nSet), group=typeF)

# create MDS plot:
pdf(file = paste0(plotDir,  "_mdslaneNormGC.pdf"), height = 15, width = 20)
plotMDS(y)
dev.off()


# calculate the deviance residuals from a first-pass GLM regression of the counts on the co-variates  of interest (p8 RUVseq manual):
# estimate dispersion:
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
# adjust values using dispersion:
fit <- glmFit(y, design)

saveRDS(fit, file=paste0(RobjectDir, "/",  "DEfitGC.rds"))

con <- c(-1, 1)
comp <- "FT_vs_prPT"

# perform likelihood ratio test:
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

saveRDS(lrt, paste0(RobjectDir, "FT_vs_prPT_lrtGC.rds"))


### 5. Calculate differential expression values ###

# PCA and RLE plots of RUVseq RUVr-normalised data looked best with etween lane normalisation =   'full', will go with this #

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

