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
library(reshape2)
library(Rmisc)

# define starting variables:
project <- "hgsoc_repeats"
expName <- "exp5"
STypes <- c("arPT5", "rcAF6")
Type <- "c"

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
#homeDir <- "/Users/jamestorpy/Documents/Garvan/phd"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/", expName, "/")
plotDir <- paste0(resultsDir, "/R/", expName, "/plots/DEplots/arPT5_vs_rcAF6/")

system(paste0("mkdir -p ", plotDir))


### 1. Load in all counts ###

cCounts <- readRDS(paste0(RobjectDir, "/", Type, "_allcounts.htseq.rds"))
gcCounts <- readRDS(paste0(RobjectDir, "/gc_allcounts.htseq.rds"))

# append GCcountsDF to each GRanges object of custom1Counts:
Counts <- rbind(cCounts, gcCounts)

# remove additional info rows:
Counts <- Counts[grep("__", Counts$gene_id, invert=T),]

# make rownames gene_id, get rid of latter column and change
# storage mode from factor to integer:
Rnames <- Counts$gene_id
Counts <- sapply(subset(Counts, select=-gene_id), unlist)
storage.mode(Counts) <- "integer"
Counts <- as.data.frame(Counts)
rownames(Counts) <- Rnames

# save Counts:
#saveRDS(Counts, file=paste0(RobjectDir, "/", Type, "_counts.RData"))
saveRDS(Counts, file=paste0(RobjectDir, "/", Type, "_GCcounts.RData"))


### 3. Perform pre-normalisation PCA and RLE plots ###

# eliminate lowly expressed genes (rows where there are less than 3 counts where df > 4):
print(paste0("No. rows before filtering is: ", nrow(Counts)))
Counts <- Counts %>%
  rownames_to_column('gene_id') %>%
  dplyr::filter(rowSums(Counts > 5) >= 1) %>%
  column_to_rownames('gene_id')
print(paste0("No. rows after  filtering: ", nrow(Counts)))

# create pre-normalised PCA plot from counts and plot:
Counts <- apply(Counts, 2, unlist)
if (ncol(Counts) > nrow(Counts)) {
  pca <- prcomp(Counts)
} else {
  pca <- princomp(Counts)	  
}

if (file.exists(paste0(plotDir, "htseq_pcaCompsPrenormGC.pdf"))) {
  print(paste0(plotDir, "htseq_pcaCompsPrenormGC.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", plotDir, "htseq_pcaCompsPrenormGC.pdf"))
  pdf(file = paste0(plotDir, "GChtseq_pcaCompsPrenormGC.pdf"))
  #pdf(file = paste0(plotDir, "htseq_pcaCompsPrenormGC.pdf"))
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
levels(typeF) <- c("arPT5", "rcAF6")

# convert matrix into SeqExpressionSet:
set <- newSeqExpressionSet(Counts, phenoData = data.frame(typeF, row.names=colnames(Counts)))

# create pre-norm RLE plot:
if (file.exists(paste0(plotDir, "htseq_RLEPrenormGC.pdf"))) {
  print(paste0(plotDir, "htseq_RLEPrenormGC.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", plotDir, "htseq_RLEPrenormGC.pdf"))
  par(mar=c(1,1,1,1))
  pdf(file = paste0(plotDir, "GChtseq_RLEPrenormGC.pdf"))
  #pdf(file = paste0(plotDir, "htseq_RLEPrenormGC.pdf"))
  plotRLE(set)
  dev.off()
}

# create RUVseq pre-norm PCA:
if (file.exists(paste0(plotDir, "htseq_pcaPrenormGC.pdf"))) {
  print(paste0(plotDir, "htseq_pcaPrenormGC.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", plotDir, "htseq_pcaPrenormGC.pdf"))
  pdf(file = paste0(plotDir, "GChtseq_pcaPrenormGC.pdf"), height = 15, width = 20)
  #pdf(file = paste0(plotDir, "htseq_pcaPrenormGC.pdf"), height = 15, width = 20)
  plotPCA(set, cex=0.7)
  dev.off()
}


### 4. perform normalisation on counts using RUVseq:

# perform between lane full normalisation:
nSet <- betweenLaneNormalization(set, which="full")
pdf(file = paste0(plotDir, "GChtseq_c_RLElaneNormGC.pdf"))
#pdf(file = paste0(plotDir, "htseq_c_RLElaneNormGC.pdf"))
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

pdf(file = paste0(plotDir, "GChtseq_c_pcalaneNormGC.pdf"), height = 15, width = 20)
#pdf(file = paste0(plotDir, "htseq_c_pcalaneNormGC.pdf"), height = 15, width = 20)
plotPCA(nSet, cex=0.7)
dev.off()




### 5. Perform differential expression comparing normalised FT controls to cancer samples ###

genes <- rownames(Counts)
# normalise using upper-quartile normalisation (http://inaykmittal.blogspot.com.au/2013/10/	pkmrpkm-normalization-caveat-and-upper.html)
# design matrix specifying all samples as the thing to be compared to:
design <- model.matrix(~0+typeF, data=pData(nSet))
# convert set into a DGElist format which specifies library size and normalisation factors:
y <- DGEList(counts=counts(nSet), group=typeF)

# calculate CPMs:
CPM <- as.data.frame(t(t(Counts)/lSizes)*1000000)

# calculate total repeat count size:
rSizes <- apply(Counts, 2, sum)
# calculate CPMRs
CPMR <- as.data.frame(t(t(Counts)/rSizes)*1000000)

# prepare each df for plotting:
df <- CPMR
sampleNames <- colnames(df)
df <- sapply(df, as.numeric)
rownames(df) <- rownames(CPM)
pDF <- melt(df, varnames = c("repeat_id", "sample"), value.name = "CPMR")

# sort scpCounts data frames according to repeat_id:
sort_rID <- function(x) {
  if (colnames(x)[1] == "repeat_id") {
    return(x[with(x, order(repeat_id)),])
  } else {
    return(x)
  }
}
pDF <- sort_rID(pDF)

# save df as list element:
#if (typeNames[n]=="CPM") {
#  if (j==1) {
#    CPMcount <- list(pDF)
#  } else {
#    CPMcount[[j]] <- pDF
#  }
#} else {
#  if (j==1) {
#    CPMRcount <- list(pDF)
#  } else {
#    CPMRcount[[j]] <- pDF
#  }
#}

# remove duplicate numbers:
#pDF$sample <- gsub("AOCS_[0-9][0-9][0-9]_", "", pDF$sample)

# order levels of sample factor, putting controls first:
orderS <- function(x) {
  #x$sample <- factor(x$sample, levels=c("FT", "erPT", "mrPT", "arPT", "prPT", "rfPT", "lrcT", "msST", "pAF", "rcAF"))
  x$sample <- factor(x$sample, levels=c("arPT", "rcAF"))
  return(x)
}
pDF <- orderS(pDF)

# create stats summary of the data:
statDF <- summarySE(pDF, measurevar = "CPMR", groupvar = c("repeat_id", "sample"))

# add 1 to all values so they are loggable:
add2zero <- function(x) {
  if (ncol(x) == 3) {
    x[,3] <- x[,3]+1
    return(x)
  } else {
    x[,2] <- x[,2]+1
    return(x)
  }
}
pDF4Log <- add2zero(pDF)

# create stats summary of the data:
statDF4Log <- summarySE(pDF, measurevar = "CPMR", groupvar = c("repeat_id", "sample"))

plot_itdf <- function(x, logNo="nope") {
  if (ncol(x) == 7) {
    p <- ggplot(x, aes(x=sample, y=x[,3], group=repeat_id, colour=repeat_id))
    p <- p + geom_errorbar(aes(ymin=eval(parse(text="CPMR"))-se, ymax=eval(parse(text="CPMR"))+se), width=0.1)
    p <- p + geom_line()
    p <- p + theme(axis.text.x = element_text(angle = 90), legend.position = "bottom")
    p <- p + ylab("CPMR")
    if (logNo=="log10") {
      p <- p + scale_y_log10()
      return(p)
    } else {
      return(p)
    }
  } else {
    p <- ggplot(x, aes(x=sample, y=x[,4], group = 1))
    #p <- p + geom_errorbar(aes(ymin=eval(parse(text="CPMR"))-se-se, ymax=eval(parse(text="CPMR"))-se+se), width=0.1)
    p <- p + geom_line()
    p <- p + theme(axis.text.x = element_text(angle = 90), legend.position = "bottom")
    p <- p + ylab("CPMR")
    if (logNo=="log10") {
      p <- p + scale_y_log10()
      return(p)
    } else {
      return(p)
    }
  }
}

pDFPlot <- plot_itdf(statDF)
pDFPlotLog10 <- plot_itdf(statDF4Log, logNo="log10")

if (file.exists(paste0(plotDir, "/", Type, "_", names(countDFs)[j], "_", "CPMR", ".pdf"))) {
  print(paste0(plotDir, "/", Type, "_", names(countDFs)[j], "_", "CPMR", ".pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", plotDir, "/", Type, "_", names(countDFs)[j], "_", "CPMR", ".pdf"))
  pdf(file = paste0(plotDir, "/", Type, "_", "CPMR", ".pdf"), width = 10, height=10)
  print(pDFPlot)
  dev.off()
}

if (file.exists(paste0(plotDir, "/", Type, "_", names(countDFs)[j], "_", "CPMR", "Log10.pdf"))) {
  print(paste0(plotDir, "/", Type, "_", names(countDFs)[j], "_", "CPMR", "Log10.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", plotDir, "/", Type, "_", names(countDFs)[j], "_", "CPMR", "Log10.pdf"))
  pdf(file = paste0(plotDir, "/", Type, "_", names(countDFs)[j], "_", "CPMR", "Log10.pdf"), width = 10, height=10)
  print(pDFPlotLog10)
  dev.off()
}
n=n+1
}
j=j+1
}

names(CPMRcount) <- names(countDFs)
saveRDS(CPMRcount, file=paste0(RobjectDir, "/", expName, "/", Type, "CPMR.rds"))

names(CPMcount) <- names(countDFs)
saveRDS(CPMcount, file=paste0(RobjectDir, "/", expName, "/", Type, "CPM.rds"))
