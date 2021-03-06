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
library(plyr)

# define starting variables:
project <- "hgsoc_repeats"
expName <- "exp9"
sTypes <- c("AOCS_091_arPT", "AOCS_091_rcAF")
fTypes <- c("AOCS_091_acquired_resistance_primary_tumour", "AOCS_091_relapse_tumour_ascites_fluid")

######
# test using other samples:
#sTypes <- c("AOCS_172_FT", "AOCS_075_prPT")
#fTypes <- c("AOCS_172_fallopian_tube_control", "AOCS_075_primary_resistant_tumour")
######

Type <- "gc"
ctls <- c("ENSG00000111640", "ENSG00000012048",
	"ENSG00000156414", "ENSG00000187605",
	"ENSG00000077800", "ENSG00000173809")
names(ctls) <- c("GAPDH", "BRCA1", "TDRD9", "TET3",
	"FKBP6", "TDRD12")

# define directories:
#homeDir <- "/Users/jamestorpy/clusterHome/"
homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
	expName, "/")
plotDir <- paste0(resultsDir, "/R/", expName,
	"/plots/conf/arPT5_vs_rcAF6/")

system(paste0("mkdir -p ", plotDir))


### 1. Load in gc counts ###

Counts <- readRDS(paste0(RobjectDir,
	"/gc_allcounts.htseq.rds"))

# remove additional info rows:
Counts <- Counts[grep("__", Counts$gene_id, invert=T),]

# add ensembl ids to rownames:
rownames(Counts) <- Counts$gene_id

# select only relevant samples:
Counts <- subset(Counts, select=c(sTypes))


### 2. Calculate CPMs ###

# calculate library sizes:
lSizes <- apply(Counts, 2, sum)
saveRDS(lSizes, file=paste0(RobjectDir, "/", Type,
	"_libsizes.RData"))

# calculate CPMs:
CPM <- as.data.frame(t(t(Counts)/lSizes)*1000000)

# subset to include only control genes and repeats:
ctlCPM <- CPM[rownames(CPM) %in% ctls,]
ctlCPM$id <- names(ctls)[match(rownames(ctlCPM),
	ctls)]

# prepare ctlCPM for plotting:
ctlCPM <- melt(ctlCPM, varnames = c("id", "sample"),
	value.name = "CPM")
colnames(ctlCPM)[2] <- "sample"

ctlCPM <- ctlCPM[c(match(names(ctls), ctlCPM$id),
	(match(names(ctls), ctlCPM$id))+length(ctls)),]

# relevel factors to put negative ctls first:
ctlCPM$id <- factor(ctlCPM$id)

Types <- gsub("^.*_", "", sTypes)

ctlCPM$sample <- as.character(ctlCPM$sample)

# rename samples:
for (i in 1:length(ctlCPM$sample)) {
  print(i)
  if (ctlCPM$sample[i] == sTypes[1]) {
    ctlCPM$sample[i] <- fTypes[1]
  } else {
    ctlCPM$sample[i] <- fTypes[2]
  }
}

# make sample column a factor column and relevel:
ctlCPM$sample <- factor(ctlCPM$sample, fTypes)
ctlCPM$id <- factor(ctlCPM$id, levels = names(ctls))

# add 1 to all values so they are loggable:
ctlCPM4Log <- data.frame(ctlCPM[,1:2], ctlCPM$CPM+1)
colnames(ctlCPM4Log)[3] <- "CPM"

# create barplot of control gene CPMs:
p <- ggplot(ctlCPM, aes(x=sample, y=CPM))
p <- p + geom_bar(stat = "identity", aes(fill = id), position = "dodge")
p <- p + xlab("samples") + ylab("CPM")
pdf(file = paste0(plotDir, "/ctlCPM_barplot.pdf"), width = 10, height=10)
p
dev.off()

p <- ggplot(ctlCPM4Log, aes(x=sample, y=CPM))
p <- p + geom_bar(stat = "identity", aes(fill = id), position = "dodge")
p <- p + scale_y_log10()
p <- p + xlab("samples") + ylab("log10_CPM")
pdf(file = paste0(plotDir, "/ctlCPM_barplot_log10.pdf"), width = 10, height=10)
p
dev.off()

