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
expName <- "exp7"
#sTypes <- c("AOCS_091_arPT", "AOCS_091_rcAF")

######
# test using other samples:
sTypes <- c("AOCS_172_FT", "AOCS_075_prPT")
fTypes <- factor(c("AOCS_172_fallopian_tube_control", "AOCS_075_primary_resistant_tumour"))
######

Type <- "gc"
ctls <- c("ENSG00000111640", "ENSG00000012048",
	"ENSG00000156414", "ENSG00000187605",
	"ENSG00000077800", "ENSG00000173809",
	"ENSG00000105173")
names(ctls) <- c("GAPDH", "BRCA1", "TDRD9", "TET3",
	"FKBP6", "TDRD12", "CCNE1")

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
	expName, "/")
plotDir <- paste0(resultsDir, "/R/", expName,
	"/plots/arPT5_vs_rcAF6/test/")

system(paste0("mkdir -p ", plotDir))


### 1. Load in all counts ###

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
ctlCPM$id <- rownames(ctlCPM)

# prepare ctlCPM for plotting:
ctlCPM <- melt(ctlCPM, varnames=c("id", "sample"), value.name = "CPM")
colnames(ctlCPM)[2] <- "sample"

ctlCPM <- ctlCPM[c(match(names(ctls), ctlCPM$id),
	(match(names(ctls), ctlCPM$id))+7),]

# relevel factors to put negative ctls first:
ctlCPM$id <- factor(ctlCPM$id)

Types <- gsub("^.*_", "", sTypes)

# rename samples:
for (i in length(ctlCPM$sample)) {
  if (ctlCPM$sample[i] == sTypes[1]) {
    ctlCPM$sample[i] <- as.factor(fTypes[1])
  } else {
    ctlCPM$sample[i] <- fTypes[2]
  }
}


ctlCPM$sample <- gsub(
	Types[1], exp[1], gsub(
    	Types[2], exp[2],
    	ctlCPM$sample
	)
)

#ctlCPM$sample <- factor(ctlCPM$sample,
#	levels=c("AOCS_091_acquired_resistant_primary_tumour",
#		"AOCS_091_recurrent_ascites_fluid"))

ctlCPM$sample <- factor(ctlCPM$sample,
	levels=c("AOCS_172_FT",
		"AOCS_091_recurrent_ascites_fluid"))

# add 1 to all values so they are loggable:
ctlCPM4Log <- data.frame(ctlCPM[,1:2], ctlCPM$CPM+1)
colnames(ctlCPM4Log)[3] <- "CPM"

# create barplot of control gene CPMs:
p <- ggplot(ctlCPM, aes(x=sample, y=CPM))
p <- p + geom_bar(stat = "identity", aes(fill = id), position = "dodge")
p <- p + xlab("samples") + ylab("CPM")
pdf(file = paste0(plotDir, "/ctlCPM_barplot.pdf"), width = 10, height=10)
print(p)
dev.off()

p <- ggplot(ctlCPM4Log, aes(x=sample, y=CPM))
p <- p + geom_bar(stat = "identity", aes(fill = id), position = "dodge")
p <- p + scale_y_log10()
p <- p + xlab("samples") + ylab("log10_CPM")
pdf(file = paste0(plotDir, "/ctlCPM_barplot_log10.pdf"), width = 10, height=10)
print(p)
dev.off()

