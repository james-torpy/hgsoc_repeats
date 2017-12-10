### 5.HGSOC_CPMRs.R ###

# This script takes a list of dfs of different classes and types of repeat counts for
# primary resistant and FT control RNA-seq data sets and calculates CPMRs:


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
expName <- "exp9"

Type <- "custom3"

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- "/Users/jamestorpy/clusterHome//projects/hgsoc_repeats/RNA-seq/Robjects/exp9/"
plotDir <- paste0(resultsDir, "/R/", expName,
                  "/plots/DEplots/")

system(paste0("mkdir -p ", plotDir))


### 1. Load in all counts and subset for genes of interest ###

Counts <- readRDS(paste0(RobjectDir, "/", Type,
                                "_allcounts.htseq.rds"))

# remove pAF sample:
rownames(Counts) <- Counts$gene_id
Counts <- subset(Counts, select=-c(AOCS_093_pAF, gene_id))

STypes <- unique(
  grep(
    "id", gsub(
      "^.*\\_", "", colnames(Counts)
    ), value=T, invert = T
  )
)

#sigReps <- readRDS(paste0(RobjectDir, "/sigReps_pvalue_0.01_FC_1.rds"))
sigReps <- c("(CATTC)n", "HSAT5", "Helitron1Na_Mam", "Helitron1Nb_Mam", "L1P5", "GSATX", "GSAT","HSATII", "ACRO1")

Counts <- Counts[sigReps,]


### 2. Calculate CPMRs ###

# calculate total repeat count size:
rSizes <- apply(Counts, 2, sum)

# calculate CPMRs
CPMR <- as.data.frame(t(t(Counts)/rSizes)*1000000)

# log CPMRs:
logCPMR = log10(CPMR+1)

### 3. Create logFC and FDR heatmaps:

pheatmap(logCPMR, fontsize = 7, cluster_cols = F)

#colnames(fcDF) <- c("acquired_resistance", "extreme_response", "multiple_response", "metastatic", "primary_ascites", "primary_resistant", "refractory_ascites", "refractory")
Names <- gsub(
  "FT_vs_", "", gsub(
    "prPT", "primary_resistant", gsub(
      "rfPT", "primary_refractory", gsub("typeF", "", colnames(fcDF))
    )
  )
)

colnames(fcDF) <- Names
fcDF <- cbind(fcDF[,grep("resistant", colnames(fcDF))], fcDF[,grep("refractory", colnames(fcDF))])

colnames(fdrDF) <- Names
fdrDF <- cbind(fdrDF[,grep("resistant", colnames(fdrDF))], fdrDF[,grep("refractory", colnames(fdrDF))])


######

fake_cat <- c("AOCS_055_primary_resistant", "AOCS_079_primary_resistant", "AOCS_081_primary_resistant")
ind <- colnames(fcDF) %in% fake_cat

# arrange dfs to group annotation together:
fcDF <- cbind(fcDF[,which(ind)], fcDF[,which(!ind)])
fdrDF <- cbind(fdrDF[,which(ind)], fdrDF[,which(!ind)])

# create annotation around fake cat:
Annotation <- data.frame(fakeCat = factor(ind))
rownames(Annotation) <- colnames(fcDF)
Annotation$fakeCat <- as.factor(
  gsub(
    "TRUE", "mutated", gsub(
      "FALSE", "wt", Annotation$fakeCat
    )
  )
)

# change the colors of annotation:
fakeCat <- c("navy", "darkgreen")
names(fakeCat) <- c("wt", "mutated")
anno_cols <- list(fakeCat = fakeCat)

# create heatmap with annotation bar:
pheatmap(fcDF, color = colorRampPalette(c("#08519C", "white", "firebrick3"))(50),
         display_numbers = as.matrix(ifelse(fdrDF < 0.1, "*", "")), fontsize = 8, scale="column",
         cluster_cols = F, annotation = Annotation, annotation_colors = anno_cols)



