### CPMRheatmap.R ###

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
library(pheatmap)

# define starting variables:
project <- "hgsoc_repeats"
expName <- "exp9"

Type <- "custom3"

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/", expName, "/hrd_ccne_unknown/")
plotDir <- paste0(resultsDir, "/R/", expName, "/plots/DEplots/hrd_ccne_unknown/")

system(paste0("mkdir -p ", plotDir))


### 1. Load in all counts and subset for genes of interest ###

Counts <- readRDS(file=paste0(RobjectDir, "/", Type, "_counts.RData"))

sTypes <- unique(
  grep(
    "id", gsub(
      "^.*\\_", "", colnames(Counts)
    ), value=T, invert = T
  )
)

sTypes <- sTypes[order(sTypes)]

# define repeat species of heatmap as those with FDR<0.05 and FC>1 in any sample comparison with controls
sigReps <- readRDS(file=paste0(RobjectDir, "/", Type, "_DEsigReps.rds"))
#sigReps <- c("(CATTC)n", "Helitron1Na_Mam", "Helitron1Nb_Mam", "L1HS", "L1MD3", "L1P2", "L1P4d", "L1PA2", 
#             "AluYi6", "GSATX", "GSAT","HSATII", "HSAT4", "HSAT5", "ACRO1", "FAM", "REP522", "CER", "SATR2")
sigRepNames <- rownames(sigReps[[1]])
Counts <- Counts[sigRepNames,]

Counts <- Counts[order(gsub("^.*_", "", colnames(Counts)))]


### 2. Calculate CPMRs ###

# calculate total repeat count size:
rSizes <- apply(Counts, 2, sum)

# calculate CPMRs
CPMR <- as.data.frame(t(t(Counts)/rSizes)*1000000)

# log CPMRs:
logCPMR = log10(CPMR+1)


### 3. Create logFC and FDR heatmaps:
p <- pheatmap(logCPMR, fontsize = 7, cluster_cols = F)

pdf(paste0(plotDir, "/repeatLogCPMRhmap.pdf"))
p
dev.off()

# fetch clusters:
logCPMRclust <- cbind(logCPMR, 
                   cluster = cutree(p$tree_row, 
                                    k = 4))
for (i in 1:4) {
  if (i==1) {
    clusters <- list(logCPMRclust[logCPMRclust$cluster==i,])
    df <- subset(clusters[[i]], select=-cluster)
    pdf(paste0(plotDir, "/repeatLogCPMRhmap", as.character(i), ".pdf"))
    pheatmap(df, fontsize = 7, cluster_cols = F)
    dev.off()
  } else {
    clusters[[i]] <- logCPMRclust[logCPMRclust$cluster==i,]
    df <- subset(clusters[[i]], select=-cluster)
    pdf(paste0(plotDir, "/repeatLogCPMRhmap", as.character(i), ".pdf"))
    pheatmap(df, fontsize = 7, cluster_cols = F)
    dev.off()
  }
}

######


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



