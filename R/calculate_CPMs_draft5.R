### calcuate_CPMs.R ###
# This script takes a data frame of counts and a list of library sizes, and calculates CPMs
# from these

### 0. Set up variables and directories ###

# load packages needed:
library(ggplot2)
library(reshape2)

# define starting variables:
project <- "hgsoc_repeats"
methodName <- "method1"

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/")
inDir <- paste0(resultsDir, "/R/", methodName)


### 1. Load in inputs ###

# load in counts:
counts <- read.csv(file = paste0(inDir, "/counts.csv"))
# remove first column of counts:
counts <- counts[,-1]
# rename repeats column of counts:
colnames(counts)[1] <- "repeat_id"

# load in library sizes:
lSizes <- unlist(readRDS(file = paste0(RobjectDir, "/libSizes.rds")))

save.image(file = paste0("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/results/R/", methodName, "_CPM_inputs.RData"))
methodName = "method1"
#load(file = paste0("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/results/R/", methodName, "_CPM_inputs.RData"))


### 2. Calculate CPMs ###

CPMs <- as.data.frame(t(apply(counts[,-1], 1, "/", lSizes)) * 1000000)
CPMs$repeat_id <- counts$repeat_id


### 3. Create line plots of CPMs of each repeat type for each sample type ###

# create repeat_type column with NAs:
CPMs$repeat_type <- rep("NA", nrow(CPMs))

# annotate L1s in repeat_type column:
CPMs$repeat_type[grep("L1", CPMs$repeat_id)] <- "L1"
CPMs$repeat_type[grep("Alu", CPMs$repeat_id)] <- "Alu"
CPMs$repeat_type[grep("centromere", CPMs$repeat_id)] <- "centromere"
CPMs$repeat_type[grep("TTAGGG", CPMs$repeat_id)] <- "telomere"
CPMs$repeat_type[grep(")n", CPMs$repeat_id)] <- "simple"
CPMs$repeat_type[grep("NA", CPMs$repeat_type)] <- "other"

# order CPMs by repeat type:
CPMs <- CPMs[with(CPMs, order(repeat_type)),]

subset_repeats <- function(x, type) {
  return(x[grep(type, x$repeat_id), ])
}

L1s <- subset_repeats(CPMs, "L1")
Alus <- subset_repeats(CPMs, "Alu")
cents <- subset_repeats(CPMs, "centromere")
telos <- subset_repeats(CPMs, "telomere")
simples <- subset_repeats(CPMs, ")n")

oInd <- grep("L1", CPMs$repeat_id, invert = T)
oInd <- append(oInd, grep("Alu", CPMs$repeat_id, invert = T))
oInd <- append(oInd, grep("centromere", CPMs$repeat_id, invert = T))
oInd <- append(oInd, grep("telomere", CPMs$repeat_id, invert = T))
oInd <- append(oInd, grep(")n", CPMs$repeat_id, invert = T))
others <- CPMs[oInd, ]

# sum CPMs of each repeat type and put into data frame:
typeCPMs <- apply(L1s[ ,-c(ncol(L1s)-1, ncol(L1s))], 2, sum)
typeCPMs <- rbind(typeCPMs, apply(Alus[ ,-c(ncol(Alus)-1, ncol(Alus))], 2, sum))
typeCPMs <- rbind(typeCPMs, apply(cents[ ,-c(ncol(cents)-1, ncol(cents))], 2, sum))
typeCPMs <- rbind(typeCPMs, apply(telos[ ,-c(ncol(telos)-1, ncol(telos))], 2, sum))
typeCPMs <- rbind(typeCPMs, apply(simples[ ,-c(ncol(simples)-1, ncol(simples))], 2, sum))
typeCPMs <- as.data.frame(rbind(typeCPMs, apply(others[ ,-c(ncol(others)-1, ncol(others))], 2, sum)), row.names = F)
typeCPMs$repeat_type <- c("L1", "Alu", "centromere", "telomere", "simple", "other")
colnames(typeCPMs) <- gsub(".counts", "", colnames(typeCPMs))

# replace zeros with tiny number so they are loggable:
typeCPMs[typeCPMs==0] <- 0.0000001

###### order columns - will need to automate for diff sample types at some point #####
gtxInd <- grep("gtx", colnames(typeCPMs))
ftInd <- grep("FT", colnames(typeCPMs))
endoInd <- grep("endo[1-3]", colnames(typeCPMs))
endosisInd <- grep("endosis", colnames(typeCPMs))
kurInd <- grep("kur", colnames(typeCPMs))
prInd <- grep("PR", colnames(typeCPMs))
allInd <- c(gtxInd, ftInd, endoInd, endosisInd, kurInd, prInd)

typeCPMs <- typeCPMs[ , allInd]
colnames(typeCPMs) <- gsub("[0-9]", "", colnames(typeCPMs))

save.image(file = paste0("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/results/R/", methodName, "_CPM_temp.RData"))
methodName = "method1"
#load(file = paste0("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/results/R/", methodName, "_CPM_temp.RData"))

# melt data for ggplot:
plotCPMs <- melt(typeCPMs, id.vars = "repeat_type", value.name = "CPM", variable.name = "sample_type")

#p <- ggplot(plotCPMs, aes(x=sample_type, y=CPM, group=repeat_type, colour = repeat_type))
p <- ggplot(plotCPMs[grep("telomere", plotCPMs$repeat_type, invert = T),], aes(x=sample_type, y=CPM, group=repeat_type, colour = repeat_type))
p <- p + geom_line()
p <- p + scale_y_log10()
if (file.exists(paste0(outDir, "/bowtellFT1_repeat_overlaps.pdf"))) {
  print(paste0(outDir, "/bowtellFT1_repeat_overlaps.pdf already exists. No need to create."))
} else {
  print(paste0("Creating ", outDir, "/bowtellFT1_repeat_overlaps.pdf"))
  pdf(file = paste0(outDir, "/bowtellFT1_repeat_overlaps.pdf"))
  p
  dev.off()
}

save.image(file = paste0("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/results/R/", methodName, "_CPM_triplicates.RData"))
methodName = "method1"
#load(file = paste0("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/results/R/", methodName, "_CPM_triplicates.RData"))


