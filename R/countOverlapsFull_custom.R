### loadBams.R ###

# This R script takes a references to a bam file from call_loadBams.bash,
# loads the bam in and saves as an RDS file

### 0. Set up variables and directories ###

library(GenomicRanges)
library(ShortRead)
library(rtracklayer)
library("BSgenome.Hsapiens.UCSC.hg38")
library(ggplot2)

# define starting variables:
project <- "hgsoc_repeats"
methodName <- "method1"

# define directories:
#homeDir <- "/Users/jamestorpy/clusterHome/"
homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project)
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/")
bamDir <- paste0(projectDir, "/RNA-seq/results/star/", methodName)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
plotDir <- paste0(resultsDir, "/R/plots")


### 1. Load in arguments ###

args = commandArgs(trailingOnly = TRUE)
for (n in 1:2) {
  print(args[n])
}

if (!is.null(args[1])) {
  inFile <- args[1]
}

if (!is.null(args[2])) {
  uID <- args[2]
}

if (!is.null(args[3])) {
  annot <- args[3]
}

### 2. Load in repeats annotations ###

#custom1Genes <- readRDS(file = paste0(RobjectDir, "/custom1_RepeatGenes.rds"))


### 3. Load in bam file and convert to GenomicRanges object:

# index bam if needed:
if (file.exists(paste0(bamDir, "/", uID, "/Aligned.sortedByCoord.out.bam.bai"))) {
  print(paste0("No need to index, ", inFile, " has already been indexed"))
} else {
  print(paste0("Indexing ", inFile))
  indexBam(inFile)
}

# load in bam:
if (file.exists(paste0(RobjectDir, "/", uID, "_bamFile.rds"))) {
	print(paste0("Loading ", uID, "_bamFile.rds"))
	bam <- readRDS(paste0(RobjectDir, "/", uID, "_bamFile.rds"))
} else {
	print(paste0("Loading ", inFile))
 bam <- scanBam(inFile ,param=param)
	saveRDS(bam, paste0(RobjectDir, "/", uID, "_bamFile.rds"))
}

# assign column names for dataframe:
what <- c("qname","rname","strand","pos","qwidth")
# flag unmapped sequences to leave out:
flag <- scanBamFlag(isUnmappedQuery=FALSE)
# define parameters of bam scan:
param <- ScanBamParam(what=what, flag=flag)

# create data frame of human chromosome lengths:
seq_lengths <- seqlengths(Hsapiens)

bamGR <- GRanges(
    seqnames = bam[[1]]$rname,
    ranges = IRanges(start = bam[[1]]$pos, width = bam[[1]]$qwidth),
    strand = bam[[1]]$strand,
    seqlengths = seq_lengths,
    qnames = bam[[1]]$qname)

saveRDS(bamGR, file = paste0(RobjectDir, uID, "_bamGR.rds"))

### 4. Fetch library size ###

# determine library size of each data set:
libSize <- length(seqnames(bamGR))

if (file.exists(paste0(RobjectDir, "/", uID, "_libSize.rds"))) {
  print(paste0(RobjectDir, "/", uID, "_libSize.rds already exists, no need to create"))
} else {
  print(paste0("Creating", RobjectDir, "/", uID, "_libSize.rds"))
  saveRDS(libSize, file = paste0(RobjectDir, "/", uID, "_libSize.rds"))
}


### 5. Count overlaps of bams with all repeat annotation GRanges objects:

# count overlaps of bams with annotation ranges:
count_it <- function(x, bam=bamGR) {
  writeLines("\n")
  print(paste0("Counting overlaps..."))
  counts <- as.data.frame(countOverlaps(x, bam))
  return(sum(counts))
}

# for each list of cGenes, count each element of the list with each element of bamGRs:
# for each element of scGenes, count with each element of bamGRs
custom1Counts <- lapply(custom1Genes, function(x) {
  return(count_it(x))
})
custom1CountsDF <- as.data.frame(do.call("rbind", custom1Counts))

# save the counts as RDS files:
if (file.exists(paste0(RobjectDir, "/custom1_", uID, "RepeatCountsDF.rds already exists, no need to create"))) {
  print(paste0(RobjectDir, "/custom1_", uID, "RepeatCountsDF.rds"))
} else {
  print(paste0("Creating", RobjectDir, "/custom1_", uID, "RepeatCountsDF.rds"))
  saveRDS(custom1CountsDF, file = paste0(RobjectDir, "/custom1_", uID, "RepeatCountsDF.rds"))
}


### 6. Count overlaps of bams with gencode ###

# load in gencode as GRanges object:

if (file.exists(paste0(RobjectDir, "GCgenes.rds"))) {
  GCgenes <- readRDS(paste0(RobjectDir, "GCgenes.rds"))
} else {
  GCgenes <- import(paste0(genomeDir, gcName))
  saveRDS(GCgenes, file = paste0(RobjectDir, "GCgenes.rds"))
}

GCgenes <- GCgenes[grep("\\.|M", seqnames(GCgenes), invert = T),]
GCcounts <- as.data.frame(countOverlaps(GCgenes, bamGR))
GCcounts$gene_id <- gsub("\\.*", "", GCgenes$gene_id)
GCcountsDF <- aggregate(.~gene_id, GCcounts, sum)

saveRDS(GCcountsDF, file = paste0(RobjectDir, uID, "_GCcountsDF.rds"))
