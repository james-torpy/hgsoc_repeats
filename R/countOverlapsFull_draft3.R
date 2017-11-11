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
Types <- c("sc", "c", "t")

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


### 2. Load in repeats annotations ###

scGenes <- readRDS(file = paste0(RobjectDir, "/sc_RepeatGenes.rds"))
cGenes <- readRDS(file = paste0(RobjectDir, "/c_RepeatGenes.rds"))
tGenes <- readRDS(file = paste0(RobjectDir, "/t_RepeatGenes.rds"))


### 3. Load in bam file and convert to GenomicRanges object:

# assign column names for dataframe:
what <- c("qname","rname","strand","pos","qwidth")
# flag unmapped sequences to leave out:
flag <- scanBamFlag(isUnmappedQuery=FALSE)
# define parameters of bam scan:
param <- ScanBamParam(what=what, flag=flag)

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


### 4. Fetch library size ###

# determine library size of each data set:
libSize <- length(seqnames(bamGR))


### 5. Count overlaps of bams with all repeat annotation GRanges objects:

# count overlaps of bams with annotation ranges:
count_it <- function(x, bam=bamGR) {
  writeLines("\n")
  print(paste0("Counting overlaps..."))
  counts <- as.data.frame(countOverlaps(x, bam))
  return(sum(counts))
}

# for each element of scGenes, count with each element of bamGRs
scCounts <- lapply(scGenes, function(x) {
  return(count_it(x))
})
scCountsDF <- as.data.frame(do.call("rbind", scCounts))

# for each list of cGenes, count each element of the list with each element of bamGRs:
cCounts <- list()
i=1
for (e in cGenes) {
  cCounts[[i]] <- lapply(e, function(x) { 
    return(count_it(x))
  })
  writeLines("\n")
  print(names(cGenes)[i])
  i=i+1
}

names(cCounts) <- names(cGenes)
  
# concatenate list into df:
cCountsDF <- as.data.frame(unlist(cCounts))

# for each list of tGenes, count each element of the list with each element of bamGRs:
tCounts <- list()
i=1
for (e in tGenes) {
  tCounts[[i]] <- lapply(e, function(x) {
      return(count_it(x))
  })
  writeLines("\n")
  print(names(tGenes)[i])
  i=i+1
}

names(tCounts) <- names(tGenes)

# concatenate list into df:
tCountsDF <- as.data.frame(unlist(tCounts))


# save the counts as RDS files:
if (file.exists(paste0(RobjectDir, "/sc_", uID, "RepeatCountsDF.rds"))) {
  print(paste0(RobjectDir, "/sc_", uID, "RepeatCountsDF.rds already exists, no need to create"))
} else {
  print(paste0("Creating", RobjectDir, "/sc_", uID, "RepeatCountsDF.rds"))
  saveRDS(scCountsDF, file = paste0(RobjectDir, "/sc_", uID, "RepeatCountsDF.rds"))
}

if (file.exists(paste0(RobjectDir, "/c_", uID, "RepeatCountsDF.rds already exists, no need to create"))) {
  print(paste0(RobjectDir, "/c_", uID, "RepeatCountsDF.rds"))
} else {
  print(paste0("Creating", RobjectDir, "/c_", uID, "RepeatCountsDF.rds"))
  saveRDS(cCountsDF, file = paste0(RobjectDir, "/c_", uID, "RepeatCountsDF.rds"))
}

if (file.exists(paste0(RobjectDir, "/t_", uID, "RepeatCountsDF.rds already exists, no need to create"))) {
  print(paste0(RobjectDir, "/t_", uID, "RepeatCountsDF.rds"))
} else {
  print(paste0("Creating", RobjectDir, "/t_", uID, "RepeatCountsDF.rds"))
  saveRDS(tCountsDF, file = paste0(RobjectDir, "/t_", uID, "RepeatCountsDF.rds"))
}
