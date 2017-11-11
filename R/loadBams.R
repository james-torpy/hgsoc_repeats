### loadBams.R ###

### 0. Set up variables and directories ###

# load packages needed:
library(GenomicRanges)
library(ShortRead)
library("BSgenome.Hsapiens.UCSC.hg38")

# define starting variables:
project <- "hgsoc_repeats"
methodName <- "method1"

# define directories:
#homeDir <- "/Users/jamestorpy/clusterHome/"
homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "projects/", project, "/")
resultsDir <- paste0(projectDir, "RNA-seq/results/")
bamDir <- paste0(resultsDir, "star/", methodName)

RobjectDir <- paste0(projectDir, "RNA-seq/Robjects/")


### 1. Load bam file ###

# load in the command line arguments:
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

if (file.exists(paste0(RobjectDir, "/", uID, "_bamGR.rds"))) {
  print(paste0("Loading ", RobjectDir, "/", uID, "_bamGR.rds"))
  bamGRs <- readRDS(file=paste0(RobjectDir, "/", uID, "_bamGR.rds"))
} else {
  # assign column names for dataframe:
  what <- c("qname","rname","strand","pos","qwidth")
  # flag unmapped sequences to leave out:
  flag <- scanBamFlag(isUnmappedQuery=FALSE)
  # define parameters of bam scan:
  param <- ScanBamParam(what=what, flag=flag)
  
  # index bam if necessary:
  writeLines("\n")
  if (file.exists(paste0(bamDir, "/", uID, "/Aligned.sortedByCoord.out.bam.bai"))) {
    print(paste0("No need to index, ", bamDir, "/", uID, "/Aligned.sortedByCoord.out.bam has already been indexed"))
  } else {
    print(paste0("Indexing ", bamDir, "/", uID, "/Aligned.sortedByCoord.out.bam"))
    indexBam(paste0(bamDir, "/", uID, "/Aligned.sortedByCoord.out.bam"))
  }
  
  # load in bam:
  #if (file.exists(paste0(RobjectDir, "/", uID, "_bamGR.rds"))) {
  #  bamFile <- readRDS(file = paste0(RobjectDir, "/", uID, "_bamGR.rds"))
  #} else {
  print(paste0("Loading ", inFile))
  bamFile <- scanBam(inFile, param=param)
    
  # save the bamFile as RDS file:
  saveRDS(bamFile, file = paste0(RobjectDir, "/", uID, "_bamGR.rds"))
  #}
  
  # create data frame of human chromosome lengths:
  seq_lengths <- seqlengths(Hsapiens)
  
  print(str(bamFile))
  
  # convert the bams to GRanges objects:
  bamGR <- GRanges(
    seqnames = bamFile[[1]]$rname,
    ranges = IRanges(start = bamFile[[1]]$pos, width = bamFile[[1]]$qwidth),
    strand = bamFile[[1]]$strand,
    seqlengths = seq_lengths,
    qnames = bamFile[[1]]$qname)
  # save the bamGRs as RDS file:
  saveRDS(bamGR, file = paste0(RobjectDir, "/", uID, "_bamGR.rds"))
}