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
resultsDir <- paste0(projectDir, "/RNA-seq/results")


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

Genes <- readRDS(file = paste0(RobjectDir, "/", annot, "_RepeatGenes.rds"))


### 3. Load in bam file and convert to GenomicRanges object if needed:

if (file.exists(paste0(RobjectDir, uID, "_bamGR.rds"))) {
  print(paste0("Loading ", RobjectDir, uID, "_bamGR.rds"))
  bamGR <- readRDS(file = paste0(RobjectDir, uID, "_bamGR.rds"))
} else {
  # index bam if needed:
  if (file.exists(paste0(bamDir, "/", uID, "/Aligned.sortedByCoord.out.bam.bai"))) {
    print(paste0("No need to index, ", inFile, " has already been indexed"))
  } else {
    print(paste0("Indexing ", inFile))
    indexBam(inFile)
  }
  
  # assign column names for dataframe:
  what <- c("qname","rname","strand","pos","qwidth")
  # flag unmapped sequences to leave out:
  flag <- scanBamFlag(isUnmappedQuery=FALSE)
  # define parameters of bam scan:
  param <- ScanBamParam(what=what, flag=flag)
  
  # load in bam:
  if (file.exists(paste0(RobjectDir, "/", uID, "_bamFile.rds"))) {
    print(paste0("Loading ", uID, "_bamFile.rds"))
    bam <- readRDS(paste0(RobjectDir, "/", uID, "_bamFile.rds"))
  } else {
    print(paste0("Loading ", inFile))
    bam <- scanBam(inFile ,param=param)
    saveRDS(bam, paste0(RobjectDir, "/", uID, "_bamFile.rds"))
  }
  
  # create data frame of human chromosome lengths:
  seq_lengths <- seqlengths(Hsapiens)
  
  bamGR <- GRanges(
    seqnames = bam[[1]]$rname,
    ranges = IRanges(start = bam[[1]]$pos, width = bam[[1]]$qwidth),
    strand = bam[[1]]$strand,
    seqlengths = seq_lengths,
    qnames = bam[[1]]$qname)
  saveRDS(bamGR, file = paste0(RobjectDir, uID, "_bamGR.rds"))
}


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

# for each list of cGenes, count each element of the list with each element of bamGRs and return as df:
for (i in 1:length(Genes)) {
  if (i==1) {
    Counts <- list(as.data.frame(unlist(lapply(Genes[[i]], function(x) {
      return(count_it(x))
    }))))
    names(Counts)[i] <- names(Genes)[i]
    print(i)
  } else {
    Counts[[i]] <- as.data.frame(unlist(lapply(Genes[[i]], function(x) {
      return(count_it(x))
    })))
    names(Counts)[i] <- names(Genes)[i]
  }
  print(i)
}

# name column of each df:
Counts <- lapply(Counts, function(x) {
  colnames(x) <- "Counts"
  return(x)
})

# merge all dataframes with rows < 3 into 'other' df:
j=1
for (i in 1:length(Counts)) {
  if (nrow(Counts[[i]]) < 3) {
    if (j==1) {
      merged <- Counts[[i]]
      rmInd <- c(i)
    } else {
      merged <- rbind(merged, Counts[[i]])
      rmInd[j] <- i
    }
    j=j+1  
  }
}

# check if any dfs did have rows < 3:
if (exists("merged")) {
  # subset Counts to include only dfs with rows > 3:
  `%notin%` <- function(x,y) !(x %in% y) 
  ind <- seq(1:length(Counts))[seq(1:length(Counts)) %notin% rmInd]
  
  # add merged dfs form above to Counts:
  Counts <- c(Counts[ind], list(merged))
  names(Counts) <- c(names(Counts)[1:(length(Counts)-1)], "other")
}

# split dfs with rows > 10:
j=1
n=1
for (i in 1:length(Counts)) {
  
  if (nrow(Counts[[i]]) > 10) {
    if (j==1) {
      spl <- split(Counts[[i]], c( rep(1, round(nrow(Counts[[i]])/2)), rep(2, (nrow(Counts[[i]])-round(nrow(Counts[[i]])/2))) ))
      names(spl)[1:(n+1)] <- c(paste0(names(Counts)[i], 1), paste0(names(Counts)[i], 2))
      rmInd <- c(j)
    } else {
      spl <- c(spl, split(Counts[[i]], c( rep(1, round(nrow(Counts[[i]])/2)), rep(2, (nrow(Counts[[i]])-round(nrow(Counts[[i]])/2))) )))
      names(spl)[n:(n+1)] <- c(paste0(names(Counts)[i], 1), paste0(names(Counts)[i], 2))
      rmInd[j] <- i
    }
    j=j+1
    n=n+2
  }
}

# check if any dfs did have rows > 10
if (exists("spl")) {
  # subset Counts to include only dfs with rows < 10:
  `%notin%` <- function(x,y) !(x %in% y) 
  ind <- seq(1:length(Counts))[seq(1:length(Counts)) %notin% rmInd]
  
  # add merged dfs form above to Counts:
  Counts <- c(Counts[ind], spl)
}


# create output directory:
outDir <- paste0(RobjectDir, "/", annot, "_RepeatCounts/")
system(paste0("mkdir -p ", outDir))

# save the counts as RDS files:
if (file.exists(paste0(outDir, "/", uID, "_", annot, "RepeatCounts.rds already exists, no need to create"))) {
  print(paste0(outDir, "/", uID, "_", annot, "RepeatCounts.rds"))
} else {
  print(paste0("Creating", outDir, "/", uID, "_", annot, "RepeatCounts.rds"))
  saveRDS(Counts, file = paste0(outDir, "/", uID, "_", annot, "RepeatCounts.rds"))
}


### 6. Count overlaps of bams with gencode ###

# load in gencode as GRanges object:
if (!file.exists(paste0(outDir, "/", uID, "_GCcountsDF.rds"))) {
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
  
  saveRDS(GCcountsDF, file = paste0(outDir, "/", uID, "_GCcountsDF.rds"))
}
