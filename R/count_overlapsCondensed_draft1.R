# This R script counts overlaps of multiple RNA-seq data sets with a GTF file of transcripts
# and outputs a counts .csv file for each sample

### 0. Set up variables and directories ###

# load packages needed:
library(GenomicRanges)
library(ShortRead)
library(rtracklayer)
library("BSgenome.Hsapiens.UCSC.hg38")
library(ggplot2)

# define starting variables:
project <- "hgsoc_repeats"
methodName <- "method1"
sampleName <- ""

# define directories:
#homeDir <- "/Users/jamestorpy/clusterHome/"
homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project)
refDir <- paste0(projectDir, "/RNA-seq/refs/")
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/")


### 1. Load in repeats annotations ###

scGenes <- readRDS(file = paste0(RobjectDir, "/sc_RepeatGenes.rds"))
cGenes <- readRDS(file = paste0(RobjectDir, "/c_RepeatGenes.rds"))
tGenes <- readRDS(file = paste0(RobjectDir, "/t_RepeatGenes.rds"))


### 2. Load in bam files ###

# load in the bamFiles as list of GRanges objects:
# fetch vector of bam files:
inFiles <- list.files(RobjectDir, recursive = T, pattern = "bamFile.rds", full.names = T)

# define the sampleNames:
sampleNames <- strsplit(inFiles, "/")
sampleNames <- lapply(sampleNames, function(x) {
  return(gsub("_bamFile.rds", "", x[11]))
})
sampleNames <- unlist(sampleNames)

print("The samples to be processed are: ")
print(sampleNames)

# assign column names for dataframe:
what <- c("qname","rname","strand","pos","qwidth")
# flag unmapped sequences to leave out:
flag <- scanBamFlag(isUnmappedQuery=FALSE)
# define parameters of bam scan:
param <- ScanBamParam(what=what, flag=flag)

# create empty list for bamFiles:
bamFiles <- list()
i=1
for (s in sampleNames) {
  # save complete filename of bam RDS for loading:
  file <- paste0(RobjectDir, "/", s, "_bamFile.rds")
  
  # load in bams:
  print(paste0("Loading ", file))
  bamFiles[[i]] <- readRDS(file = paste0(file))
  
  i <<- i+1
}

# create data frame of human chromosome lengths:
seq_lengths <- seqlengths(Hsapiens)

# convert the bams to GRanges objects:
bamGRs <- lapply(bamFiles, function(x) {
  bamGR <- GRanges(
    seqnames = x[[1]]$rname,
    ranges = IRanges(start = x[[1]]$pos, width = x[[1]]$qwidth),
    strand = x[[1]]$strand,
    seqlengths = seq_lengths,
    qnames = x[[1]]$qname)
  return(bamGR)
})
names(bamGRs) <- sampleNames


### 3. Fetch library sizes ###

# determine library size of each data set:
if (file.exists(paste0(RobjectDir, "/lib", sampleName, "Sizes.rds"))) {
  lSizes <- readRDS(file = paste0(RobjectDir, "/lib", sampleName, "Sizes.rds"))
} else {
  libSizes <- lapply(bamGRs, function(x) {
    return(length(seqnames(x)))
  })
  # save the libSizes as RDS file:
  saveRDS(libSizes, file = paste0(RobjectDir, "/lib", sampleName, "Sizes.rds"))
}

### 4. Count overlaps of bams with all repeat annotation GRanges objects:

# count overlaps of bams with annotation ranges:
count_it <- function(x, annot = scGenes[[1]]) {
  writeLines("\n")
  print("Counting overlaps...")
  counts <- as.data.frame(countOverlaps(annot, x))
  return(sum(counts))
}

# for each element of scGenes, count with each element of bamGRs
scCounts <- lapply(scGenes, function(y) {
  return(lapply(bamGRs, count_it, annot = y))
})
scCountsDF <- as.data.frame(do.call("rbind", scCounts))

# define names of bamGRs to add as colnames for following dfs:
Names <- names(bamGRs)

# for each list of cGenes, count each element of the list with each element of bamGRs:
cCounts <- list()
i=1
for (e in cGenes) {
  cCounts[[i]] <- lapply(e, function(x) {
    result <- lapply(bamGRs, count_it, annot = x)
    return(do.call("cbind", result))
  })
  writeLines("\n")
  print(names(cGenes)[i])
  i=i+1
}

names(cCounts) <- names(cGenes)
  
# concatenate list into df:
cCountsTemp <- lapply(cCounts, function(x) {
  result <- do.call("rbind", x)
  rownames(result) <- names(x)
  return(result)
})

cCountsDF <- do.call("rbind", cCountsTemp)

# for each list of tGenes, count each element of the list with each element of bamGRs:
tCounts <- list()
i=1
for (e in tGenes) {
  tCounts[[i]] <- lapply(e, function(x) {
    result <- lapply(bamGRs, count_it, annot = x)
    return(do.call("cbind", result))
  })
  writeLines("\n")
  print(names(tGenes)[i])
  i=i+1
}

names(tCounts) <- names(tGenes)

# concatenate list into df:
tCountsTemp <- lapply(tCounts, function(x) {
  result <- do.call("rbind", x)
  rownames(result) <- names(x)
  return(result)
})

tCountsDF <- do.call("rbind", tCountsTemp)

# save the counts as RDS files:
if (file.exists(paste0(RobjectDir, "/sc_", sampleName, "RepeatCountsDF.rds"))) {
  print(paste0(RobjectDir, "/sc_", sampleName, "RepeatCountsDF.rds already exists, no need to create"))
} else {
  print(paste0("Creating", RobjectDir, "/sc_", sampleName, "RepeatCountsDF.rds"))
  saveRDS(scCountsDF, file = paste0(RobjectDir, "/sc_", sampleName, "RepeatCountsDF.rds"))
}

if (file.exists(paste0(RobjectDir, "/c_", sampleName, "RepeatCountsDF.rds already exists, no need to create"))) {
  print(paste0(RobjectDir, "/c_", sampleName, "RepeatCountsDF.rds"))
} else {
  print(paste0("Creating", RobjectDir, "/c_", sampleName, "RepeatCountsDF.rds"))
  saveRDS(cCountsDF, file = paste0(RobjectDir, "/c_", sampleName, "RepeatCountsDF.rds"))
}

if (file.exists(paste0(RobjectDir, "/t_", sampleName, "RepeatCountsDF.rds already exists, no need to create"))) {
  print(paste0(RobjectDir, "/t_", sampleName, "RepeatCountsDF.rds"))
} else {
  print(paste0("Creating", RobjectDir, "/t_", sampleName, "RepeatCountsDF.rds"))
  saveRDS(tCountsDF, file = paste0(RobjectDir, "/t_", sampleName, "RepeatCountsDF.rds"))
}


### 6. Count overlaps of bams with gencode ###

# load in gencode as GRanges object:
if (!file.exists(paste0(RobjectDir, "GCcountsDF.rds"))) {
  if (file.exists(paste0(RobjectDir, "GCgenes.rds"))) {
    GCgenes <- readRDS(paste0(RobjectDir, "GCgenes.rds"))
  } else {
    GC <- import(paste0(genomeDir, gcName))
    saveRDS(GC, file = paste0(RobjectDir, "GCgenes.rds"))
  }
  
  GCgenes <- GC[grep("\\.|M", seqnames(GC), invert = T),]
  
  GCcounts <- lapply(bamGRs, function(x) {
    result <- countOverlaps(GCgenes, x)
  })
  
  GCtemp <- as.data.frame(do.call("cbind", GCcounts))
  GCtemp$gene_id <- gsub("\\.*", "", GCgenes$gene_id)
  GCcountsDF <- aggregate(.~gene_id, GCtemp, sum)
  
  saveRDS(GCcountsDF, file = paste0(RobjectDir, "GCcountsDF.rds"))
}
