### loadBams.R ###

# This R script takes a references to a bam file from call_loadBams.bash,
# loads the bam in and saves as an RDS file

### 0. Load in arguments ###
args = commandArgs(trailingOnly = TRUE)
for (n in 1:2) {
  print(args[n])
}

if (!is.null(args[1])) {
  uID <- args[1]
}

if (!is.null(args[2])) {
  annot <- args[2]
}

print(uID)
print(annot)

### 1. Set up variables and directories ###

library(GenomicRanges)
library(ShortRead)
library(rtracklayer)
library("BSgenome.Hsapiens.UCSC.hg38")
library(ggplot2)

# define starting variables:
project <- "hgsoc_repeats"
gcName <- "gencode_v24_hg38_annotation.gtf"
expName <- "exp2"
#uID <- "bowtell_FT3_subset"
#annot <- "c"

# define directories:
#homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ClusterShare/thingamajigs/jamtor/"
homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project, "/RNA-seq/")
RobjectDir <- paste0(projectDir, "/Robjects/", expName, "/")
bamDir <- paste0(projectDir, "/results/star/GC", uID, "/")
genomeDir <- "/share/ScratchGeneral/jamtor/genomes/"

inFile <- paste0(projectDir, "/results/star/GC/", expName, "/", uID, "/Aligned.sortedByCoord.out.bam")


### 2. Load in repeats annotations ###

Genes <- readRDS(file = paste0(RobjectDir, "/", annot, "_RepeatGenes.rds"))

# split Genes into 2 and use second half:
splitInd <- c( rep(1, round(length(Genes)/2)), rep(2, (length(Genes) - round(length(Genes)/2))) )
Genes <- split(Genes, splitInd)[[2]]
Genes <- lapply(Genes, function(x) {
  strand(x) <- "*"
  return(reduce(x))
})


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

# for each list of cGenes, count each element of the list with each element of bamGRs and return as df:
for (i in 1:length(Genes)) {
  if (i==1) {
    Counts <- list(as.data.frame(unlist(lapply(Genes[[i]], function(x) {
      return(as.data.frame(countOverlaps(x, bamGR)))
    }))))
    names(Counts)[i] <- names(Genes)[i]

    tCounts <- list(as.data.frame(unlist(lapply(Genes[[i]], function(x) {
      return(sum(as.data.frame(countOverlaps(x, bamGR))))
    }))))    
    names(tCounts)[i] <- names(Genes)[i]

    print(i)
  } else {
    Counts[[i]] <- as.data.frame(unlist(lapply(Genes[[i]], function(x) {
      return(as.data.frame(countOverlaps(x, bamGR)))
    })))
    names(Counts)[i] <- names(Genes)[i]

    tCounts[[i]] <- as.data.frame(unlist(lapply(Genes[[i]], function(x) {
      return(sum(as.data.frame(countOverlaps(x, bamGR))))
    })))
    names(tCounts)[i] <- names(Genes)[i]
    print(i)
  }
}

# name column of each df:
Counts <- lapply(Counts, function(x) {
  colnames(x) <- "Counts"
  return(x)
})

tCounts <- lapply(tCounts, function(x) {
  colnames(x) <- "Counts"
  return(x)
})

# merge all dataframes with rows < 3 into 'other' df:
j=1
for (i in 1:length(tCounts)) {
  if (nrow(tCounts[[i]]) < 3) {
    if (j==1) {
      merged <- tCounts[[i]]
      rmInd <- c(i)
    } else {
      merged <- rbind(merged, tCounts[[i]])
      rmInd[j] <- i
    }
    j=j+1  
  }
}

# check if any dfs did have rows < 3:
if (exists("merged")) {
  # subset tCounts to include only dfs with rows > 3:
  `%notin%` <- function(x,y) !(x %in% y) 
  ind <- seq(1:length(tCounts))[seq(1:length(tCounts)) %notin% rmInd]
  
  # add merged dfs form above to tCounts:
  tCounts <- c(tCounts[ind], list(merged))
  names(tCounts) <- c(names(tCounts)[1:(length(tCounts)-1)], "other")
}

# split dfs with rows > 10:
j=1
n=1
for (i in 1:length(tCounts)) {
  if (nrow(tCounts[[i]]) > 10) {
    if (j==1) {
      spl <- split(tCounts[[i]], c( rep(1, round(nrow(tCounts[[i]])/2)), rep(2, (nrow(tCounts[[i]])-round(nrow(tCounts[[i]])/2))) ))
      names(spl)[1:(n+1)] <- c(paste0(names(tCounts)[i], 1), paste0(names(tCounts)[i], 2))
      rmInd <- c(j)
    } else {
      spl <- c(spl, split(tCounts[[i]], c( rep(1, round(nrow(tCounts[[i]])/2)), rep(2, (nrow(tCounts[[i]])-round(nrow(tCounts[[i]])/2))) )))
      names(spl)[n:(n+1)] <- c(paste0(names(tCounts)[i], 1), paste0(names(tCounts)[i], 2))
      rmInd[j] <- i
    }
    j=j+1
    n=n+2
  }
}

# check if any dfs did have rows > 10
if (exists("spl")) {
  # subset tCounts to include only dfs with rows < 10:
  `%notin%` <- function(x,y) !(x %in% y) 
  ind <- seq(1:length(tCounts))[seq(1:length(tCounts)) %notin% rmInd]
  
  # add merged dfs form above to tCounts:
  tCounts <- c(tCounts[ind], spl)
}

# create output directory:
outDir <- paste0(RobjectDir, "/", annot, "_RepeatCounts/")
system(paste0("mkdir -p ", outDir))

# save the counts as RDS files:
if (file.exists(paste0(outDir, "/", uID, "_", annot, "RepeatCountsExt_half2a.rds already exists, no need to create"))) {
  print(paste0(outDir, "/", uID, "_", annot, "RepeatCountsExt_half2a.rds"))
} else {
  print(paste0("Creating ", outDir, "/", uID, "_", annot, "RepeatCountsExt_half2a.rds"))
  saveRDS(Counts, file = paste0(outDir, "/", uID, "_", annot, "RepeatCountsExt_half2a.rds"))
}

if (file.exists(paste0(outDir, "/", uID, "_", annot, "RepeatCounts_half2a.rds already exists, no need to create"))) {
  print(paste0(outDir, "/", uID, "_", annot, "RepeatCounts_half2a.rds"))
} else {
  print(paste0("Creating ", outDir, "/", uID, "_", annot, "RepeatCounts_half2a.rds"))
  saveRDS(tCounts, file = paste0(outDir, "/", uID, "_", annot, "RepeatCounts_half2a.rds"))
}
