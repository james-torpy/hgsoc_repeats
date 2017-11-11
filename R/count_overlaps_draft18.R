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
refName <- "human-89.repeats.tab"
sampleName <- "subset"
gcName <- "gencode_v24_hg38_annotation.gtf"

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project)
bamDir <- paste0(projectDir, "/RNA-seq/results/star/", methodName)
refDir <- paste0(projectDir, "/RNA-seq/refs/")
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/")
outDir <- paste0(resultsDir, "/R/", methodName, "/")
bedDir <- paste0(outDir, "bed")
genomeDir <- paste0(homeDir, "genomes/hg38_ercc/")

# create outDir:
#system(paste0("mkdir -p ", outDir))
#print(paste0("The outDir is: ", outDir))
#system(paste0("mkdir -p ", bedDir))


### 1. Load in repeats annotation file ###

if (file.exists(paste0(RobjectDir, "/sc_RepeatGenes.rds"))) {
  print("Loading genes annotations")
  scGenes <- readRDS(file = paste0(RobjectDir, "/sc_RepeatGenes.rds"))
  cGenes <- readRDS(file = paste0(RobjectDir, "/c_RepeatGenes.rds"))
  tGenes <- readRDS(file = paste0(RobjectDir, "/t_RepeatGenes.rds"))
} else {
  # define gtf variable:
  ref <- paste0(refDir, refName)
  print(paste0("The ref file is: ", ref))

  # load in the gene coordinates from the tab file:
  genes <- read.table(file = ref, sep = "\t", as.is = T)
  saveRDS(genes, file = paste0(RobjectDir, "genesTemp.rds"))
  #genes <- readRDS(file = paste0(RobjectDir, "genesTemp.rds"))

  # format and convert to GRanges object:
  genes <- genes[-1,]
  # format the annotation df so it only includes chromosomes 1-22, X, Y:
  #genes <- genes[grep("\\.|MT", genes$V1, invert = T),]
  genes$V2 <- as.numeric(genes$V2)
  genes$V3 <- as.numeric(genes$V3)
  genesGR <- GRanges(
    seqnames = genes$V1,
    ranges = IRanges(start=genes$V2, end=genes$V3),
    strand = genes$V4,
    score = genes$V7
  )
  values(genesGR) <- data.frame(genes$V5, genes$V6, genes$V7, genes$V8, genes$V9)
  colnames(values(genesGR)) <- c("type", "analysis", "score", "class", "class_desc")
  
  # save annotation as RDF file:
  saveRDS(genesGR, file=paste0(RobjectDir, "repeatsGR.rds"))
  #genesGR <- readRDS(file=paste0(RobjectDir, "repeatsGR.rds"))
  
  
  ### 2. Set up lists ###
  
  # Each GRangesList corresponds to a plot #
  # Each GRanges object corresponds to a line on a plot #
  
  
  ### 2a. Set up superclass GRangesList ###
  
  # modify class names to group as superclasses:
  genesGR$superclass <- gsub("\\?.*", "", gsub("/.*", "", genesGR$class))
  # Split genesGR by superclasses and reduce to avoid overlaps between superclasses:
  scGenes <- endoapply(split(genesGR, genesGR$superclass), reduce)
  
  
  ### 2b. Set up class list of GRangesLists ###
  
  # convert class and type columns of genesGR from factors to characters
  genesGR$class <- as.character(genesGR$class)
  genesGR$type <- as.character(genesGR$type)
  
  # split genesGR by superclass:
  scTemp <- split(genesGR, genesGR$superclass)
  # make a list of GRanges superclass objects:
  scTempList <- list()
  i=1
  scTempList <- lapply(scTemp, function(x) {
    scTempList[[i]] <- x
    return(scTempList)
    i <<- i+1
  })
  
  # split each GRanges element of scTempList into GRanges lists by class:
  cGenes <- lapply(scTempList, function(x) {
    spl <- split(x[[1]], x[[1]]$class)
    return(reduce(spl))
  })
  
  
  ### 2b. Set up type list of GRangesLists ###
  
  # split genesGR by class:
  cTemp <- split(genesGR, genesGR$class)
  # make a list of GRanges superclass objects:
  cTempList <- list()
  i=1
  cTempList <- lapply(cTemp, function(x) {
    cTempList[[i]] <- x
    return(cTempList)
    i <<- i+1
  })
  
  
  
  # split each GRanges element of cTempList into GRanges lists by type:
  tGenes <- lapply(cTempList, function(x) {
    spl <- split(x[[1]], x[[1]]$type)
    return(reduce(spl))
  })
  
  create_rds <- function(x, rds) {
    if (file.exists(rds)) {
      print(paste0(rds, " exists, no need to create"))
    } else {
      print(paste0("Creating ", rds))
      saveRDS(x, file = rds)
    }
  }

  readRDS(scGenes, file = paste0(RobjectDir, "/sc_RepeatGenes.rds"))
  readRDS(cGenes, file = paste0(RobjectDir, "/c_RepeatGenes.rds"))
  readRDS(tGenes, file = paste0(RobjectDir, "/t_RepeatGenes.rds"))
}

  if (file.exists(paste0(RobjectDir, "/custom1_RepeatGenes.rds"))) {
    readRDS(file = paste0(RobjectDir, "/custom1_RepeatGenes.rds"))
  } else {
    # load in genes of interest from subtype DE and create annotation from
  # this:
  GOI <- readRDS(file = "/share/ScratchGeneral/jamtor//projects/hgsoc_repeats/RNA-seq/Robjects/GOIsubsets.rds")
  # customise to break into interesting categories:
  GOI <- list(GOI[[2]], GOI[[3]])
  GOI[[1]][1] <- "RC\\?/Helitron\\?"
  GOI[[1]][4] <- "RC/Helitron"
  GOI[[2]] <- c("L1MD", "MER", "MST", "PABL")
  
  # find these classes and categories and append into GRangesList:
  classInd <- list()
  i=1
  for (class in GOI[[1]]) {
    print(paste0("class is ", class))
    classInd[[i]] <- lapply(cGenes, function(x) {
      return(grep(class, names(x)))
    })
  i <<- i+1
  }
  
  typeInd <- list()
  i=1
  for (type in GOI[[2]]) {
    print(paste0("type is ", type))
    typeInd[[i]] <- lapply(tGenes, function(x) {
      return(grep(type, names(x)))
    })
  i <<- i+1
  }
  
  
  custom1 <- c(GRangesList(cGenes[[5]][[2]]),
    GRangesList(cGenes[[3]][[7]]),
    GRangesList(cGenes[[1]][[21]]),
    window(tGenes[[27]], start=61, end=66),
    window(tGenes[[22]], start=1, end=17),
    window(tGenes[[39]], start=61, end=75),
    window(tGenes[[35]], start=290, end=293)
  )
  
  names(custom1)[1:3] <- c("RC?/Helitron?", "LTR/ERVL?", "DNA/hAT?")
  saveRDS(custom1, file = paste0(RobjectDir, "/custom1_RepeatGenes.rds"))
}


### 3. Load in bam files ###

if (file.exists) {

  bamGRs <- readRDS(paste0(RobjectDir, "/bam",  sampleName, "GRs.rds"))

} else {

    # load in the bamFiles as list of GRanges objects:
  # fetch vector of bam files:
  inFiles <- list.files(bamDir, recursive = T, pattern = "Aligned.sortedByCoord.out.bam", full.names = T)
  inFiles <- grep(".bai", inFiles, invert = T, value = T)
  if (sampleName == "subset") {
    inFiles <- grep("subset", inFiles, value = T)
  } else {
    inFiles <- grep("subset", inFiles, invert = T, value = T)
  }
  inFiles <- grep("ribosome", inFiles, invert = T, value = T)
  
  # define the sampleNames:
  sampleNames <- strsplit(inFiles, "/")
  sampleNames <- lapply(sampleNames, function(x) {
    return(x[12])
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
    # index bams if necessary:
    writeLines("\n")
    if (file.exists(paste0(bamDir, "/", s, "/Aligned.sortedByCoord.out.bam.bai"))) {
      print(paste0("No need to index, ", bamDir, "/", s, "/Aligned.sortedByCoord.out.bam has already been indexed"))
    } else {
      print(paste0("bowIndexing ", bamDir, "/", s, "/Aligned.sortedByCoord.out.bam"))
      indexBam(paste0(bamDir, "/", s, "/Aligned.sortedByCoord.out.bam"))
    }
    # save complete filename of bam for loading:
    file <- paste0(bamDir, "/", s, "/Aligned.sortedByCoord.out.bam")
    
    # load in bams:
    print(paste0("Loading ", file, " into bams list."))
    bamFiles[[i]] <- scanBam(file ,param=param)
    
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
  saveRDS(bamGRs, file=paste0(RobjectDir, "/bam",  sampleName, "GRs.rds"))
}



#create_rds(bamGRs, paste0(RobjectDir, "/bamGRs.rds"))


# determine library size of each data set:
#libSizes <- lapply(bamGRs, function(x) {
#  return(length(seqnames(x)))
#3})

# save the libSizes as RDS file:
#saveRDS(libSizes, file = paste0(RobjectDir, "/lib", sampleName, "Sizes.rds"))
lSizes <- readRDS(file = paste0(RobjectDir, "/lib", sampleName, "Sizes.rds"))


### 4. Count overlaps of bams with all repeat annotation GRanges objects:

# count overlaps of bams with annotation ranges:
count_it <- function(x, annot = scGenes[[1]]) {
  writeLines("\n")
  print("Counting overlaps...")
  counts <- as.data.frame(countOverlaps(annot, x))
  return(sum(counts))
}

scCounts <- lapply(scGenes, function(y) {
  return(lapply(bamGRs, count_it, annot = y))
})
scCountsDF <- as.data.frame(do.call("rbind", scCounts))

Names <- names(bamGRs)

rbindList <- function(x) {
  counts <- as.data.frame(do.call("rbind", x))
  colnames(counts) <- Names
  return(counts)
}

cCounts <- lapply(cGenes, function(x) {
  return(lapply(bamGRs, count_it, annot = x))
})
cCountsDF <- rbind_list(lapply(cCounts, function(x) {
  counts <- as.data.frame(do.call("cbind", x))
  colnames(counts) <- Names
  return(counts)
}))

tCounts <- lapply(tGenes, function(x) {
  return(lapply(bamGRs, count_it, annot = x))
})
tCounts <- rbind_list(lapply(tCounts, rbind_list))

custom1Counts <- lapply(custom1, function(x) {
  return(lapply(bamGRs, count_it, annot = x))
})
custom1CountsDF <- rbindList(custom1Counts)


# save the counts and countsDFs as RDS files:
if (file.exists(paste0(RobjectDir, "/sc_", sampleName, "RepeatCounts.rds"))) {
  print(paste0(RobjectDir, "/sc_", sampleName, "RepeatCounts.rds already exists, no need to create"))
} else {
  print(paste0("Creating", RobjectDir, "/sc_", sampleName, "RepeatCounts.rds"))
  saveRDS(scCounts, file = paste0(RobjectDir, "/sc_", sampleName, "RepeatCounts.rds"))
}

if (file.exists(paste0(RobjectDir, "/c_", sampleName, "RepeatCounts.rds already exists, no need to create"))) {
  print(paste0(RobjectDir, "/c_", sampleName, "RepeatCounts.rds"))
} else {
  print(paste0("Creating", RobjectDir, "/c_", sampleName, "RepeatCounts.rds"))
  saveRDS(cCounts, file = paste0(RobjectDir, "/c_", sampleName, "RepeatCounts.rds"))
}

if (file.exists(paste0(RobjectDir, "/t_", sampleName, "RepeatCounts.rds already exists, no need to create"))) {
  print(paste0(RobjectDir, "/t_", sampleName, "RepeatCounts.rds"))
} else {
  print(paste0("Creating", RobjectDir, "/t_", sampleName, "RepeatCounts.rds"))
  saveRDS(tCounts, file = paste0(RobjectDir, "/t_", sampleName, "RepeatCounts.rds"))
}

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
  saveRDS(tCounts, file = paste0(RobjectDir, "/t_", sampleName, "RepeatCountsDF.rds"))
}




### 6. Count overlaps of bams with gencode ###

# load in gencode as GRanges object:

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