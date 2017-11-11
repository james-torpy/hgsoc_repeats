# This R script counts overlaps of multiple RNA-seq data sets with a GTF file of transcripts

# load packages needed:
library(GenomicRanges)
library(ShortRead)
library(rtracklayer)
library("BSgenome.Hsapiens.UCSC.hg38")

# define starting variables:
project <- "hgsoc_repeats"
sampleName <- "L1subset5"
gtfName <- "L1_consensus4.gtf"

# define directories:
#homeDir <- "/Users/jamestorpy/clusterHome"
homeDir <- "/Users/jamestorpy/Documents/Garvan/phd"
projectDir <- paste0(homeDir, "/projects/", project)
inDir <- paste0(projectDir, "/RNA-seq/results/star/", sampleName)
refDir <- paste0(projectDir, "/RNA-seq/refs/")
resultsDir <- paste0(projectDir, "/RNA-seq/results")
outDir <- paste0(resultsDir, "/R/", sampleName)

# create outDir:
system(paste0("mkdir -p ", outDir))
print(paste0("The outDir file is: ", outDir))

# define gtf variable:
gtf <- paste0(refDir, gtfName)
print(paste0("The gtf file is: ", gtf))

# load in the gene coordinates from the gtf file:
genes <- import(gtf)
# change chromosome name to chr1 for bam loading purposes:
genes@seqnames <- gsub("chrL1", "chr1", genes@seqnames)

# save RData file:
#save.image(file="/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/repGTF_loaded.RData")
#load(file="/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/repGTF_loaded.RData")
save.image(file="/Users/jamestorpy/Documents/Garvan/phd/projects/hgsoc_repeats/RNA-seq/Robjects/repGTF_loaded.RData")
#load(file="/Users/jamestorpy/Documents/Garvan/phd/projects/hgsoc_repeats/RNA-seq/Robjects/repGTF_loaded.RData")

# define the inFiles:
inFiles <- list.files(path = inDir, pattern = "Aligned.sortedByCoord.out.sorted.bam", recursive = T)
inFiles <- grep(".bai", inFiles, invert = T, value = T)
# or manually define inFiles:
#inFile <- paste0(inDir, "/kur_v1_subset/Aligned.sortedByCoord.out.sorted.bam")
# define the samplenames:
samplenames <- gsub("/Aligned.sortedByCoord.out.sorted.bam", "", inFiles)
print("The samples to be processed are: ")
print(samplenames)

# load the bamfile as GRanges object:
# assign column names for dataframe:
what <- c("qname","rname","strand","pos","qwidth")
# flag unmapped sequences to leave out:
flag <- scanBamFlag(isUnmappedQuery=FALSE)
# define parameters of bam scan:
param <- ScanBamParam(what=what, flag=flag)

# create empty list for bams:
bams <- list()

# set counter i to 1:
i=1
# index and load bams into list:
for (file in inFiles) {
  # define unique IDs for each file:
  uID <- strsplit(file, "/")[[1]][1]
  
  # index bams if necessary:
  if (file.exists(paste0(inDir, "/", uID, "/Aligned.sortedByCoord.out.sorted.bam.bai"))) {
    print("No need to index, bam has already been indexed.")
  } else {
    print(paste0("Indexing ", file, "."))
    indexBam(paste0(inDir, "/", file))
  }
  # load in bams:
  print(paste0("Loading ", file, " into bams list."))
  bams[i] <- scanBam(paste0(inDir, "/", file) ,param=param)
  
  i <<- i+1
}

#save.image(file="/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/repbam_loaded2.RData")
#load(file="/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/repbam_loaded2.RData")
save.image(file="/Users/jamestorpy/Documents/Garvan/phd/projects/hgsoc_repeats/RNA-seq/Robjects/repbam_loaded2.RData")
#load(file="/Users/jamestorpy/Documents/Garvan/phd/projects/hgsoc_repeats/RNA-seq/Robjects/repbam_loaded2.RData")

# create data frame of human chromosome lengths
seq_lengths <- seqlengths(Hsapiens)

# convert the bams to GRanges objects (chr1 actually represents chrL1):
bamGRs <- lapply(bams, function(x) {
  bamGR <- GRanges(
    seqnames = c(rep("chr1", length(x[[1]]))),
    ranges = IRanges(start = x$pos, width = x$qwidth),
    strand = x$strand,
    seqlengths = seq_lengths)
  return(bamGR)
})

# count overlaps of bams with annotation ranges:
counts <- lapply(bamGRs, function(x) {
  writeLines("\n")
  print("Counting overlaps...")
  overlaps <- countOverlaps(genes, x)
  counts <- as.data.frame(overlaps)
  counts$counts <- counts[ ,1]
  counts[ ,1] <- genes$gene_id
  colnames(counts) <- c("genes", "counts")
  print(paste0("Start of counts for sample:"))
  print(head(counts))
  return(counts)
  i <<- i+1
})
# name the count dataframes in the list:
names(counts) <- samplenames

# filter only counts above 0:
nonZero <- lapply(counts, function(x) {
  result <- x[grep(0, x$counts, invert = T),]
  return(result)
})
# name the nonZero dataframes in the list:
names(nonZero) <- samplenames

# filter only counts above 1:
aboveOne <- lapply(nonZero, function(x) {
  result <- x[grep(1, x$counts, invert = T),]
  return(result)
})
# name the aboveOne dataframes in the list:
names(aboveOne) <- samplenames

# determine the count instances to assess the range of counts:
cRange <- lapply(counts, function(x) {
  result <- unique(x$counts)
  return(result)
})
# name the count instances dataframes in the list:
names(cRange) <- samplenames

# save the counts as csv files:
for (n in 1:2) {
  if (file.exists(paste0(outDir, "/", samplenames[n], "/overlaps.csv"))) {
      print(paste0(outDir, "/", samplenames[n], "/overlaps.csv already exists"))
    } else {
      print(paste0("creating", outDir, "/", samplenames[n], "/overlaps.csv"))
      #write.csv(counts, file = paste0(outDir, "/", samplenames[n], "/overlaps.csv"))
    }
}

#save.image(file="/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/countsOverlapped-2.RData")
#load(file="/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/countsOverlapped-2.RData")
save.image(file="/Users/jamestorpy/Documents/Garvan/phd/projects/hgsoc_repeats/RNA-seq/Robjects/countsOverlapped-2.RData")
#load(file="/Users/jamestorpy/Documents/Garvan/phd/projects/hgsoc_repeats/RNA-seq/Robjects/countsOverlapped-2.RData")

