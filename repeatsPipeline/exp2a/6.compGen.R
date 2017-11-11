# This script takes the following inputs:
# 1. gencode gtf
# 2. counts bam
# 3. repeat counts RDS
# 4. STAR Log.final.out from ribosomal mapping
# and creates barplots with the composition of protein-coding, non-coding, other, repeats and ribosomal counts
# with one bar per sample

### 0. Set up variables and directories ###

args = commandArgs(trailingOnly = TRUE)
for (n in 1:2) {
  print(args[n])
}

if (!is.null(args[1])) {
  uID <- args[1]
}

if (!is.null(args[2])) {
  projectDir <- args[2]
}

library(GenomicRanges)
library(ShortRead)
library(rtracklayer)
library("BSgenome.Hsapiens.UCSC.hg38")
library(ggplot2)

# define starting variables:
project <- "hgsoc_repeats"
gencodeName <- "gencode_v24_hg38_annotation.gtf"
ribName <- "human_rRNA.fa"
refName <- "human-89.repeats.tab"

# define directories:
homeDir <- "/share/ClusterShare/thingamajigs/jamtor/"
RobjectDir <- paste0(projectDir, "/Robjects/")
bamDir <- paste0(projectDir, "/results/star/GC", uID)
genomeDir <- paste0("share/ScratchGeneral/jamtor/genomes/")
gencodeDir <- paste0(genomeDir, "hg38_ercc/")
ribDir <- paste0(projectDir, "/results/star/ribo")
refDir <- paste0(projectDir, "/refs")

inFile <- paste0(projectDir, "/results/star/GC/", uID, "/Aligned.sortedByCoord.out.bam")

### 1. Load gencode gtf and split into protein_coding and non_coding/linc gene entries ###

# define gencode variable:
gencodeFile <- paste0(gencodeDir, gencodeName)
if (file.exists(paste0(RobjectDir, "GCgenes.rds"))) {
  print("Loading gencode...")
  gencode <- readRDS(file=paste0(RobjectDir, "GCgenes.rds"))
} else {
  print(paste0("The gencode file is: ", gencodeFile))
  # load in the gene coordinates from the gencode file:
  gencode <- import(gencodeFile)
  # save gencode as RDF file:
  saveRDS(gencode, file=paste0(RobjectDir, "GCgenes.rds"))
}

# select only gene entries from gencode:
temp <- (gencode[gencode$type %in% "exon"])

# select relevant entries for non-coding or protein coding annotations and reduce to aggregate overlapping annotations:
pcGencode <- reduce(temp[temp$gene_type %in% "protein_coding"])
ncGencode <- reduce(temp[temp$gene_type %in% c("lincRNA", "non_coding")])
otherGencode <- reduce(temp[!temp$gene_type %in% c("protein_coding", "lincRNA", "non_coding")])


### 2. Load in bam file and convert to GenomicRanges object if needed:

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


### 3. Fetch library size ###

# determine library size of each data set:
libSize <- length(seqnames(bamGR))

print(paste0("The libSize is: ", libSize))


### 4. Count overlaps between bams and gencode annotation:

# count overlaps of bams with annotation ranges:
count_it <- function(x, bam=bamGR) {
  writeLines("\n")
  print(paste0("Counting overlaps..."))
  counts <- as.data.frame(countOverlaps(x, bam))
  return(sum(counts))
}

# count overlaps and aggregate gene types protein_coding, non_coding and other:
if (file.exists(paste0(RobjectDir, "/", uID, "_gcCounts.rds"))) {
  Counts <- readRDS(file=paste0(RobjectDir, "/", uID, "_gcCounts.rds"))
} else {
  writeLines("\n")
  print("Counting overlaps...")
  Counts <- list(sum(countOverlaps(pcGencode, bamGR)), sum(countOverlaps(ncGencode, bamGR)), sum(countOverlaps(otherGencode, bamGR)))
  names(Counts) <- c("protein_coding", "non_coding", "other")
  saveRDS(Counts, file=paste0(RobjectDir, "/", uID, "_gcCounts.rds"))
}


### 4. Count all repeats in samples ####

# load in repeats annotation and reduce:
if (file.exists(paste0(RobjectDir, "/", uID, "_TotalRepeatCounts.rds"))) {
  print(paste0("Loading ", RobjectDir, "/", uID, "_TotalRepeatCounts.rds"))
  repeatCount <- readRDS(file=paste0(RobjectDir, "/", uID, "_TotalRepeatCounts.rds"))
} else {
  print(paste0("Loading", RobjectDir, "/repeatsGR.rds"))
  rGenes <- reduce(readRDS(file=paste0(RobjectDir, "/repeatsGR.rds")))
  print("Counting overlaps with total repeats annotation")
  repeatCount <- count_it(bamGR, rGenes)
  print(paste0("Generating ", RobjectDir, "/", uID, "_TotalRepeatCounts.rds"))
  saveRDS(repeatCount, file=paste0(RobjectDir, "/", uID, "_TotalRepeatCounts.rds"))
}

# add and label repeats count to 'Counts':
Counts[[4]] <- repeatCount
names(Counts)[4] <- "repeats"

# remove the bam file and bamGR file:
print(paste0("Removing ", inFile, "..."))
system(paste0("rm ", inFile))

print(paste0("Removing ", RobjectDir, uID, "_bamGR.rds..."))
system(paste0("rm ", RobjectDir, uID, "_bamGR.rds"))


### 5. Fetch number of reads mapped to ribosomal RNA ###

riboLog <- paste0(ribDir, "/", uID, "/Log.final.out")

tab <- read.table(riboLog, sep = "\t", fill = T, as.is = T)
Counts[5] <- as.numeric(tab[8,2]) + as.numeric(tab[23,2])
names(Counts)[5] <- "ribosome"


### 6. Load in report and add Counts to it:
saveRDS(Counts, file=paste0(RobjectDir, "/", uID, "_compCounts.rds"))



