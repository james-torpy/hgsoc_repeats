# This script takes the following inputs:
# 1. gencode gtf
# 2. counts bam
# 3. repeat counts RDS
# 4. STAR Log.final.out from ribosomal mapping
# and creates barplots with the composition of protein-coding, non-coding, other, repeats and ribosomal counts
# with one bar per sample

### 0. Set up variables and directories ###

# load packages needed:
library(GenomicRanges)
library(ShortRead)
library("BSgenome.Hsapiens.UCSC.hg38")
library(reshape2)
library(ggplot2)

# define starting variables:
project <- "hgsoc_repeats"
methodName <- "method1"
gencodeName <- "gencode_v24_hg38_annotation.gtf"
ribName <- "human_rRNA.fa"
refName <- "human-89.repeats.tab"

# define directories:
#homeDir <- "/Users/jamestorpy/clusterHome/"
homeDir <- "/share/ScratchGeneral/jamtor/"
genomeDir <- paste0(homeDir, "genomes/")
gencodeDir <- paste0(genomeDir, "hg38_ercc/")
projectDir <- paste0(homeDir, "projects/", project, "/")
resultsDir <- paste0(projectDir, "RNA-seq/results/")
RobjectDir <- paste0(projectDir, "RNA-seq/Robjects/")
ribDir <- paste0(resultsDir, "/star/", methodName, "/ribosome_12/")
refDir <- paste0(projectDir, "/RNA-seq/refs")

### 1. Load in bamGRs and put into a list ###
# load in the command line arguments:
args = commandArgs(trailingOnly = TRUE)
for (n in 1:length(args)) {
  print(args[n])
}

if (!is.null(args[1])) {
  bamGRfile <- args[1]
}
print(paste0("bamGRfile is: ", bamGRfile))

if (!is.null(args[2])) {
  uID <- args[2]
}
print(paste0("uID is: ", uID))

if (!is.null(args[3])) {
  sampleType <- args[3]
}
print(paste0("sampleType is: ", sampleType))


### 1. Load bamGR ###
print(paste0("Reading in ", bamGRfile))
bamGR <- readRDS(file = paste0(RobjectDir, "/", bamGRfile))


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
temp <- (gencode[gencode$type %in% "gene"])

# select relevant entries for non-coding or protein coding annotations and reduce to aggregate overlapping annotations:
pcGencode <- reduce(temp[temp$gene_type %in% "protein_coding"])
ncGencode <- reduce(temp[temp$gene_type %in% c("lincRNA", "non_coding")])
otherGencode <- reduce(temp[!temp$gene_type %in% c("protein_coding", "lincRNA", "non_coding")])


### 2. Determine library size of data set ###
if (file.exists(paste0(RobjectDir, "/", uID, "_libSize.rds"))) {
  print(paste0("Loading ", RobjectDir, "/", uID, "_libSize.rds"))
  libSize <- readRDS(file = paste0(RobjectDir, "/", uID, "_libSize.rds"))
} else {
  libSize <- length(bamGR@seqnames)
  # save the libSize as RDS file:
  saveRDS(libSize, file = paste0(RobjectDir, "/", uID, "_libSize.rds"))
}

print(paste0("The libSize is: ", libSize))


### 3. Count overlaps between bams and gencode annotation:

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

# create function to count overlaps of bams with annotation ranges:
count_it <- function(x, annot = rGenes) {
  writeLines("\n")
  print("Counting overlaps...")
  counts <- as.data.frame(countOverlaps(annot, x))
  return(sum(counts))
}

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


### 5. Fetch number of reads mapped to ribosomal RNA ###

riboLog <- paste0(ribDir, "/", uID, "/Log.final.out")

tab <- read.table(riboLog, sep = "\t", fill = T, as.is = T)
Counts[5] <- as.numeric(tab[8,2]) + as.numeric(tab[23,2])
names(Counts)[5] <- "ribosome"


### 6. Save Counts as RDS object:
saveRDS(Counts, file=paste0(RobjectDir, "/", uID, "_compCounts.rds"))
