# This R script counts overlaps of multiple RNA-seq data sets with a GTF file of transcripts
# and outputs a counts .csv file for each sample

### 0. Set up variables and directories ###
methodName <- "method1"
load(file=paste0("/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/Robjects/", methodName, "repGTF_loaded.RData"))

# load packages needed:
library(GenomicRanges)
library(ShortRead)
library(rtracklayer)
library("BSgenome.Hsapiens.UCSC.hg38")
library(ggplot2)

# define starting variables:
project <- "hgsoc_repeats"
methodName <- "method1"
gtfName <- "human-89.repeats4.gtf"

# define directories:
#homeDir <- "/Users/jamestorpy/clusterHome/"
homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project)
bamDir <- paste0(projectDir, "/RNA-seq/results/star/", methodName)
refDir <- paste0(projectDir, "/RNA-seq/refs/")
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/")
outDir <- paste0(resultsDir, "/R/", methodName)


load(file=paste0("/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/Robjects/", methodName, "_repbam_loaded.RData"))

# create data frame of human chromosome lengths:
seq_lengths <- seqlengths(Hsapiens)

# convert the bams to GRanges objects:
bamGRs <- lapply(bamFiles, function(x) {
  bamGR <- GRanges(
    seqnames = x$rname,
    ranges = IRanges(start = x$pos, width = x$qwidth),
    strand = x$strand,
    seqlengths = seq_lengths,
    qnames = x$qname)
  return(bamGR)
})
names(bamGRs) <- sampleNames


# determine library size of each data set:
libSizes <- lapply(bamGRs, function(x) {
  return(length(x@seqnames))
})

# save the libSizes as RDS file:
saveRDS(libSizes, file = paste0(RobjectDir, "/libSizes.rds"))


### 4. Count overlaps of bams with repeat annotation gtf:

# count overlaps of bams with annotation ranges:
counts <- lapply(bamGRs, function(x) {
  writeLines("\n")
  print("Counting overlaps...")
  overlaps <- countOverlaps(fGenes, x)
  counts <- as.data.frame(overlaps)
  counts$counts <- counts[ ,1]
  counts[ ,1] <- fGenes$gene_id
  colnames(counts) <- c("genes", "counts")
  print(paste0("Start of counts for sample:"))
  print(head(counts))
  return(counts)
  i <<- i+1
})
# name the count dataframes in the list:
names(counts) <- sampleNames

# convert counts list into dataframe with samples as rows
counts <- do.call("cbind", counts)
colnames(counts)[1] <- "repeat"
nameInd <- grep("genes", colnames(counts))
counts <- counts[ , -nameInd]

# select only first of triplicate for development purposes:
#sCounts <- counts[ , c(-3, -4, -6, -7, -9, -10, -12, -13, -15, -16, -18, -19)]

# save the counts as csv files:
if (file.exists(paste0(outDir, "/counts.csv"))) {
  print(paste0(outDir, "/counts.csv"))
} else {
  print(paste0("creating", outDir, "/counts.csv"))
  write.csv(counts, file = paste0(outDir, "/counts.csv"))
}

save.image(file=paste0("/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/Robjects/", methodName, "_countsOverlapped.RData"))

