# This R script outputs a report of the STAR mapping method chosen to handle repeats

### 0. Set up variables and directories ###

# load packages needed:
library(GenomicRanges)
library(ShortRead)
library(rtracklayer)
library("BSgenome.Hsapiens.UCSC.hg38")
library(ggplot2)

#read arguments from the command line:
args = commandArgs(trailingOnly = TRUE)
print(args)

# define report of sample:
if (!is.null(args[1])) {
  methodName = args[1]
}

if (!is.null(args[2])) {
  rFile = args[2]
}

# define starting variables:
project <- "hgsoc_repeats"
gtfName <- "human-89.repeats.gtf"
sampleName <- "fulltest"

# define directories:
homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
bamDir <- paste0(resultsDir, "/star/", methodName)
refDir <- paste0(projectDir, "/RNA-seq/refs/")
fastqDir <- paste0(projectDir, "/RNA-seq/raw_files/", sampleName)
rawDir <- paste0(projectDir, "/RNA-seq/raw_files")
genomeDir <- paste0(homeDir, "/genomes/L1")
reportDir <- paste0(resultsDir, "/method_reports/", methodName)


### 1. Load in the input files ###

# load in the report file:
Report <- read.table(file = rFile, sep = "\t")
Report[,1] <- as.character(Report[,1])

# define gtf variable:
gtf <- paste0(refDir, gtfName)
print(paste0("The gtf file is: ", gtf))
writeLines("\n")
# load in the gene coordinates from the gtf file:
genes <- import(gtf)

# define the bamFile:
dirName <- gsub("_report.txt", "", basename(rFile))
bamFile <- paste0(bamDir, "/", dirName, "/Aligned.sortedByCoord.out.bam")

print("The bam file to be processed is: ")
print(bamFile)
writeLines("\n")

# index the bam file if necessary:
if (file.exists(paste0(bamDir, "/Aligned.sortedByCoord.out.bam.bai"))) {
    print(paste0("No need to index, ", bamFile, ", has already been indexed"))
  } else {
    print(paste0("Indexing ", bamFile))
    indexBam(bamFile)
  }
writeLines("\n")

# assign column names for bam dataframe:
what <- c("qname","rname","strand","pos","qwidth")
# flag unmapped sequences to leave out:
flag <- scanBamFlag(isUnmappedQuery=FALSE)
# define parameters of bam scan:
param <- ScanBamParam(what=what, flag=flag)

# load in the bam file:
print(paste0("Loading ", bamFile))
Bam <- scanBam(bamFile, param=param)

# create data frame of human chromosome lengths:
seq_lengths <- seqlengths(Hsapiens)

# convert the bams to GRanges objects:
bamGR <- GRanges(
  seqnames = Bam[[1]]$rname,
  ranges = IRanges(start = Bam[[1]]$pos, width = Bam[[1]]$qwidth),
  strand = Bam[[1]]$strand,
  seqlengths = seq_lengths,
  qnames = Bam[[1]]$qname
)

#methodName <- "method1"
#save.image(file=paste0("/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/Robjects/", methodName, "_reportInputs_loaded.RData"))
#load(file=paste0("/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/Robjects/", methodName, "_reportInputs_loaded.RData"))


### 2. Find counts and qnames of overlaps with L1 annotations ###

# subset genes to only include L1s:
L1genes <- genes[grep("L1", genes$gene_id), ]

# fetch total number of counts of L1 overlaps:
print("Counting overlaps...")
oCounts <- countOverlaps(L1genes, bamGR)
L1counts <- sum(oCounts)

# enter L1counts into Report:
Report[8,1] <- "Total number of L1 counts:"
Report[grep("L1 counts", Report[,1]),2] <- L1counts

# write Report as tab-delimited table:
write.table(Report, file = paste0(reportDir, "/", dirName, "_final_report.txt", sep = "\t"), quote = F, col.names = F, row.names = F)
