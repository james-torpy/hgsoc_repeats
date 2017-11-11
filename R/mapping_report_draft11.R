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

if (!is.null(args[3])) {
  uID = args[3]
}

# define starting variables:
project <- "hgsoc_repeats"
sampleName <- "subset"

# define directories:
homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
reportDir <- paste0(resultsDir, "/method_reports/", methodName)
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/")


### 1. Load in the input files ###

# load in the report file:
Report <- read.table(file = rFile, sep = "\t")
Report[,1] <- as.character(Report[,1])

# load in the scCounts file:
scCounts <- readRDS(file = paste0(RobjectDir, "sc_", sampleName, "RepeatCounts.rds"))
######
sampleDF <- scCounts[[grepl(uID, names(scCounts))]]
LINEtypeI <- sampleDF[grepl("LINE Type I transposons", sampleDF$repeat_id), 1]
######


### 2. Enter LINEtypeI counts into Report: ###
Report[8,1] <- "Total number of LINEtypeI counts:"
Report[grep("LINEtypeI counts", Report[,1]),2] <- LINEtypeI

# write Report as tab-delimited table:
write.table(Report, file = paste0(reportDir, "/", uID, "_final_report.txt", sep = "\t"), quote = F, col.names = F, row.names = F)
