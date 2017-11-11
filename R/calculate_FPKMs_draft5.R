### calcuate_FPKMs.R ###
# This script takes a data frame of counts and a list of library sizes, and calculates FPKMs
# from these

### 0. Set up variables and directories ###

# load packages needed:
library(ggplot2)
library(reshape2)
library(edgeR)

# define starting variables:
project <- "hgsoc_repeats"
methodName <- "method1"

# define directories:
homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/")
refDir <- paste0(projectDir, "/RNA-seq/refs")
inDir <- paste0(resultsDir, "/R/", methodName)


### 1. Load in inputs ###

# load in counts:
counts <- read.csv(file = paste0(inDir, "/counts.csv"))
# remove first column of counts:
counts <- counts[,-1]
# rename repeats column of counts:
colnames(counts)[1] <- "repeat_id"

# subset counts:
#counts <- head(counts, 4000000)

# load in gene widths
widths <- readRDS(paste0(RobjectDir, "widths.rds"))

# load in library sizes:
lSizes <- unlist(readRDS(file = paste0(RobjectDir, "/libSizes.rds")))

# load in repeat classes and types:
fullClasses <- read.table(paste0(refDir, "/fullClasses.txt"))
fullTypes <- read.table(paste0(refDir, "/fullTypes.txt"))

# add classes and types columns to counts:
counts$class <- fullClasses$V1
counts$type <- fullTypes$V1
counts <- counts[,c(20, 21, 2:19)]

classes <- sort(unique(fullClasses))
types <- sort(unique(fullTypes))

# define superclasses and split classes into groups accordingly:
sClasses <- unique(gsub("/.*|\\?.*", "", Classes))
classGroups <- split(Classes, gsub("/.*|\\?.*", "", Classes))







