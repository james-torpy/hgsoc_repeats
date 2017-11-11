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

# load in repeat classes and types:
fullClasses <- read.table(paste0(refDir, "/fullClasses.txt"))
fullTypes <- read.table(paste0(refDir, "/fullTypes.txt"))

###### fix the following in countOverlaps.R ######
# delete erroneous telomere rows:
counts <- counts[-grep("telomere", counts$repeat_id),]
# add duplicate of first row:
counts <- rbind(counts[1,], counts)

# add classes and types columns to counts:
counts$type <- as.character(fullTypes$V1)
counts$class <- as.character(fullClasses$V1)
counts <- counts[,c(20, 21, 2:19)]

# aggregate counts so there is only one row for each type:
counts <- aggregate(. ~ type+class, counts, sum)

# save counts as RDS object:
saveRDS(counts, file = paste0(RobjectDir, "/counts.RDS"))

# define superclasses and group types
Classes <- sort(unique(fullClasses$V1))
Types <- sort(unique(fullTypes$V1))

# define superclasses and split classes into groups accordingly:
sClasses <- unique(gsub("/.*|\\?.*", "", Classes))
classGroups <- split(Classes, gsub("/.*|\\?.*", "", Classes))

# make a list with a data frame for each superclass with relevant classes:

scCounts <- lapply(classGroups, function(x) {
	# create empty data frame:
	m <- data.frame(matrix(ncol = ncol(counts)-2, nrow = 0))
	result <- cbind(data.frame(column1 = character, column2 = character), m)
	for (class in x) {
		result <- rbind(result, counts[counts$class %in% class,])
	}
	return(result)
})

# make a list with a data frame for each class with relevant types:
cCounts <- list()
i=1
for (c in Classes) {
	cCounts[[i]] <- counts[counts$class %in% c,]
	i <<- i+1
}
names(cCounts) <- Classes

# save scCounts and cCounts:
saveRDS(scCounts, file = paste0(RobjectDir, "/scCounts.RDS"))
saveRDS(cCounts, file = paste0(RobjectDir, "/cCounts.RDS"))



