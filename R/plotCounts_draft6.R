# load packages needed:
library(ggplot2)
library(reshape2)

# define starting variables:
project <- "hgsoc_repeats"
methodName <- "method1"
sampleType <- "subset"
Types <- c("sc", "c", "t")

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/")
inDir <- paste0(resultsDir, "/R/", methodName)
plotDir <- paste0(inDir, "/plots")

system(paste0("mkdir -p ", plotDir))


### 1. Load in inputs ###

# load in inputs:
scCounts <- readRDS(file = paste0(RobjectDir, "/sc_", sampleType, "RepeatCountsDF.rds"))
cCounts <- readRDS(file = paste0(RobjectDir, "/c_", sampleType, "RepeatCountsDF.rds"))
tCounts <- readRDS(file = paste0(RobjectDir, "/t_", sampleType, "RepeatCountsDF.rds"))

lSizes <- unlist(readRDS(file = paste0(RobjectDir, "/", "lib", sampleType, "Sizes.rds")))


### 2. Plot sc counts ###

# prepare data frame for plotting:
melt_it <- function(x) {
  return(melt(x, value.name = "counts", variable.name = "sample"))
}

scpCounts <- lapply(scCounts, melt_it)

# remove replicate numbers from sample IDs:
remove_rep <- function(x) {
  temp <- gsub("[1-4]", "", x[,2])
  x[,2] <- as.factor(gsub(".counts", "", temp))
  return(x)
}
scpCounts <- lapply(scpCounts, remove_rep)

###### debug following ######
# sort scpCounts data frames according to class:
scpCounts <- lapply(scpCounts, function(x) {
  return(x[with(x, order(class)),])
})

# aggregate duplicates and all of each class by mean:
scpCounts <- lapply(scpCounts, function(x) {
  temp <- aggregate(counts~sample+class, x, mean)
})


# order levels of sample factor, putting controls first:
orderS <- function(x) {
  x$sample <- factor(x$sample,levels(x$sample)[c(1,5,3,4,6,2)])
  return(x)
}

scpCounts <- lapply(scpCounts, orderS)

# replace all zeros with tiny number so they are loggable:
add2zero <- function(x) {
  counts <- x[,3]
  counts[counts==0] <- 1
  x[,3] <- counts
  return(x)
}
scpCounts4Log <- lapply(scpCounts, add2zero)


### 3. Plot other counts ###

i=2
for (o in c(cCounts, tCounts)) {
	#  prepare data for plotting:
	pCounts <- lapply(Counts, melt_it)
	# remove replicate numbers from sample IDs:
	pCounts <- lapply(pCounts, remove_rep)

	# sort cpCounts data frames according to type:
	pCounts <- lapply(pCounts, function(x) {
  		return(x[with(x, order(type)),])
	})

	pCounts <- lapply(pCounts, function(x) {
  		temp <- aggregate(counts~sample+type, x, mean)
	})

	pCounts <- lapply(pCounts, orderS)

	# remove slash from DNA/hATCounts for cpCounts as this causes problems when saving:
	names(pCounts) <- gsub("/", "", names(pCounts))

	assign(pCounts, paste0(Types[i], "pCounts"))

	assign(lapply(pCounts, add2zero), paste0(Types[i], "pCounts4Log"))

}



scpPlots <- lapply(scpCounts, plot_it)
scpPlotsLog10 <- lapply(scpCounts4Log, plot_it, cat="class", logNo="log10")

cpPlots <- lapply(cpCounts, plot_it, cat="type")
cpPlotsLog10 <- lapply(cpCounts4Log, plot_it, cat="type", logNo="log10")

pdf_it <- function(x, log="nope") {
  for (e in names(x)) {
    if (log=="nope") {
      print(paste0("Creating ", plotDir, "/", e, "Counts.pdf"))
      pdf(file = paste0(plotDir, "/", e, "Counts.pdf"), width = 10, height=10)
      print(x[[i]])
      dev.off()
      i <<- i+1
    } else {
      print(paste0("Creating ", plotDir, "/", e, "CountsLog10.pdf"))
      pdf(file = paste0(plotDir, "/", e, "CountsLog10.pdf"), width = 10, height=10)
      print(x[[i]])
      dev.off()
      i <<- i+1
    }
  }
}

i=1
pdf_it(scpPlots)
i=1
pdf_it(scpPlotsLog10, "log10")
i=1
pdf_it(cpPlots)
i=1
pdf_it(cpPlotsLog10, "log10")