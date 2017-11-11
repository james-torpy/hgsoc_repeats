### rpCompBarPlot_prepare_cCounts.R ###

# breaks cCounts and tCounts up into dataframe for each class and saves each a RData objects #

# define starting variables:
project <- "hgsoc_repeats"
methodName <- "method1"
sampleName <- ""
types <- c("c", "t")

# define directories:
#homeDir <- "/share/ScratchGeneral/jamtor/"
homeDir <- "/Users/jamestorpy/clusterHome"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results/")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/")
inDir <- paste0(resultsDir, "/R/", methodName)
plotDir <- paste0(inDir, "/plots/linePlots/")

system(paste0("mkdir -p ", plotDir))


### 1. Load in inputs ###

# load in repeat counts into list:

for (Type in types) {
  countFiles <- grep(paste0("subset|_custom|all"), list.files(paste0(RobjectDir, "/", Type, "_RepeatCountDFs/"), pattern="RepeatCountsDF.rds", full.names = T), value=T, invert=T)
  countFiles <- grep(paste0("_", Type), countFiles, value=T)
  for (i in 1:length(countFiles)) {
    sampleName <- basename(gsub(paste0("_", Type, "RepeatCountsDF.rds"), "", countFiles[i]))
    if (i==1) {
      countL <- list(readRDS(file=countFiles[i]))
      names(countL)[i] <- sampleName
    } else {
      countL[[i]] <- readRDS(file=countFiles[i])
      names(countL)[i] <- sampleName
    }
  }
  
  
  ### 2. Prepare and save counts dfs ###
  
  # convert lists into df with all samples:
  countsDF <- do.call("cbind", countL)
  colnames(countsDF) <- names(countL)
  
  saveRDS(countsDF, file=paste0(RobjectDir, "/", Type, "_RepeatCountDFs/all_", Type, "RepeatCountsDF.rds"))
}

