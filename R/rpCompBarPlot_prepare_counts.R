### rpCompBarPlot_prepare_cCounts.R ###

# breaks cCounts and tCounts up into dataframe for each class and saves each a RData objects #

# define starting variables:
project <- "hgsoc_repeats"
methodName <- "method1"
sampleName <- "subset"

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

# load in inputs:
cCounts <- readRDS(file = paste0(RobjectDir, "/c_", sampleName, "RepeatCounts.rds"))
tCounts <- readRDS(file = paste0(RobjectDir, "/t_", sampleName, "RepeatCounts.rds"))


### 2. Prepare and save counts ###

# for each class, cbind sample columns:
cCountsDFs <- lapply(cCounts, function(x) {
  result <- do.call("cbind", x)
  colnames(result) <- names(x)
  return(result)
})

for (i in 1:length(cCountsDFs)) {
  df <- cCountsDFs[[i]]
  saveRDS(df, file=paste0(RobjectDir, "/", names(cCountsDFs)[i], "_", sampleName, "_classCountsDF.rds"))
  print(i)
}

# for each type, cbind sample columns:
tCountsDFs <- lapply(tCounts, function(x) {
  result <- do.call("cbind", x)
  colnames(result) <- names(x)
  return(result)
})

for (i in 1:length(tCountsDFs)) {
  df <- tCountsDFs[[i]]
  saveRDS(df, file=paste0(RobjectDir, "/", names(tCountsDFs)[i], "_", sampleName, "_CountsDF.rds"))
  print(i)
}

