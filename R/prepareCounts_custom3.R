### prepareCounts.R ###

# breaks cCounts and tCounts up into dataframe for each class and saves each a RData objects #

# define starting variables:
project <- "hgsoc_repeats"
methodName <- "method1"
expName <- "exp1"
Type <- "custom3"

# define directories:
#homeDir <- "/share/ScratchGeneral/jamtor/"
homeDir <- "/Users/jamestorpy/clusterHome"
#homeDir <- "/Users/jamestorpy/Documents/Garvan/phd/"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results/")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/")
inDir <- paste0(resultsDir, "/R/", methodName)
plotDir <- paste0(inDir, "/plots/linePlots/")

system(paste0("mkdir -p ", plotDir))


### 1. Load in inputs ###

# load in repeat counts, append each half for eah sample and add into list:
countFiles <- grep(paste0("subset|all"), list.files(paste0(RobjectDir, "/", expName, "/", Type, "_RepeatCounts/"), pattern="RepeatCounts", full.names = T), value=T, invert=T)
#countFiles <- grep(paste0("_", Type), countFiles, value=T)
sampleNames <- unique(gsub(paste0("_", Type, "RepeatCounts.rds"), "", basename(countFiles)))
uIDs <- unique(gsub("_custom3RepeatCounts_half[1-2].rds", "", sampleNames))

j=1
for (i in seq(1, length(countFiles), 2)) {
  # append half 1 of sample counts to half 2 and add to list:
  if (i==1) {
    Counts <- list(c(readRDS(file=countFiles[i]), readRDS(file=countFiles[i+1])))
  } else {
    Counts[[j]] <- c(readRDS(file=countFiles[i]), readRDS(file=countFiles[i+1]))
  }
j=j+1
}
names(Counts) <- uIDs

for (i in 1:length(Counts)) {
  j=1
  if (i==1) {
    for (j in 1:length(Counts[[i]])) {
      if (j==1) {
        newCounts <- list(data.frame(Counts[[i]][[j]]))
        j=j+1    
      } else {
        newCounts[[j]] <- Counts[[i]][[j]]
        j=j+1
      }
    }
  } else {
    j=1
    for (j in 1:length(Counts[[i]])) {
      newCounts[[j]][,i] <- Counts[[i]][[j]]
    }  
  }
}

newCounts <- lapply(newCounts, function(x) {
  colnames(x) <- uIDs
  return(x)
})
names(newCounts) <- names(Counts[[1]])
  
if (!file.exists(paste0(RobjectDir, "/", expName, "/", Type, "_RepeatCounts/all_", Type, "RepeatCountDFs.rds"))) {
  saveRDS(newCounts, file=paste0(RobjectDir, "/", expName, "/", Type, "_RepeatCounts/all_", Type, "RepeatCountDFs.rds"))
}
