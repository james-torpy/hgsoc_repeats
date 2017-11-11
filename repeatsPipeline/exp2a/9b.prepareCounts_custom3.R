### prepareCounts.R ###

# breaks cCounts and tCounts up into dataframe for each class and saves each a RData objects #

# define starting variables:
project <- "hgsoc_repeats"
expName <- "exp2"
Type <- "custom3"

# define directories:
homeDir <- "/share/ScratchGeneral/jamtor/"
#homeDir <- "/Users/jamestorpy/clusterHome"
#homeDir <- "/Users/jamestorpy/Documents/Garvan/phd/"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results/")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/")
inDir <- paste0(resultsDir, "/R/", expName)
plotDir <- paste0(inDir, "/plots/linePlots/")
rawDir <- paste0(projectDir, "/RNA-seq/raw_files/fullsamples/bowtell_primary/")

system(paste0("mkdir -p ", plotDir))


### 1. Load in inputs ###

# load in repeat counts, append each half for eah sample and add into list:
countFiles <- grep("Ext", grep(paste0("subset|all"), list.files(paste0(RobjectDir,
  "/", expName, "/", Type, "_RepeatCounts/"), pattern="RepeatCounts", full.names = T), 
  value=T, invert=T), value=T, invert=T)
#countFiles <- grep(paste0("_", Type), countFiles, value=T)
uIDs <- unique(gsub(paste0("_", Type, "RepeatCounts_half[1-2]a.rds"), "", basename(countFiles)))

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
  
# load in patient ids:
Key <- read.table(file=paste0(rawDir, "/sampleKey.txt"), sep=" ", fill=T)[,3:4]
Key <- Key[3:nrow(Key),]
Key$V4 <- gsub("_[A-Z].*$", "", gsub("^.*AOCS_", "AOCS_", Key$V4))
Key <- Key[with(Key, order(V3)), ]
# only include samples in newCounts:
Key <- dplyr::filter(Key, Key$V3 %in% colnames(newCounts[[1]]))

# add patient ids to sample names insteat of arbitrary numbers:
newCounts <- lapply(newCounts, function(x) {
  colnames(x) <- paste0(Key$V4, "_", gsub("[0-9]", "", colnames(x)))
  return(x)
})

if (!file.exists(paste0(RobjectDir, "/", expName, "/", Type, "_RepeatCounts/all_", Type, "RepeatCountDFs_a.rds"))) {
  saveRDS(newCounts, file=paste0(RobjectDir, "/", expName, "/", Type, "_RepeatCounts/all_", Type, "RepeatCountDFs_a.rds"))
}
