
library(reshape2)
library(ggplot2)
library(pheatmap)


# define starting variables:
project <- "hgsoc_repeats"
methodName <- ""
refName <- "human-89.repeats.tab"
expName <- "exp1"
STypes <- c("FT", "erPT", "mrPT", "arPT", "prPT", "rfPT", "lrcT", "msST", "pAF", "rcAF")
CTypes <- c("arPT", "erPT", "FT", "lrcT", "mrPT", "msST", "pAF", "prPT", "rcAF", "rfPT")
annot <- "c"

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
#homeDir <- "/Users/jamestorpy/Documents/Garvan/phd"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/", expName, "/")
plotDir <- paste0(resultsDir, "/R/", expName, "/plots/DEplots/")

system(paste0("mkdir -p ", plotDir))


### 1. Load in DE results ###

# topRep is top 15 with FDR<0.1, strRep is stricter - top 15 with FDR<0.01
repGeneL <- readRDS(file=paste0(RobjectDir, "/", annot, "_DEreps.rds"))
topRep <- readRDS(file=paste0(RobjectDir, "/", annot, "_topRep.rds"))
strRep <- readRDS(file=paste0(RobjectDir, "/", annot, "_strRep.rds"))

allRep <- list(repGeneL, topRep, strRep)
names(allRep) <- c("repGene", "topRep", "strRep")

for (i in 1:length(allRep)) {
  # remove lrcT as only has 1 replicate:
  rep <- c(allRep[[i]][1:2], allRep[[i]][4:length(allRep[[i]])])
  j=1
  rep <- lapply(rep, function(x) {
    x$repeat_id <- rownames(x)
    x$sample <- rep(names(rep)[j], nrow(x))
    j <<- j+1
    return(x)
  })

  repDF <- do.call("rbind", rep)
  df <- subset(repDF, select=c(repeat_id, sample, logFC))

  pDF <- dcast(df, repeat_id ~ sample)
  rownames(pDF) <- pDF$repeat_id
  pDF <- subset(pDF, select=-repeat_id)

  #  log10DF <- apply(pDF, 2, function(x) {
  #    return(log10(2^x))
  #  })

  if (i==1) {
    pDFs <- list(pDF)
  } else {
    pDFs[[i]] <- pDF
  }
}

# need to add missing cols to all the samples in DE script!

  