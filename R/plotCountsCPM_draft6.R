# load packages needed:
library(ggplot2)
library(reshape2)
library(dplyr)

# define starting variables:
project <- "hgsoc_repeats"
methodName <- "method1"
sampleName <- ""
Types <- c("sc", "c", "t")

# define directories:
#homeDir <- "/Users/jamestorpy/clusterHome/"
homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/old_240817/")
inDir <- paste0(resultsDir, "/R/", methodName)
plotDir <- paste0(inDir, "/plots/linePlots")

system(paste0("mkdir -p ", plotDir))

### 1. Load in inputs ###

# load in inputs:
scCountsDF <- readRDS(file = paste0(RobjectDir, "/sc_", sampleName, "RepeatCountsDF.rds"))
cCountsDF <- readRDS(file = paste0(RobjectDir, "c_RepeatCountDFs/", "/c_", sampleName, "RepeatCountsDF.rds"))
tCountsDF <- readRDS(file = paste0(RobjectDir, "/t_", sampleName, "RepeatCountDFs/t_", sampleName, "RepeatCountsDF.rds"))

if (sampleName == "subset") {
  lSizes <- unlist(readRDS(file = paste0(RobjectDir, "/lib", sampleName, "Sizes.rds")))
} else {
  if (file.exists(paste0(RobjectDir, "/lib", sampleName, "Sizes.rds"))) {
    lSizes <- readRDS(paste0(RobjectDir, "/lib", sampleName, "Sizes.rds"))
  } else {
    # load in library sizes:
    libFiles <- grep("subset", list.files(RobjectDir, pattern="libSize", full.names = T), value = T, invert = T)
    for (i in 1:length(libFiles)) {
      if (i==1) {
        lSizes <- c(readRDS(file=libFiles[i]))
      } else {
        lSizes[i] <- readRDS(file=libFiles[i])
      }
      print(i)
      i=i+1
    }
    saveRDS(lSizes, file=paste0(RobjectDir, "/lib", sampleName, "Sizes.rds"))
  }
}


### 2. Plot counts ###
j=1
for (Counts in c(scCountsDF, cCountsDF, tCountsDF)) {
  Type <- Types[j]
  
 
  CPM <- Counts
  i=1
  for (i in 1:length(lSizes)) {
    CPM[,i] <- (CPM[,i]/lSizes[i])*1000000
    print(i)
  }
  
  # calculate total repeat count size:
  rSizes <- apply(Counts, 2, sum)
  CPMR <- as.data.frame(t(t(Counts)/rSizes)*1000000)
  
  for (df in c(CPM, CPMR)) {
    CPM <- sapply(CPM, as.numeric)
    rownames(CPM) <- rownames(Counts)
    
    # order by row means:
    CPM <- as.data.frame(CPM[order(rowMeans(CPM)), ])
    CPM$repeat_id <- rownames(CPM)
    
    # split into 8-row dfs:
    CPM1 <- CPM[1:8,]
    CPM2 <- CPM[8:16,]
    CPM3 <- CPM[17:24,]
    CPM4 <- CPM[25:32,]
    CPM5 <- CPM[33:40,]
    CPM6 <- CPM[41:48,]
    CPM7 <- CPM[49:56,]
    CPM8 <- CPM[57:64,]
    
    i=1
    for (CPM in list(CPM1, CPM2, CPM3, CPM4, CPM5, CPM6, CPM7, CPM8)) {
      pCPM <- melt(CPM, variable.name = "sample", value.name = "CPM")
      
      # sort scpCounts data frames according to repeat_id:
      sort_rID <- function(x) {
        if (colnames(x)[1] == "repeat_id") {
          return(x[with(x, order(repeat_id)),])
        } else {
          return(x)
        }
      }
      
      pCPM <- sort_rID(pCPM)
      
      # remove duplicate numbers:
      pCPM$sample <- gsub("[1-4]", "", pCPM$sample)
      
      # aggregate duplicates and all of each class by mean:
      agg_dupesCPM <- function(x) {
        if (colnames(x)[1] == "repeat_id") {
          return(aggregate(CPM~sample+repeat_id, x, mean))
        } else {
          return(aggregate(CPM~sample, x, mean))
        }
      }
      
      pCPM <- agg_dupesCPM(pCPM)
      
      # order levels of sample factor, putting controls first:
      orderS <- function(x) {
        x$sample <- factor(x$sample, levels=x$sample[c(1,5,3,4,6,2)]) 
        return(x)
      }
      
      pCPM <- orderS(pCPM)
      
      # add 1 to all values so they are loggable:
      add2zero <- function(x) {
        if (ncol(x) == 3) {
          x[,3] <- x[,3]+1
          return(x)
        } else {
          x[,2] <- x[,2]+1
          return(x)
        }
      }
      pCounts4Log <- add2zero(pCounts)
      pCPM4Log <- add2zero(pCPM)
      
      plot_itCPM <- function(x, logNo="nope") {
        if (ncol(x) == 3) {
          p <- ggplot(x, aes(x=sample, y=CPM, group=repeat_id, colour=repeat_id))
          p <- p + geom_line()
          p <- p + theme(axis.text.x = element_text(angle = 90), legend.position = "bottom")
          if (logNo=="log10") {
            p <- p + scale_y_log10()
            return(p)
          } else {
            return(p)
          }
        } else {
          p <- ggplot(x, aes(x=sample, y=CPM, group = 1))
          p <- p + geom_line()
          p <- p + theme(axis.text.x = element_text(angle = 90), legend.position = "bottom")
          if (logNo=="log10") {
            p <- p + scale_y_log10()
            return(p)
          } else {
            return(p)
          }
        }
      }
      
      pCPMPlot <- plot_itCPM(pCPM)
      pCPMPlotLog10 <- plot_itCPM(pCPM4Log, logNo="log10")
      
      if (file.exists(paste0(plotDir, "/tCPM", i, ".pdf"))) {
        print(paste0(plotDir, "/tCPM", i, ".pdf already exists, no need to create"))
      } else {
        print(paste0("Creating ", plotDir, "/tCPM", i, ".pdf"))
        pdf(file = paste0(plotDir, "/tCPM", i, ".pdf"), width = 10, height=10)
        print(pCPMPlot)
        dev.off()
      }
      
      if (file.exists(paste0(plotDir, "/tCPM", i, "Log10.pdf"))) {
        print(paste0(plotDir, "/tCPM", i, "Log10.pdf already exists, no need to create"))
      } else {
        print(paste0("Creating ", plotDir, "/tCPM", i, "Log10.pdf"))
        pdf(file = paste0(plotDir, "/tCPM", i, "Log10.pdf"), width = 10, height=10)
        print( pCPMPlotLog10)
        dev.off()
      }
      i=i+1  
    }
  }
  
  }

  