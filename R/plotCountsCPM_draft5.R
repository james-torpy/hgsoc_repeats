# load packages needed:
library(ggplot2)
library(reshape2)

# define starting variables:
project <- "hgsoc_repeats"
methodName <- "method1"
sampleName <- "subset"
Types <- c("sc", "c", "t")

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/")
inDir <- paste0(resultsDir, "/R/", methodName)
plotDir <- paste0(inDir, "/plots/linePlots")

system(paste0("mkdir -p ", plotDir))

### 1. Load in inputs ###

# load in inputs:
scCountsDF <- readRDS(file = paste0(RobjectDir, "/sc_", sampleName, "RepeatCountsDF.rds"))
cCountsDF <- readRDS(file = paste0(RobjectDir, "/c_", sampleName, "RepeatCountsDF.rds"))
tCountsDF <- readRDS(file = paste0(RobjectDir, "/t_", sampleName, "RepeatCountsDF.rds"))
lSizes <- unlist(readRDS(file = paste0(RobjectDir, "/lib", sampleName, "Sizes.rds")))


### 2. Plot counts ###
j=1
for (Counts in c(scCountsDF, cCountsDF, tCountsDF)) {
  Type <- Types[j]
  
  CPM <- Counts
  for (s in lSizes) {
    while (i < length(lSizes)) {
      Counts[,i] <- (Counts[,i]/s)*1000000
      i=i+1
    }
  }
  pCounts <- melt(df, varnames = c("repeat_id", "sample"), value.name = "counts")
  pCPM <- melt(CPM, varnames = c("repeat_id", "sample"), value.name = "CPM")
  
  # sort scpCounts data frames according to repeat_id:
  sort_rID <- function(x) {
    if (colnames(x)[1] == "repeat_id") {
      return(x[with(x, order(repeat_id)),])
    } else {
      return(x)
    }
  }
  
  pCounts <- sort_rID(pCounts)
  pCPM <- sort_rID(pCPM)
  
  # remove duplicate numbers:
  pCounts$sample <- gsub("[1-4]", "", pCounts$sample)
  pCPM$sample <- gsub("[1-4]", "", pCPM$sample)
  
  # aggregate duplicates and all of each class by mean:
  if (colnames(x)[1] == "repeat_id") {
    pCounts <- aggregate(counts~sample+repeat_id, pCounts, mean)
    } else {
      aggregate(counts~sample, pCounts, mean)
    }
  
  agg_dupesCPM <- function(x) {
    if (colnames(x)[1] == "repeat_id") {
      return(aggregate(CPM~sample+repeat_id, x, mean))
    } else {
      return(aggregate(CPM~sample, x, mean))
    }
  }
  
  scpCounts <- agg_dupes(scpCounts)
  scpCPM <- agg_dupesCPM(scpCPM)
  
  cpCounts <- lapply(cpCounts, agg_dupes)
  cpCPM <- lapply(cpCPM, agg_dupesCPM)
  
  # order levels of sample factor, putting controls first:
  orderS <- function(x) {
    x$sample <- factor(x$sample, levels=x$sample[c(1,5,3,4,6,2)]) 
    return(x)
  }
  
  scpCounts <- orderS(scpCounts)
  scpCPM <- orderS(scpCPM)
  
  cpCounts <- lapply(cpCounts, orderS)
  cpCPM <- lapply(cpCPM, orderS)
  
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
  scpCounts4Log <- add2zero(scpCounts)
  scpCPM4Log <- add2zero(scpCPM)
  
  cpCounts4Log <- lapply(cpCounts, add2zero)
  cpCPM4Log <- lapply(cpCPM, add2zero)
  
  plot_it <- function(x, logNo="nope") {
    if (ncol(x) == 3) {
      p <- ggplot(x, aes(x=sample, y=counts, group=repeat_id, colour=repeat_id))
      p <- p + geom_line()
      p <- p + theme(axis.text.x = element_text(angle = 90), legend.position = "bottom")
      if (logNo=="log10") {
        p <- p + scale_y_log10()
        return(p)
      } else {
        return(p)
      }
    } else {
      p <- ggplot(x, aes(x=sample, y=counts, group = 1))
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
  
  scpPlot <- plot_it(scpCounts)
  scpPlotLog10 <- plot_it(scpCounts4Log, logNo="log10")
  scpCPMPlot <- plot_itCPM(scpCPM)
  scpCPMPlotLog10 <- plot_itCPM(scpCPM4Log, logNo="log10")
  
  
  cpPlots <- lapply(cpCounts, plot_it)
  cpPlotsLog10 <- lapply(cpCounts4Log, plot_it, logNo="log10")
  cpCPMPlots <- lapply(cpCPM, plot_itCPM)
  cpCPMPlotsLog10 <- lapply(cpCPM4Log, plot_itCPM, logNo="log10")
  
  # remove slash from DNA/hATCounts for tpCounts as this causes problems when saving:
  #names(tpCounts) <- gsub("/", "", names(tpCounts))
  
  if (file.exists(paste0(plotDir, "/scCounts.pdf"))) {
    print(paste0(plotDir, "/scCounts.pdf already exists, no need to create"))
  } else {
    print(paste0("Creating ", plotDir, "/scCounts.pdf"))
    pdf(file = paste0(plotDir, "/scCounts.pdf"), width = 10, height=10)
    print(scpPlot)
    dev.off()
    i <<- i+1
  }
  
  if (file.exists(paste0(plotDir, "/scCountsLog10.pdf"))) {
    print(paste0(plotDir, "/scCountsLog10.pdf already exists, no need to create"))
  } else {
    print(paste0("Creating ", plotDir, "/scCountsLog10.pdf"))
    pdf(file = paste0(plotDir, "/scCountsLog10.pdf"), width = 10, height=10)
    print(scpPlot)
    dev.off()
    i <<- i+1
  }
  
  pdf_it <- function(x, log="nope") {
    for (e in names(x)) {
      if (log=="nope") {
        if (file.exists(paste0(plotDir, "/", e, "Counts.pdf"))) {
          print(paste0(plotDir, "/", e, "Counts.pdf already exists, no need to create"))
        } else {
          print(paste0("Creating ", plotDir, "/", e, "Counts.pdf"))
          pdf(file = paste0(plotDir, "/", e, "Counts.pdf"), width = 10, height=10)
          print(x[[i]])
          dev.off()
          i <<- i+1
        }
      } else {
        if (file.exists(paste0(plotDir, "/", e, "CountsLog10.pdf"))) {
          print(paste0(plotDir, "/", e, "CountsLog10.pdf already exists, no need to create"))
        } else {
          print(paste0("Creating ", plotDir, "/", e, "CountsLog10.pdf"))
          pdf(file = paste0(plotDir, "/", e, "CountsLog10.pdf"), width = 10, height=10)
          print(x[[i]])
          dev.off()
          i <<- i+1  
        }
      }
    }
  }
  
  
  i=1
  pdf_it(cpPlots)
  i=1
  pdf_it(cpPlotsLog10, "log10")
  
  
  if (file.exists(paste0(plotDir, "/scCPM.pdf"))) {
    print(paste0(plotDir, "/scCPM.pdf already exists, no need to create"))
  } else {
    print(paste0("Creating ", plotDir, "/scCPM.pdf"))
    pdf(file = paste0(plotDir, "/scCPM.pdf"), width = 10, height=10)
    print(scpCPMPlot)
    dev.off()
    i <<- i+1
  }
  
  if (file.exists(paste0(plotDir, "/scCPMLog10.pdf"))) {
    print(paste0(plotDir, "/scCPMLog10.pdf already exists, no need to create"))
  } else {
    print(paste0("Creating ", plotDir, "/scCPMLog10.pdf"))
    pdf(file = paste0(plotDir, "/scCPMLog10.pdf"), width = 10, height=10)
    print(scpPlot)
    dev.off()
    i <<- i+1
  }
  
  pdf_it <- function(x, log="nope") {
    for (e in names(x)) {
      if (log=="nope") {
        if (file.exists(paste0(plotDir, "/", e, "CPM.pdf"))) {
          print(paste0(plotDir, "/", e, "CPM.pdf already exists, no need to create"))
        } else {
          print(paste0("Creating ", plotDir, "/", e, "CPM.pdf"))
          pdf(file = paste0(plotDir, "/", e, "CPM.pdf"), width = 10, height=10)
          print(x[[i]])
          dev.off()
          i <<- i+1
        }
      } else {
        if (file.exists(paste0(plotDir, "/", e, "CPMLog10.pdf"))) {
          print(paste0(plotDir, "/", e, "CPMLog10.pdf already exists, no need to create"))
        } else {
          print(paste0("Creating ", plotDir, "/", e, "CPMLog10.pdf"))
          pdf(file = paste0(plotDir, "/", e, "CPMLog10.pdf"), width = 10, height=10)
          print(x[[i]])
          dev.off()
          i <<- i+1  
        }
      }
    }
  }
  
  
  i=1
  pdf_it(cpCPMPlots)
  i=1
  pdf_it(cpCPMPlotsLog10, "log10")
  
  j <<- j+1
}
