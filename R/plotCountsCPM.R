# load packages needed:
library(ggplot2)
library(reshape2)
library(dplyr)
library(Rmisc)

# define starting variables:
project <- "hgsoc_repeats"
methodName <- "method1"
expName <- "exp1"
Type <- "custom3"

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
#homeDir <- "/Users/jamestorpy/Documents/Garvan/phd/"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/")
inDir <- paste0(resultsDir, "/R/", methodName)
plotDir <- paste0(inDir, "/plots/linePlots")

system(paste0("mkdir -p ", plotDir))

### 1. Load in inputs ###

# load in inputs:
countDFs <- readRDS(file = paste0(RobjectDir, "/", expName, "/", Type, "_RepeatCounts/all_",  Type, "RepeatCountDFs.rds"))

if (file.exists(paste0(RobjectDir, "/", expName, "/", "/lib", "Sizes.rds"))) {
  lSizes <- readRDS(paste0(RobjectDir, "/", expName, "/", "/lib", "Sizes.rds"))
} else {
  # load in library sizes:
  libFiles <- grep("subset", list.files(paste0(RobjectDir, "/", expName, "/libSizes"), pattern="\\_libSize", full.names = T), value = T, invert = T)
  for (i in 1:length(libFiles)) {
    if (i==1) {
      lSizes <- c(readRDS(file=libFiles[i]))
    } else {
      lSizes[i] <- readRDS(file=libFiles[i])
    }
    print(i)
    i=i+1
  }
  saveRDS(lSizes, file=paste0(RobjectDir, "/", expName, "/", "/libSizes.rds"))
}


# remove second element from countDFs as it is a duplicate:
countDFs <- c(countDFs[1], countDFs[3:length(countDFs)])

# split dfs with rows > 8:
c=1
n=1
for (i in 1:length(countDFs)) {
  if (nrow(countDFs[[i]]) > 8) {
    if (c==1) {
      spl <- split(countDFs[[i]], c( rep(1, round(nrow(countDFs[[i]])/3)), rep(2, round(nrow(countDFs[[i]])/3)), rep(3, (nrow(countDFs[[i]])-2*(round(nrow(countDFs[[i]])/3)))) ))
      names(spl)[1:(n+2)] <- c(paste0(names(countDFs)[i], 1), paste0(names(countDFs)[i], 2), paste0(names(countDFs)[i], 3))
      rmInd <- c(i)
    } else {
      spl <- c(spl,split(countDFs[[i]], c( rep(1, round(nrow(countDFs[[i]])/3)), rep(2, round(nrow(countDFs[[i]])/3)), rep(3, (nrow(countDFs[[i]])-2*(round(nrow(countDFs[[i]])/3)))) )))
      names(spl)[n:(n+2)] <- c(paste0(names(countDFs)[i], 1), paste0(names(countDFs)[i], 2), paste0(names(countDFs)[i], 3))
      rmInd[c] <- i
    }
    c=c+1
    n=n+3
  }
}

# check if any dfs did have rows > 10
if (exists("spl")) {
  # subset Counts to include only dfs with rows < 10:
  `%notin%` <- function(x,y) !(x %in% y) 
  ind <- seq(1:length(countDFs))[seq(1:length(countDFs)) %notin% rmInd]
  
  # add merged dfs from above to Counts:
  countDFs <- c(countDFs[ind], spl)
}

### 2. Calculate CPM and counts per million repeat reads, prepare for plotting ###
j=1
for (Counts in countDFs) {
  
  # calculate CPMs:
  CPM <- as.data.frame(t(t(Counts)/lSizes)*1000000)
  
  # calculate total repeat count size:
  rSizes <- apply(Counts, 2, sum)
  # calculate CPMRs
  CPMR <- as.data.frame(t(t(Counts)/rSizes)*1000000)
  
  # prepare each df for plotting:
  typeNames <- c("CPM", "CPMR")
  n=1
  for (df in list(CPM, CPMR)) {
    sampleNames <- colnames(df)
    df <- sapply(df, as.numeric)
    rownames(df) <- rownames(CPM)
    pDF <- melt(df, varnames = c("repeat_id", "sample"), value.name = typeNames[n])
      
    # sort scpCounts data frames according to repeat_id:
    sort_rID <- function(x) {
        if (colnames(x)[1] == "repeat_id") {
          return(x[with(x, order(repeat_id)),])
        } else {
          return(x)
        }
      }
    pDF <- sort_rID(pDF)
      
    # remove duplicate numbers:
    pDF$sample <- gsub("[1-4]", "", pDF$sample)
    
    # order levels of sample factor, putting controls first:
    orderS <- function(x) {
      x$sample <- factor(x$sample, levels=c("bowtell_FT", "gtx_ft", "grant_endo", "grant_endosis", "kur_v", "bowtell_PR"))
      return(x)
    }
    pDF <- orderS(pDF)
      
    # create stats summary of the data:
    statDF <- summarySE(pDF, measurevar = typeNames[n], groupvar = c("repeat_id", "sample"))
    
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
    pDF4Log <- add2zero(pDF)
    
    # create stats summary of the data:
    statDF4Log <- summarySE(pDF, measurevar = typeNames[n], groupvar = c("repeat_id", "sample"))
      
    plot_itdf <- function(x, logNo="nope") {
      if (ncol(x) == 7) {
        p <- ggplot(x, aes(x=sample, y=x[,4], group=repeat_id, colour=repeat_id))
        p <- p + geom_errorbar(aes(ymin=eval(parse(text=typeNames[n]))-se, ymax=eval(parse(text=typeNames[n]))+se), width=0.1)
        p <- p + geom_line()
        p <- p + theme(axis.text.x = element_text(angle = 90), legend.position = "bottom")
        p <- p + ylab(typeNames[n])
        if (logNo=="log10") {
          p <- p + scale_y_log10()
          return(p)
        } else {
          return(p)
        }
      } else {
        p <- ggplot(x, aes(x=sample, y=x[,4], group = 1))
        p <- p + geom_errorbar(aes(ymin=eval(parse(text=typeNames[n]))-se-se, ymax=eval(parse(text=typeNames[n]))-se+se), width=0.1)
        p <- p + geom_line()
        p <- p + theme(axis.text.x = element_text(angle = 90), legend.position = "bottom")
        p <- p + ylab(typeNames[n])
        if (logNo=="log10") {
          p <- p + scale_y_log10()
            return(p)
          } else {
            return(p)
          }
        }
      }
      
      pDFPlot <- plot_itdf(statDF)
      pDFPlotLog10 <- plot_itdf(statDF4Log, logNo="log10")
      
      if (file.exists(paste0(plotDir, "/", Type, "_", names(countDFs)[j], "_", typeNames[n], ".pdf"))) {
        print(paste0(plotDir, "/", Type, "_", names(countDFs)[j], "_", typeNames[n], ".pdf already exists, no need to create"))
      } else {
        print(paste0("Creating ", plotDir, "/", Type, "_", names(countDFs)[j], "_", typeNames[n], ".pdf"))
        pdf(file = paste0(plotDir, "/", Type, "_", names(countDFs)[j], "_", typeNames[n], ".pdf"), width = 10, height=10)
        print(pDFPlot)
        dev.off()
      }
      
      if (file.exists(plotDir, "/", Type, "_", names(countDFs)[j], "_", typeNames[n], "Log10.pdf")) {
        print(paste0(plotDir, "/", Type, "_", names(countDFs)[j], "_", typeNames[n], "Log10.pdf already exists, no need to create"))
      } else {
        print(paste0("Creating ", plotDir, "/", Type, "_", names(countDFs)[j], "_", typeNames[n], "Log10.pdf"))
        pdf(file = paste0(plotDir, "/", Type, "_", names(countDFs)[j], "_", typeNames[n], "Log10.pdf"), width = 10, height=10)
        print(pDFPlotLog10)
        dev.off()
      }
    n=n+1
  }
  j=j+1
}

  