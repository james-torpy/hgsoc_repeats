# load packages needed:
library(ggplot2)
library(reshape2)
library(RColorBrewer)

# define starting variables:
project <- "hgsoc_repeats"
methodName <- "method1"
sampleType <- ""
types <- c("c", "t")

# define directories:
#homeDir <- "/share/ScratchGeneral/jamtor/"
homeDir <- "/Users/jamestorpy/clusterHome"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results/")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/")
inDir <- paste0(resultsDir, "/R/", methodName)
plotDir <- paste0(inDir, "/plots/compBarplot/")

system(paste0("mkdir -p ", plotDir))


### 1. Load in inputs and QC check ###

# load cCounts into list:
Counts <- list(readRDS(file=paste0(RobjectDir, "/c_RepeatCountDFs/all_cRepeatCountsDF.rds")))
names(Counts)[1] <- "cCounts"

# load tCounts and split into groups:
tCounts <- readRDS(file=paste0(RobjectDir, "/t_RepeatCountDFs/all_", Type, "RepeatCountsDF.rds"))
groups <- unique(gsub("\\/.*|\\?", "", rownames(tCounts)))
i=2
for (g in groups) {
  Counts[[i]] <- tCounts[grep(g, rownames(tCounts)),]
  names(Counts)[i] <- paste0(g, "Counts")
  i=i+1
}

# select elements of counts elements with 1 or less repeat categories:
j=1
n=1
for (i in 1:length(Counts)) {
  if (nrow(Counts[[i]]) < 3) {
    if (j==1) {
      tempCounts <- list(Counts[[i]])
      j=j+1
    } else {
      tempCounts[[j]] <- Counts[[i]]
      j=j+1
    }
  } else {
    if (n==1) {
      newCounts <- list(Counts[[i]])
      names(newCounts)[n] <- names(Counts)[i]
      n=n+1
    } else {
      newCounts[[n]] <- Counts[[i]]
      names(newCounts)[n] <- names(Counts)[i]
      n=n+1
    }
  }
}

# convert nCounts to df and add as Counts list element:
newCounts[[length(newCounts)+1]] <- do.call("rbind", nCounts)
names(newCounts)[length(newCounts)] <- "other"


### 2. Calculate percentages and plot both dfs ###

for (i in 1:length(newCounts)) {
  # calculate percentages as new data frame:
  perCounts <- apply(newCounts[[i]], 2, function(x) {
    return((x/sum(x))*100)
  })
    
  # convert df to long format for plotting:
  plotCounts <- melt(perCounts, varnames = c("type", "sample"), value.name = "percentage")
    
  # define plot colour scheme:
  #tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
  #pal(tol21rainbow)
    
  # plot data as barplot:
  p <- ggplot(plotCounts, aes(x=sample, y=percentage))
  p <- p + geom_bar(stat="identity", aes(fill=type))
  p <- p + theme(axis.text.x = element_text(angle = 90), legend.position = "bottom")
  pdf(file = paste0(plotDir, "/", names(newCounts)[i], "_", sampleType, "_repeat_percent_compBarplot", ".pdf"), height = 10, width = 10)
  print(p)
  dev.off()
}
  


######
# remove null values from list of count elements:
#oneCounts <- Filter(Negate(is.null), oneCounts)

# select elements of counts elements with 2 or less repeat categories:
#twoCounts <- list()
#i=1
#twoCounts <- lapply(Counts, function(x) {
#  twoCounts[[i]] <- Filter(function(y) length(y)==2, x)
#  return(twoCounts)
#  i <<- i+1
#})