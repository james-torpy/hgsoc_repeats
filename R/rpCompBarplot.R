# load packages needed:
library(ggplot2)
library(reshape2)

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
plotDir <- paste0(inDir, "/plots/compBarplots/")

system(paste0("mkdir -p ", plotDir))


### 1. Load in inputs ###

# load in inputs:
inFiles <- list.files(RobjectDir, pattern=paste0("_subset_typeCountsDF.rds|_subset_classCountsDF.rds|sc_", sampleName, "RepeatCountsDF.rds"), full.names=T)

for (File in inFiles) {
  print(paste0("Loading ", File))
  df <- readRDS(file=File)
  
  if (nrow(df) < 2) {
    print("Not enough rows")
  } else {
    print("Generating barplot...")
    
    Type <- gsub(paste0("_", sampleName, ".*"), "", basename(File))
    
    # If df values are lists, unlist them:
    if (class(df[,1]) == "list") {
      print("Unlisting dataframe values")
      df <- sapply(df, unlist)
    }
    
    ### 2. Calculate percentages and plot ###
    
    # calculate percentages as new data frame:
    perDF <- apply(df, 2, function(x) {
      return((x/sum(x))*100)
    })
    
    plotDF <- melt(perDF, varnames = c("type", "sample"), value.name = "percentage")
    
    # plot data as barplot:
    p <- ggplot(plotDF, aes(x=sample, y=percentage))
    p <- p + geom_bar(stat="identity", aes(fill=type))
    p <- p + theme(axis.text.x = element_text(angle = 90), legend.position = "bottom")
    pdf(file = paste0(plotDir, "/", Type, "_", sampleName, "_repeat_percent_compBarplot", ".pdf"), height = 7, width = 7)
    print(p)
    dev.off()
  }
}
 