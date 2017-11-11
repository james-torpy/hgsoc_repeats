# This script takes the following inputs:
# 1. gencode gtf
# 2. counts bam
# 3. repeat counts RDS
# 4. STAR Log.final.out from ribosomal mapping
# and creates barplots with the composition of protein-coding, non-coding, other, repeats and ribosomal counts
# with one bar per sample

### 0. Set up variables and directories ###

# load packages needed:
library(GenomicRanges)
library(ShortRead)
library("BSgenome.Hsapiens.UCSC.hg38")
library(reshape2)
library(ggplot2)

# define starting variables:
project <- "hgsoc_repeats"
methodName <- "method1"
sampleType <- ""

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "projects/", project, "/")
RobjectDir <- paste0(projectDir, "RNA-seq/Robjects/")
resultsDir <- paste0(projectDir, "/RNA-seq/results")
plotDir <- paste0(resultsDir, "/R/", methodName, "/plots/compBarplot/")

# create plotDir:
system(paste0("mkdir -p ", plotDir))
print(paste0("The outDir is: ", plotDir))


### 1. Load in counts lists for each sample ###

allCounts <- list.files(RobjectDir, pattern="_compCounts.rds")

for (i in 1:length(allCounts)) {
  sampleName <- gsub("_compCounts.rds", "", allCounts[i])
  if (i==1) {
    Counts <- list(readRDS(file=paste0(RobjectDir, allCounts[i])))
    names(Counts)[i] <- sampleName
  } else {
    Counts[[i]] <- readRDS(file=paste0(RobjectDir, allCounts[i]))
    names(Counts)[i] <- sampleName
  }
}


### 2. Create barplot ###

# convert Counts to a dataframe, convert each individual value from a list to an integer:
countsDF <- as.data.frame(sapply(as.data.frame(do.call("cbind", Counts)), unlist))
# add rownames as column and melt dataframe:
countsDF$gene_type <- rownames(countsDF)

# calculate percentages as new data frame:
perCountsDF <- as.data.frame(apply(countsDF[,1:ncol(countsDF)-1], 2, function(x) {
  return(as.integer(x)/as.integer(sum(x))*100)
}))
perCountsDF$gene_type <- countsDF$gene_type

# create composition barplots of CountsDF and perCountsDF:
cDFs <- list(countsDF, perCountsDF)
Plots <- list()
for (i in 1:2) {
  pCounts <- melt(cDFs[i], variable.name = "sample")
  pCounts$gene_type <- factor(pCounts$gene_type, levels = c("non_coding", "ribosome", "other", "repeats", "protein_coding"))
  
  # plot data as barplot:
  p <- ggplot(pCounts, aes(x=sample, y=value))
  p <- p + geom_bar(stat="identity", aes(fill=gene_type))
  Plots[[i]] <- p + theme(axis.text.x = element_text(angle = 90))
  i=i+1
}


pdf(file = paste0(plotDir, sampleType, "compBarplotCounts.pdf"))
Plots[[1]]
dev.off()

pdf(file = paste0(plotDir, sampleType, "compBarplotPercent.pdf"))
Plots[[2]]
dev.off()

# convert counts to CPMs:
CPMdf <- countsDF
i=1
for (s in libSizes) {
  CPMdf[,i] <- (countsDF[,i]/libSizes[[i]])*1000000
  i=i+1
}

# melt dataframe:
pCPM <- melt(CPMdf, variable.name = "sample")
pCPM$gene_type <- factor(pCPM$gene_type, levels = c("non_coding", "ribosome", "other", "repeats", "protein_coding"))

# plot data as barplot:
p <- ggplot(pCPM, aes(x=sample, y=value))
p <- p + geom_bar(stat="identity", aes(fill=gene_type))
p <- p + theme(axis.text.x = element_text(angle = 90))
pdf(file = paste0(plotDir, sampleType, "compBarplotCPM.pdf"))
p
dev.off()

