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
expName <- "exp2"

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "projects/", project, "/")
RobjectDir <- paste0(projectDir, "RNA-seq/Robjects/")
resultsDir <- paste0(projectDir, "/RNA-seq/results")
plotDir <- paste0(resultsDir, "/R/", expName, "/plots/compBarplot/")

# create plotDir:
system(paste0("mkdir -p ", plotDir))
print(paste0("The outDir is: ", plotDir))


### 1. Load in counts lists for each sample ###

allCounts <- grep("subset", list.files(paste0(RobjectDir, "/", expName, "/compCounts"), pattern="_compCounts.rds"), value=T, invert=T)

for (i in 1:length(allCounts)) {
  sampleName <- gsub("_compCounts.rds", "", allCounts[i])
  if (i==1) {
    Counts <- list(readRDS(file=paste0(RobjectDir, expName, "/compCounts/", allCounts[i])))
    names(Counts)[i] <- sampleName
  } else {
    Counts[[i]] <- readRDS(file=paste0(RobjectDir, expName, "/compCounts/", allCounts[i]))
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
  Plots[[i]] <- p + theme(text = element_text(size=20), axis.text.x = element_text(angle = 90))
  i=i+1
}


pdf(file = paste0(plotDir, "compBarplotCounts.pdf"), height=20, width=35)
Plots[[1]]
dev.off()

pdf(file = paste0(plotDir, "compBarplotPercent.pdf"), height=20, width=35)
Plots[[2]]
dev.off()

# load or create libSize:
if (file.exists(paste0(RobjectDir, "/", expName, "/", "/lib", "Sizes.rds"))) {
  lSizes <- readRDS(paste0(RobjectDir, "/", expName, "/", "/lib", "Sizes.rds"))
} else {
  # load in library sizes:
  libFiles <- grep("subset", list.files(paste0(RobjectDir, "/", expName, "/libSize"), pattern="\\_libSize", full.names = T), value = T, invert = T)
  for (i in 1:length(libFiles)) {
    if (i==1) {
      lSizes <- c(readRDS(file=libFiles[i]))
    } else {
      lSizes[i] <- readRDS(file=libFiles[i])
    }
    print(i)
    i=i+1
  }
  saveRDS(lSizes, file=paste0(RobjectDir, "/", expName, "/", "/libSize.rds"))
}

# convert counts to CPMs:
CPMdf <- countsDF
i=1
for (s in libSize) {
  CPMdf[,i] <- (countsDF[,i]/libSize[[i]])*1000000
  i=i+1
}

# melt dataframe:
pCPM <- melt(CPMdf, variable.name = "sample")
pCPM$gene_type <- factor(pCPM$gene_type, levels = c("non_coding", "ribosome", "other", "repeats", "protein_coding"))

# plot data as barplot:
p <- ggplot(pCPM, aes(x=sample, y=value))
p <- p + geom_bar(stat="identity", aes(fill=gene_type))
p <- p + theme(text = element_text(size=20), axis.text.x = element_text(angle = 90))
pdf(file = paste0(plotDir, "compBarplotCPM.pdf"), height=20, width=35)
p
dev.off()

