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
gencodeName <- "gencode_v24_hg38_annotation.gtf"
ribName <- "human_rRNA.fa"
refName <- "human-89.repeats.tab"

# define directories:
#homeDir <- "/Users/jamestorpy/clusterHome/"
homeDir <- "/share/ScratchGeneral/jamtor/"
genomeDir <- paste0(homeDir, "genomes/")
gencodeDir <- paste0(genomeDir, "hg38_ercc/")
projectDir <- paste0(homeDir, "projects/", project, "/")
resultsDir <- paste0(projectDir, "RNA-seq/results/")
bamDir <- paste0(resultsDir, "star/", methodName)

RobjectDir <- paste0(projectDir, "RNA-seq/Robjects/")
ribDir <- paste0(resultsDir, "/star/", methodName, "/ribosome_12/")
plotDir <- paste0(resultsDir, "/R/", methodName, "/plots/compBarplot/")
refDir <- paste0(projectDir, "/RNA-seq/refs")

# create plotDir:
system(paste0("mkdir -p ", plotDir))
print(paste0("The outDir is: ", plotDir))


### 1. Load in bamGRs and put into a list ###
# load in the command line arguments:
args = commandArgs(trailingOnly = TRUE)
for (n in 1:length(args)) {
  print(args[n])
}

if (!is.null(args[1])) {
  sampleType <- args[1]
}
if (sampleType == "full") {
  sampleType <- ""
}
print(paste0("The sampleType is: ", sampleType))

bamGRfiles <- c()
j=1
for (i in 2:length(args)) {
  if (!is.null(args[i])) {
    bamGRfiles[j] <- args[i]
  }
  j=j+1
}
print(paste0("bamGRfile is: ", bamGRfiles))

sNames <- gsub("_bamGR.rds", "", bamGRfiles)

for (i in 1:length(bamGRfiles)) {
  if (i==1) {
    bamGRs <- list(readRDS(file=paste0(RobjectDir, "/", bamGRfiles[i])))
  } else {
    bamGRs[[i]] <- readRDS(file=paste0(RobjectDir, "/", bamGRfiles[i]))
  }
}


### 1. Load gencode gtf and split into protein_coding and non_coding/linc gene entries ###

# define gencode variable:
gencodeFile <- paste0(gencodeDir, gencodeName)
if (file.exists(paste0(RobjectDir, "GCgenes.rds"))) {
  print("Loading gencode...")
  gencode <- readRDS(file=paste0(RobjectDir, "gencode.RDS"))
} else {
  print(paste0("The gencode file is: ", gencodeFile))
  # load in the gene coordinates from the gencode file:
  gencode <- import(gencodeFile)
  # save gencode as RDF file:
  saveRDS(gencode, file=paste0(RobjectDir, "GCgenes.rds"))
}

# select only gene entries from gencode:
temp <- (gencode[gencode$type %in% "gene"])

# select relevant entries for non-coding or protein coding annotations and reduce to aggregate overlapping annotations:
pcGencode <- reduce(temp[temp$gene_type %in% "protein_coding"])
ncGencode <- reduce(temp[temp$gene_type %in% c("lincRNA", "non_coding")])
otherGencode <- reduce(temp[!temp$gene_type %in% c("protein_coding", "lincRNA", "non_coding")])



# determine library size of each data set:
if (file.exists(paste0(RobjectDir, "/lib", sampleType, "Sizes.rds"))) {
  print(paste0("Loading ", RobjectDir, "/lib", sampleType, "Sizes.rds"))
  libSizes <- readRDS(file = paste0(RobjectDir, "/lib", sampleType, "Sizes.rds"))
} else {
  libSizes <- lapply(bamGRs, function(x) {
    return(length(x@seqnames))
  })
  # save the libSizes as RDS file:
  saveRDS(libSizes, file = paste0(RobjectDir, "/lib", sampleType, "Sizes.rds"))
}

print(paste0("The libSizes are: ", libSizes))


### 3. Count overlaps between bams and gencode annotation:

# count overlaps and aggregate gene types protein_coding, non_coding and other:
if (file.exists(paste0(RobjectDir, "/CGtype", sampleType, "Counts.rds"))) {
  Counts <- readRDS(file=paste0(RobjectDir, "/CGtype", sampleType, "Counts.rds"))
} else {
  Counts <- lapply(bamGRs, function(x) {
    writeLines("\n")
    print("Counting overlaps...")
    result <- list(sum(countOverlaps(pcGencode, x)), sum(countOverlaps(ncGencode, x)), sum(countOverlaps(otherGencode, x)))
    names(result) <- c("protein_coding", "non_coding", "other")
    return(result)
  })
  saveRDS(Counts, file=paste0(RobjectDir, "/CGtype", sampleType, "Counts.rds"))
}


### 4. Count all repeats in samples ####

# create function to count overlaps of bams with annotation ranges:
count_it <- function(x, annot = scGenes[[1]]) {
  writeLines("\n")
  print("Counting overlaps...")
  counts <- as.data.frame(countOverlaps(annot, x))
  return(sum(counts))
}

# load in repeats annotation and reduce:
if (file.exists(paste0(RobjectDir, "/", sampleType, "TotalRepeatCounts.rds"))) {
  print(paste0("Loading ", RobjectDir, "/", sampleType, "TotalRepeatCounts.rds"))
  repeatCounts <- readRDS(file=paste0(RobjectDir, "/", sampleType, "TotalRepeatCounts.rds"))
} else {
  print(paste0("Loading", RobjectDir, "/repeatsGR.rds"))
  rGenes <- reduce(readRDS(file=paste0(RobjectDir, "/repeatsGR.rds")))
  print("Counting overlaps with total repeats annotation")
  repeatCounts <- lapply(bamGRs, count_it, rGenes)
  print(paste0("Generating ", RobjectDir, "/", sampleType, "TotalRepeatCounts.rds"))
  saveRDS(repeatCounts, file=paste0(RobjectDir, "/", sampleType, "TotalRepeatCounts.rds"))
}

i=1
for (e in repeatCounts) {
  Counts[[i]][4] <- e
  names(Counts[[i]])[4] <- "repeats"
  i=i+1
}


### 5. Fetch number of reads mapped to ribosomal RNA ###

riboLogs <- list.files(path = ribDir, pattern = "Log.final.out", recursive = T, full.names = T)
if (sampleType == "subset") {
  riboLogs <- grep("subset", riboLogs, value=T)
} else {
  riboLogs <- grep("subset", riboLogs, value=T, invert=T)
}

i=1
for (f in riboLogs) {
  tab <- read.table(f, sep = "\t", fill = T, as.is = T)
  Counts[[i]][5] <- as.numeric(tab[8,2]) + as.numeric(tab[23,2])
  names(Counts[[i]])[5] <- "ribosome"
  i=i+1
}


### 7. Create barplot ###

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

#melt dataframe:
pCPM <- melt(CPMdf, variable.name = "sample")
pCPM$gene_type <- factor(pCPM$gene_type, levels = c("non_coding", "ribosome", "other", "repeats", "protein_coding"))

# plot data as barplot:
p <- ggplot(pCPM, aes(x=sample, y=value))
p <- p + geom_bar(stat="identity", aes(fill=gene_type))
p <- p + theme(axis.text.x = element_text(angle = 90))
pdf(file = paste0(plotDir, sampleType, "compBarplotCPM.pdf"))
p
dev.off()

