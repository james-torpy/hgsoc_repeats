### loadBams.R ###

# This R script takes a references to a bam file from call_loadBams.bash,
# loads the bam in and saves as an RDS file

### 0. Load in arguments ###
args = commandArgs(trailingOnly = TRUE)
for (n in 1:2) {
  print(args[n])
}

if (!is.null(args[1])) {
  uID <- args[1]
}
uID <- "FT2"

if (!is.null(args[2])) {
  annot <- args[2]
}
annot <- "c"

print(uID)
print(annot)

### 1. Set up variables and directories ###

library(GenomicRanges)
library(ShortRead)
library(rtracklayer)
library(Rsubread)

# define starting variables:
project <- "hgsoc_repeats"
gcName <- "gencode_v24_hg38_annotation.gtf"
expName <- "exp3"

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ClusterShare/thingamajigs/jamtor/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project, "/RNA-seq/")
RobjectDir <- paste0(projectDir, "/Robjects/", expName, "/")
bamDir <- paste0(projectDir, "/results/star/GC/comparison/", uID, "/")
genomeDir <- "/share/ScratchGeneral/jamtor/genomes/"

inFile <- paste0(bamDir, "/Aligned.sortedByCoord.out.bam")


### 2. Load in repeats annotations ###

Genes <- readRDS(file = paste0(RobjectDir, "/", annot, "_RepeatGenes.rds"))

# split Genes into 2 and use first half:
splitInd <- c( rep(1, round(length(Genes)/2)), rep(2, (length(Genes) - round(length(Genes)/2))) )
Genes <- split(Genes, splitInd)[[1]]
Genes <- lapply(Genes, function(x) {
  strand(x) <- "*"
  return(reduce(x))
})


### 3. Count overlaps of bams with all repeat annotation GRanges objects using featureCounts:

for (i in 1:length(Genes)) {
  j=1
  if (i==1) {
    SAF <- list(lapply(Genes[[i]], function(x) {
      df <- data.frame(rep(names(Genes)[j], length(x)), seqnames(x), start(x), end(x), strand(x))
      names(df) <- c("GeneID", "Chr", "Start", "End", "Strand")
      j=j+1
      return(df)
    }))
  } else {
    SAF[[i]] <- lapply(Genes[[i]], function(x) {
      df <- data.frame(rep(names(Genes)[j], length(x)), seqnames(x), start(x), end(x), strand(x))
      names(df) <- c("GeneID", "Chr", "Start", "End", "Strand")
      j=j+1
      return(df)
    })
  }
}
names(SAF) <- names(Genes)

# split lists with length > 8:
#SAF[[length(SAF)+1]] <- SAF$DNA[5:8]
#SAF[[length(SAF)+1]] <- SAF$DNA[9:12]
#SAF[[length(SAF)+1]] <- SAF$DNA[13:16]
#SAF[[length(SAF)+1]] <- SAF$DNA[17:20]
#SAF[[length(SAF)+1]] <- SAF$DNA[21:22]
#SAF[[1]] <- SAF$DNA[1:4]
#names(SAF)[(length(SAF)-6):length(SAF)] <- c("DNA2", "DNA3", "DNA4")

#SAF[[length(SAF)+1]] <- SAF$LTR[6:10]
#SAF[[3]] <- SAF$LTR[1:5]
#names(SAF)[length(SAF)] <- "LTR2"

#i=1
#system.time(DNA_counts <- lapply(SAF[[1]], function(x) {
#  print(i)
#  print(paste0("Counting overlaps with repeat: ", names(SAF[[1]])[i]))
#  i <<- i+1
#  return(featureCounts(inFile, annot.ext=x))
#}))
#saveRDS(DNA_counts, file=paste0(RobjectDir, "/", uID, "_DNAfeatureCounts.rds"))


system.time(allResults <- for (i in 1:length(SAF)) {
  if (i==1) {
    j <<- 1

    SAFcounts <- list(lapply(SAF[[i]], function(x) {
      print(paste0("Counting overlaps with repeat: ", names(SAF[[i]])[j]))
      print(j)
      j <<- j+1
      return(featureCounts(inFile, annot.ext=x))
    }))
  } else {
    j <<- 1
    SAFcounts[[i]] <- lapply(SAF[[i]], function(x) {
      print(paste0("Counting repeat: ", names(SAF[[i]])[j]))
      print(j)
      j <<- j+1
      return(featureCounts(inFile, annot.ext=x))
    })
  }
  print(i)
})
saveRDS(allResults, file=paste0(RobjectDir, "/", uID, "_allFeatureResults.rds"))

allCounts <- lapply(Counts, function(x) {
  stat <- x$stat[1,]
  colnames(stat) <- c("Type", "Counts")
  return(stat)
})
saveRDS(allCounts, file=paste0(RobjectDir, "/", uID, "_allFeatureCounts.rds"))



  
  #assign(paste0(names(SAF)[i], "_counts"), lapply(SAF[[i]], function(x) {
  #  print(i)
  #  print(j)
  #  print(paste0("Variable is: ", names(SAF[i]), "_counts"))
  #  j=j+1
  #  return(featureCounts(inFile, annot.ext=x))
  #}))
#}


### 6. Count overlaps of bams with gencode ###

# load in GCgenes:
if (file.exists(paste0(RobjectDir, "GCgenes.rds"))) {
  GCgenes <- readRDS(paste0(RobjectDir, "GCgenes.rds"))
} else {
  GCgenes <- import(paste0(genomeDir, gcName))
  saveRDS(GCgenes, file = paste0(RobjectDir, "GCgenes.rds"))
}

GCgenes <- GCgenes[grep("\\.|M", seqnames(GCgenes), invert = T),]
GCgenes <- GCgenes[grep("exon", GCgenes$type),]

SAFgc <- data.frame(GCgenes$gene_id, seqnames(GCgenes), start(GCgenes), end(GCgenes), strand(GCgenes))
names(SAFgc) <- c("GeneID", "Chr", "Start", "End", "Strand")

system.time(GCcounts <- featureCounts(inFile, annot.ext=SAFgc))

saveRDS(GCcounts, file=paste0(RobjectDir, "/", uID, "_GCfeatureCounts.rds"))