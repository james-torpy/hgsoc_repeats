### loadBams.R ###

# This R script takes a references to a bam file from call_loadBams.bash,
# loads the bam in and saves as an RDS file

### 0. Load in arguments ###
#args = commandArgs(trailingOnly = TRUE)
#for (n in 1:2) {
#  print(args[n])
#}

#if (!is.null(args[1])) {
#  uID <- args[1]
#}
uID <- "bowtell_FT3_subset"

#if (!is.null(args[2])) {
#  annot <- args[2]
#}
#annot <- "c"

print(uID)
#print(annot)

# to test: qsub -q short.q -N paircoGC_$uID -b y -wd /share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/exp6a/logs -j y -R y -P GenomeInformatics -pe smp 12 -V "/home/jamtor/local/lib/r/R-3.2.2/bin/R --vanilla < /share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/exp6a/1a.paircountOverlaps_GC.R"

### 1. Set up variables and directories ###

library(GenomicRanges)
library(ShortRead)
library(rtracklayer)
library("BSgenome.Hsapiens.UCSC.hg38")
library(ggplot2)

# define starting variables:
project <- "hgsoc_repeats"
gcName <- "gencode_v24_hg38_annotation.gtf"
expName <- "exp6a"

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ClusterShare/thingamajigs/jamtor/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project, "/RNA-seq/")
RobjectDir <- paste0(projectDir, "/Robjects/", expName, "/")
resultsDir <- paste0(projectDir, "/results/")
bamDir <- paste0(resultsDir, "/star/GC/comparison/", uID, "/")
genomeDir <- "/share/ScratchGeneral/jamtor/genomes/"
refDir <- paste0(projectDir, "/refs/")

inFile <- paste0(bamDir, "/Aligned.sortedByCoord.out.bam")

system(paste0("mkdir ", RobjectDir))


### 2. Load in bam file and convert to GenomicRanges object if needed:

if (file.exists(paste0(RobjectDir, uID, "_bamGR.rds"))) {
  print(paste0("Loading ", RobjectDir, uID, "_bamGR.rds"))
  bamGR <- readRDS(file = paste0(RobjectDir, uID, "_bamGR.rds"))
} else {
  # index bam if needed:
  if (file.exists(paste0(bamDir, "/", uID, "/Aligned.sortedByCoord.out.bam.bai"))) {
    print(paste0("No need to index, ", inFile, " has already been indexed"))
  } else {
    print(paste0("Indexing ", inFile))
    indexBam(inFile)
  }
  
  # assign column names for dataframe:
  what <- c("qname","rname","strand","pos","qwidth", "flag")
  # flag unmapped sequences to leave out:
  flag <- scanBamFlag(isUnmappedQuery=FALSE)
  # define parameters of bam scan:
  param <- ScanBamParam(what=what, flag=flag)
  
  # load in bam:
  if (file.exists(paste0(RobjectDir, "/", uID, "_bamFile.rds"))) {
    print(paste0("Loading ", uID, "_bamFile.rds"))
    bam <- readRDS(paste0(RobjectDir, "/", uID, "_bamFile.rds"))
  } else {
    print(paste0("Loading ", inFile))
    bam <- scanBam(inFile ,param=param)
  }
  
  # create data frame of human chromosome lengths:
  seq_lengths <- seqlengths(Hsapiens)
  
  bamGR <- GRanges(
    seqnames = bam[[1]]$rname,
    ranges = IRanges(start = bam[[1]]$pos, width = bam[[1]]$qwidth),
    strand = bam[[1]]$strand,
    seqlengths = seq_lengths,
    qnames = bam[[1]]$qname,
    flags = bam[[1]]$flag)
  saveRDS(bamGR, file = paste0(RobjectDir, uID, "_bamGR.rds"))
}

# load in list of samflags specifying 1st read mates:
samflag <- read.table(file=paste0(refDir, "/first_mate_flags.txt"))[,1]

######

testGR <- bamGR[1:20]
testPairGR <- GRangesList
chek <- c()
i=1

testPairGR <- sapply(testGR, function(x) {
  pair <- c(x, testGR[grep(x$qnames, testGR$qnames),])
  print(paste0("For qname ", x$qnames))
  j=1
  for (j in 1:length(pair)) {
    if (pair[j]$flags %in% samflag) {
      print(paste0("Entry ", j, " is first mate - changing to ", strand(pair[j]), " strand"))
      strand(pair) <- strand(pair[j])
    }
    j=j+1
  }
  chek[i] <- x$qnames
  i <<- i+1
  return(pair)
})




######

testGR <- bamGR[1:20]

# convert bam GR to GR list:
i=1
testGRlist <- sapply(testGR, function(x) {
  if (i==1) {
    testGRlist <- GRangesList(x)
  } else {
    testGRlist[[i]] <- x
  }
  return(unlist(testGRlist))
  i <<- i+1
})

# grep for common read ids to group read pairs into individual GR objects, convert strand info to that of the first mate
# and create list from these:
i=1
testPairList <- lapply(testGRlist, function(x) {
  if (!x$qnames %in% chek) {
    pair <- testGR[grep(x$qnames, testGR$qnames),]
    print(paste0("For qname ", x$qnames))
    j=1
    for (j in 1:length(pair)) {
      if (pair[j]$flags %in% samflag) {
        print(paste0("Entry ", j, " is first mate - changing to ", strand(pair[j]), " strand"))
        strand(pair) <- strand(pair[j])
      }
    j=j+1
    }
  } else {
    pair <- NULL
  }
  if (i==1) {
    chek <- c(x$qnames)
  } else {
    chek[i] <<- x$qnames
  }
  i <<- i+1
  return(pair)
})

# remove NULL elements and convert into GR list:
testPairGR <- GRangesList(testPairList[!sapply(testPairList, is.null)])


######

bamPairGR <- GRangesList()
chek <- c()
i=1
for (l in bamGR$qnames) {
	if (!l %in% chek) {
		ind <- grep (l, bamGR$qnames)
		pair <- bamGR[ind]
	if (pair$flags[1] %in% samflag) {
		print("First entry is 1st mate - using this strand info")
		strand(pair) <- strand(pair)[1]
	} else if (pair$flags[2] %in% samflag) {
		print("Second entry is 1st mate - using this strand info")
		strand(pair) <- strand(pair)[2]
	} else {
		print("First mate not present")
	}
	bamPairGR[[i]] <- pair

	chek[i] <- l
	i=i+1
	}
}
saveRDS(bamPairGR, file=paste0(RobjectDir, "/bamPairGR.rds"))


### 3. Count overlaps of bams with gencode ###

# load in GCgenes:
if (file.exists(paste0(RobjectDir, "GCgenes.rds"))) {
  GCgenes <- readRDS(paste0(RobjectDir, "GCgenes.rds"))
} else {
  GCgenes <- import(paste0(genomeDir, gcName))
  saveRDS(GCgenes, file = paste0(RobjectDir, "GCgenes.rds"))
}
  
GCgenes <- GCgenes[grep("\\.|M", seqnames(GCgenes), invert = T),]
GCgenes <- GCgenes[grep("exon", GCgenes$type),]

# load in gencode as GRanges object:
if (!file.exists(paste0(outDir, "/", uID, "_GCcountsDFa.rds"))) {
  GCcounts <- as.data.frame(countOverlaps(GCgenes, bamPairGR))
  GCcounts$gene_id <- gsub("\\.*", "", GCgenes$gene_id)
  GCcountsDF <- aggregate(.~gene_id, GCcounts, sum)
  
  saveRDS(GCcountsDF, file = paste0(outDir, "/", uID, "_GCcountsDFa.rds"))
}

if (!file.exists(paste0(outDir, "/", uID, "_GCcountsDFext_a.rds"))) {
  GCcountsDFext <- as.data.frame(countOverlaps(GCgenes, bamPairGR))
  saveRDS(GCcountsDFext, file = paste0(outDir, "/", uID, "_GCcountsDFext_a.rds"))
}

save.image(file=paste0(RobjectDir, uID, "_paircountOverlaps.rds"))