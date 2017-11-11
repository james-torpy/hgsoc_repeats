# This R script counts overlaps of multiple RNA-seq data sets with a GTF file of transcripts
# and outputs a counts .csv file for each sample

### 0. Set up variables and directories ###

# load packages needed:
library(GenomicRanges)
library(ShortRead)
library(rtracklayer)
library("BSgenome.Hsapiens.UCSC.hg38")
library(ggplot2)

# define starting variables:
project <- "hgsoc_repeats"
methodName <- "method1"
refName <- "human-89.repeats.tab"
sampleName <- "subset"
gcName <- "gencode_v24_hg38_annotation.gtf"

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project)
bamDir <- paste0(projectDir, "/RNA-seq/results/star/", methodName)
refDir <- paste0(projectDir, "/RNA-seq/refs/")
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/")
outDir <- paste0(resultsDir, "/R/", methodName, "/")
bedDir <- paste0(outDir, "bed")
genomeDir <- paste0(homeDir, "genomes/hg38_ercc/")

# create outDir:
#system(paste0("mkdir -p ", outDir))
#print(paste0("The outDir is: ", outDir))
#system(paste0("mkdir -p ", bedDir))


### 1. Load in repeats annotation file ###

# define gtf variable:
ref <- paste0(refDir, refName)
print(paste0("The ref file is: ", ref))

# load in the gene coordinates from the tab file:
genes <- read.table(file = ref, sep = "\t", as.is = T)
saveRDS(genes, file = paste0(RobjectDir, "genesTemp.rds"))
#genes <- loadRDS(file = paste0(RobjectDir, "genesTemp.rds"))

# format and convert to GRanges object:
genes <- genes[-1,]
# format the annotation df so it only includes chromosomes 1-22, X, Y:
genes <- genes[grep("\\.|MT", genes$V1, invert = T),]
genes$V2 <- as.numeric(genes$V2)
genes$V3 <- as.numeric(genes$V3)
genesGR <- GRanges(
  seqnames = genes$V1,
  ranges = IRanges(start=genes$V2, end=genes$V3),
  strand = genes$V4,
  score = genes$V7
)
values(genesGR) <- data.frame(genes$V5, genes$V6, genes$V7, genes$V8, genes$V9)
colnames(values(genesGR)) <- c("type", "analysis", "score", "class", "class_desc")

# save annotation as RDF file:
saveRDS(genesGR, file=paste0(RobjectDir, "repeatsGR.rds"))
#genesGR <- readRDS(file=paste0(RobjectDir, "repeatsGR.rds"))

# convert GR lists into bed files to check overlaps using IGV:
#i=1
#GR2bed <- function(x, reduced = "reduced", Names = gsub(" ", "_", gsub("/", "", names(scGenes))), Class = "sc") {
#  df <- data.frame(seqnames=seqnames(x),
#    starts=as.integer(start(x)-1),
#    ends=end(x),
#    names=c(rep(".", length(x))),
#    score=c(rep(".", length(x))),
#    strands=strand(x))
#  if (reduced == "reduced") {
#    write.table(df, file=paste0(bedDir, "/", Names[i], "_reduced_", Class, ".bed"), quote=F, sep="\t", row.names=F, col.names=F)
#  } else {
#    write.table(df, file=paste0(bedDir, "/", Names[i], "_", Class, ".bed"), quote=F, sep="\t", row.names=F, col.names=F)
#  }
#  i <<- i+1
#  return(df)
#}

#i=1
#scGenesBed <- lapply(scGenes, GR2bed, reduced="no")
#i=1
#scGenesReducedBed <- lapply(scGenesReduced, GR2bed, reduced="reduced")
#i=1
#cGenesBed <- lapply(scGenes, GR2bed, reduced="no", Names = gsub("/", "", names(cGenes)), Class = "c")
#i=1
#cGenesReducedBed <- lapply(scGenesReduced, GR2bed, reduced="reduced", Names = gsub("/", "", names(cGenes)), Class = "c")


### 2. Set up lists ###

# Each GRangesList corresponds to a plot #
# Each GRanges object corresponds to a line on a plot #


### 2a. Set up superclass GRangesList ###

# modify class names to group as superclasses:
genesGR$superclass <- gsub("\\?.*", "", gsub("/.*", "", genesGR$class))
# Split genesGR by superclasses and reduce to avoid overlaps between superclasses:
scGenes <- endoapply(split(genesGR, genesGR$superclass), reduce)


### 2b. Set up class list of GRangesLists ###

# convert class column of genesGR from factors to characters
genesGR$class <- as.character(genesGR$class)

# split genesGR by superclass:
scTemp <- split(genesGR, genesGR$superclass)
# make a list of GRanges superclass objects:
scTempList <- list()
i=1
scTempList <- lapply(scTemp, function(x) {
  scTempList[[i]] <- x
  return(scTempList)
  i <<- i+1
})

# split each GRanges element of scTempList into GRanges lists by class:
cGenes <- lapply(scTempList, function(x) {
  spl <- split(x[[1]], x[[1]]$class)
  return(reduce(spl))
})


### 2b. Set up type list of GRangesLists ###

# split genesGR by class:
cTemp <- split(genesGR, genesGR$class)
# make a list of GRanges superclass objects:
cTempList <- list()
i=1
cTempList <- lapply(cTemp, function(x) {
  cTempList[[i]] <- x
  return(cTempList)
  i <<- i+1
})

# split each GRanges element of cTempList into GRanges lists by type:
tGenes <- lapply(cTempList, function(x) {
  spl <- split(x[[1]], x[[1]]$type)
  return(reduce(spl))
})

create_rds <- function(x, rds) {
  if (file.exists(rds)) {
    print(paste0(rds, " exists, no need to create"))
  } else {
    print(paste0("Creating ", rds))
    saveRDS(x, file = rds)
  }
}

for (e in list(scGenes, cGenes, tGenes)) {
  create_rds(e, paste0(RobjectDir, "/", deparse(substitute(e)), "_RepeatGenes.rds"))
}

scGenes <- readRDS(file = paste0(RobjectDir, "/sc_RepeatGenes.rds"))
cGenes <- readRDS(file = paste0(RobjectDir, "/c_RepeatGenes.rds"))
tGenes <- readRDS(file = paste0(RobjectDir, "/t_RepeatGenes.rds"))


### 2b. Append gencode annotation to each GRanges object:

# load in gencode as GRanges object:
GC <- import(paste0(genomeDir, gcName))
GCtemp <- GC[grep("\\.|M", seqnames(GC), invert = T),]
GCtemp2 <- GCtemp[grep("gene", GCtemp$type),]
GCgenes <- GRanges(seqnames = seqnames(GCtemp2),
                   ranges = ranges(GCtemp2),
                   strand = strand(GCtemp2),
                   type = rep("gencode", length(GCtemp2))
)

# append gencode annotation to each GRanges object in repeat annotation lists:
scGCgenes <- endoapply(scGenes, function(x) {
  x$type <- "repeat"
  return(c(x, GCgenes))
})
cGCgenes <- lapply(cGenes, function(x) {
  return(endoapply(x, function(x) {
    x$type <- "repeat"
    return(c(x, GCgenes))
  }))
})


### 3. Calculate gene widths and save as RDS file for FPKM generation:
# sum the width of all entires for each gene_id, providing the total length of instances of the gene in the genome:
#width_it <- function(x) {
#  return(as.data.frame(sum(width(reduce(x)))))
#}

#scWidths <- sapply(scGenes, width_it)

#cWidths <- lapply(cGenes, function(x) {
#  sapply(x, width_it)
#})

#tWidths <- lapply(tGenes, function(x) {
#  sapply(x, width_it)
#})
# save the widths as RDS file:
#saveRDS(scWidths, file = paste0(RobjectDir, "/sc_", sampleName, "RepeatWidths.rds"))
#saveRDS(cWidths, file = paste0(RobjectDir, "/c_", sampleName, "RepeatWidths.rds"))
#saveRDS(cWidths, file = paste0(RobjectDir, "/t_", sampleName, "RepeatWidths.rds"))
#cWidths <- readRDS(file = paste0(RobjectDir, "/sc_RepeatWidths.rds"))
#cWidths <- readRDS(file = paste0(RobjectDir, "/c_RepeatWidths.rds"))
#tWidths <- readRDS(file = paste0(RobjectDir, "/t_RepeatWidths.rds"))


### 4. Load in bam files ###

# load in the bamFiles as list of GRanges objects:
# fetch vector of bam files:
inFiles <- list.files(bamDir, recursive = T, pattern = "Aligned.sortedByCoord.out.bam", full.names = T)
inFiles <- grep(".bai", inFiles, invert = T, value = T)
if (sampleName == "subset") {
  inFiles <- grep("subset", inFiles, value = T)
} else {
  inFiles <- grep("subset", inFiles, invert = T, value = T)
}
inFiles <- grep("ribosome", inFiles, invert = T, value = T)

# define the sampleNames:
sampleNames <- strsplit(inFiles, "/")
sampleNames <- lapply(sampleNames, function(x) {
  return(x[12])
})
sampleNames <- unlist(sampleNames)

print("The samples to be processed are: ")
print(sampleNames)

# assign column names for dataframe:
what <- c("qname","rname","strand","pos","qwidth")
# flag unmapped sequences to leave out:
flag <- scanBamFlag(isUnmappedQuery=FALSE)
# define parameters of bam scan:
param <- ScanBamParam(what=what, flag=flag)

# create empty list for bamFiles:
bamFiles <- list()
i=1
for (s in sampleNames) {
  # index bams if necessary:
  writeLines("\n")
  if (file.exists(paste0(bamDir, "/", s, "/Aligned.sortedByCoord.out.bam.bai"))) {
    print(paste0("No need to index, ", bamDir, "/", s, "/Aligned.sortedByCoord.out.bam has already been indexed"))
  } else {
    print(paste0("bowIndexing ", bamDir, "/", s, "/Aligned.sortedByCoord.out.bam"))
    indexBam(paste0(bamDir, "/", s, "/Aligned.sortedByCoord.out.bam"))
  }
  # save complete filename of bam for loading:
  file <- paste0(bamDir, "/", s, "/Aligned.sortedByCoord.out.bam")
  
  # load in bams:
  print(paste0("Loading ", file, " into bams list."))
  bamFiles[[i]] <- scanBam(file ,param=param)
  
  i <<- i+1
}

# create data frame of human chromosome lengths:
seq_lengths <- seqlengths(Hsapiens)

# convert the bams to GRanges objects:
bamGRs <- lapply(bamFiles, function(x) {
  bamGR <- GRanges(
    seqnames = x[[1]]$rname,
    ranges = IRanges(start = x[[1]]$pos, width = x[[1]]$qwidth),
    strand = x[[1]]$strand,
    seqlengths = seq_lengths,
    qnames = x[[1]]$qname)
  return(bamGR)
})
names(bamGRs) <- sampleNames


#create_rds(bamGRs, paste0(RobjectDir, "/bam", sampleName, "GRs.rds"))
if (!file.exists(paste0(RobjectDir, "/bam", sampleName, "GRs.rds"))) {
  readRDS(file = paste0(RobjectDir, "/bam", sampleName, "GRs.rds"))
}


# determine library size of each data set:
#libSizes <- lapply(bamGRs, function(x) {
#  return(length(seqnames(x)))
#3})

# save the libSizes as RDS file:
#saveRDS(libSizes, file = paste0(RobjectDir, "/lib", sampleName, "Sizes.rds"))
lSizes <- readRDS(file = paste0(RobjectDir, "/lib", sampleName, "Sizes.rds"))


# This R script counts overlaps of multiple RNA-seq data sets with a GTF file of transcripts
# and outputs a counts .csv file for each sample

### 0. Set up variables and directories ###

# load packages needed:
library(GenomicRanges)
library(ShortRead)
library(rtracklayer)
library("BSgenome.Hsapiens.UCSC.hg38")
library(ggplot2)

# define starting variables:
project <- "hgsoc_repeats"
methodName <- "method1"
refName <- "human-89.repeats.tab"
sampleName <- "subset"
gcName <- "gencode_v24_hg38_annotation.gtf"

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project)
bamDir <- paste0(projectDir, "/RNA-seq/results/star/", methodName)
refDir <- paste0(projectDir, "/RNA-seq/refs/")
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/")
outDir <- paste0(resultsDir, "/R/", methodName, "/")
bedDir <- paste0(outDir, "bed")
genomeDir <- paste0(homeDir, "genomes/hg38_ercc/")

# create outDir:
#system(paste0("mkdir -p ", outDir))
#print(paste0("The outDir is: ", outDir))
#system(paste0("mkdir -p ", bedDir))


### 1. Load in repeats annotation file ###

# define gtf variable:
ref <- paste0(refDir, refName)
print(paste0("The ref file is: ", ref))

# load in the gene coordinates from the tab file:
genes <- read.table(file = ref, sep = "\t", as.is = T)
saveRDS(genes, file = paste0(RobjectDir, "genesTemp.rds"))
#genes <- loadRDS(file = paste0(RobjectDir, "genesTemp.rds"))

# format and convert to GRanges object:
genes <- genes[-1,]
# format the annotation df so it only includes chromosomes 1-22, X, Y:
genes <- genes[grep("\\.|MT", genes$V1, invert = T),]
genes$V2 <- as.numeric(genes$V2)
genes$V3 <- as.numeric(genes$V3)
genesGR <- GRanges(
  seqnames = genes$V1,
  ranges = IRanges(start=genes$V2, end=genes$V3),
  strand = genes$V4,
  score = genes$V7
)
values(genesGR) <- data.frame(genes$V5, genes$V6, genes$V7, genes$V8, genes$V9)
colnames(values(genesGR)) <- c("type", "analysis", "score", "class", "class_desc")

# save annotation as RDF file:
saveRDS(genesGR, file=paste0(RobjectDir, "repeatsGR.rds"))
#genesGR <- readRDS(file=paste0(RobjectDir, "repeatsGR.rds"))

# convert GR lists into bed files to check overlaps using IGV:
#i=1
#GR2bed <- function(x, reduced = "reduced", Names = gsub(" ", "_", gsub("/", "", names(scGenes))), Class = "sc") {
#  df <- data.frame(seqnames=seqnames(x),
#    starts=as.integer(start(x)-1),
#    ends=end(x),
#    names=c(rep(".", length(x))),
#    score=c(rep(".", length(x))),
#    strands=strand(x))
#  if (reduced == "reduced") {
#    write.table(df, file=paste0(bedDir, "/", Names[i], "_reduced_", Class, ".bed"), quote=F, sep="\t", row.names=F, col.names=F)
#  } else {
#    write.table(df, file=paste0(bedDir, "/", Names[i], "_", Class, ".bed"), quote=F, sep="\t", row.names=F, col.names=F)
#  }
#  i <<- i+1
#  return(df)
#}

#i=1
#scGenesBed <- lapply(scGenes, GR2bed, reduced="no")
#i=1
#scGenesReducedBed <- lapply(scGenesReduced, GR2bed, reduced="reduced")
#i=1
#cGenesBed <- lapply(scGenes, GR2bed, reduced="no", Names = gsub("/", "", names(cGenes)), Class = "c")
#i=1
#cGenesReducedBed <- lapply(scGenesReduced, GR2bed, reduced="reduced", Names = gsub("/", "", names(cGenes)), Class = "c")


### 2. Set up lists ###

# Each GRangesList corresponds to a plot #
# Each GRanges object corresponds to a line on a plot #


### 2a. Set up superclass GRangesList ###

# modify class names to group as superclasses:
genesGR$superclass <- gsub("\\?.*", "", gsub("/.*", "", genesGR$class))
# Split genesGR by superclasses and reduce to avoid overlaps between superclasses:
scGenes <- endoapply(split(genesGR, genesGR$superclass), reduce)


### 2b. Set up class list of GRangesLists ###

# convert class column of genesGR from factors to characters
genesGR$class <- as.character(genesGR$class)

# split genesGR by superclass:
scTemp <- split(genesGR, genesGR$superclass)
# make a list of GRanges superclass objects:
scTempList <- list()
i=1
scTempList <- lapply(scTemp, function(x) {
  scTempList[[i]] <- x
  return(scTempList)
  i <<- i+1
})

# split each GRanges element of scTempList into GRanges lists by class:
cGenes <- lapply(scTempList, function(x) {
  spl <- split(x[[1]], x[[1]]$class)
  return(reduce(spl))
})


### 2b. Set up type list of GRangesLists ###

# split genesGR by class:
cTemp <- split(genesGR, genesGR$class)
# make a list of GRanges superclass objects:
cTempList <- list()
i=1
cTempList <- lapply(cTemp, function(x) {
  cTempList[[i]] <- x
  return(cTempList)
  i <<- i+1
})

# split each GRanges element of cTempList into GRanges lists by type:
tGenes <- lapply(cTempList, function(x) {
  spl <- split(x[[1]], x[[1]]$type)
  return(reduce(spl))
})

createORload_rds <- function(x, rds) {
  if (file.exists(rds)) {
    print(paste0(rds, " exists, loading"))
    return(readRDS(rds))
  } else {
    print(paste0("Creating ", rds))
    saveRDS(x, file = rds)
  }
}

for (e in list(scGenes, cGenes, tGenes)) {
  create_rds(e, paste0(RobjectDir, "/", deparse(substitute(e)), "_RepeatGenes.rds"))
}

scGenes <- readRDS(file = paste0(RobjectDir, "/sc_RepeatGenes.rds"))
cGenes <- readRDS(file = paste0(RobjectDir, "/c_RepeatGenes.rds"))
tGenes <- readRDS(file = paste0(RobjectDir, "/t_RepeatGenes.rds"))


### 2b. Append gencode annotation to each GRanges object:

# load in gencode as GRanges object:
GC <- import(paste0(genomeDir, gcName))
GCtemp <- GC[grep("\\.|M", seqnames(GC), invert = T),]
GCtemp2 <- GCtemp[grep("gene", GCtemp$type),]
GCgenes <- GRanges(seqnames = seqnames(GCtemp2),
                   ranges = ranges(GCtemp2),
                   strand = strand(GCtemp2),
                   type = rep("gencode", length(GCtemp2))
)

if (file.exists(paste0(RobjectDir, "GCgenes.rds"))) {
  GCgenes <- readRDS(GCgenes)
} else {
  saveRDS(GCgenes, file = paste0(RobjectDir, "GCgenes.rds"))
}

# append gencode annotation to each GRanges object in repeat annotation lists:
scGCgenes <- endoapply(scGenes, function(x) {
  x$type <- "repeat"
  return(c(x, GCgenes))
})

cGCgenes <- lapply(cGenes, function(x) {
  return(endoapply(x, function(x) {
    x$type <- "repeat"
    return(c(x, GCgenes))
  }))
})


### 3. Calculate gene widths and save as RDS file for FPKM generation:
# sum the width of all entires for each gene_id, providing the total length of instances of the gene in the genome:
#width_it <- function(x) {
#  return(as.data.frame(sum(width(reduce(x)))))
#}

#scWidths <- sapply(scGenes, width_it)

#cWidths <- lapply(cGenes, function(x) {
#  sapply(x, width_it)
#})

#tWidths <- lapply(tGenes, function(x) {
#  sapply(x, width_it)
#})
# save the widths as RDS file:
#saveRDS(scWidths, file = paste0(RobjectDir, "/sc_", sampleName, "RepeatWidths.rds"))
#saveRDS(cWidths, file = paste0(RobjectDir, "/c_", sampleName, "RepeatWidths.rds"))
#saveRDS(cWidths, file = paste0(RobjectDir, "/t_", sampleName, "RepeatWidths.rds"))
#cWidths <- readRDS(file = paste0(RobjectDir, "/sc_RepeatWidths.rds"))
#cWidths <- readRDS(file = paste0(RobjectDir, "/c_RepeatWidths.rds"))
#tWidths <- readRDS(file = paste0(RobjectDir, "/t_RepeatWidths.rds"))


### 4. Load in bam files ###

# load in the bamFiles as list of GRanges objects:
# fetch vector of bam files:
inFiles <- list.files(bamDir, recursive = T, pattern = "Aligned.sortedByCoord.out.bam", full.names = T)
inFiles <- grep(".bai", inFiles, invert = T, value = T)
if (sampleName == "subset") {
  inFiles <- grep("subset", inFiles, value = T)
} else {
  inFiles <- grep("subset", inFiles, invert = T, value = T)
}
inFiles <- grep("ribosome", inFiles, invert = T, value = T)

# define the sampleNames:
sampleNames <- strsplit(inFiles, "/")
sampleNames <- lapply(sampleNames, function(x) {
  return(x[12])
})
sampleNames <- unlist(sampleNames)

print("The samples to be processed are: ")
print(sampleNames)

# assign column names for dataframe:
what <- c("qname","rname","strand","pos","qwidth")
# flag unmapped sequences to leave out:
flag <- scanBamFlag(isUnmappedQuery=FALSE)
# define parameters of bam scan:
param <- ScanBamParam(what=what, flag=flag)

# create empty list for bamFiles:
bamFiles <- list()
i=1
for (s in sampleNames) {
  # index bams if necessary:
  writeLines("\n")
  if (file.exists(paste0(bamDir, "/", s, "/Aligned.sortedByCoord.out.bam.bai"))) {
    print(paste0("No need to index, ", bamDir, "/", s, "/Aligned.sortedByCoord.out.bam has already been indexed"))
  } else {
    print(paste0("bowIndexing ", bamDir, "/", s, "/Aligned.sortedByCoord.out.bam"))
    indexBam(paste0(bamDir, "/", s, "/Aligned.sortedByCoord.out.bam"))
  }
  # save complete filename of bam for loading:
  file <- paste0(bamDir, "/", s, "/Aligned.sortedByCoord.out.bam")
  
  # load in bams:
  print(paste0("Loading ", file, " into bams list."))
  bamFiles[[i]] <- scanBam(file ,param=param)
  
  i <<- i+1
}

# create data frame of human chromosome lengths:
seq_lengths <- seqlengths(Hsapiens)

# convert the bams to GRanges objects:
bamGRs <- lapply(bamFiles, function(x) {
  bamGR <- GRanges(
    seqnames = x[[1]]$rname,
    ranges = IRanges(start = x[[1]]$pos, width = x[[1]]$qwidth),
    strand = x[[1]]$strand,
    seqlengths = seq_lengths,
    qnames = x[[1]]$qname)
  return(bamGR)
})
names(bamGRs) <- sampleNames

#createORload_rds(bamGRs, paste0(RobjectDir, "/bam", sampleName, "GRs.rds"))
saveRDS(bamGRs, paste0(RobjectDir, "/bam", sampleName, "GRs.rds"))

# determine library size of each data set:
libSizes <- lapply(bamGRs, function(x) {
  return(length(seqnames(x)))
})

# save or load the libSizes as RDS file:
#saveRDS(libSizes, file = paste0(RobjectDir, "/lib", sampleName, "Sizes.rds"))
lSizes <- readRDS(file = paste0(RobjectDir, "/lib", sampleName, "Sizes.rds"))
#lSizes <- createORload_rds(lSizes, paste0(RobjectDir, "/lib", sampleName, "Sizes.rds"))


### 4. Count overlaps of bams with all repeat annotation GRanges objects:

# count overlaps of bams with annotation ranges:
count_it <- function(x, annot = scGCgenes[[1]]) {
  writeLines("\n")
  print("Counting overlaps...")
  counts <- as.data.frame(countOverlaps(annot, x))
  return(sum(counts))
}

scGCcounts <- lapply(scGCgenes, function(y) {
  return(lapply(bamGRs, count_it, annot = y))
})
scGCcountsDF <- as.data.frame(do.call("rbind", scGCcounts))

Names <- names(bamGRs)
cGCcounts <- lapply(cGCgenes, function(x) {
  return(lapply(bamGRs, count_it, annot = x))
})
cGCcountsDF <- lapply(cGCcounts, function(x) {
  counts <- as.data.frame(do.call("cbind", x))
  colnames(counts) <- Names
  return(counts)
})

tCounts <- lapply(tGenes, function(x) {
  return(lapply(bamGRs, count_it, annot = x))
})
tCountsDF <- lapply(tCounts, function(x) {
  counts <- as.data.frame(do.call("cbind", x))
  colnames(counts) <- Names
  return(counts)
})

# save the counts as RDS files:
if (file.exists(paste0(RobjectDir, "/scGC_", sampleName, "RepeatCountsDF.rds"))) {
  print(paste0(RobjectDir, "/scGC_", sampleName, "RepeatCountsDF.rds already exists, no need to create"))
} else {
  print(paste0("Creating", RobjectDir, "/scGC_", sampleName, "RepeatCountsDF.rds"))
  saveRDS(scGCcountsDF, file = paste0(RobjectDir, "/scGC_", sampleName, "RepeatCountsDF.rds"))
}

if (file.exists(paste0(RobjectDir, "/cGC_", sampleName, "RepeatCountsDF.rds already exists, no need to create"))) {
  print(paste0(RobjectDir, "/cGC_", sampleName, "RepeatCountsDF.rds"))
} else {
  print(paste0("Creating", RobjectDir, "/cGC_", sampleName, "RepeatCountsDF.rds"))
  saveRDS(cGCcountsDF, file = paste0(RobjectDir, "/cGC_", sampleName, "RepeatCountsDF.rds"))
}

if (file.exists(paste0(RobjectDir, "/t_", sampleName, "RepeatCounts.rds already exists, no need to create"))) {
  print(paste0(RobjectDir, "/t_", sampleName, "RepeatCounts.rds"))
} else {
  print(paste0("Creating", RobjectDir, "/t_", sampleName, "RepeatCounts.rds"))
  saveRDS(tCounts, file = paste0(RobjectDir, "/t_", sampleName, "RepeatCounts.rds"))
}
