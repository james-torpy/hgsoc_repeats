library(GenomicRanges)
library(ShortRead)
library(rtracklayer)
library("BSgenome.Hsapiens.UCSC.hg38")

project <- "hgsoc_repeats"
gcName <- "gencode_v24_hg38_annotation.gtf"
annot <- "c"

### 1. Define directories/variables ###

homeDir <- "/Users/jamestorpy/clusterHome/"
projectDir <- paste0(homeDir, "/projects/hgsoc_repeats/RNA-seq/")
RobjectDir <- paste0(projectDir, "/Robjects/exp1/")
bamDir <- paste0(projectDir, "/results/star/prelim/cluster/method1/")

inFile <- paste0(projectDir, "/results/star/prelim/cluster/method1/bowtell_FT1_subset/Aligned.sortedByCoord.out.bam")
uID <- "bowtell_FT1_subset"


### 2. Load in repeats annotations ###

Genes <- readRDS(file = paste0(RobjectDir, "/", annot, "_RepeatGenes.rds"))

# split Genes into 2 and use first half:
splitInd <- c( rep(1, round(length(Genes)/2)), rep(2, (length(Genes) - round(length(Genes)/2))) )
Genes <- split(Genes, splitInd)[[1]]


### 3. Load in GC annotation ###

if (file.exists(paste0(RobjectDir, "GCgenes.rds"))) {
  GCgenes <- readRDS(paste0(RobjectDir, "GCgenes.rds"))
} else {
  GCgenes <- import(paste0(genomeDir, gcName))
  saveRDS(GCgenes, file = paste0(RobjectDir, "GCgenes.rds"))
}

GCgenes <- GCgenes[grep("\\.|M", seqnames(GCgenes), invert = T),]
GCgenes <- GCgenes[grep("exon", GCgenes$type),]
#GCgenes <- GCgenes[grep("\\+", as.character(strand(GCgenes)))]

# load in gencode as GRanges object:
if (!file.exists(paste0(outDir, "/", uID, "_GCcountsDF.rds"))) {
  GCcounts <- as.data.frame(countOverlaps(GCgenes, bamGR))
  GCcounts$gene_id <- gsub("\\.*", "", GCgenes$gene_id)
  GCcountsDF <- aggregate(.~gene_id, GCcounts, sum)
  
  saveRDS(GCcountsDF, file = paste0(outDir, "/", uID, "_GCcountsDF.rds"))
}


### 4. Load in bam file ###

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
  what <- c("qname","rname","strand","pos","qwidth")
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
    qnames = bam[[1]]$qname)
  saveRDS(bamGR, file = paste0(RobjectDir, uID, "_bamGR.rds"))
}


### 5. Subset LTR Gypsy from repeats annotation and the first 500 ranges
# of gencode chr 14 annotation ###

LTRGypsy <- unlist(Genes[[3]][[8]])
LTRGypsy14 <- LTRGypsy[seqnames(LTRGypsy)=="chr14"]
strand(LTRGypsy14) <- "*"
LTRGypsy14 <- reduce(LTRGypsy14)

LTRcounts <- countOverlaps(LTRGypsy14, bamGR)


GCgenesR <- reduce(GCgenes)
GC14 <- GCgenesR[seqnames(GCgenesR)=="chr14"]
GC14 <- GC14[1:500]

GCcounts <- countOverlaps(GC14, bamGR)





Helitron2 <- unlist(Genes[[5]][[2]])
strand(Helitron2) <- "*"
Helitron2 <- reduce(Helitron2)

Helitron2counts <- countOverlaps(Helitron2, bamGR)

GCgenesR <- reduce(GCgenes)

GCYR <- GCgenesR[seqnames(GCgenesR)=="chrY"]
GCYR <- GCYR[1:500]
GCcountsR <- countOverlaps(GCYR, bamGR)
findGCR <- findOverlaps(GCYR[77], bamGR)

GCY <- GCgenes[seqnames(GCgenes)=="chrY"]
GCY <- GCY[1:500]
GCcounts <- countOverlaps(GCY, bamGR)


GCselfoverlaps <- countOverlaps(GCgenes, GCgenes)

GCnstr <- GCgenesR
strand(GCnstr) <- "*"
countOverlaps(GCnstr, GCnstr)
