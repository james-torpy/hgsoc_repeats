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
gtfName <- "human-89.repeats4.gtf"

# read arguments from the command line:
args = commandArgs(trailingOnly = TRUE)
for (n in 1:6) {
  print(args(n))
}

if (!is.null(args[1])) {
  methodName = args[1]
}

if (!is.null(args[2])) {
  bamDir = args[2]
}

if (!is.null(args[3])) {
  RobjectDir = args[3]
  }

if (!is.null(args[4])) {
  outDir = args[4]
  }

if (!is.null(args[5])) {
  gtf = args[5]
  }

if (!is.null(args[4])) {
  inFiles = args[6]
  }

### 1. Load in repeats annotation file ###

# load in the gene coordinates from the gtf file:
genes <- import(gtf)


### 2. Add telomeres annotation to repeats gtf ###

# fetch telomere annotations from UCSC using rtracklayer:
# open a UCSC session for hg38:
session <- browserSession("UCSC")
genome(session) <- "hg38"
# define track to query - can fetch list of possible tracks using trackNames(session):
gapTrack <- ucscTableQuery(session, "Gap", GRangesForUCSCGenome("hg38"))
# define table to query and fetch it - can fetch list of possible tables using tableNames(gapTrack):
tableName(gapTrack) <- "gap"
gapTrack <- getTable(gapTrack)
# filter gapTrack for telomere entries using grep:
teloInd <- grep("telomere", gapTrack$type)
telomeres <- gapTrack[teloInd, ]
# add telomeres info to GRanges object and append to Granges object 'genes':
teloGR <- GRanges(
  seqnames = telomeres$chrom,
  ranges= IRanges(start = telomeres$chromStart, end = telomeres$chromEnd),
  strand = "*",
  source = rep("UCSC_Gaps", length(telomeres$chrom)),
  type = rep("genomic_loci", length(telomeres$chrom)),
  score = rep("<NA>", length(telomeres$chrom)),
  phase = rep("<NA>", length(telomeres$chrom)),
  frame = rep("<NA>", length(telomeres$chrom)),
  gene_id = rep("telomere", length(telomeres$chrom)),
  transcript_id = rep("<NA>", length(telomeres$chrom)),
  gene_type = rep("Simple repeats", length(telomeres$chrom)),
  gene_status = rep("KNOWN", length(telomeres$chrom)),
  gene_name = rep("telomere", length(telomeres$chrom)),
  transcript_type = rep("<NA>", length(telomeres$chrom)),
  transcript_status = rep("<NA>", length(telomeres$chrom)),
  transcript_name = rep("<NA>", length(telomeres$chrom)),
  exon_number = rep("<NA>", length(telomeres$chrom)),
  exon_id = rep("<NA>", length(telomeres$chrom)),
  level = rep(3, length(telomeres$chrom))
) 

fGenes <- append(genes, teloGR)

# save RData file:
save.image(file=paste0("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/", methodName, "repGTF_loaded.RData"))
methodName <- "method1"
#load(file=paste0("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/", methodName, "repGTF_loaded.RData"))


### 3. Load in bam files ###

# load in the bamFiles as list of GRanges objects:

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
    print(paste0("Indexing ", bamDir, "/", s, "/Aligned.sortedByCoord.out.bam"))
    indexBam(paste0(bamDir, "/", s, "/Aligned.sortedByCoord.out.bam"))
  }
  # save complete filename of bam for loading:
  file <- paste0(bamDir, "/", s, "/Aligned.sortedByCoord.out.bam")
  
  # load in bams:
  print(paste0("Loading ", file, " into bams list."))
  bamFiles[i] <- scanBam(file ,param=param)

  i <<- i+1
}
  
save.image(file=paste0("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/", methodName, "_repbam_loaded.RData"))
methodName <- "method1"
#load(file=paste0("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/", methodName, "_repbam_loaded.RData"))

# create data frame of human chromosome lengths:
seq_lengths <- seqlengths(Hsapiens)

# convert the bams to GRanges objects:
bamGRs <- lapply(bamFiles, function(x) {
  bamGR <- GRanges(
    seqnames = x$rname,
    ranges = IRanges(start = x$pos, width = x$qwidth),
    strand = x$strand,
    seqlengths = seq_lengths,
    qnames = x$qname)
  return(bamGR)
})
names(bamGRs) <- sampleNames


# determine library size of each data set:
libSizes <- lapply(bamGRs, function(x) {
  return(length(x@seqnames))
})

# save the libSizes as RDS file:
saveRDS(libSizes, file = paste0(RobjectDir, "/libSizes.rds"))


### 4. Count overlaps of bams with repeat annotation gtf:

# count overlaps of bams with annotation ranges:
counts <- lapply(bamGRs, function(x) {
  writeLines("\n")
  print("Counting overlaps...")
  overlaps <- countOverlaps(fGenes, x)
  counts <- as.data.frame(overlaps)
  counts$counts <- counts[ ,1]
  counts[ ,1] <- fGenes$gene_id
  colnames(counts) <- c("genes", "counts")
  print(paste0("Start of counts for sample:"))
  print(head(counts))
  return(counts)
  i <<- i+1
})
# name the count dataframes in the list:
names(counts) <- sampleNames

# convert counts list into dataframe with samples as rows
counts <- do.call("cbind", counts)
colnames(counts)[1] <- "repeat"
nameInd <- grep("genes", colnames(counts))
counts <- counts[ , -nameInd]

# save the counts as csv files:
if (file.exists(paste0(outDir, "/counts.csv"))) {
  print(paste0(outDir, "/counts.csv"))
} else {
  print(paste0("creating", outDir, "/counts.csv"))
  write.csv(counts, file = paste0(outDir, "/counts.csv"))
}

save.image(file=paste0("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/", methodName, "_countsOverlapped.RData"))
methodName <- "method1"
#load(file=paste0("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/", methodName, "_countsOverlapped.RData"))