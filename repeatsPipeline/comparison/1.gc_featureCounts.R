### 0. Load in arguments ###
args = commandArgs(trailingOnly = TRUE)
for (n in 1:2) {
  print(args[n])
}

if (!is.null(args[1])) {
  uID <- args[1]
}
uID <- "FT1"

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

expName <- "comparison"

# define directories:
#homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ClusterShare/thingamajigs/jamtor/"
homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project, "/RNA-seq/")
RobjectDir <- paste0(projectDir, "/Robjects/", expName, "/")
bamDir <- paste0(projectDir, "/results/star/GC", uID, "/")
genomeDir <- "/share/ScratchGeneral/jamtor/genomes/"

inFile <- paste0(projectDir, "/results/star/GC/", expName, "/", uID, "/Aligned.sortedByCoord.out.bam")


### 2. Count overlaps of bams with gencode using featureCounts ###

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

saveRDS(GCcounts, file=paste0(RobjectDir, "/", uID, "_", GCcounts.rds"))
