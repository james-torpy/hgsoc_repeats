# This R script counts overlaps of multiple RNA-seq data sets with a GTF file of transcripts

# load packages needed:
library(GenomicRanges)
library(ShortRead)
library(rtracklayer)
library("BSgenome.Hsapiens.UCSC.hg38")

# define starting variables:
project <- "hgsoc_repeats"
sampleName <- "L1subset3"
gtfName <- "human-89.repeats4.gtf"

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome"
projectDir <- paste0(homeDir, "/projects/", project)
inDir <- paste0(projectDir, "/RNA-seq/results/star/", sampleName)
refDir <- paste0(projectDir, "/RNA-seq/refs/")
resultsDir <- paste0(projectDir, "/RNA-seq/results")
outDir <- paste0(resultsDir, "/overlaps")

# create outDir:
system(paste0("mkdir -p ", outDir))
print(paste0("The outDir file is: ", outDir))

# define gtf variable:
gtf <- paste0(refDir, gtfName)
print(paste0("The gtf file is: ", gtf))

# load in the gene coordinates from the gtf file:
genes <- import(gtf)

# save RData file:
save.image(file="/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/repGTF_loaded.RData")
#load(file="/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/repGTF_loaded.RData")

# define the inFiles:
inFiles <- list.files(path = inDir, pattern = "Aligned.sortedByCoord.out.sorted.bam", recursive = T)
# or manually define inFiles:
inFile <- paste0(inDir, "/bowtell_FT1_subset/Aligned.sortedByCoord.out.sorted.bam")

# load the bamfile as GRanges object:
# assign column names for dataframe:
what <- c("qname","rname","strand","pos","qwidth")
# flag unmapped sequences to leave out:
flag <- scanBamFlag(isUnmappedQuery=FALSE)
# define parameters of bam scan:
param <- ScanBamParam(what=what, flag=flag)

# index bam:
indexBam(inFile)

# load bam:
bam <- scanBam(inFile ,param=param)

save.image(file="/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/repbam_loaded.RData")
#load(file="/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/repbam_loaded.RData")

# assign the imported bam file to the variable name 'results':
results <- bam[[1]]

# assign list of transcripts, sorted by reference chromosome name, to the variable 'id':
id <- as.vector(results$qname)
chrs <- unique(results$rname)

# create vector of chromosomes present in the annotation:
goodChrs <- unique(as.vector(genes@seqnames))

# filter out reads with chromosomes not contained in goodChrs:
# change rname and strand from factors to characters:
results$rname <- as.character(results$rname)
results$strand <- as.character(results$strand)
# convert to dataframe:
df <- as.data.frame(do.call("cbind", results))
# change rname and strand from factors to characters:
df$rname <- as.character(temp$rname)
df$strand <- as.character(temp$strand)

# filter out reads not in goodChars:
m <- match(goodChrs, df$rname)
result <- temp[temp$rname %in% goodChrs, ]

# convert result$pos to numeric vector:
result$rname <- as.factor(result$rname)
result$pos <- as.integer(result$pos)
result$qwidth <- as.integer(result$qwidth)

save.image(file="/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/temp1.RData")
#load(file="/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/temp1.RData")


# convert the bam to GRanges object:
bamGR <- GRanges(
  seqnames = result$rname,
  ranges = IRanges(start = result$pos, width = results$qwidth)
)

ranges = IRanges(start = result$pos, width = results$qwidth)







######
#assign length of all chromosomes without "_" to the variable name 'seq_lengths':
m <- match(chrs, names(seqlengths(Hsapiens)))
seq_lengths <- seqlengths(Hsapiens)[m]

#convert the inFile to GRanges object:
bam_gr=GRanges(
  seqnames = results$rname[correct_chromosomes],
  ranges = IRanges(start=bam[[1]]$pos[correct_chromosomes], width=bam[[1]]$qwidth[correct_chromosomes]),
  strand = bam[[1]]$strand[correct_chromosomes],
  seqlengths = seq_lengths
)

gene_list = split(genes, genes$gene_id)



