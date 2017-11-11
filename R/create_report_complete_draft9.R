# This R script counts overlaps of multiple RNA-seq data sets with a GTF file of transcripts

### 0. Set up variables and directories ###

# load packages needed:
library(GenomicRanges)
library(ShortRead)
library(rtracklayer)
library("BSgenome.Hsapiens.UCSC.hg38")
library(seqinr)
library(ggplot2)

# define starting variables:
project <- "hgsoc_repeats"
methodName <- "method1"
gtfName <- "human-89.repeats4.gtf"

# define directories:
#homeDir <- "/Users/jamestorpy/clusterHome/Users/jamestorpy/clusterHome"
homeDir <- "/Users/jamestorpy/clusterHome/"
projectDir <- paste0(homeDir, "/projects/", project)
reportDir <- paste0(projectDir, "/RNA-seq/results/method_reports/")
bamDir <- paste0(projectDir, "/RNA-seq/results/star/", methodName)
refDir <- paste0(projectDir, "/RNA-seq/refs/")
resultsDir <- paste0(projectDir, "/RNA-seq/results")
L1dir <- paste0(refDir, "star/L1")
fastqDir <- paste0(projectDir, "/RNA-seq/results/trimgalore/subset/")
outDir <- paste0(resultsDir, "/R/", methodName)

# create outDir:
system(paste0("mkdir -p ", outDir))
print(paste0("The outDir file is: ", outDir))

# fetch vector of report files:
rFiles <- list.files(reportDir, full.names = T)

# define the sampleNames:
sampleNames <- gsub(paste0(methodName, "_"), "", basename(rFiles))
sampleNames <- gsub("_report.txt", "", sampleNames)

print("The samples to be processed are: ")
print(sampleNames)


### 1. Load in fastq, STAR output and bam information ###

# load in reports for fastq, STAR output and bam information as list:
reports <- list()
i=1
for (f in rFiles) {
  writeLines("\n")
  print(paste0("Reading in ", f))
  r <- read.table(file = f, sep = "\t")
  reports[[i]] <- r
  colnames(reports[[i]]) <- c("entry", "value")
  reports[[i]]$entry <- gsub(":", "", reports[[i]]$entry)
  i <<- i+1
}

# name each element of reports list:
names(reports) <- sampleNames


### 2. Load in repeats annotation file ###

# define gtf variable:
gtf <- paste0(refDir, gtfName)
print(paste0("The gtf file is: ", gtf))

# load in the gene coordinates from the gtf file:
genes <- import(gtf)

# save RData file:
save.image(file=paste0("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/", methodName, "repGTF_loaded.RData"))
methodName <- "method1"
#load(file=paste0("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/", methodName, "repGTF_loaded.RData"))


### 3. Add telomeres annotation to repeats gtf ###

# fetch telomere annotations from UCSC using rtracklayer:
# open a UCSC session for hg38:
session <- browserSession("UCSC")
genome(session) <- "hg38"
# define track to query - can fetch list of possible tracks using trackNames(session):
gapTrack <- ucscTableQuery(session, "Gap", GRangesForUCSCGenome("hg38"))
# define table to query and fetch it - can fetch list of possible tables using tableNames(gapTrack):
tableName(gapTrack) <- "gap"
gapTrack <- getTable(teloAnn)
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


### 4. Load in bam files ###

# load in the bamFiles as list of GRanges objects:
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
  if (file.exists(paste0(bamDir, "/", s, "/Aligned.sortedByCoord.out.sorted.bam.bai"))) {
    print(paste0("No need to index, ", bamDir, "/", s, "/Aligned.sortedByCoord.out.sorted.bam has already been indexed"))
  } else {
    print(paste0("Indexing ", bamDir, "/", s, "/Aligned.sortedByCoord.out.sorted.bam"))
    indexBam(paste0(bamDir, "/", s, "/Aligned.sortedByCoord.out.sorted.bam"))
  }
  # save complete filename of bam for loading:
  file <- paste0(bamDir, "/", s, "/Aligned.sortedByCoord.out.sorted.bam")
  
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


### 5. Count overlaps of bams with repeat annotation gtf:

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

# filter only counts above 0:
nonZero <- lapply(counts, function(x) {
  result <- x[grep(0, x$counts, invert = T),]
  return(result)
})
# name the nonZero dataframes in the list:
names(nonZero) <- sampleNames

# count total number of reads overlapping with repeats:
cSum <- lapply(nonZero, function(x) {
  result <- sum(x$counts)
  return(result)
})

# filter only counts above 1:
aboveOne <- lapply(nonZero, function(x) {
  result <- x[grep(1, x$counts, invert = T),]
  return(result)
})
# name the aboveOne dataframes in the list:
names(aboveOne) <- sampleNames

# count total number of reads with at least 2 overlaps with repeats:
over2Sum <- lapply(aboveOne, function(x) {
  result <- sum(x$counts)
  return(result)
})

# determine the count instances to assess the range of counts:
cRange <- lapply(counts, function(x) {
  result <- unique(x$counts)
  return(result)
})
# name the count instances dataframes in the list:
names(cRange) <- sampleNames

# save the counts as csv files:
i=1
for (s in sampleNames) {
  if (file.exists(paste0(outDir, "/", s, "_overlaps.csv"))) {
      print(paste0(outDir, "/", s, "_overlaps.csv already exists"))
    } else {
      print(paste0("creating", outDir, "/", s, "_overlaps.csv"))
      write.csv(counts[[i]], file = paste0(outDir, "/", s, "_overlaps.csv"))
    }
  
  i <<- i+1
}

save.image(file=paste0("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/", methodName, "_countsOverlapped.RData"))
methodName <- "method1"
#load(file=paste0("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/", methodName, "_countsOverlapped.RData"))

L1s <- lapply(nonZero, function(x) {
  ind <- grep("L1", x$fGenes)
  result <- x[ind, ]
  return(result)
})

Alus <- lapply(nonZero, function(x) {
  ind <- grep("Alu", x$fGenes)
  result <- x[ind, ]
  return(result)
})

centSats <- lapply(nonZero, function(x) {
  ind <- grep("centromere", x$fGenes)
  result <- x[ind, ]
  return(result)
})

simpleRs <- lapply(nonZero, function(x) {
  ind <- grep(")n", x$fGenes)
  result <- x[ind, ]
  return(result)
})

others <- lapply(nonZero, function(x) {
  ind <- grep("L1", x$fGenes, invert = T)
  ind <- grep("Alu", x$fGenes, invert = T)  
  ind <- grep("centromere", x$fGenes, invert = T)
  ind <- grep(")n", x$fGenes, invert = T)
  result <- x[ind, ]
  return(result)
})

bowtellFT1 <- list(L1s[[1]], Alus[[1]], centSats[[1]], simpleRs[[1]], others[[1]])
names(bowtellFT1) <- c("L1s", "Alu", "centromere", "simple", "other")

kur_v1 <- list(L1s[[2]], Alus[[2]], centSats[[2]], simpleRs[[2]], others[[2]])
names(kur_v1) <- c("L1s", "Alu", "centromere", "simple", "other")

bowtellCounts <- lapply(bowtellFT1, function(x) {
  counts <- sum(x$counts)
  return(counts)
})

bowtellCountDF <- data.frame(column1 = names(bowtellCounts), column2 = integer(5))

for (i in 1:5) {
  bowtellCountDF[i, 2] <- bowtellCounts[[i]]
}
colnames(bowtellCountDF) <- c("repeats", "overlaps")

bp <- ggplot(bowtellCountDF, aes(x=repeats, y=overlaps, fill=repeats))
bp <- bp + geom_bar(stat="identity", position = "dodge")
bp <- bp + scale_y_log10()
pdf(file = paste0(outDir, "/bowtellFT1_repeat_overlaps.pdf"))
bp
dev.off()



kurCounts <- lapply(kur_v1, function(x) {
  counts <- sum(x$counts)
  return(counts)
})

kurCountDF <- data.frame(column1 = names(kurCounts), column2 = integer(5))

for (i in 1:5) {
  kurCountDF[i, 2] <- kurCounts[[i]]
}
colnames(kurCountDF) <- c("repeats", "overlaps")

kp <- ggplot(kurCountDF, aes(x=repeats, y=overlaps, fill=repeats)
kp <- kp + geom_bar(stat="identity", position = "dodge")
kp <- kp + scale_y_log10()
pdf(file = paste0(outDir, "/kur_v1_repeat_overlaps.pdf"))
kp
dev.off()

###### repeats mapped to gencode_v24 and hg38_ercc with 4% mismatches allowed, multimappers were mapped to only one of the best matched loci chosen at random ######





### 6. Create bar plots of counts of each repeat type ###





### 7. Fetch L1 qnames from bam so sequences can be identified within  ###

# subset genes to only include L1s:
L1index <- grep("L1", fGenes$gene_id)
L1genes <- fGenes[L1index, ]

# find overlaps of bams with L1 annotation ranges to identify bam entries:
# create empty list to put qnames of each bam into:
L1qnames <- lapply(bamGRs, function(x) {
  writeLines("\n")
  print("Finding overlaps...")
  olaps <- findOverlaps(L1genes, x)
  # subset the indices of the bamGR which overlaps with L1s:
  #L1ind <- olaps@to
  # for commandline R:
  L1ind <- olaps@subjectHits
  # match indices to rows of bamGR to find those overlapping with L1s
  L1s <- x[L1ind,]
  # fetch vector of qnames of the L1s and put into list element:
  qnames <- L1s$qnames
  # convert qnames to dataframe:
  qnames <- as.data.frame(qnames)
  colnames(qnames) <- NULL
  return(qnames)
})
# name the qname vectors in the list:
names(L1qnames) <- sampleNames

# save the L1qnames as csv files:
i=1
for (s in sampleNames) {
  if (file.exists(paste0(outDir, "/", s, "_L1qnames.csv"))) {
    print(paste0(outDir, "/", s, "_L1qnames.csv already exists"))
  } else {
    print(paste0("creating ", outDir, "/", s, "_L1qnames.csv"))
    write.csv(L1qnames[[i]], file = paste0(outDir, "/", s, "_L1qnames.csv"), row.names = F, quote = F)
  }
  
  i <<- i+1
}


### 8. Map L1 overlaps to L1 consensus using STAR ###

# check for existence of STAR L1 consensus reference:
if (file.exists(paste0(L1dir, "/SA"))) {
  print("L1 consensus STAR reference exists")
} else {
  print("L1 consensus STAR reference does not exist and must be built")
}

# Fetch fastq filenames in fastq list:
fastqs <- list.files(fastqDir, recursive = T, pattern = "val_[1-2].fq.gz")
fastqs <- grep("grant", fastqs, invert = T, value = T)
fastqs <- grep("gtx", fastqs, invert = T, value = T)
fastqs <- grep("PR", fastqs, invert = T, value = T)
fastqs <- grep("FT[2-3]", fastqs, invert = T, value = T)
fastqs <- grep("v[2-3]", fastqs, invert = T, value = T)

save.image(file=paste0("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/", methodName, "_temp1.RData"))
methodName <- "method1"
load(file=paste0("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/", methodName, "_temp1.RData"))


for (f in fastqs) {
  uID <- gsub("_val_[1-2].fq.gz", "", basename(f))
  refID <- gsub("_R1", "", uID)
  system(paste0("gunzip -c ", fastqDir, "/", f, " | grep -A3 -Ff ", outDir, "/", refID, "_L1qnames.csv | sed '/^--$/d' > ", outDir, "/", uID, "_L1s.fq"))
}




# Fetch L1qname entries from original fastqs, put them into seperate fastqs and map them to L1 consensus using STAR



