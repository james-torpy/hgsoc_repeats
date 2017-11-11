
library(GenomicRanges)
library(rtracklayer)

# define directories:
inDir <- "/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/annotations"

# read in repeat annotations file:
repGTF <- read.table(file=paste0(inDir, "/gencode_v24_hg38_repeats_final.gtf"), sep = "\t", header = F)

# add colnames:
colnames(repGTF) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

# change "*" to "+" in strand field:
repGTF[,7] <- gsub("\\*", "\\+", repGTF[,7])

save.image("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/annotations/GTF_loaded.RData")
#load("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/annotations/GTF_loaded.RData")

# convert repGTF data frame to GRanges object
repGR <- makeGRangesFromDataFrame(repGTF)

# add extra columns from repGTF:
repGR$source <- repGTF$source
repGR$score <- repGTF$score
repGR$gene_id <- strsplit(" ", repGTF$attribute)

test <- strsplit(" ", repGTF$attribute)

# export as GTF:
export(repGR, "/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/annotations/gencode_v24_hg38_repeats_final.gtf")

save.image("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/annotations/GTF_made.RData")
#load("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/annotations/GTF_made.RData")

######

repGR <- GRanges(
  seqnames = Rle(repGTF[ ,1]),
  source = repGTF[ ,2],
  feature = repGTF[ ,3],
  ranges = IRanges()
)

test <- repGTF[16842200:16842242,]
test[,7] <- as.character(test[,7])

test[,7]
gsub("\\*", "\\+", test[,7])

test[,7] <- gsub("\\*", "\\+", test[,7])
######