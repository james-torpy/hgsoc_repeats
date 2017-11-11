
library(GenomicRanges)
library(rtracklayer)

# define directories:
inDir <- "/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/annotations"

# read in repeat annotations file:
repGR <- import(paste0(inDir, "/repeats.gtf"))

# remove phase, feature columns from repGR:
repGR$phase <- NULL
repGR$feature <- NULL

# split attribute column into gene_id, gene_type:
spl <- strsplit(repGR$attribute, ";")

gene_id <- lapply(spl, function(x) {
  entry <- x[2]
  return(entry)
})

gene_type <- lapply(spl, function(x) {
  entry <- x[1]
  return(entry)
})

# convert from lists to dataframes:
gene_id <- as.data.frame(do.call("rbind", gene_id))
gene_type <- as.data.frame(do.call("rbind", gene_type))

# define gene_status, gene_name:
gene_status <- as.data.frame(rep("NOVEL", nrow(gene_id)))
colnames(gene_status) <- "V1"
gene_name <- gene_id

# add attribute columns to repGR and remove original column:
repGR$gene_id <- gene_id$V1
repGR$transcript_id <- gene_id$V1
repGR$gene_type <- gene_type$V1
repGR$gene_status <- gene_status$V1
repGR$gene_name <- gene_name$V1
repGR$transcript_type <- gene_type$V1
repGR$transcript_status <- gene_status$V1
repGR$transcript_name <- gene_name$V1
repGR$exon_number <- "1"

repGR$attribute <- NULL

# change 'type' column to 'transcript' to assume all repeats are transcribed:
repGR$type <- rep("transcript", length(repGR$type))

# add 'level' column with value '3' for all - represents entries that have been automatically annotated:
repGR$level <- rep("3", length(repGR$source))

# export as GTF:
export(repGR, "/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/refs/human-89.repeats3.gtf")

save.image(file="/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/rep_annot_gen_done.RData")
#load(file="/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/rep_annot_gen_done.RData")

repGR <- repGR[, c(1:7, 15, 8:14)]
