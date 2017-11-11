library(GenomicRanges)
library(rtracklayer)

grL1 <- GRanges(
        seqnames = rep("chrL1:1-558", 2),
        strand = rep("+", 2),
        source = rep("EMBL-EBI", 2),
        type = c("transcript", "exon"),
        score = rep(NA, 2),
        frame = rep(".", 2),
        gene_id = rep("L1", 2),
        transcript_id = rep("L1", 2),
        gene_type = rep("type_i_transposons_LINE", 2),
        gene_status = rep("KNOWN", 2),
        gene_name = rep("L1", 2),
        transcript_type = rep("type_i_transposons_LINE", 2),
        transcript_status = rep("KNOWN", 2),
        transcript_name = rep("L1", 2),
        exon_number = c(NA, "1"),
        exon_id = c(NA, "L1"),
        level = rep(2, 2)
)


# export as GTF:
export(grL1, "/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/refs/L1_consensus3.gtf")

save.image(file="/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/L1_annot_gen_done.RData")
#load(file="/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/rep_annot_gen_done.RData")
