library(rtracklayer)
library(GenomicRanges)


rCounts1 <- readRDS("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/exp2/c_RepeatCounts/FT1_cRepeatCounts_half1a.rds")

rCounts2 <- readRDS("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/exp3/FT1_DNAfeatureCounts.rds")

rCounts2 <- lapply(rCounts2, function(x) {
  stat <- x$stat[1,]
  colnames(stat) <- c("Type", "featureCounts")
  return(stat)
})

bothCountsDF <- cbind(rbind(rCounts1$DNA1, 
                            rCounts1$DNA2), do.call("rbind", rCounts2)[,2])

colnames(bothCountsDF) <- c("countOverlaps_counts", 
                            "featureCounts_counts")

write.table(bothCountsDF, file="/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/Robjects/FT1repeatCountsComparison.tab", quote=F)


ids <- c("ENSG00000091831", "ENSG00000196776", 
         "ENSG00000105173", "ENSG00000111640", 
         "ENSG00000075624", "ENSG00000204574", 
         "ENSG00000104904", "ENSG00000089157")

gcCounts1 <- readRDS("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/exp2/GCcountsDF/FT1_GCcountsDFa.rds")
for (i in 1:length(ids)) {
  if (i==1) {
    ctlCounts1 <- as.data.frame(gcCounts1[grep(ids[i],
                                               gcCounts1$gene_id),])
  } else {
    ctlCounts1[i,] <- gcCounts1[grep(ids[i], 
                                     gcCounts1$gene_id),]
  }
}

gcCounts2 <- readRDS(file="/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/exp3/GCcounts1.rds")
gcCounts2 <- as.data.frame(gcCounts2$counts)
gcCounts2$gene_id <- rownames(gcCounts2)
colnames(gcCounts2) <- c("counts", "gene_id")
gcCounts2$gene_id <- gsub("\\..*$", "", gcCounts2$gene_id)
gcCounts2 <- aggregate(.~gene_id, gcCounts2, sum)
rownames(gcCounts2) <- gcCounts2$gene_id

for (i in 1:length(ids)) {
  if (i==1) {
    ctlCounts2 <- as.data.frame(gcCounts2[grep(ids[i],
                                               gcCounts2$gene_id),])
  } else {
    ctlCounts2[i,] <- gcCounts2[grep(ids[i], 
                                     gcCounts2$gene_id),]
  }
}


ctlCounts <- cbind(ctlCounts1[,2], ctlCounts2[,2])
colnames(ctlCounts) <- c("countOverlaps", "featureCounts")

saveRDS(ctlCounts, file="/Users/jamestorpy/clusterHome/projects/gsoc_repeats/RNA-seq/Robjects/FT1ctlCountsComparison.tab")

pdf("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/FT1ctlCountsComparison.pdf")
plot(log10(gcCounts1[,2]+1), log10(gcCounts2[,2]+1), xlab="countOverlaps_log10_counts", ylab="featureCounts_log10_counts")
dev.off()


# check individual ranges from gcCounts 2 for IGV viewing:
# import gencode annotation:
GCgenes <- readRDS("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/exp3/GCgenes.rds")

gcCounts2[4,]
GCgenes[grep("45713", GCgenes$gene_id),]




