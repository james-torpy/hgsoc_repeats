exp3_lrt <- readRDS("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/exp3/FT_vs_prPT_lrtGC.rds")
exp3Genes <- as.data.frame(topTags(exp3_lrt, n=Inf))
exp3Genes$threshold <- exp3Genes$FDR < 0.05

exp4_lrt <- readRDS("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/exp4/FT_vs_prPT_lrtGC.rds")
exp4Genes <- as.data.frame(topTags(exp4_lrt, n=Inf))
exp4Genes$threshold <- exp4Genes$FDR < 0.05

exp5_lrt <- readRDS("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/exp5/FT_vs_prPT_lrtGC.rds")
exp5Genes <- as.data.frame(topTags(exp5_lrt, n=Inf))
exp5Genes$threshold <- exp5Genes$FDR < 0.05

df3v4 <- subset(merge(exp3Genes, exp4Genes, by="row.names"), select=c(Row.names, logFC.x, FDR.x, threshold.x, logFC.y, FDR.y, threshold.y))
colnames(df3v4) <- c("rownames", "exp3logFC", "exp3FDR", "exp3threshold", "exp4logFC", "exp4FDR", "exp4threshold")

i=1
for (i in 1:nrow(df3v4)) {
  if (i==1) {
    if (df3v4[i,]$exp3threshold == TRUE & df3v4[i,]$exp4threshold == TRUE) {
      all_threshold <- c("Both FDR < 0.05")
    } else if (df3v4[i,]$exp3threshold == TRUE & df3v4[i,]$exp4threshold == FALSE) {
      all_threshold <- c("featureCounts FDR < 0.05")
    } else if (df3v4[i,]$exp3threshold == FALSE & df3v4[i,]$exp4threshold == TRUE) {
      all_threshold <- c("RSEM FDR < 0.05")
    } else {
      all_threshold <- c("None significant")
    }
  } else {
    if (df3v4[i,]$exp3threshold == TRUE & df3v4[i,]$exp4threshold == TRUE) {
      all_threshold[i] <- "Both FDR < 0.05"
    } else if (df3v4[i,]$exp3threshold == TRUE & df3v4[i,]$exp4threshold == FALSE) {
      all_threshold[i] <- "featureCounts FDR < 0.05"
    } else if (df3v4[i,]$exp3threshold == FALSE & df3v4[i,]$exp4threshold == TRUE) {
      all_threshold[i] <- "RSEM FDR < 0.05"
    } else {
      all_threshold[i] <- "None significant"
    }
  }
}

df3v4 <- cbind(subset(df3v4, select=-c(exp3threshold, exp4threshold)), all_threshold)

p <- ggplot(data=df3v4, aes(x=exp4logFC, y=exp3logFC, color=all_threshold))
p <- p + geom_point(data=df3v4)
p <- p + labs(x="RSEM log2FC", y="featureCounts log2FC")
pdf("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/results/R/comparison/featureCounts_vs_RSEM_GClog2FC.pdf")
p
dev.off()



df5v4 <- subset(merge(exp5Genes, exp4Genes, by="row.names"), select=c(Row.names, logFC.x, FDR.x, threshold.x, logFC.y, FDR.y, threshold.y))
colnames(df5v4) <- c("rownames", "exp5logFC", "exp5FDR", "exp5threshold", "exp4logFC", "exp4FDR", "exp4threshold")

i=1
for (i in 1:nrow(df5v4)) {
  if (i==1) {
    if (df5v4[i,]$exp5threshold == TRUE & df5v4[i,]$exp4threshold == TRUE) {
      all_threshold <- c("Both FDR < 0.05")
    } else if (df5v4[i,]$exp5threshold == TRUE & df5v4[i,]$exp4threshold == FALSE) {
      all_threshold <- c("HTSeq FDR < 0.05")
    } else if (df5v4[i,]$exp5threshold == FALSE & df5v4[i,]$exp4threshold == TRUE) {
      all_threshold <- c("RSEM FDR < 0.05")
    } else {
      all_threshold <- c("None significant")
    }
  } else {
    if (df5v4[i,]$exp5threshold == TRUE & df5v4[i,]$exp4threshold == TRUE) {
      all_threshold[i] <- "Both FDR < 0.05"
    } else if (df5v4[i,]$exp5threshold == TRUE & df5v4[i,]$exp4threshold == FALSE) {
      all_threshold[i] <- "HTSeq FDR < 0.05"
    } else if (df5v4[i,]$exp5threshold == FALSE & df5v4[i,]$exp4threshold == TRUE) {
      all_threshold[i] <- "RSEM FDR < 0.05"
    } else {
      all_threshold[i] <- "None significant"
    }
  }
}

df5v4 <- cbind(subset(df5v4, select=-c(exp5threshold, exp4threshold)), all_threshold)

p <- ggplot(data=df5v4, aes(x=exp4logFC, y=exp5logFC, color=all_threshold))
p <- p + geom_point(data=df5v4)
p <- p + labs(x="RSEM log2FC", y="HTSeq log2FC")
pdf("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/results/R/comparison/HTSeq_vs_RSEM_GClog2FC.pdf")
p
dev.off()


