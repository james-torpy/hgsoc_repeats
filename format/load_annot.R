inDir <- "/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/annotations"

reps <- read.table(file = paste0(inDir, "/human-89.repeats.tab"), sep = "\t", na.strings = T, header = T)

head(reps)

gtf <- data.frame(reps$chr, reps$analysis, reps$class, reps$start, reps$end, reps$score, reps$strand, rep(".", nrow(reps)), paste0(reps$class_desc, "; ", reps$type))

write.table(gtf, file=paste0(inDir, "/human-89.repeats_formatted.tab"), sep="\t", col.names = F, row.names = F)
