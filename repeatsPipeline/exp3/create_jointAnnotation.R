> library("rtracklayer")
gc <- import("/share/ScratchGeneral/jamtor/genomes/gc_and_repeat_noOL/gencode_v24_hg38_annotation.gtf")
gc <- gc[grep("exon", gc$type),]

rp <- import("/share/ScratchGeneral/jamtor/genomes/gc_and_repeat_noOL/human-89.repeats2.gtf")

ol <- findOverlaps(rp, gc)

gcAdj <- gc[-subjectHits(ol),]

rp <- GRanges(seqnames = seqnames(rp),
	ranges = ranges(rp),
	strand = strand(rp),
	source = rp$source,
	type = rp$type,
	score = rp$score,
	phase = rp$phase,
	gene_id = rp$gene_id,
	gene_type = rp$gene_type,
	gene_status = rp$gene_status,
	gene_name = rp$gene_name,
	transcript_id = rp$transcript_id,
	transcript_type = rp$transcript_type,
	transcript_status = rp$transcript_status,
	transcript_name = rp$transcript_name,
	exon_number = rp$exon_number,
	exon_id = rp$exon_id,
	protein_id = rep("NA", length(rp))
	)

gcAdj <- subset(gcAdj, select=-c(tag, transcript_support_level, havana_transcript))