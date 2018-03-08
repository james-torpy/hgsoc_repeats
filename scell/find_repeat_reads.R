### find_repeat_reads.R ###

# Finds overlaps of repeat annotation gtf and single-cell sam

# run on cluster with:
#briR
#qsub -N scellRP -b y -wd \
#/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/single-cell/logs \
#-j y -R y -pe smp 4 -V "Rscript /share/ScratchGeneral/jamtor/projects/hgsoc_repeats/single-cell/scripts/find_repeat_reads.R"


library(rtracklayer)
library(GenomicRanges)

project <- "hgsoc_repeats"
Type <- "single-cell"

# define directories:
#homeDir <- "/Users/jamestorpy/clusterHome/"
homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project, "/")
resultsDir <- paste0(projectDir, "/", Type, "/results/")
RobjectDir <- paste0(projectDir, "/", Type, "/RobjectDir/")
genomeDir <- "/home/jamtor/genomes/scell/GRCh38_rp_w_artificial/"

system(paste0("mkdir -p ", RobjectDir))


### 0. Load in files ###

if ( !file.exists(paste0(RobjectDir, "/samAlign.rds")) ) {
  
  print("Loading repeats annotation...")
	annot <- import(paste0(genomeDir, "/repeats_CD3G.gtf"))
	
	print("Loading sam file...")
	sam <- read.table(file=paste0(resultsDir, "/4404_primary_rp_w_artificial/outs/possorted_genome_sam.sam"), header=F, fill=T)
  
	print("Saving image...")
	save.image(paste0(RobjectDir, "/annot_and_sam_read.RData"))

	print("Creating sam alignment object...")
	samAlign <- readGAlignmentPairs(sam)

	print("Saving sam alignment onbject...")
	saveRDS(samAlign, paste0(RobjectDir, "/samAlign.rds"))
} else {
	samAlign <- load(paste0(RobjectDir, "/samAlign.rds"))
}

print("Creating sam GR object...")
samGR <- GRanges(
	seqnames=seqnames(samAlign),
	ranges=IRanges(start(samAlign), end(samAlign)),
	strand=strand(samAlign),
	cigar=cigar(samAlign)
)

print("Saving sam GR object...")
saveRDS(samGR, paste0(RobjectDir, "/samGR.rds"))

olaps <- findOverlaps(samGR, annot)

print("Saving overlaps...")
saveRDS(olaps, paste0(RobjectDir, "/overlaps.rds"))

save.image(paste0(RobjectDir, "/overlaps_done.RData"))


