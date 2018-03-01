### 3.repeat_fusions_circos.R ###

# This script takes filtered output from find_repeat_fusions.R 
# and creates circos plots of repeat fusions with at least 6 supporting
# reads

# run on cluster with:
#briR
#qsub -N rpcirc -b y -wd \
#/share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/logs/repeat_fusions/ \
#-j y -R y -pe smp 2 -V "Rscript /share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/fusion_repeats/3.repeat_fusions_circos.R"


### 0. Define variables/paths ###

library(rtracklayer)
library(GenomicRanges)
library("BiocParallel")
library(RCircos)

project <- "hgsoc_repeats"
expName <- "exp9"

Type <- "fusion_repeats"
descrip <- "fusion_repeats_6_supporting_reads"


# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project, "/")
resultsDir <- paste0(projectDir, "/RNA-seq/results")
starDir <- paste0(resultsDir, "/star/", descrip, "/")
rpAnnotDir <- paste0(projectDir, "/RNA-seq/refs/repeats/")
gcAnnotDir <- paste0(projectDir, "/RNA-seq/refs/")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
                     expName, "/fusion_repeats/")
plotDir <- paste0(resultsDir, "/R/", expName,
                  "/plots/", descrip)
newRobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
                        expName, "/", descrip, "/")

system(paste0("mkdir -p ", plotDir))
system(paste0("mkdir -p ", newRobjectDir))


### 1. Load data ###

fusions <- readRDS(paste0(RobjectDir, "/moderate_fusions.rds"))

# reduce to fusions with at least 6 supporting reads:
fusions <- fusions[fusions$readNo >5]


### 2. Prepare circos plot inputs ###

# make vector of genes to label:
genes2label <- na.omit(unique(c(fusions$repeatHits, fusions$acceptorHits)))

# load in annotations and fetch co-ordinates of genes to label:
repeatAnnot <- import(paste0(rpAnnotDir, "/human-89.repeats.gtf"))
# reduce annotations to seqnames, ranges, strand and ID, then concatenate:
values(repeatAnnot) <- c(repeatAnnot$type)
colnames(values(repeatAnnot)) <- c("gene_id")

gcAnnot <- import(paste0(gcAnnotDir, "/gencode_v24_hg38_annotation.gtf"))
values(gcAnnot) <- c(gcAnnot$gene_id)
colnames(values(gcAnnot)) <- c("gene_id", "symbol")

allAnnot <- c(repeatAnnot[repeatAnnot$gene_id %in% genes2label], gcAnnot[gcAnnot$gene_id %in% genes2label])
allAnnot <- allAnnot[grep("[0-9].[0-9]", seqnames(allAnnot), invert=T),]

# only keep regions where fusions occurred:
oLaps <- findOverlaps(allAnnot, fusions)
fusAnnot <- unique(allAnnot[queryHits(oLaps)])

# replace ensembl_ids with symbols, and remove symbol column:
for ( i in 1:length(fusAnnot) ) {
  if ( (length(grep("ENSG", fusAnnot$gene_id[i])) > 0) ) {
    sym <- fusions$symbol[grep(fusAnnot$gene_id[i], fusions$ensembl_id)]
    if ( !(is.na(sym)) ) {
      levels(fusAnnot$gene_id) <- c(levels(fusAnnot$gene_id), sym)
      fusAnnot$gene_id[i] <- sym
    }
  }
}

# create gene labels dataframe:
geneLabels <- data.frame( seqnames(fusAnnot), start(ranges(fusAnnot)), end(ranges(fusAnnot)), fusAnnot$gene_id )
colnames(geneLabels) <- c("Chromosome", "chromStart", "chromEnd", "Gene")

# create links dataframe:
links <- data.frame(seqnames(fusions), start(ranges(fusions)), end(ranges(fusions)), fusions$acceptorChr, 
                    fusions$acceptorBP1, fusions$acceptorBP1)
colnames(links) <- c("Chromosome", "chromStart", "chromEnd", "Chromosome.1", "chromStart.1", "chromEnd.1")


### 3. Create circos plot ###

# load in hg38 genome:
data(UCSC.HG38.Human.CytoBandIdeogram)
hg38 <- UCSC.HG38.Human.CytoBandIdeogram

# set circos plot basic parameters:
chrExclude <- NULL
cytoInfo <- hg38
tracksInside <- 3
tracksOutside <- 0
RCircos.Set.Core.Components(cytoInfo, chrExclude, tracksInside, tracksOutside)

# create circos plot:
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

# create gene label connectors:
nameCol <- 4
side <- "in"
trackNo <- 1
RCircos.Gene.Connector.Plot(geneLabels, trackNo, side)

# create gene labels:
trackNo <- 2
RCircos.Gene.Name.Plot(gene.data = geneLabels, name.col = 4, track.num = 2, side = "in")

# create links plot:
RCircos.Link.Plot(link.data = links, track.num = 3, by.chromosome = TRUE)


### 4. Create circos plot pdf ###

pdf(file=paste0(plotDir, "/fusionCircos.pdf"), height=8, width=8, compress=TRUE)
# set circos plot basic parameters:
chrExclude <- NULL
cytoInfo <- hg38
tracksInside <- 3
tracksOutside <- 0
RCircos.Set.Core.Components(cytoInfo, chrExclude, tracksInside, tracksOutside)

# create circos plot:
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

# create gene label connectors:
nameCol <- 4
side <- "in"
trackNo <- 1
RCircos.Gene.Connector.Plot(geneLabels, trackNo, side)

# create gene labels:
trackNo <- 2
RCircos.Gene.Name.Plot(gene.data = geneLabels, name.col = 4, track.num = 2, side = "in")

# create links plot:
RCircos.Link.Plot(link.data = links, track.num = 3, by.chromosome = TRUE)
dev.off()


