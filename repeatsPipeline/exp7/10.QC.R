# This script takes the following inputs:
# 1. STAR gc log
# 2. STAR ribo log
# 3. htseq gc counts
# 4. htseq total repeat counts
# 5. htseq custom3 repeat counts
# and creates:
# a) scatter plot of total mapped gc reads vs total gc counts
# b) barplots with the composition of protein-coding, non-coding, other, repeats and ribosomal counts
# with one bar per sample


### 0. Set up variables and directories ###

library(GenomicRanges)
library(ShortRead)
library(rtracklayer)
library("BSgenome.Hsapiens.UCSC.hg38")
library(ggplot2)
library(reshape2)

# define starting variables:
project <- "hgsoc_repeats"
expName <- "exp7"

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "projects/", project, "/RNA-seq/")
htseqDir <- paste0(projectDir, "/results/htseq/", expName, "/")
STARgcDir <- paste0(projectDir, "/results/star/GC/", expName, "/")
STARriboDir <- paste0(projectDir, "/results/star/ribo/", expName, "/")
RobjectDir <- paste0(projectDir, "/Robjects/", expName, "/")
plotDir <- paste0(projectDir, "/results/R/", expName, "/plots/QC/")
rawDir <- paste0(projectDir, "/raw_files/fullsamples/bowtell_primary/")
gencodeDir <- "/share/ScratchGeneral/jamtor/genomes/hg38_ercc/"
gencodeName <- "gencode_v24_hg38_annotation.gtf"

# create plotDir:
system(paste0("mkdir -p ", plotDir))
print(paste0("The outDir is: ", plotDir))


### 1. Fetch total gc counts ###

gcCounts <- readRDS(file=paste0(RobjectDir, "/gc_allcounts.htseq.rds"))
rownames(gcCounts) <- gcCounts$gene_id
gcCounts <- subset(gcCounts, select=-gene_id)

# define gencode variable:
gencodeFile <- paste0(gencodeDir, gencodeName)
if (file.exists(paste0(RobjectDir, "GCgenes.rds"))) {
  print("Loading gencode...")
  gencode <- readRDS(file=paste0(RobjectDir, "GCgenes.rds"))
} else {
  print(paste0("The gencode file is: ", gencodeFile))
  # load in the gene coordinates from the gencode file:
  gencode <- import(gencodeFile)
  # save gencode as RDF file:
  saveRDS(gencode, file=paste0(RobjectDir, "GCgenes.rds"))
}

# select only gene entries from gencode:
temp <- (gencode[gencode$type %in% "exon"])

# select relevant entries for non-coding or protein coding annotations:
pcGencode <- unique(
  gsub(
    "\\..*$", "", temp[temp$gene_type %in% "protein_coding"]$gene_id
  )
)
            
ncGencode <- unique(
  gsub(
    "\\..*$", "", temp[temp$gene_type %in% c("lincRNA", "non_coding")]$gene_id
  )
)              
                
otherGencode <- unique(
  gsub(
    "\\..*$", "", temp[!temp$gene_type %in% c("protein_coding", "lincRNA", "non_coding")]$gene_id
  )
)

# group the total counts per sample of gcCounts according to above annotations:
Counts <- list(
  apply(
    gcCounts[rownames(gcCounts) %in% pcGencode,], 2, sum
  )
)
    
Counts[[2]] <- apply(
    gcCounts[rownames(gcCounts) %in% ncGencode,], 2, sum
)

Counts[[3]] <- apply(
  gcCounts[rownames(gcCounts) %in% otherGencode,], 2, sum
)


### 2. Fetch the sum of total repeat counts for each sample ###
Counts[[4]] <- readRDS(file=paste0(RobjectDir, "/all_allcounts.htseq.rds"))


### 3. Load STAR gc and ribo logs and fetch total reads mapped ###

# load in patient ids:
Key <- read.table(file=paste0(rawDir, "/sampleKey.txt"), sep=" ", fill=T)[,3:4]
Key <- Key[3:nrow(Key),]
Key$V4 <- gsub("_[A-Z].*$", "", gsub("^.*AOCS_", "AOCS_", Key$V4))
Key <- Key[with(Key, order(V3)), ]

# only include samples in counts above and add patient codes to samplenames:
Key <- dplyr::filter(Key, Key$V4 %in% gsub("\\_[a-zA-Z][a-zA-Z]*$", "", colnames(gcCounts)))

sampleNames <- list.files(STARgcDir)[list.files(STARgcDir) %in% as.character(Key$V3)]
for (i in 1:length(sampleNames)) {
  if (i==1) {
    gcSTARfiles <- c(list.files(paste0(STARgcDir, "/", sampleNames)[i], pattern="Log.final.out", full.names=T))
  } else {
    gcSTARfiles[i] <- list.files(paste0(STARgcDir, "/", sampleNames)[i], pattern="Log.final.out", full.names=T)
  }
}

for (i in 1:length(gcSTARfiles)) {
  Log <- read.table(file=gcSTARfiles[i], header=T, sep="\t", fill=T)
  Log[,2] <- as.character(Log[,2])
  num <- as.numeric(Log[7,2]) + as.numeric(Log[22,2])
  if (i==1) {
    gcMapped <- c(num)
  } else {
    gcMapped[i] <- num
  }
}
names(gcMapped) <- Key$V4

for (i in 1:length(sampleNames)) {
  if (i==1) {
    riboSTARfiles <- c(list.files(paste0(STARriboDir, "/", sampleNames)[i], pattern="Log.final.out", full.names=T))
  } else {
    riboSTARfiles[i] <- list.files(paste0(STARriboDir, "/", sampleNames)[i], pattern="Log.final.out", full.names=T)
  }
}

for (i in 1:length(riboSTARfiles)) {
  Log <- read.table(file=riboSTARfiles[i], header=T, sep="|", fill=T)
  Log[,2] <- as.character(Log[,2])
  num <- as.numeric(gsub("\t", "", Log[7,2])) + as.numeric(gsub("\t", "", Log[22,2]))
  if (i==1) {
    riboMapped <- c(num)
  } else {
    riboMapped[i] <- num
  }
}
names(riboMapped) <- Key$V4

Counts[[5]] <- riboMapped
names(Counts) <- c("protein-coding", "non-coding", "other", "repeat", "ribosome")


### 4. Create total mapped vs total counts scatter plot ###

total_gcCounts <- apply(gcCounts, 2, sum)
mapVcount <- data.frame(gcMapped, total_gcCounts)
mapVcount$samples <- rownames(mapVcount)

p <- ggplot(mapVcount, aes(x=gcMapped, y=total_gcCounts))
p <- p + geom_point()
p

pdf(paste0(plotDir, "/mapped_vs_counts.pdf"))
p
dev.off()


### 5. Create composition barplots ###

# convert Counts to a dataframe:
countsDF <- as.data.frame(t(do.call("cbind", Counts)))
# add rownames as column and melt dataframe:
countsDF$gene_type <- rownames(countsDF)

# calculate percentages as new data frame:
perCountsDF <- as.data.frame(apply(countsDF[,1:ncol(countsDF)-1], 2, function(x) {
  return(as.integer(x)/as.integer(sum(x))*100)
}))
perCountsDF$gene_type <- countsDF$gene_type
 

# create composition barplots of CountsDF and perCountsDF:
cDFs <- list(countsDF, perCountsDF)
Plots <- list()
for (i in 1:2) {
  pCounts <- melt(cDFs[i], variable.name = "sample")
  pCounts$gene_type <- factor(pCounts$gene_type, levels = c("non-coding", "ribosome", "other", "repeat", "protein-coding"))
  
  # plot data as barplot:
  p <- ggplot(pCounts, aes(x=sample, y=value))
  p <- p + geom_bar(stat="identity", aes(fill=gene_type))
  p <- p + scale_fill_manual(values = c("#00BF7D", "#AA4499", "#00B0F6", "#DDCC77", "#F8766D"))
  Plots[[i]] <- p + theme(text = element_text(size=20), axis.text.x = element_text(angle = 90))
  i=i+1
}



pdf(file = paste0(plotDir, "compBarplotCounts.pdf"), height=20, width=35)
Plots[[1]]
dev.off()

pdf(file = paste0(plotDir, "compBarplotPercent.pdf"), height=20, width=35)
Plots[[2]]
dev.off()
  
  
save.image(file=paste0(RobjectDir, "/QCimage.RData"))







