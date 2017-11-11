### rpCompBarPlot_prepare_cCounts.R ###

# breaks cCounts and tCounts up into dataframe for each class and saves each a RData objects #

# define starting variables:
project <- "hgsoc_repeats"
expName <- "exp5"
Type <- "custom3"

# define directories:
#homeDir <- "/share/ScratchGeneral/jamtor/"
homeDir <- "/Users/jamestorpy/clusterHome"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results/")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/")
inDir <- paste0(resultsDir, "/htseq/", expName)
plotDir <- paste0(inDir, "/plots/linePlots/")
rawDir <- paste0(projectDir, "/RNA-seq/raw_files/fullsamples/bowtell_primary/")

system(paste0("mkdir -p ", plotDir))


### 1. Load in repeat counts ###

countFiles <- grep(
  ".sam", list.files(inDir, pattern="custom3rep", recursive = T, full.names = T), value=T, invert=T
)
countFiles <- grep("arPT|rcAF", countFiles, invert=T, value=T)
sampleNames <- gsub("\\..*$", "", basename(countFiles))

######
countFiles <- c(countFiles[2:4], countFiles[9:10], countFiles[21])
######

for (i in 1:length(countFiles)) {
  if (i==1) {
    custom3Counts <- read.table(file=countFiles[i])
  } else {
    custom3Counts <- cbind(custom3Counts, read.table(file=countFiles[i])[,2])
  }
}

colnames(custom3Counts) <- c("gene_id", sampleNames)


### 2. Load in gc counts ###

GCcountFiles <- grep(
  ".sam", list.files(inDir, pattern="gc_", recursive = T, full.names = T), value=T, invert=T
)

for (i in 1:length(GCcountFiles)) {
  if (i==1) {
    gcCounts <- read.table(file=GCcountFiles[i])
  } else {
    gcCounts <- cbind(gcCounts, read.table(file=GCcountFiles[i])[,2])
  }
}
colnames(gcCounts) <- c("gene_id", sampleNames)

# aggregate multiple types of same gene:
gcCounts$gene_id <- gsub("\\..*$", "", gcCounts$gene_id)
gcCounts <- aggregate(.~gene_id, gcCounts, mean)

# remove stat rows:
gcCounts <- gcCounts[-grep("__", gcCounts$gene_id),]

# load in patient ids:
Key <- read.table(file=paste0(rawDir, "/sampleKey.txt"), sep=" ", fill=T)[,3:4]
Key <- Key[3:nrow(Key),]
Key$V4 <- gsub("_[A-Z].*$", "", gsub("^.*AOCS_", "AOCS_", Key$V4))
Key <- Key[with(Key, order(V3)), ]

# only include samples in cCounts:
Key <- dplyr::filter(Key, Key$V3 %in% colnames(custom3Counts))
newNames <- paste0(Key$V4, "_", gsub("[0-9]", "", colnames(subset(custom3Counts, select=-gene_id))))
colnames(custom3Counts)[2:ncol(custom3Counts)] <- newNames
colnames(gcCounts)[2:ncol(gcCounts)] <- newNames


# save the counts:
if (!file.exists(paste0(RobjectDir, "/",  expName, "/", Type, "_allcounts.htseq.rds"))) {
  saveRDS(custom3Counts, file=paste0(RobjectDir, "/",  expName, "/custom3_allcounts.htseq.rds"))
}

if (!file.exists(paste0(RobjectDir, "/",  expName, "/gc_allcounts.htseq.rds"))) {
  saveRDS(gcCounts, file=paste0(RobjectDir, "/",  expName, "/gc_allcounts.htseq.rds"))
}
