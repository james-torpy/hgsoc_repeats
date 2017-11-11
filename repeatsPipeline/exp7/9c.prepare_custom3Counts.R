### 4c.prepare_gcCounts.R ###

# breaks cCounts up into dataframe for each class and saves each as RData objects #

# define starting variables:
project <- "hgsoc_repeats"
expName <- "exp7"
Type <- "custom3"

# define directories:
#homeDir <- "/share/ScratchGeneral/jamtor/"
homeDir <- "/Users/jamestorpy/clusterHome"
projectDir <- paste0(homeDir, "/projects/", project, "/")
resultsDir <- paste0(projectDir, "/RNA-seq/results/")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/", expName, "/")
inDir <- paste0(resultsDir, "/htseq/", expName, "/")
plotDir <- paste0(inDir, "/plots/linePlots/")
rawDir <- paste0(projectDir, "/RNA-seq/raw_files/fullsamples/bowtell_primary/")

system(paste0("mkdir -p ", RobjectDir))
system(paste0("mkdir -p ", plotDir))


### 1. Load in repeat counts ###

# load in custom3 counts for samples run through the whole pipeline:
compFiles <- as.character(read.table(paste0(inDir, "/completed_files.txt"))[,1])

# load in custom3_countFiles into a df:
for (i in 1:length(compFiles)) {
  if (i==1) {
    print(compFiles[i])
    custom3Counts <- data.frame(read.table(file=paste0(inDir, "/", compFiles[i], "/", compFiles[i], ".custom3.htseq.txt")))
  } else {
    print(compFiles[i])
    custom3Counts <- cbind(custom3Counts, read.table(file=paste0(inDir, "/", compFiles[i], "/", compFiles[i], ".custom3.htseq.txt"))[,2])
  }
}
colnames(custom3Counts) <- c("gene_id", compFiles)

# remove duplicate rows:
custom3Counts <- unique(custom3Counts)
# remove specs lines:
custom3Counts <- custom3Counts[grep("__", custom3Counts$gene_id, invert=T),]


### 2. Replace uIDs with patient ids:

# load in patient ids:
Key <- read.table(file=paste0(rawDir, "/sampleKey.txt"), sep=" ", fill=T)[,3:4]
Key <- Key[3:nrow(Key),]
Key$V4 <- gsub("_[A-Z].*$", "", gsub("^.*AOCS_", "AOCS_", Key$V4))
Key <- Key[with(Key, order(V3)), ]

# only include samples in cCounts:
Key <- dplyr::filter(Key, Key$V3 %in% colnames(custom3Counts))
newNames <- paste0(Key$V4, "_", gsub("[0-9]", "", colnames(subset(custom3Counts, select=-gene_id))))
colnames(custom3Counts)[2:ncol(custom3Counts)] <- newNames

# save the counts:
if (!file.exists(paste0(RobjectDir, "/", Type, "_allcounts.htseq.rds"))) {
  saveRDS(custom3Counts, file=paste0(RobjectDir, "/custom3_allcounts.htseq.rds"))
}