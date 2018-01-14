### 10a.prepare_gcCounts.R ###

# breaks gcCounts up into dataframe for each class and saves each as RData objects #

# define starting variables:
project <- "hgsoc_repeats"
expName <- "exp9"

# define directories:
#homeDir <- "/share/ScratchGeneral/jamtor/"
homeDir <- "/Users/jamestorpy/clusterHome"
projectDir <- paste0(homeDir, "/projects/", project, "/")
resultsDir <- paste0(projectDir, "/RNA-seq/results/")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/", expName, "/")
inDir <- paste0(resultsDir, "/htseq/", expName)
plotDir <- paste0(inDir, "/plots/linePlots/")
rawDir <- paste0(projectDir, "/RNA-seq/raw_files/fullsamples/bowtell_primary/")

system(paste0("mkdir -p ", RobjectDir))
system(paste0("mkdir -p ", plotDir))


### 1. Load in inputs ###

# load in gc counts, append each half for each sample and add into list:
compFiles <- as.character(read.table(paste0(inDir, "/completed_files.txt"))[,1])

# check files have the same number of lines:
for (i in 1:length(compFiles)) {
  if (i==1) {
    lineNo1 <- strsplit(system(paste0("wc -l ", inDir, "/", compFiles[i], "/", compFiles[i], ".gc.htseq.txt"), intern=T), " ")[[1]][4]
    print(paste0(compFiles[i], " has ", lineNo1, " lines"))
    print(paste0("Appending ", compFiles[i], " to goodFilesGC"))
    goodFilesGC <- c(compFiles[i])
    badFilesGC <- c()
  } else {
    lineNo <- strsplit(system(paste0("wc -l ", inDir, "/", compFiles[i], "/", compFiles[i], ".gc.htseq.txt"), intern=T), " ")[[1]][4]
    print(paste0(compFiles[i], " has ", lineNo, " lines"))
    if (lineNo == lineNo1) {
      print(paste0("Appending ", compFiles[i], " to goodFilesGC"))
      goodFilesGC <- append(goodFilesGC, compFiles[i])
    } else {
      print(paste0("Appending ", compFiles[i], " to badFilesGC"))
      badFilesGC <- append(badFilesGC, compFiles[i])
    }
  }
}

# load in gc_countFiles into a df:
for (i in 1:length(compFiles)) {
  if (i==1) {
    print(compFiles[i])
    Counts <- data.frame(read.table(file=paste0(inDir, "/", compFiles[i], "/", compFiles[i], ".gc.htseq.txt")))
    lineNo <- nrow(Counts)
    print(paste0("The number of lines is: ", lineNo))
  } else {
    print(compFiles[i])
    Counts <- cbind(Counts, read.table(file=paste0(inDir, "/", compFiles[i], "/", compFiles[i], ".gc.htseq.txt"))[,2])
  }
}
colnames(Counts) <- c("gene_id", compFiles)

######
saveRDS(Counts, paste0(RobjectDir, "gcCounts_temp.rds"))
######

# remove duplicate rows:
Counts <- unique(Counts)
# remove specs lines:
Counts <- Counts[grep("__", Counts$gene_id, invert=T),]


### 2. Replace uIDs with patient ids:

# load in patient ids:
Key <- read.table(file=paste0(rawDir, "/sampleKey.txt"), sep=" ", fill=T)[,3:4]
Key <- Key[3:nrow(Key),]
Key$V4 <- gsub("_[A-Z].*$", "", gsub("^.*AOCS_", "AOCS_", Key$V4))
Key <- Key[with(Key, order(V3)), ]

# renames duplicate sample names:
Key <- dplyr::filter(Key, Key$V3 %in% colnames(Counts))
newNames <- paste0(Key$V4, "_", gsub("[0-9]", "", colnames(subset(Counts, select=-gene_id))))
newNames[duplicated(newNames)] <- gsub("AOCS", "AOCS2", newNames[duplicated(newNames)])
newNames[duplicated(newNames)] <- gsub("AOCS2", "AOCS3", newNames[duplicated(newNames)])
Key$V4 <- gsub("\\_[a-zA-Z].*$", "", newNames)

# change colnames of Counts:
m <- match(colnames(Counts), Key$V3)
colnames(Counts) <- gsub("NA_", "", 
                         gsub("[0-9]$", "",
                              gsub(
                                "[0-9][0-9]$", "", paste0(Key$V4[m], "_", colnames(Counts))
                              )
                         )
)


## load in patient ids:
#Key <- read.table(file=paste0(rawDir, "/sampleKey.txt"), sep=" ", fill=T)[,3:4]
#Key <- Key[3:nrow(Key),]
#Key$V4 <- gsub("_[A-Z].*$", "", gsub("^.*AOCS_", "AOCS_", Key$V4))
#Key <- Key[with(Key, order(V3)), ]
#
## only include samples in cCounts:
#Key <- dplyr::filter(Key, Key$V3 %in% colnames(Counts))
#newNames <- paste0(Key$V4, "_", gsub("[0-9]", "", colnames(subset(Counts, select=-gene_id))))
#newNames[duplicated(newNames)] <- gsub("AOCS", "AOCS2", newNames[duplicated(newNames)])
#newNames[duplicated(newNames)] <- gsub("AOCS2", "AOCS3", newNames[duplicated(newNames)])
#colnames(Counts)[2:ncol(Counts)] <- newNames

Counts <- apply(Counts[,2:ncol(Counts)], 2, sum)


# save the counts:
if (!file.exists(paste0(RobjectDir, "/gc_allcounts.htseq.rds"))) {
  saveRDS(Counts, file=paste0(RobjectDir, "/gc_allcounts.htseq.rds"))
}
