### 10.prepareCounts.R ###

# breaks gcCounts up into dataframe for each class and
# saves each as RData objects #

# define starting variables:
project <- "hgsoc_repeats"
expName <- "rsem_no_mmappers"

# define directories:
#homeDir <- "/share/ScratchGeneral/jamtor/"
homeDir <- "/Users/jamestorpy/clusterHome"
projectDir <- paste0(homeDir, "/projects/", 
  project, "/")
resultsDir <- paste0(projectDir, "/RNA-seq/results/")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/", 
  expName, "/")
inDir <- paste0(resultsDir, "/htseq/", expName)
plotDir <- paste0(inDir, "/plots/linePlots/")
rawDir <- paste0(projectDir, 
  "/RNA-seq/raw_files/fullsamples/bowtell_primary/")
starGC_dir <- paste0(resultsDir, "/star/GC/", expName)
starRibo_dir <- paste0(resultsDir, "/star/ribo/", 
  expName)

system(paste0("mkdir -p ", RobjectDir))
system(paste0("mkdir -p ", plotDir))


### 1. Check for and discard bad inputs ###

# fetch counts file ids:
compFiles <- as.character(read.table(paste0(inDir, 
  "/completed_files.txt"))[,1])

# check for presence of STAR log file for GC mapping:
for (i in 1:length(compFiles)) {
  logFile <- paste0(starGC_dir, "/", compFiles[i],
    "/Log.final.out")
  if (i==1) {
    if (file.exists(logFile)) {
      print(paste0(logFile, " exists"))
      print(paste0("Appending ", logFile,
        " to goodFiles_GClog"))
      goodFiles_GClog <- c(compFiles[i])
      badFiles_GClog <- c()
    } else {
      print(paste0(logFile, " does not exist"))
      print(paste0("Appending ", logFile,
        " to badFiles_GClog"))
      badFiles_GClog <- c(compFiles[i])
      goodFiles_GClog <- c()
    }
  } else {
    if (file.exists(logFile)) {
      print(paste0(logFile, " exists"))
      print(paste0("Appending ", logFile, 
        " to goodFiles_GClog"))
      goodFiles_GClog <- append(goodFiles_GClog, 
        compFiles[i])
    } else {
      print(paste0(logFile, " does not exist"))
      print(paste0("Appending ", logFile, 
        " to badFiles_GClog"))
      badFiles_GClog <- append(badFiles_GClog, 
        compFiles[i])
    }
  }
}

# check for presence of STAR log file for ribo mapping:
for (i in 1:length(goodFiles_GClog)) {
  logFile <- paste0(starRibo_dir, "/", goodFiles_GClog[i],
    "/Log.final.out")
  if (i==1) {
    if (file.exists(logFile)) {
      print(paste0(logFile, " exists"))
      print(paste0("Appending ", logFile,
        " to goodFiles_riboLog"))
      goodFiles_riboLog <- c(goodFiles_GClog[i])
      badFiles_riboLog <- c()
    } else {
      print(paste0(logFile, " does not exist"))
      print(paste0("Appending ", logFile,
        " to badFiles_riboLog"))
      badFiles_riboLog <- c(goodFiles_GClog[i])
      goodFiles_riboLog <- c()
    }
  } else {
    if (file.exists(logFile)) {
      print(paste0(logFile, " exists"))
      print(paste0("Appending ", logFile, 
        " to goodFiles_riboLog"))
      goodFiles_riboLog <- append(goodFiles_riboLog, 
        goodFiles_GClog[i])
    } else {
      print(paste0(logFile, " does not exist"))
      print(paste0("Appending ", logFile, 
        " to badFiles_riboLog"))
      badFiles_riboLog <- append(badFiles_riboLog, 
        goodFiles_GClog[i])
    }
  }
}

# check gc files have the same number of lines and put
# these in a vector:
for (i in 1:length(goodFiles_riboLog)) {
  if (i==1) {
    lineNo1 <- strsplit(system(paste0("wc -l ", inDir, 
      "/", goodFiles_riboLog[i], "/", 
      goodFiles_riboLog[i], ".gc.htseq.txt"), intern=T), 
    " ")[[1]][1]
    print(paste0(goodFiles_riboLog[i], " has ", lineNo1, 
      " lines"))
    print(paste0("Appending ", goodFiles_riboLog[i], 
      " to goodFilesGC"))
    goodFilesGC <- c(goodFiles_riboLog[i])
    badFilesGC <- c()
  } else {
    lineNo <- strsplit(system(paste0("wc -l ", inDir, 
      "/", goodFiles_riboLog[i], "/", 
      goodFiles_riboLog[i], 
      ".gc.htseq.txt"), intern=T), " ")[[1]][1]
    print(paste0(goodFiles_riboLog[i], " has ", lineNo, 
      " lines"))
    if (lineNo == lineNo1) {
      print(paste0("Appending ", goodFiles_riboLog[i], 
        " to goodFilesGC"))
      goodFilesGC <- append(goodFilesGC, 
        goodFiles_riboLog[i])
    } else {
      print(paste0("Appending ", goodFiles_riboLog[i], 
        " to \
        badFilesGC"))
      badFilesGC <- append(badFilesGC, 
        goodFiles_riboLog[i])
    }
  }
}

# check allfiles corresponding with ids of goodFilesGC have the same number of lines and put this into a vector:
for (i in 1:length(goodFilesGC)) {
  if (i==1) {
    lineNo1 <- strsplit(system(paste0("wc -l ", inDir, "/", goodFilesGC[i], "/", goodFilesGC[i], ".all.htseq.txt"), intern=T), " ")[[1]][1]
    print(paste0(goodFilesGC[i], " has ", lineNo1, " lines"))
    print(paste0("Appending ", goodFilesGC[i], " to goodFilesAll"))
    goodFilesAll <- c(goodFilesGC[i])
    badFilesAll <- c()
  } else {
    lineNo <- strsplit(system(paste0("wc -l ", inDir, "/", goodFilesGC[i], "/", goodFilesGC[i], ".all.htseq.txt"), intern=T), " ")[[1]][1]
    print(paste0(goodFilesGC[i], " has ", lineNo, " lines"))
    if (lineNo == lineNo1) {
      print(paste0("Appending ", goodFilesGC[i], " to goodFilesAll"))
      goodFilesAll <- append(goodFilesAll, goodFilesGC[i])
    } else {
      print(paste0("Appending ", goodFilesGC[i], " to badFilesAll"))
      badFilesAll <- append(badFilesAll, goodFilesGC[i])
    }
  }
}

### remove FT1 from the list as custom3 results have too few lines ###
#goodFilesAll <- goodFilesAll[2:length(goodFilesAll)]

# check custom3files corresponding with ids of goodFilesAll have the same number of lines and put this into a vector:
for (i in 1:length(goodFilesAll)) {
  if (i==1) {
    lineNo1 <- strsplit(system(paste0("wc -l ", inDir, 
      "/", goodFilesAll[i], "/", goodFilesAll[i], 
      ".custom3.htseq.txt"), intern=T), " ")[[1]][1]
    print(paste0(goodFilesAll[i], " has ", lineNo1, 
      " lines"))
    print(paste0("Appending ", goodFilesAll[i], 
      " to goodFilesCustom3"))
    goodFilesCustom3 <- c(goodFilesAll[i])
    badFilesCustom3 <- c()
  } else {
    lineNo <- strsplit(system(paste0("wc -l ", inDir, 
      "/", goodFilesAll[i], "/", goodFilesAll[i], 
      ".custom3.htseq.txt"), intern=T), " ")[[1]][1]
    print(paste0(goodFilesAll[i], " has ", lineNo, 
      " lines"))
    if (lineNo == lineNo1) {
      print(paste0("Appending ", goodFilesAll[i], 
        " to goodFilesCustom3"))
      goodFilesCustom3 <- append(goodFilesCustom3, 
        goodFilesAll[i])
    } else {
      print(paste0("Appending ", goodFilesAll[i], 
        " to badFilesCustom3"))
      badFilesCustom3 <- append(badFilesCustom3, 
        goodFilesAll[i])
    }
  }
}

goodFiles <- goodFilesCustom3

saveRDS(goodFiles, paste0(RobjectDir, "/goodFiles.rds"))


### 2. Load in GC inputs ###

# load in gc_countFiles into a df:
for (i in 1:length(goodFiles)) {
  if (i==1) {
    print(goodFiles[i])
    CountsGC <- data.frame(read.table(file=paste0(inDir, "/", goodFiles[i], "/", goodFiles[i], ".gc.htseq.txt")))
  } else {
    print(goodFiles[i])
    CountsGC <- cbind(CountsGC, read.table(file=paste0(inDir, "/", goodFiles[i], "/", goodFiles[i], ".gc.htseq.txt"))[,2])
  }
}
colnames(CountsGC) <- c("gene_id", goodFiles)

######
saveRDS(CountsGC, paste0(RobjectDir, "gcCounts_temp.rds"))
######

# aggregate multiple types of same gene:
CountsGC$gene_id <- gsub("\\..*$", "", CountsGC$gene_id)
CountsGC <- aggregate(.~gene_id, CountsGC, mean)
# remove duplicate rows:
CountsGC <- unique(CountsGC)
# remove specs lines:
CountsGC <- CountsGC[grep("__", CountsGC$gene_id, invert=T),]


# replace uIDs with patient ids:

# load in patient ids:
Key <- read.table(file=paste0(rawDir, "/sampleKey.txt"), sep=" ", fill=T)[,3:4]
Key <- Key[3:nrow(Key),]
Key$V4 <- gsub("_[A-Z].*$", "", gsub("^.*AOCS_", "AOCS_", Key$V4))
Key <- Key[with(Key, order(V3)), ]

# renames duplicate sample names:
Key <- dplyr::filter(Key, Key$V3 %in% colnames(CountsGC))
newNames <- paste0(Key$V4, "_", gsub("[0-9]", "", colnames(subset(CountsGC, select=-gene_id))))
newNames[duplicated(newNames)] <- gsub("AOCS", "AOCS2", newNames[duplicated(newNames)])
newNames[duplicated(newNames)] <- gsub("AOCS2", "AOCS3", newNames[duplicated(newNames)])
Key$V4 <- gsub("\\_[a-zA-Z].*$", "", newNames)

# change colnames of CountsGC:
m <- match(colnames(CountsGC), Key$V3)
colnames(CountsGC) <- gsub("NA_", "", 
                         gsub("[0-9]$", "",
                              gsub(
                                "[0-9][0-9]$", "", paste0(Key$V4[m], "_", colnames(CountsGC))
                              )
                         )
)


# save the counts:
if (!file.exists(paste0(RobjectDir, "/gc_allcounts.htseq.rds"))) {
  saveRDS(CountsGC, file=paste0(RobjectDir, "/gc_allcounts.htseq.rds"))
}


### 3. Load in 'all' inputs ###

# load in all_countFiles into a df:
for (i in 1:length(goodFiles)) {
  if (i==1) {
    print(goodFiles[i])
    CountsAll <- data.frame(read.table(file=paste0(inDir, "/", goodFiles[i], "/", goodFiles[i], ".all.htseq.txt")))
  } else {
    print(goodFiles[i])
    CountsAll <- cbind(CountsAll, read.table(file=paste0(inDir, "/", goodFiles[i], "/", goodFiles[i], ".all.htseq.txt"))[,2])
  }
}
colnames(CountsAll) <- c("gene_id", goodFiles)

######
saveRDS(CountsAll, paste0(RobjectDir, "allCounts_temp.rds"))
######

# remove duplicate rows:
CountsAll <- unique(CountsAll)
# remove specs lines:
CountsAll <- CountsAll[grep("__", CountsAll$gene_id, invert=T),]

# replace uIDs with patient ids:

# load in patient ids:
Key <- read.table(file=paste0(rawDir, "/sampleKey.txt"), sep=" ", fill=T)[,3:4]
Key <- Key[3:nrow(Key),]
Key$V4 <- gsub("_[A-Z].*$", "", gsub("^.*AOCS_", "AOCS_", Key$V4))
Key <- Key[with(Key, order(V3)), ]

# renames duplicate sample names:
Key <- dplyr::filter(Key, Key$V3 %in% colnames(CountsAll))
newNames <- paste0(Key$V4, "_", gsub("[0-9]", "", colnames(subset(CountsAll, select=-gene_id))))
newNames[duplicated(newNames)] <- gsub("AOCS", "AOCS2", newNames[duplicated(newNames)])
newNames[duplicated(newNames)] <- gsub("AOCS2", "AOCS3", newNames[duplicated(newNames)])
Key$V4 <- gsub("\\_[a-zA-Z].*$", "", newNames)

# change colnames of CountsAll:
m <- match(colnames(CountsAll), Key$V3)
colnames(CountsAll) <- gsub("NA_", "", 
  gsub("[0-9]$", "",
    gsub(
      "[0-9][0-9]$", "", paste0(Key$V4[m], "_", colnames(CountsAll))
    )
  )
)

CountsAll <- apply(CountsAll[,2:ncol(CountsAll)], 2, sum)

# save the counts:
if (!file.exists(paste0(RobjectDir, "/all_allcounts.htseq.rds"))) {
  saveRDS(CountsAll, file=paste0(RobjectDir, "/all_allcounts.htseq.rds"))
}


### 4. Load in 'custom3' inputs ###

# load in custom3_countFiles into a df:
for (i in 1:length(goodFiles)) {
  if (i==1) {
    print(goodFiles[i])
    CountsCustom3 <- data.frame(read.table(file=paste0(inDir, "/", goodFiles[i], "/", goodFiles[i], ".custom3.htseq.txt")))
  } else {
    print(goodFiles[i])
    CountsCustom3 <- cbind(CountsCustom3, read.table(file=paste0(inDir, "/", goodFiles[i], "/", goodFiles[i], ".custom3.htseq.txt"))[,2])
  }
}
colnames(CountsCustom3) <- c("gene_id", goodFiles)

# remove duplicate rows:
CountsCustom3 <- unique(CountsCustom3)
# remove specs lines:
CountsCustom3 <- CountsCustom3[grep("__", CountsCustom3$gene_id, invert=T),]


### replace uIDs with patient ids:

# load in patient ids:
Key <- read.table(file=paste0(rawDir, "/sampleKey.txt"), sep=" ", fill=T)[,3:4]
Key <- Key[3:nrow(Key),]
Key$V4 <- gsub("_[A-Z].*$", "", gsub("^.*AOCS_", "AOCS_", Key$V4))
Key <- Key[with(Key, order(V3)), ]

# renames duplicate sample names:
Key <- dplyr::filter(Key, Key$V3 %in% colnames(CountsCustom3))
newNames <- paste0(Key$V4, "_", gsub("[0-9]", "", colnames(subset(CountsCustom3, select=-gene_id))))
newNames[duplicated(newNames)] <- gsub("AOCS", "AOCS2", newNames[duplicated(newNames)])
newNames[duplicated(newNames)] <- gsub("AOCS2", "AOCS3", newNames[duplicated(newNames)])
Key$V4 <- gsub("\\_[a-zA-Z].*$", "", newNames)

# change colnames of CountsCustom3:
m <- match(colnames(CountsCustom3), Key$V3)
colnames(CountsCustom3) <- gsub("NA_", "", 
  gsub("[0-9]$", "",
    gsub(
      "[0-9][0-9]$", "", paste0(Key$V4[m], "_", colnames(CountsCustom3))
    )
  )
)

# save the counts:
if (!file.exists(paste0(RobjectDir, "/custom3_allcounts.htseq.rds"))) {
  saveRDS(CountsCustom3, file=paste0(RobjectDir, "/custom3_allcounts.htseq.rds"))
}