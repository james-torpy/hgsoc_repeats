### rpCompBarPlot_prepare_cCounts.R ###

# breaks cCounts and tCounts up into dataframe for each class and saves each a RData objects #

# define starting variables:
project <- "hgsoc_repeats"
expName <- "exp5"
sampleName <- c("arPT5", "rcAF6")

# define directories:
homeDir <- "/share/ScratchGeneral/jamtor/"
#homeDir <- "/Users/jamestorpy/clusterHome"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results/")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/", expName, "/")
inDir <- paste0(resultsDir, "/htseq/", expName)
plotDir <- paste0(inDir, "/plots/linePlots/")
rawDir <- paste0(projectDir, "/RNA-seq/raw_files/fullsamples/bowtell_primary/")

system(paste0("mkdir -p ", plotDir))


### 1. Load in inputs ###

for (Type in list("c", "custom3"))
    # load in repeat counts, append each half for eah sample and add into list:
    countFiles <- grep(
      paste0("\\.", Type, "\\."), 
      list.files(inDir, pattern="htseq.txt", 
      full.names = T, recursive=T),
    value=T)

    countFiles <- grep(paste0("subset|all"), countFiles, value=T,
      invert=T)

    countFiles <- countFiles[gsub("\\.c..*$", "",
      basename(countFiles)) %in% sampleName]
    
    countFiles <- grep("gc_", countFiles, value=T, invert=T)
    gc_countFiles <- grep("gc_", countFiles, value=T)
    
    rep_ids <- unique(
      gsub(
        "\\.htseq.txt", "",
        gsub(
          paste0("^.*", Type, "\\."), "", basename(countFiles)
        )
      )
    )

    countFiles <- split(countFiles, c(rep(1, length(grep(sampleName[1], countFiles))), rep(2, length(grep(sampleName[2], countFiles)))))

    # calculate the sums of each repeat type for samples and add to
    # dataframe:
    countDF <- lapply(countFiles, function(x) {
      for (i in 1:length(x)) {
        if (i==1) {
          tab <- read.table(x[i])
          Counts <- data.frame(sum(tab[grep("__", tab$V1, invert=T),]$V2))
          rownames(Counts)[i] <- rep_ids[i]
        } else {
          tab <- read.table(x[i])
          Counts[i,] <- sum(tab[grep("__", tab$V1, invert=T),]$V2)
          rownames(Counts)[i] <- rep_ids[i]
        }
      }
      return(Counts)
    })





    # load in c_countFiles into a df:
    for (i in 1:length(c_countFiles)) {
      # define sample name:
      Name <- sampleName[i]
      if (i==1) {
      	cCounts <- data.frame(read.table(file=c_countFiles[i]))
      } else {
      	cCounts <- cbind(cCounts, read.table(file=c_countFiles[i])[,2])
      }
    }
    colnames(cCounts) <- c("gene_id", sampleName)
    
    
    # load in gc_countFiles into a df:
    for (i in 1:length(gc_countFiles)) {
      # define sampleName:
      sampleName <- sampleNames[i]
      if (i==1) {
      	gcCounts <- data.frame(read.table(file=gc_countFiles[i]))
      } else {
      	gcCounts <- cbind(gcCounts, read.table(file=gc_countFiles[i])[,2])
      }
    }
    colnames(gcCounts) <- c("gene_id", sampleNames)
    
    # aggregate multiple types of same gene:
    gcCounts$gene_id <- gsub("\\..*$", "", gcCounts$gene_id)
    gcCounts <- aggregate(.~gene_id, gcCounts, mean)
    
    # load in patient ids:
    Key <- read.table(file=paste0(rawDir, "/sampleKey.txt"), sep=" ", fill=T)[,3:4]
    Key <- Key[3:nrow(Key),]
    Key$V4 <- gsub("_[A-Z].*$", "", gsub("^.*AOCS_", "AOCS_", Key$V4))
    Key <- Key[with(Key, order(V3)), ]
    
    # only include samples in cCounts:
    Key <- dplyr::filter(Key, Key$V3 %in% colnames(cCounts))
    newNames <- paste0(Key$V4, "_", gsub("[0-9]", "", colnames(subset(cCounts, select=-gene_id))))
    colnames(cCounts)[2:ncol(cCounts)] <- newNames
    colnames(gcCounts)[2:ncol(gcCounts)] <- newNames
    
    
    # save the counts:
    if (!file.exists(paste0(RobjectDir, "/",  expName, "/", Type, "_allcounts.htseq.rds"))) {
      saveRDS(cCounts, file=paste0(RobjectDir, "/",  expName, "/", Type, "_allcounts.htseq.rds"))
    }
    
    if (!file.exists(paste0(RobjectDir, "/",  expName, "/gc_allcounts.htseq.rds"))) {
      saveRDS(gcCounts, file=paste0(RobjectDir, "/",  expName, "/gc_allcounts.htseq.rds"))
    }
}