### 13.DE_FT_vs_HGSOC.R ###

# This script takes a list of dfs of different classes and types of repeat counts for
# HGSOC and FT control RNA-seq data sets and performs DE analysis:

### 0. Define variables/paths ###

# load packages needed:
library(tibble)
library(dplyr)
library(RColorBrewer)
library(RUVSeq)
library(ggplot2)
library(ggrepel)
library(preprocessCore)

# define starting variables:
project <- "hgsoc_repeats"
expName <- "exp9"
descrip <- "CPMR_vs_structural_rearrangements"
Type <- "custom3"

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project)
rawDir <- paste0("/Users/jamestorpy/clusterHome2/projects/hgsoc_repeats/RNA-seq/raw_files")
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
                     expName, "/")
interesting_fcRobjectDir <- paste0(RobjectDir, "/htseq_hgsoc_split_by_driver_mutation_no_ascites_or_metastatic_vs_FT")
newRobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
                        expName, "/", descrip, "/")
plotDir <- paste0(resultsDir, "/R/", expName,
                  "/plots/DEplots/", descrip, "/")

system(paste0("mkdir -p ", plotDir))
system(paste0("mkdir -p ", newRobjectDir))


### 1. Load in all counts ###

custom3Counts <- readRDS(paste0(RobjectDir, "/", Type,
                                "_allcounts.htseq.rds"))
gcCounts <- readRDS(paste0(RobjectDir,
                           "/gc_allcounts.htseq.rds"))

# append gcCounts to custom3Counts:
Counts <- rbind(custom3Counts, gcCounts)

# make rownames gene_id, get rid of latter column and change
# storage mode from factor to integer:
rownames(Counts) <- Counts$gene_id
Counts <- subset(Counts, select=-gene_id)
Counts <- Counts[,grep("FT", colnames(Counts), invert=T)]


### 2. Calculate CPMR ###

# calculate total repeat count size:
rSizes <- apply(Counts, 2, sum)

# calculate CPMRs
CPMR <- as.data.frame(t(t(Counts)/rSizes)*1000000)


### 3. Isolate interesting repeats and add structural rearrangement info to each sample ###

intRepeats <- readRDS(paste0(interesting_fcRobjectDir, "/interesting_fc.rds"))
CPMR <- CPMR[intRepeats,]

srNo <- read.csv(file=paste0(rawDir, "/somatic mutation counts.csv"), header=T)[,2:6]
srNo$Case.ID <- as.character(srNo$Case.ID)

for (j in 1:nrow(srNo)) {
  if (srNo$Sample.time.point[j] == "primary") {
    if (srNo$Chemotherapy.response[j] == "refractory") {
      srNo$Case.ID[j] <- paste0(srNo$Case.ID[j], "_rfPT")
    } else if (srNo$Chemotherapy.response[j] == "resistant") {
      srNo$Case.ID[j] <- paste0(srNo$Case.ID[j], "_prPT")
    } else if (length(grep(srNo$Case.ID[j], names(CATTC_df))) < 2 & length(grep(srNo$Case.ID[j], names(CATTC_df))) > 0) {
      srNo$Case.ID[j] <- names(CATTC_df)[grep(srNo$Case.ID[j], names(CATTC_df))]
    } else if (srNo$Chemotherapy.response[j] == "sensitive") {
      srNo$Case.ID[j] <- paste0(srNo$Case.ID[j], "_arPT")
    }
  } else if (srNo$Sample.time.point[j] == "relapse") {
    if (srNo$Chemotherapy.response[j] == "resistant") {
      srNo$Case.ID[j] <- paste0(srNo$Case.ID[j], "_rcAF")
    } else if (length(grep(srNo$Case.ID[j], names(CATTC_df))) < 2 & length(grep(srNo$Case.ID[j], names(CATTC_df))) > 0) {
        srNo$Case.ID[j] <- names(CATTC_df)[grep(srNo$Case.ID[j], names(CATTC_df))]
    }
  } else if (srNo$Sample.time.point[j] == "metastasis") {
    srNo$Case.ID[j] <- paste0(srNo$Case.ID[j], "_msST")
  } else if (srNo$Sample.time.point[j] == "Primary") {
    if (srNo$Chemotherapy.response[j] == "sensitive") {
      srNo$Case.ID[j] <- paste0(srNo$Case.ID[j], "_arPT")
    }
  }
}

srNo <- srNo[srNo$Case.ID %in% names(CATTC_df),]



CATTCdf <- CPMR[1,]

for (j in 1:length(CATTC_df)) {
  if (j==1) {
    CATTCsrNo <- grep(
      gsub(
      "\\_[a-z].*[A-Z][A-Z]", "", names(CATTC_df)[j]
      ), srNo$Case.ID
    )
    
    #CATTCsrNo <- c(srNo[grep(srNo$Case.ID[j], names(CATTC_row)), 2])
  }
  CATTCsrNo <- srNo[grep(id, names(CATTC_row)), 2]
}




for (i in nrow(CPMR)) {
  
  
  if (i==1) {
    repeatDF <- list()
  }
}





