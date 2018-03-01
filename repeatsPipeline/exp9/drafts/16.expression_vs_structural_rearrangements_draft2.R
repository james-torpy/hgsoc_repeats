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
                  "/plots/", descrip, "/")

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
srNo$Total.number.of.Structural.variant.events.detected <- as.numeric(srNo$Total.number.of.Structural.variant.events.detected)

for (j in 1:nrow(srNo)) {
  if (srNo$Sample.time.point[j] == "primary") {
    if (srNo$Chemotherapy.response[j] == "refractory") {
      srNo$Case.ID[j] <- paste0(srNo$Case.ID[j], "_rfPT")
    } else if (srNo$Chemotherapy.response[j] == "resistant") {
      srNo$Case.ID[j] <- paste0(srNo$Case.ID[j], "_prPT")
    } else if (length(grep(srNo$Case.ID[j], CPMR[j,])) < 2 & length(grep(srNo$Case.ID[j], CPMR[j,])) > 0) {
      srNo$Case.ID[j] <- CPMR[j,][grep(srNo$Case.ID[j], CPMR[j,])]
    } else if (srNo$Chemotherapy.response[j] == "sensitive") {
      srNo$Case.ID[j] <- paste0(srNo$Case.ID[j], "_arPT")
    }
  } else if (srNo$Sample.time.point[j] == "relapse") {
    if (srNo$Chemotherapy.response[j] == "resistant") {
      srNo$Case.ID[j] <- paste0(srNo$Case.ID[j], "_rcAF")
    } else if (length(grep(srNo$Case.ID[j], CPMR[j,])) < 2 & length(grep(srNo$Case.ID[j], CPMR[j,])) > 0) {
        srNo$Case.ID[j] <- CPMR[j,][grep(srNo$Case.ID[j], CPMR[j,])]
    }
  } else if (srNo$Sample.time.point[j] == "metastasis") {
    srNo$Case.ID[j] <- paste0(srNo$Case.ID[j], "_msST")
  } else if (srNo$Sample.time.point[j] == "Primary") {
    if (srNo$Chemotherapy.response[j] == "sensitive") {
      srNo$Case.ID[j] <- paste0(srNo$Case.ID[j], "_arPT")
    }
  }
}

srNo <- subset(srNo, select=c("Case.ID", "Total.number.of.Structural.variant.events.detected"))[srNo$Case.ID %in% names(CPMR[1,]),]
colnames(srNo) <- c("sample", "structural_rearrangement_no")

write.csv(srNo, file=paste0(RobjectDir, "/somatic_mutation_counts_formatted.csv"))

repCPMR <- as.data.frame(t(CPMR[1,]))
repCPMR$sample <- rownames(repCPMR)
colnames(repCPMR)[1] <- "CPMR"

str_vs_rep <- merge(repCPMR, srNo)
str_vs_rep <- str_vs_rep[order(str_vs_rep$structural_rearrangement_no),]
#str_vs_rep$sample <- paste0(str_vs_rep$sample, " (", str_vs_rep$structural_rearrangement_no, " structural rearrangements)")

vec <- str_vs_rep$CPMR
names(vec) <- str_vs_rep$sample
barplot(vec, cex.names=0.5, las=2)



#p <- ggplot(data=str_vs_rep, aes(x=sample, y=CPMR))
#p <- p + geom_bar(stat = "identity")

#pdf(file = paste0(plotDir, "CATTC.pdf"), width=10, height=30)
#p
#dev.off()
