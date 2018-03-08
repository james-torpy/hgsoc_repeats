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
descrip <- "expression_vs_structural_rearrangements_indivHGSOC"
Type <- "custom3"
metric <- "HGSOC_vs_FT"
repeatList <- c("(CATTC)n", "HSATII", "L1M3a", "L1M3c", "L1M4c", "L1M5", "L1M8", "L1MA5", "L1MC")
repeatPattern <- ""
FDRthresh <- 0.05
removeSamples <- c("pAF", "rcAF")

# define plot colours:
no <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
rawDir <- paste0("/Users/jamestorpy/clusterHome2/projects/hgsoc_repeats/RNA-seq/raw_files")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
       expName, "/")
newRobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
                        expName, "/", descrip, "/")
plotDir <- paste0(resultsDir, "/R/", expName,
                  "/plots/", descrip, "/")

system(paste0("mkdir -p ", plotDir))
system(paste0("mkdir -p ", newRobjectDir))


### 1. Load in DESeq results and create logFC data frame ###

DElist <- readRDS(file=paste0(RobjectDir, "/htseq_single_HGSOCs_vs_FT_DElist.rds"))

# Concatenate all logFC values into data frame
i=1
logFC_list <- lapply(DElist, function(x) {
  df <- data.frame(x$log2FoldChange)
  rownames(df) <- x$id
  colnames(df) <- names(DElist)[i]
  df[df[,1]=="-Inf",] <- 0
  i <<- i+1
  return(df)
})

logFC <- do.call("cbind", logFC_list)

# Concatenate all padj values into data frame
i=1
padj_list <- lapply(DElist, function(x) {
  df <- data.frame(x$padj)
  rownames(df) <- x$id
  colnames(df) <- names(DElist)[i]
  i <<- i+1
  return(df)
})

padj <- do.call("cbind", padj_list)


# Concatenate all FC values into data frame
i=1
FC_list <- lapply(DElist, function(x) {
  df <- data.frame(x$foldChange)
  rownames(df) <- x$id
  colnames(df) <- names(DElist)[i]
  df[df[,1]=="-Inf",] <- 0
  i <<- i+1
  return(df)
})

FC <- do.call("cbind", FC_list)


### 2. Order samples by structural rearrangement number ###

srNo <- read.csv(file=paste0(rawDir, "/somatic mutation counts.csv"), header=T)[,2:6]
srNo$Case.ID <- as.character(srNo$Case.ID)
srNo$Total.number.of.Structural.variant.events.detected <- as.numeric(srNo$Total.number.of.Structural.variant.events.detected)

for (j in 1:nrow(srNo)) {
  if (srNo$Sample.time.point[j] == "primary") {
    if (srNo$Chemotherapy.response[j] == "refractory") {
      srNo$Case.ID[j] <- paste0(srNo$Case.ID[j], "_rfPT")
    } else if (srNo$Chemotherapy.response[j] == "resistant") {
      srNo$Case.ID[j] <- paste0(srNo$Case.ID[j], "_prPT")
    } else if (length(grep(srNo$Case.ID[j], colnames(Counts[j,]))) < 2 & length(grep(srNo$Case.ID[j], colnames(Counts[j,]))) > 0) {
      srNo$Case.ID[j] <- colnames(Counts[j,])[grep(srNo$Case.ID[j], colnames(Counts[j,]))]
    } else if (srNo$Chemotherapy.response[j] == "sensitive") {
      srNo$Case.ID[j] <- paste0(srNo$Case.ID[j], "_arPT")
    }
  } else if (srNo$Sample.time.point[j] == "relapse") {
    if (srNo$Chemotherapy.response[j] == "resistant") {
      srNo$Case.ID[j] <- paste0(srNo$Case.ID[j], "_rcAF")
    } else if (length(grep(srNo$Case.ID[j], colnames(Counts[j,]))) < 2 & length(grep(srNo$Case.ID[j], colnames(Counts[j,]))) > 0) {
      srNo$Case.ID[j] <- colnames(Counts[j,])[grep(srNo$Case.ID[j], colnames(Counts[j,]))]
    }
  } else if (srNo$Sample.time.point[j] == "metastasis") {
    srNo$Case.ID[j] <- paste0(srNo$Case.ID[j], "_msST")
  } else if (srNo$Sample.time.point[j] == "Primary") {
    if (srNo$Chemotherapy.response[j] == "sensitive") {
      srNo$Case.ID[j] <- paste0(srNo$Case.ID[j], "_arPT")
    }
  }
}

srNo <- subset(srNo, select=c("Case.ID", "Total.number.of.Structural.variant.events.detected"))[srNo$Case.ID %in% colnames(Counts),]
colnames(srNo) <- c("sample", "structural_rearrangement_no")

# remove additional column:
srNo <- subset(srNo, select=-X)

# order srNo by number of rearrangements:
srNo <- srNo[order(srNo$structural_rearrangement_no),]

write.csv(srNo, file=paste0(RobjectDir, "/somatic_mutation_counts_formatted.csv"))
#srNo <- read.csv(file=paste0(RobjectDir, "/somatic_mutation_counts_formatted.csv"))

srOrder <- as.character(unique(srNo$sample))

### 3. Order logFC, padj and FC data frames by structural rearrangement
# and perform Pearson correlations for each repeat:

for (k in 1:length(srOrder)) {
  if (k==1) {
    col <- which(grepl(strsplit(srOrder[k], "_")[[1]][2], colnames(logFC)) & 
                   grepl(strsplit(srOrder[k], "_")[[1]][3], colnames(logFC)))
    logFC_ordered <- data.frame(logFC[,col])
    rownames(logFC_ordered) <- rownames(logFC)
    if (length(col)==1) {
      colnames(logFC_ordered)[ncol(logFC_ordered)] <- as.character(paste0(srOrder[k], "_vs_FT"))
    }
  } else {
    col <- which(grepl(strsplit(srOrder[k], "_")[[1]][2], colnames(logFC)) & 
                   grepl(strsplit(srOrder[k], "_")[[1]][3], colnames(logFC)))
    logFC_ordered <- cbind(logFC_ordered, logFC[,col])
    if (length(col)==1) {
      colnames(logFC_ordered)[ncol(logFC_ordered)] <- as.character(paste0(srOrder[k], "_vs_FT"))
    }
  }
}

padj_ordered <- padj[,colnames(logFC_ordered)]
FC_ordered <- FC[,colnames(logFC_ordered)]

for (j in 1:nrow(logFC_ordered)) {
  if (length(grep("ENSG", rpName)) < 1) {
    rp <- as.data.frame(t(logFC_ordered[j,]))
    rp$sample <- factor(rownames(rp), levels = rownames(rp))
    rp$padj <- as.numeric(t(padj_ordered[j,]))
    rp$pThresh <- as.logical((rp$padj < 0.05))
    rp$sig <- NA
    rp$sig[rp$pThresh] <- "< 0.05"
    rp$sig[!(rp$pThresh)] <- "non-significant"

    rp$FC <- as.numeric(t(FC_ordered[j,]))
    rp$str_no <- srNo$structural_rearrangement_no
    
    colnames(rp) <- c("log2FC", "sample", "padj", "pThresh", "sig", "FC", "str_no")
    if (j==1) {
      rpList <- list(rp)
    } else {
      rpList[[j]] <- rp
    }
  }
}
rpList <- rpList[grep("ENSG", names(rpList), invert=T)]

names(rpList) <- rownames(logFC_ordered)[1:length(rpList)]

saveRDS(rpList, file=paste0(RobjectDir, "/rpList.rds"))

# check if repeat FC distribution is normal:
#allFC <- lapply(rpList, function(x) {
#  return(x$FC)
#})
#allFC <- unlist(allFC)
#repFC <- allFC[grep("ENSG", names(allFC), invert=T)]

#dFC <- density(repFC)
#plot(dFC)

m=1
corList <- lapply(rpList, function(x) {
  df <- data.frame(x$str_no, x$FC)
  df <- df[rpList[[m]]$sig == "< 0.05",]
  m <<- m+1
  return(cor(df, method = "spearman")[1,2])
})
names(corList) <- names(rpList)
cors <- na.omit(unlist(corList))

rpInt <- rpList[names(cors)]


### 4. Plot logFC vs structural rearrangement for each repeat ###

for (k in 1:length(rpInt)) {
  rpName <- gsub("\\(|\\)", "", names(rpInt)[k])
  pdf(file = paste0(plotDir, "/", rpName, ".pdf"), width=10, height=10)
  p <- ggplot(rpInt[[k]], aes(x = factor(rpInt[[k]]$sample), y = rpInt[[k]]$log2FC))
  p <- p + geom_bar(stat="identity", position = "dodge", aes(fill = rpInt[[k]]$sig))
  p <- p + ylab("log2FC")
  p <- p + xlab("sample")
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.text = "significance")
  p <- p + annotate(cors[k])
  print(p)
  dev.off()
}






