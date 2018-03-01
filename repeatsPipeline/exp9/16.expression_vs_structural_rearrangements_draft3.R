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
descrip <- "expression_vs_structural_rearrangements"
Type <- "custom3"
metric <- "HGSOC_vs_FT"
#metric <- "CPMR"
repeatList <- c("(CATTC)n", "HSATII", "L1M3a", "L1M3c", "L1M4c", "L1M5", "L1M8", "L1MA5", "L1MC")
#repeatList <- ""
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


# create a sampleKey to check annotation of structural rearrangements df:
#sampleKey <- read.table(file = paste0(homeDir, "/projects/hgsoc_repeats/RNA-seq/raw_files/fullsamples/bowtell_primary/sampleKey.txt"), fill = T)

#sampleKey$V4 <- gsub(
#  "^.*AOCS", "AOCS", gsub(
#    "_ICGC.*$", "", gsub(
#      "_EXTERNA.*$", "", sampleKey$V4
#    )
#  )
#)

#sampleKey <- sampleKey[3:122, 3:4]

#sampleKey <- sort(paste0(sampleKey$V4, "_", sampleKey$V3))


### 1. Load in all counts and structural rearrangement counts ###

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

write.csv(srNo, file=paste0(RobjectDir, "/somatic_mutation_counts_formatted.csv"))
#srNo <- read.csv(file=paste0(RobjectDir, "/somatic_mutation_counts_formatted.csv"))

intRepeats <- readRDS(paste0(interesting_fcRobjectDir, "/interesting_fc.rds"))

if (metric == "CPMR") {
 
  ### 2. Calculate CPMR ###
  
  # calculate total repeat count size:
  rSizes <- apply(Counts, 2, sum)
  
  # calculate CPMRs
  CPMR <- as.data.frame(t(t(Counts)/rSizes)*1000000)
  
  
  ### 3. Isolate interesting repeats, add structural rearrangement info and plot ###
  
  CPMR <- CPMR[intRepeats,]
  
  for (i in 1:nrow(CPMR)) {
    print(i)
    repID <- gsub("\\(|\\)", "", rownames(CPMR)[i])
    
    repCPMR <- as.data.frame(t(CPMR[i,]))
    repCPMR$sample <- rownames(repCPMR)
    colnames(repCPMR)[1] <- "CPMR"
    
    str_vs_rep <- merge(repCPMR, srNo)
    str_vs_rep <- str_vs_rep[order(str_vs_rep$structural_rearrangement_no),]
    #str_vs_rep$sample <- paste0(str_vs_rep$sample, " (", str_vs_rep$structural_rearrangement_no, " structural rearrangements)")
    
    vec <- str_vs_rep$CPMR
    names(vec) <- str_vs_rep$sample
    
    pdf(file = paste0(plotDir, "/", repID, "_", metric, "_", descrip, ".pdf"), width=10, height=10)
    barplot(vec, cex.names=0.5, las=2)
    dev.off()
  }
  
} else if (metric == "HGSOC_vs_FT") {
  
  
  ### 2. Load in differential expression of repeats and fetch fold change and FDR values ###
  
  DEsigReps <- readRDS(file = paste0(RobjectDir, "/custom3_DEsigReps.rds"))
  
  # fetch fold change columns, and FDR columns, of each data frame for each HGSOC vs FT comparison
  # and create new data frame for each:
  for (k in 1:length(DEsigReps)) {
    if (k==1) {
      FCreps <- as.data.frame(DEsigReps[[k]]$logFC, drop = F, row.names = rownames(DEsigReps[[k]]))
      FDRreps <- as.data.frame(DEsigReps[[k]]$FDR, drop = F, row.names = rownames(DEsigReps[[k]]))
    } else {
      FCreps[,k] <- DEsigReps[[k]]$logFC
      FDRreps[,k] <- DEsigReps[[k]]$FDR
    }
  }
  
  # name the columns of the data frame according to comparison, and change "FT_vs" label to "_sv_FT" if necessary:
  for (n in 1:ncol(FCreps)) {
    if (startsWith(names(DEsigReps)[n], "FT_vs_")) {
      colnames(FCreps)[n] <- paste0(strsplit(names(DEsigReps[n]), "_")[[1]][3], "_vs_", strsplit(names(DEsigReps[n[]]), "_")[[1]][1])
    } else {
      colnames(FCreps) <- names(DEsigReps)
    }
  }
  
  for (n in 1:ncol(FDRreps)) {
    if (startsWith(names(DEsigReps)[n], "FT_vs_")) {
      colnames(FDRreps)[n] <- paste0(strsplit(names(DEsigReps[n]), "_")[[1]][3], "_vs_", strsplit(names(DEsigReps[n[]]), "_")[[1]][1])
    } else {
      colnames(FDRreps) <- names(DEsigReps)
    }
  }
 
  
  ### 3. Isolate interesting repeats, add structural rearrangement info and fetch
  # FDR for each DE value ###
  
  FCreps <- na.omit(FCreps[intRepeats,])
  FDRreps <- na.omit(FCreps[intRepeats,])

  # calculate the mean structural rearrangements for each HGSOC type:
  srNo$sample <- gsub("^.*\\_", "", srNo$sample)
  srNo <- aggregate(srNo$structural_rearrangement_no, list(srNo$sample), mean)
  colnames(srNo) <- c("sample", "mean_rearrangements")
  
  # add "_vs_FT" onto subtype names for the purpose of merging to FCreps:
  srNo$sample <- paste0(srNo$sample, "_vs_FT")
  
  colnam <- colnames(FCreps)
  
  if ((repeatList != "")[1]) {
    FCreps <- FCreps[repeatList,]
    colnames(FCreps) <- colnam
    FCreps$repeat_id <- rownames(FCreps)
    
    # remove samples not needed:
    toRemove <- paste0(removeSamples, "_vs_FT")
    FCreps <- FCreps[,!colnames(FCreps) %in% toRemove]
    srNo <- srNo[!srNo$sample %in% toRemove,]
    
    # order by increasing no structural rearrangements:
    FCreps <- FCreps[,srNo$sample[order(srNo$mean_rearrangements)]]
    
    # add repeat_id column for melting:
    FCreps$repeat_id <- rownames(FCreps)
    
    # prepare for plotting:
    FCrepsM <- melt(FCreps)
    FCrepsM$variable <- factor(FCrepsM$variable, levels = srNo$sample[order(srNo$mean_rearrangements)])
                                  
    # prepare FDRs:
    #FDRreps <- FCreps[repeatList,]
    #colnames(FDRreps) <- colnam
    
    # plot barplot with grouping according to DE comparison, one bars per repeat:
    p <- ggplot(FCrepsM, aes(factor(variable), value, fill = repeat_id))
    p <- p + geom_bar(stat="identity", position = "dodge")
    p <- p + scale_fill_manual(values=col_vector[1:nrow(FCreps)])
    p <- p + ylab("log2FC")
    pdf(file = paste0(plotDir, "/", metric, "_", descrip, "_selected_repeats.pdf"), width=10, height=10)
    p
    dev.off()
    
  } else if (repeatPattern != "") {
    
    
  } else {
    for (i in 1:nrow(FCreps)) {
      print(i)
      repID <- gsub("\\(|\\)", "", rownames(FCreps)[i])
      
      fcDE <- as.data.frame(t(FCreps[i,]))
      fcDE$sample <- rownames(fcDE)
      colnames(fcDE)[1] <- "FC"
      
      # remove pAF:
      fcDE <- fcDE[grep("pAF", rownames(fcDE), invert=T),]
      
      str_vs_rep <- merge(fcDE, srNo)
      str_vs_rep <- str_vs_rep[order(str_vs_rep$mean_rearrangements),]
      #str_vs_rep$sample <- paste0(str_vs_rep$sample, " (", str_vs_rep$structural_rearrangement_no, " structural rearrangements)")
      
      vec <- str_vs_rep$FC
      names(vec) <- str_vs_rep$sample
      
      par(mar=c(4, 4, 4, 4))
      pdf(file = paste0(plotDir, "/", repID, "_", metric, "_", descrip, ".pdf"), width=10, height=10)
      barplot(vec, cex.names=0.5, las=2)
      dev.off()
    }
  }
}





