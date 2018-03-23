### 13.DE_FT_vs_HGSOC.R ###

# This script takes a list of dfs of different classes and types of repeat counts for
# HGSOC and FT control RNA-seq data sets and performs DE analysis:

# Run on cluster with:
#briR
#qsub -N exp9EdgeR -b y -wd \
#/share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/logs/exp9/DE \
#-j y -R y -pe smp 2 -V "Rscript /share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/exp9/13.DE_master.R"

### 0. Define variables/paths ###

# load packages needed:
library(tibble)
library(dplyr)
library(RColorBrewer)
library(RUVSeq)
library(ggplot2)
library(ggrepel)
library(preprocessCore)
library(edgeR)
library(org.Hs.eg.db)

# define starting variables:
project <- "hgsoc_repeats"
expName <- "exp9"

Type <- "custom3"
sTypes <- c("FT", "primary_resistant", "acquired_resistant", "drug_responders", 
  "recurrent_ascites", "metastatic")
#sTypes <- c("HRD", "CCNEamp")
descrip <- "htseq_hgsoc_split_by_drug_response_acquired_resistant_vs_rest"

# choose analysis type:
#analysis <- "edger"
analysis <- "edger_eda"
#analysis <- "DESeq"
#analysis <- "DESeq_eda"

# define sample groups to compare:

sGroups <- list(c("prPT", "rfPT"), c("arPT"), c("mrPT", "erPT"), "rcAF", "msST")
sGroups <- list(c("primary_resistant", "acquired_resistant", "drug_responders", "recurrent_ascites", "metastatic"))
#sGroups <- c("HRD", "CCNEamp")
names(sGroups) <- sTypes
allHGSOC_vs_FT <- FALSE
cat_by_driver <- FALSE

# define sample group to use as control:
ctl <- "acquired_resistant"

# specify what combination of repeat genes (repeats), epigenetic modulators (epiMods),
# RNAi genes (RNAi) and protein-coding genes (pCoding) should contribute to the results:
#resultTypes <- c("repeats", "epiMods")
resultTypes <- c("all", "repeats")

# specify what FDR and log2 fold change thresholds to use:
FDRthresh <- 0.05
FCthresh <- 0

# specify control genes to include:
posGeneIDs <- c("ENSG00000111640", "ENSG00000196776")
posGeneNames <- c("GAPDH", "CD47")
negGeneIDs <- c("ENSG00000075624", "ENSG00000169919")
negGeneNames <- c("beta-actin", "GUSB")

# specify epigenetic modifier genes to include if necessary:
epiIDs <- c("ENSG00000130816", "ENSG00000119772", "ENSG00000088305", "ENSG00000276043", 
            "ENSG00000138336", "ENSG00000168769", "ENSG00000187605", "ENSG00000101945")

epiSym <- c("DNMT1", "DNMT3A", "DNMT3B", "UHRF1", "TET1", "TET2", "TET3", "SUV39H1")

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project)
rawDir <- paste0("/Users/jamestorpy/clusterHome2/projects/hgsoc_repeats/RNA-seq/raw_files")
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
                     expName, "/")
newRobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
                        expName, "/", descrip, "/")
plotDir <- paste0(resultsDir, "/R/", expName,
                  "/plots/DEplots/", descrip, "/")

system(paste0("mkdir -p ", plotDir))
system(paste0("mkdir -p ", newRobjectDir))


### 1. Load in all counts ###

if (!file.exists(paste0(newRobjectDir, "/", Type, "_counts.RData"))) {
  custom3Counts <- readRDS(paste0(RobjectDir, "/", Type, "_allcounts.htseq.rds"))
  gcCounts <- readRDS(paste0(RobjectDir, "/gc_allcounts.htseq.rds"))
  
  # append gcCounts to custom3Counts:
  Counts <- rbind(custom3Counts, gcCounts)
  
  # make rownames gene_id, get rid of latter column and change
  # storage mode from factor to integer:
  rownames(Counts) <- Counts$gene_id
  Counts <- subset(Counts, select=-gene_id)
  # save Counts:
  saveRDS(Counts, file=paste0(newRobjectDir, "/", Type, "_counts.RData"))
} else {
  Counts <- readRDS(file=paste0(newRobjectDir, "/", Type, "_counts.RData"))
}

# re-categorize samples as HRD, CCNE_amp, both_drivers or unknown_drivers:
if (cat_by_driver) {
  # load in sample key for categories homologous repair deficient (HRD) and cyclin E gain/amplification (CCNE):
  HRDkey <- read.table(paste0(rawDir, "/hrd_samples.txt"), header=F, sep="\t")
  HRDnos <- gsub("AOCS_", "", HRDkey[,1])
  CCNEkey <- read.table(paste0(rawDir, "/ccne_gain_or_amp_samples.txt"), header=F, sep="\t")
  CCNEnos <- gsub("AOCS_", "", CCNEkey[,1])
  Counts <- Counts[,grep("rcAF|pAF|msST", colnames(Counts), invert=T)]
  for (i in 1:length(colnames(Counts))) {
    print(i)
    if (length(grep("FT", colnames(Counts)[i]))<1) {
      no <- gsub(
        "_[a-zA-Z].*$", "", gsub("AOCS_", "", colnames(Counts)[i])
      )
      print(paste0("ID number is ", no))
      if (no %in% HRDnos & no %in% CCNEnos) {
        colnames(Counts)[i] <- paste0(colnames(Counts)[i], "_bothDrivers")
      } else if (no %in% HRDnos & !(no %in% CCNEnos)) {
        colnames(Counts)[i] <- paste0(colnames(Counts)[i], "_HRD")
      } else if (no %in% CCNEnos & !(no %in% HRDnos)) {
        colnames(Counts)[i] <- paste0(colnames(Counts)[i], "_CCNEamp")
      } else if (!(no %in% HRDnos | no %in% CCNEnos)) {
        colnames(Counts)[i] <- paste0(colnames(Counts)[i], "_unknown_drivers")
      }
      colnames(Counts)[i] <- gsub("[a-z].*[A-Z][A-Z]_", "", colnames(Counts)[i])
    }
  }
  
  # remove any samples not belonging to any group:
  Counts <- Counts[, colnames(
    Counts[, gsub(
      ".*\\_", "", colnames(Counts)
    ) %in% unlist(sGroups)]
  )]
  
} else {
  # remove any samples not belonging to any group:
  Counts <- Counts[, colnames(
    Counts[, gsub(
      ".*\\_", "", colnames(Counts)
    ) %in% unlist(sGroups)]
  )]
  
  # change sample names according to grouping:
  for (i in 1:length(sGroups)) {
    for (n in sGroups[[i]]) {
      colnames(Counts) <- gsub(n, names(sGroups)[i], colnames(Counts))
    }
  }
}

# eliminate lowly expressed genes (rows where there are less than 3 counts where df > 4):
print(paste0("No. rows before filtering is: ", nrow(Counts)))
Counts <- Counts %>%
  rownames_to_column('gene_id') %>%
  dplyr::filter(rowSums(Counts > 5) >= (ncol(Counts)/3)) %>%
  column_to_rownames('gene_id')
print(paste0("No. rows after  filtering: ", nrow(Counts)))


### 2. Perform pre-normalisation PCA and RLE plots ###

# create pre-normalised PCA plot from counts and plot:
if (ncol(Counts) > nrow(Counts)) {
  pca <- prcomp(Counts)
} else {
  pca <- princomp(Counts)	  
}

if (file.exists(paste0(plotDir, "/", Type, "_pcaCompsPrenormGC.pdf"))) {
  print(paste0(plotDir, "/", Type, "_pcaCompsPrenormGC.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", plotDir, "/", Type, "_pcaCompsPrenormGC.pdf"))
  pdf(file = paste0(plotDir, "/", Type, "_pcaCompsPrenormGC.pdf"))
  plot(pca)
  dev.off()
}

# change the order of columns of Counts to alphabetical order:
Counts <- Counts[,order(
  gsub(
    "AOCS.*_[0-9][0-9][0-9]_", "", colnames(Counts)
  )
)]

# define sample groups:
if (allHGSOC_vs_FT == TRUE) {
  splt <- unlist(
    lapply(
      split(
        colnames(Counts), gsub(
          "AOCS.*$", "HGSOC", gsub(
            "^.*FT", "FT", colnames(Counts)
          )
        )
      ), length
    )
  )
} else {
  splt <- unlist(
    lapply(
      split(
        colnames(Counts), gsub(
          "AOCS.*_[0-9][0-9][0-9]_", "", colnames(Counts)
        )
      ), length
    )
  )
}

for (i in 1:length(splt)) {
  if (i==1) {
    typeF <- c(rep(names(splt)[i], splt[i]))
  } else {
    typeF <- c(typeF, rep(names(splt)[i], splt[i]))
  }
}
levels(typeF) <- sTypes

# save number of samples in each group:
sampleNos <- unlist(
  lapply(
    split(
      colnames(Counts), gsub(
        "\\.1", "", gsub(
          "AOCS.*_[0-9][0-9][0-9]_", "", colnames(Counts)
        )
      )
    ), length
  )
)

saveRDS(sampleNos, file = paste0(newRobjectDir, "/sample_no_per_cat.rds"))

# convert Counts into SeqExpressionSet - elements need to be delisted and changed to integers first:
Counts <- apply(Counts, 2, unlist)
storage.mode(Counts) <- "integer"

if (analysis == "edger") {
  
  ### 3. perform normalisation and estimate dispersion off counts using EdgeR ###
  # create EDAseq expression set for RLE/PCA plot:
  set <- newSeqExpressionSet(Counts, phenoData = data.frame(typeF, row.names=colnames(Counts)))
  
  # create pre-norm RLE plot:
  if (file.exists(paste0(plotDir, "/", Type, "_", analysis, "_RLEPrenormGC.pdf"))) {
    print(paste0(plotDir, "/", Type, "_", analysis, "_RLEPrenormGC.pdf already exists, no need to create"))
  } else {
    print(paste0("Creating ", plotDir, "/", Type, "_", analysis, "_RLEPrenormGC.pdf"))
    par(mar=c(1,1,1,1))
    pdf(file = paste0(plotDir, "/", Type, "_", analysis, "_RLEPrenormGC.pdf"))
    plotRLE(set)
    dev.off()
  }
  
  # create RUVseq pre-norm PCA:
  if (file.exists(paste0(plotDir, "/", Type, "_", analysis, "_pcaPrenormGC.pdf"))) {
    print(paste0(plotDir, "/", Type, "_", analysis, "_pcaPrenormGC.pdf already exists, no need to create"))
  } else {
    print(paste0("Creating ", plotDir, "/", Type, "_", analysis, "_pcaPrenormGC.pdf"))
    pdf(file = paste0(plotDir, "/", Type, "_", analysis, "_pcaPrenormGC.pdf"), height = 10, width = 12)
    plotPCA(set, cex=0.7)
    dev.off()
  }
  
  # perform between lane full normalisation:
  y <- DGEList(counts = Counts, group = typeF)
  
  # normalise for library size:
  y <- calcNormFactors(y)
  
  # create an MDS plot to show relative similarities of the samples and save to plotDir:
  if (file.exists(paste0(plotDir, "/", Type, "_", analysis, "_MDS.pdf"))) {
    paste0(plotDir, "/", Type, "_", analysis, "_MDS.pdf already exists")
    pdf(paste0(plotDir, "/", Type, "_", analysis, "_MDS.pdf"),width=16,height=12)
    plotMDS(y)
  } else {
    paste0("Generating ", plotDir, "/", Type, "_", analysis, "_MDS.pdf")
    pdf(paste0(plotDir, "/", Type, "_", analysis, "_MDS.pdf"),width=16,height=12)
    plotMDS(y)
    dev.off()
  }
  
  if (!file.exists(paste0(newRobjectDir, "/edger_norms.rds"))) {
    # create EDAseq expression set for RLE/PCA plot:
    for (n in 1:nrow(Counts)) {
      print(n)
      if (n==1) {
        norms <- t(as.matrix(y$samples$norm.factors))
      } else {
        norms <- rbind(norms, norms[1,])
      }
    }
    saveRDS(norms, file=paste0(newRobjectDir, "/edger_norms.rds"))
  } else {
    norms <- readRDS(file=paste0(newRobjectDir, "/edger_norms.rds"))
  }

  set <- newSeqExpressionSet(Counts, offset = norms, phenoData = data.frame(typeF, row.names=colnames(Counts)))
  
  # create post-norm RLE plot:
  if (file.exists(paste0(plotDir, "/", Type, "_", analysis, "_RLElaneNormGC.pdf"))) {
    print(paste0(plotDir, "/", Type, "_", analysis, "_RLElaneNormGC.pdf already exists, no need to create"))
  } else {
    print(paste0("Creating ", plotDir, "/", Type, "_", analysis, "_RLElaneNormGC.pdf"))
    pdf(file = paste0(plotDir, "/", Type, "_", analysis, "_RLElaneNormGC.pdf"))
    plotRLE(set, outline=FALSE, ylim=c(-4, 4))
    dev.off()
  }
  
  # create RUVseq post-norm PCA:
  if (file.exists(paste0(plotDir, "/", Type, "_", analysis, "_pcalaneNormGC.pdf"))) {
    print(paste0(plotDir, "/", Type, "/", Type, "_", analysis, "_pcalaneNormGC.pdf already exists, no need to create"))
  } else {
    print(paste0("Creating ", plotDir, "/", Type, "_", analysis, "_pcalaneNormGC.pdf"))
    pdf(file = paste0(plotDir, "/", Type, "_", analysis, "_pcalaneNormGC.pdf"), height = 15, width = 20)
    plotPCA(set, cex=0.7)
    dev.off()
  }
  
  # design matrix labelling all sample types:
  design <- model.matrix(~0+typeF)
  
  # estimate dispersion:
  disp <- estimateDisp(y, design=design)
  
  # adjust values using dispersion:
  fit <- glmFit(disp, design=design, robust=TRUE)
 
} else if (analysis == "edger_eda") {
  
  ### 3. perform normalisation and estimate dispersion off counts using EDAseq and EdgeR ###
  
  set <- newSeqExpressionSet(Counts, phenoData = data.frame(typeF, row.names=colnames(Counts)))
  
  # create pre-norm RLE plot:
  if (file.exists(paste0(plotDir, "/", Type, "_", analysis, "_RLEPrenormGC.pdf"))) {
    print(paste0(plotDir, "/", Type, "_", analysis, "_RLEPrenormGC.pdf already exists, no need to create"))
  } else {
    print(paste0("Creating ", plotDir, "/", Type, "_", analysis, "_RLEPrenormGC.pdf"))
    par(mar=c(1,1,1,1))
    pdf(file = paste0(plotDir, "/", Type, "_", analysis, "_RLEPrenormGC.pdf"))
    plotRLE(set)
    dev.off()
  }
  
  # create RUVseq pre-norm PCA:
  if (file.exists(paste0(plotDir, "/", Type, "_", analysis, "_pcaPrenormGC.pdf"))) {
    print(paste0(plotDir, "/", Type, "_", analysis, "_pcaPrenormGC.pdf already exists, no need to create"))
  } else {
    print(paste0("Creating ", plotDir, "/", Type, "_", analysis, "_pcaPrenormGC.pdf"))
    pdf(file = paste0(plotDir, "/", Type, "_", analysis, "_pcaPrenormGC.pdf"), height = 10, width = 12)
    plotPCA(set, cex=0.7)
    dev.off()
  }
  
  # perform between lane full normalisation:
  nSet <- betweenLaneNormalization(set, which="full",offset=TRUE)
  
  # create post-norm RLE plot:
  if (file.exists(paste0(plotDir, "/", Type, "_", analysis, "_RLElaneNormGC.pdf"))) {
    print(paste0(plotDir, "/", Type, "_", analysis, "_RLElaneNormGC.pdf already exists, no need to create"))
  } else {
    print(paste0("Creating ", plotDir, "/", Type, "_", analysis, "_RLElaneNormGC.pdf"))
    pdf(file = paste0(plotDir, "/", Type, "_", analysis, "_RLElaneNormGC.pdf"))
    plotRLE(nSet, outline=FALSE, ylim=c(-4, 4))
    dev.off()
  }
  
  # create RUVseq post-norm PCA:
  if (file.exists(paste0(plotDir, "/", Type, "_", analysis, "_pcalaneNormGC.pdf"))) {
    print(paste0(plotDir, "/", Type, "_", analysis, "_pcalaneNormGC.pdf already exists, no need to create"))
  } else {
    print(paste0("Creating ", plotDir, "/", Type, "_", analysis, "_pcalaneNormGC.pdf"))
    pdf(file = paste0(plotDir, "/", Type, "_", analysis, "_pcalaneNormGC.pdf"), height = 15, width = 20)
    plotPCA(nSet, cex=0.7)
    dev.off()
  }
  
  # design matrix labelling all sample types:
  design <- model.matrix(~0+typeF, data=pData(nSet))
  
  # calculate the deviance residuals from a first-pass GLM regression of the counts on the co-variates of interest (p8 RUVseq manual):
  # estimate dispersion:
  disp <- estimateGLMCommonDisp(counts(nSet), design, offset=-offst(nSet))
  
  # adjust values using dispersion:
  fit <- glmFit(counts(nSet), design, disp, offset=-offst(nSet))
  
  saveRDS(fit, file=paste0(newRobjectDir, "/", Type, "_", analysis, "DEfit.rds"))
  save.image(paste0(newRobjectDir, "/", Type, "_", analysis, "DEdone.rds"))
  
}


### 4. Perform differential expression comparing normalised FT controls to cancer samples ###

# determine which column has FT control:
ctlInd <- grep(ctl, colnames(design))
con <- c(rep(0, (ctlInd - 1) ), -1, rep(0, (ncol(design) - ctlInd)))

# put sTypes in alphabetical order:
sTypes <- sTypes[order(sTypes)]

for (i in 1:ncol(design)) {
  print(i)
  if (i!=ctlInd) {
    comp <- paste0(sTypes[i], "_vs_", ctl)
    
    # perform likelihood ratio test:
    con[i] <- 1
    lrt <- glmLRT(fit, contrast = con)
    
    # determine the top DE genes:
    topTags(lrt)
    
    if (file.exists(paste0(newRobjectDir, sTypes[i], "_vs_", ctl, "_lrt.rds"))) {
      print(paste0(newRobjectDir, sTypes[i], "_vs_", ctl, "_lrt.rds already exists, no need to create"))
    } else {
      print(paste0("Creating ", newRobjectDir, sTypes[i], "_vs_", ctl, "_lrt.rds"))
      saveRDS(lrt, file = paste0(newRobjectDir, sTypes[i], "_vs_", ctl, "_lrt.rds"))
    }
    
    # fetch summary of differentially expressed genes (those  with FDR =< 0.05:
    DEs <- summary(result <- decideTestsDGE((lrt)))
    
    # fetch all gene DE info, 
    allGenes <- as.data.frame(topTags(lrt, n=Inf))
    
    # create volcano plots for allGenes:
    if ("all" %in% resultTypes) {
      # define sig DE df:
      if (length(FCthresh) == 0) {
        sigGenes <- filter(allGenes, FDR < FDRthresh)
        allGenes$threshold <- as.factor(allGenes$FDR < FDRthresh)
      } else {
        sigGenes <- filter(allGenes, (FDR < FDRthresh & logFC < -(FCthresh))|(FDR < FDRthresh & logFC > FCthresh))
        allGenes$threshold <- as.factor((allGenes$FDR < FDRthresh & allGenes$logFC < -(FCthresh))|(allGenes$FDR <  FDRthresh & allGenes$logFC > FCthresh))
      }
      
      sig <- subset(allGenes, threshold == T)
      
      # include the control genes for labelling:
      for (j in 1:length(posGeneIDs)) {
        if (j==1) {
          posGenes <- allGenes[ posGeneIDs[j],]
        } else {
          posGenes <- rbind(posGenes,   allGenes[posGeneIDs[j],])
        }
      }
      rownames(posGenes) <- posGeneNames
      
      for (j in 1:length(negGeneIDs)) {
        if (j==1) {
          negGenes <- allGenes[ negGeneIDs[j],]
        } else {
          negGenes <- rbind(negGenes,   allGenes[negGeneIDs[j],])
        }
      }
      rownames(negGenes) <- negGeneNames
      
      # set default threshold statuses  for control genes:
      posGenes$threshold <- "POSITIVE"
      if (nrow(posGenes[posGenes$FDR< FDRthresh,])>0) {
        posGenes[posGenes$FDR<  FDRthresh,]$threshold <- "POSSIG"
      }
      
      negGenes$threshold = "NEGATIVE"
      if (nrow(negGenes[negGenes$FDR< FDRthresh,])>0) {
        negGenes[negGenes$FDR<  FDRthresh,]$threshold <-  "NEGSIG"
      }
      
      #lab <- rbind(rbind(sig,   posGenes), negGenes)
      lab <- sig
      #allGenes <- rbind(rbind(allGenes,   posGenes), negGenes)
      lab$genes <- rownames(lab)
      
      # set max genes to label:
      lab <- head(lab, 200)
      
      # fetch symbols for labelling:
      egENSEMBL <- toTable(org.Hs.egENSEMBL)
      egSYMBOL <- toTable(org.Hs.egSYMBOL)
      m <- match(row.names(lab), egENSEMBL$ensembl_id)
      lab$gene_id <- egENSEMBL$gene_id[m]
      m <- match(lab$gene_id, egSYMBOL$gene_id)
      lab$symbol <- egSYMBOL$symbol[m]
      
      # remove NAs:
      lab <- na.omit(lab)
      
      #######
      reported <- read.table(file=paste0(plotDir, "/reportedDEs.txt"), header=T, sep="\t")
      reported <- reported[,c(1:2, 6)]
      colnames(reported) <- c("symbol", "Bowtell_logFC", "Bowtell_adj_pvalue")
      
      m <- match(row.names(allGenes), egENSEMBL$ensembl_id)
      allGenes$gene_id <- egENSEMBL$gene_id[m]
      m <- match(allGenes$gene_id, egSYMBOL$gene_id)
      allGenes$symbol <- egSYMBOL$symbol[m]
      
      matching <- allGenes$symbol %in% reported[,1]
      
      reportedDE <- subset(
        allGenes[matching,], select=c(logFC, FDR, symbol)
      )
      
      htseq_vs_Bowtell <- merge(reportedDE, reported, by.x = "symbol")
      
      # add 'type' column to reportedDE, reported and rbind together for volcano plot:
      reportedDE <- reportedDE[,c(3, 1, 2)]
      reportedDE$type <- "htseq"
      colnames(reportedDE) <- c("symbol", "logFC", "FDR_or_adj_pvalue", "type")
      
      reported$type <- "Bowtell"
      colnames(reported) <- c("symbol", "logFC", "FDR_or_adj_pvalue", "type")
      
      htseq_Bowtell <- rbind(reportedDE, reported)
      
      # add label data frame:
      lab1 <- reportedDE[order(reportedDE$FDR_or_adj_pvalue),][1:25,]
      lab2 <- reported[order(reported$FDR_or_adj_pvalue),][1:25,]
      lab <- rbind(lab1, lab2)
      
      # prepare data frame for heatmap:
      hmDF <- data.frame(htseq_vs_Bowtell$logFC, htseq_vs_Bowtell$Bowtell_logFC)
      rownames(hmDF) <- htseq_vs_Bowtell$symbol
      
      # prepare adj_pvalue/FDR data frame:
      pDF <- data.frame(htseq_vs_Bowtell$FDR, htseq_vs_Bowtell$Bowtell_adj_pvalue)
      rownames(pDF) <- htseq_vs_Bowtell$symbol
      
      # create heatmap:
      pheatmap(hmDF, color = colorRampPalette(c("#08519C", "white", "firebrick3"))(50),
               display_numbers = as.matrix(ifelse(interesting_fdr < 0.1, "*", "")), fontsize = 8,
               cluster_cols = T)
      
      
      p <- ggplot(data=htseq_Bowtell, aes( x=logFC, y=-log10(FDR_or_adj_pvalue),    color=type) )
      p <- p + geom_point(data=htseq_Bowtell)
      p <- p + geom_text_repel(data=lab, aes(label=symbol))
      p <- p + theme(legend.position =  "none")
      p <- p + labs(x="log2 fold change   vs FT control", y="-log10   FDR")
      p
      
      ######
      
      if (!(ctlInd==1)) {
        if (i==1) {
          allReps <- list(allGenes)
        } else {
          allReps[[i]] <- allGenes
        }
        
        if (i==1) {
          sigReps <- list(sig)
        } else {
          sigReps[[i]] <- sig
        }
      } else {
        if (i==2) {
          allReps <- list(allGenes)
        } else {
          allReps[[i]] <- allGenes
        }
        
        if (i==2) {
          sigReps <- list(sig)
        } else {
          sigReps[[i]] <- sig
        }
      }
      
      # plot on volcano plot:
      p <- ggplot(data=allGenes, aes( x=logFC, y=-log10(FDR),    color=threshold))
      p <- p + geom_point(data=allGenes)
      p <- p + geom_text_repel(data=lab, aes(label=symbol))
      p <- p + theme(legend.position =  "none")
      p <- p + labs(x="log2 fold change   vs FT control", y="-log10   FDR")
      #p <- p +  xlim(c(-4, 4))
      if (length(FCthresh) == 0) {
        if (file.exists(paste0(plotDir,   "/", Type,  "_", analysis, "_volcano_FDR_",   FDRthresh, "_", comp, ".pdf"))) {
          print(paste0(plotDir, "/",  Type,  "_", analysis, "_volcano_FDR_",   FDRthresh, "_", comp, ".pdf"))
          p
        } else {
          print(paste0("Creating  ",plotDir, "/", Type,    "_volcano_FDR_", FDRthresh, "_", comp, ".pdf"))
          pdf(file = paste0(plotDir, "/",   Type,  "_", analysis, "_volcano_FDR_",  FDRthresh, "_", comp, ".pdf"))
          print(p)
          dev.off()
        }
      } else {
        if (file.exists(paste0(plotDir, "/",  Type,  "_", analysis, "_volcano_FDR_",   FDRthresh, "_FC", FCthresh, "_", comp, ".pdf"))) {
          print(paste0(plotDir, "/",  Type,  "_", analysis, "_volcano_FDR_",   FDRthresh, "_FC", FCthresh, "_", comp, ".pdf already exists"))
          p
        } else {
          print(paste0("Creating  ", plotDir, "/",  Type,  "_", analysis, "_volcano_FDR_", FDRthresh, "_FC", FCthresh, "_", comp, ".pdf"))
          pdf(file = paste0(plotDir, "/",  Type,  "_", analysis, "_volcano_FDR_",   FDRthresh, "_FC", FCthresh, "_", comp, ".pdf"))
          print(p)
          dev.off()
        }
      }
    }
    
    
    ### 5. Calculate differential expression values of repeats ###
    
    if ("repeats" %in% resultTypes) {
      # define repeat and sig DE repeat dfs:
      repGenes <- allGenes[grep("ENS",  rownames(allGenes), invert = T),]
      print(repGenes)
      
      if (length(FCthresh) == 0) {
        sigGenes <- filter(repGenes, FDR < FDRthresh)
        repGenes$threshold <- as.factor(repGenes$FDR < FDRthresh)
      } else {
        sigGenes <- filter(repGenes, (FDR < FDRthresh & logFC < -(FCthresh))|(FDR < FDRthresh & logFC > FCthresh))
        repGenes$threshold <- as.factor((repGenes$FDR < FDRthresh & repGenes$logFC < -(FCthresh))|(repGenes$FDR <  FDRthresh & repGenes$logFC > FCthresh))
      }
      
      sig <- subset(repGenes, threshold == T)
      
      # include the control genes for labelling:
      for (j in 1:length(posGeneIDs)) {
        if (j==1) {
          posGenes <- allGenes[ posGeneIDs[j],]
        } else {
          posGenes <- rbind(posGenes,   allGenes[posGeneIDs[j],])
        }
      }
      rownames(posGenes) <- posGeneNames
      
      for (j in 1:length(negGeneIDs)) {
        if (j==1) {
          negGenes <- allGenes[ negGeneIDs[j],]
        } else {
          negGenes <- rbind(negGenes,   allGenes[negGeneIDs[j],])
        }
      }
      rownames(negGenes) <- negGeneNames
      
      # set default threshold statuses  for control genes:
      posGenes$threshold <- "POSITIVE"
      if (nrow(posGenes[posGenes$FDR< FDRthresh,])>0) {
        posGenes[posGenes$FDR<  FDRthresh,]$threshold <- "POSSIG"
      }
      
      negGenes$threshold = "NEGATIVE"
      if (nrow(negGenes[negGenes$FDR< FDRthresh,])>0) {
        negGenes[negGenes$FDR<  FDRthresh,]$threshold <-  "NEGSIG"
      }
      
      lab <- rbind(rbind(sig,   posGenes), negGenes)
      repGenes <- rbind(rbind(repGenes,   posGenes), negGenes)
      lab$genes <- rownames(lab)
      
      if (!(ctlInd==1)) {
        if (i==1) {
          allReps <- list(repGenes)
        } else {
          allReps[[i]] <- repGenes
        }
        
        if (i==1) {
          sigReps <- list(sig)
        } else {
          sigReps[[i]] <- sig
        }
      } else {
        if (i==2) {
          allReps <- list(repGenes)
        } else {
          allReps[[i]] <- repGenes
        }
        
        if (i==2) {
          sigReps <- list(sig)
        } else {
          sigReps[[i]] <- sig
        }
      }
      
      # plot on volcano plot:
      p <- ggplot(data=repGenes, aes( x=logFC, y=-log10(FDR),    color=threshold))
      p <- p + geom_point(data=repGenes)
      p <- p + geom_text_repel(data=lab, aes(label=genes))
      p <- p + theme(legend.position =  "none")
      p <- p + labs(x="log2 fold change   vs FT control", y="-log10   FDR")
      p <- p +  xlim(c(-4, 4))
      if (length(FCthresh) == 0) {
        if (file.exists(paste0(plotDir,   "/", Type,  "_", analysis, "_volcano_FDR_",   FDRthresh, "_", comp, ".pdf"))) {
          print(paste0(plotDir, "/",  Type,  "_", analysis, "_volcano_FDR_",   FDRthresh, "_", comp, ".pdf"))
          p
        } else {
          print(paste0("Creating  ",plotDir, "/", Type,    "_volcano_FDR_", FDRthresh, "_", comp, ".pdf"))
          pdf(file = paste0(plotDir, "/",   Type,  "_", analysis, "_volcano_FDR_",  FDRthresh, "_", comp, ".pdf"))
          print(p)
          dev.off()
        }
      } else {
        if (file.exists(paste0(plotDir, "/",  Type,  "_", analysis, "_volcano_FDR_",   FDRthresh, "_FC", FCthresh, "_", comp, ".pdf"))) {
          print(paste0(plotDir, "/",  Type,  "_", analysis, "_volcano_FDR_",   FDRthresh, "_FC", FCthresh, "_", comp, ".pdf already exists"))
          p
        } else {
          print(paste0("Creating  ", plotDir, "/",  Type,  "_", analysis, "_volcano_FDR_", FDRthresh, "_FC", FCthresh, "_", comp, ".pdf"))
          pdf(file = paste0(plotDir, "/",  Type,  "_", analysis, "_volcano_FDR_",   FDRthresh, "_FC", FCthresh, "_", comp, ".pdf"))
          print(p)
          dev.off()
        }
      }
      
      if ("epiMods" %in% resultTypes) {
        
        
        ### 6. Calculate differential expression values of epigenetic modifier genes ###
        
        # create epiGenes df:
        epiGenes <- allGenes[epiIDs,]
        rownames(epiGenes) <- epiSym
        
        # set threshold status:
        if (length(FCthresh) == 0) {
          epiGenes$threshold <- as.factor(epiGenes$FDR < FDRthresh)
        } else {
          epiGenes$threshold <- as.factor((epiGenes$FDR < FDRthresh & epiGenes$logFC < -(FCthresh))|(epiGenes$FDR < 
                                                                                                       FDRthresh & epiGenes$logFC > FCthresh))
        }
        
        # create significant epiGenes df:
        epiSig <- subset(epiGenes, threshold == T)
        epiGenes$genes <- rownames(epiGenes)
        
        if (!(ctlInd==1)) {
          if (i==1) {
            allEpi <- list(epiGenes)
          } else {
            allEpi[[i]] <- epiGenes
          }
        } else {
          if (i==2) {
            allEpi <- list(epiGenes)
          } else {
            allEpi[[i]] <- epiGenes
          }
        }
        
        # create volcano plots with repeat values in grey in background:
        p <- ggplot(data=epiGenes, aes(x=logFC, y=-log10(FDR),color=threshold))
        p <- p + geom_point(data=epiGenes)
        p <- p + geom_text_repel(data=epiGenes, aes(label=genes))
        p <- p + theme(legend.position = "none")
        p <- p + labs(x="log2 fold change vs FT control", y="-log10   FDR")
        p <- p +  xlim(c(-2, 2))
        if (length(FCthresh) == 0) {
          if (file.exists(paste0(plotDir,   "/", Type,  "_", analysis, "_volcano_FDR_",   FDRthresh, "_", comp, "_epi.pdf"))) {
            print(paste0(plotDir, "/",  Type,  "_", analysis, "_volcano_FDR_",   FDRthresh, "_", comp, "_epi.pdf"))
            p
          } else {
            print(paste0("Creating  ",plotDir, "/", Type,    "_volcano_FDR_", FDRthresh, "_", comp, "_epi.pdf"))
            pdf(file = paste0(plotDir, "/",   Type,  "_", analysis, "_volcano_FDR_",  FDRthresh, "_", comp, "_epi.pdf"))
            print(p)
            dev.off()
          }
        } else {
          if (file.exists(paste0(plotDir, "/",  Type,  "_", analysis, "_volcano_FDR_",   FDRthresh, "_FC", FCthresh, "_", comp, "_epi.pdf"))) {
            print(paste0(plotDir, "/",  Type,  "_", analysis, "_volcano_FDR_",   FDRthresh, "_FC", FCthresh, "_", comp, "_epi.pdf already exists"))
            p
          } else {
            print(paste0("Creating  ", plotDir, "/",  Type,  "_", analysis, "_volcano_FDR_", FDRthresh, "_FC", FCthresh, "_", comp, "_epi.pdf"))
            pdf(file = paste0(plotDir, "/",  Type,  "_", analysis, "_volcano_FDR_",   FDRthresh, "_FC", FCthresh, "_", comp, "_epi.pdf"))
            print(p)
            dev.off()
          }
        }
      }
    }
    con <- c(rep(0, (ctlInd - 1) ), -1, rep(0, (ncol(design) - ctlInd)))
  }
}


### 8. Save all repeat DE, and all repeat DE that is significant and log2 FC > 1 for at least one comparison:

# remove the NULL list dfs created when avoiding clt vs ctl:
if (sTypes>2) {
  allReps <- allReps[-ctlInd]
  sigReps <- sigReps[-ctlInd]
}

# name the list elements:
names(allReps) <- paste0(sTypes[-ctlInd], "_vs_", ctl)
names(sigReps) <- paste0(sTypes[-ctlInd], "_vs_", ctl)

# save each list for downstream analysis:
saveRDS(allReps, file=paste0(newRobjectDir, "/", Type, "_DEreps.rds"))

# do above for epi genes if necessary:
if ("epiMods" %in% resultTypes) {
  allEpi <- allEpi[-ctlInd]
  names(allEpi) <- paste0(sTypes[-ctlInd], "_vs_", ctl)
  
  saveRDS(allEpi, file=paste0(newRobjectDir, "/", Type, "_DEepiGenes.rds"))
}

# don't just want significant genes for each individual group, so find names of
# all the sig genes from all groups by taking the rownames and unlisting into one vector:
totalSig <- unique(unlist(lapply(sigReps, function(x) {
  return(rownames(x))
})))
names(totalSig) <- NULL

# fetch all the entries for the above significant repeat names and save:
allSig <- lapply(allReps, function(x) {
  return(x[totalSig,])
})

saveRDS(allSig, file=paste0(newRobjectDir, "/", Type, "_DEsigReps.rds"))

if (!file.exists(paste0(newRobjectDir, "DEImg_", expName, ".RData"))) {
  save.image(file = paste0(newRobjectDir, "DEImg_", expName, ".RData"))
}
