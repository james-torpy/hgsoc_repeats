
# Run on cluster with:
#briR
#qsub -N DESeq -b y -wd \
#/share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/logs/exp9/DE \
#-j y -R y -pe smp 4 -V "Rscript /share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/exp9/17a.DESeq_indivHGSOC_vs_FT.R"


### 0. Define variables/paths ###

# load packages needed:
library(DESeq2)
library(tibble)
library(dplyr)
library(RColorBrewer)
library(EDASeq)
library(ggplot2)
library(ggrepel)

# define starting variables:
project <- "hgsoc_repeats"
expName <- "exp9"
Type <- "custom3"

sTypes <- c("CCNEamp", "HRD")
descrip <- "htseq_DESeq2_primary_CCNEamp_vs_HRD"

# define sample group to use as control:
ctl <- "HRD"

# define sample groups to compare:
allHGSOC_vs_FT <- FALSE
cat_by_driver <- TRUE
EDAnormalise <- FALSE
primaryOnly <- TRUE

# specify what combination of repeat genes (repeats), epigenetic modulators (epiMods),
# RNAi genes (RNAi) and protein-coding genes (pCoding) should contribute to the results:
#resultTypes <- c("repeats", "epiMods")
resultTypes <- c("repeats", "all")

# specify what padj and log2 fold change thresholds to use:
padjthresh <- 0.05
FCthresh <- 0

# specify control genes to include:
posGeneIDs <- c("ENSG00000111640", "ENSG00000196776")
posGeneNames <- c("GAPDH", "CD47")
negGeneIDs <- c("ENSG00000075624", "ENSG00000169919")
negGeneNames <- c("beta-actin", "GUSB")


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
                  "/plots/DEplots/DEcompare/", descrip, "/")

system(paste0("mkdir -p ", plotDir))
system(paste0("mkdir -p ", newRobjectDir))

  
### 1. Load in all counts ###

if ( !file.exists(file=paste0(newRobjectDir, "/", Type, "_counts.RData")) ) {
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
  
  # if necessary, select primary samples only:
  if ( primaryOnly == TRUE) {
    Counts <- Counts[,grep("PT|FT", colnames(Counts))]
  }
  
  # re-categorize samples as HRD, CCNE_amp, both_drivers or unknown_drivers:
  if (cat_by_driver) {
    # load in sample key for categories homologous repair deficient (HRD) and cyclin E gain/mplification (CCNE):
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
  
  # save Counts:
  saveRDS(Counts, file=paste0(newRobjectDir, "/", Type, "_counts.RData"))
} else {
  Counts <- readRDS(file=paste0(newRobjectDir, "/", Type, "_counts.RData"))
}


### 2. Perform pre-normalisation PCA and RLE plots ###

# eliminate lowly expressed genes (rows where there are less than 3 counts where df > 4):
print(paste0("No. rows before filtering is: ", nrow(Counts)))
Counts <- Counts %>%
  rownames_to_column('gene_id') %>%
  dplyr::filter(rowSums(Counts > 5) >= (ncol(Counts)/3)) %>%
  column_to_rownames('gene_id')
print(paste0("No. rows after  filtering: ", nrow(Counts)))

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

# remove bothDrivers, FT, unknownDrivers:
Counts <- Counts[,grep("bothDrivers|FT|unknown_drivers", colnames(Counts), invert=T)]

# change the order of columns of Counts so control is first:
Counts <- Counts[,order(
  gsub(
    "AOCS_[0-9][0-9][0-9]_", "", colnames(Counts)
  ), decreasing = T
)]

# define sample groups:
if (allHGSOC_vs_FT == TRUE) {
  coldata <- data.frame(
    gsub(
      "AOCS.*$", "HGSOC", gsub(
        "^.*FT", "FT", colnames(Counts)
      )
    ),
    c(rep("paired-end", ncol(Counts)))
  )
  colnames(coldata) <- c("condition", "type")
  
} else {
  
  coldata <- data.frame(
    gsub("AOCS.*_[0-9][0-9][0-9]_", "", colnames(Counts)),
    c(rep("paired-end", ncol(Counts)))
  )
  colnames(coldata) <- c("condition", "type")
  rownames(coldata) <- colnames(Counts)
}

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

# convert Counts into SeqExpressionSet - elements need to be delisted and changed to integers irst:
Counts <- apply(Counts, 2, unlist)
storage.mode(Counts) <- "integer"
set <- newSeqExpressionSet(Counts, phenoData = data.frame(condition, row.names=colnames(Counts)))

# create pre-norm RLE plot:
if (file.exists(paste0(plotDir, "/", Type, "_RLEPrenormGC.pdf"))) {
  print(paste0(plotDir, "/", Type, "_RLEPrenormGC.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", plotDir, "/", Type, "_RLEPrenormGC.pdf"))
  par(mar=c(1,1,1,1))
  pdf(file = paste0(plotDir, "/", Type, "_RLEPrenormGC.pdf"))
  plotRLE(set)
  dev.off()
}

# create RUVseq pre-norm PCA:
if (file.exists(paste0(plotDir, "/", Type, "_pcaPrenormGC.pdf"))) {
  print(paste0(plotDir, "/", Type, "_pcaPrenormGC.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", plotDir, "/", Type, "_pcaPrenormGC.pdf"))
  pdf(file = paste0(plotDir, "/", Type, "_pcaPrenormGC.pdf"), height = 10, width = 12)
  plotPCA(set)
  dev.off()
}


### 3. perform normalisation on counts:
if ( EDAnormalise == TRUE ) {
  
  nSet <- betweenLaneNormalization(set, which="full")
  
  # create post-norm RLE plot:
  if (file.exists(paste0(plotDir, "/", Type, "_", sName, "_vs_FT__RLElaneNormGC.pdf"))) {
    print(paste0(plotDir, "/", Type, "_", sName, "_vs_FT__RLElaneNormGC.pdf already exists, no need to create"))
  } else {
    print(paste0("Creating ", plotDir, "/", Type, "_", sName, "_vs_FT__RLElaneNormGC.pdf"))
    pdf(file = paste0(plotDir, "/", Type, "_", sName, "_vs_FT__RLElaneNormGC.pdf"))
    plotRLE(nSet, outline=FALSE, ylim=c(-4, 4))
    dev.off()
  }
  
  # create RUVseq post-norm PCA:
  if (file.exists(paste0(plotDir, "/", Type, "_", sName, "_vs_FT__pcalaneNormGC.pdf"))) {
    print(paste0(plotDir, "/", Type, "_", sName, "_vs_FT__pcalaneNormGC.pdf already exists, no need to create"))
  } else {
    print(paste0("Creating ", plotDir, "/", Type, "_", sName, "_vs_FT__pcalaneNormGC.pdf"))
    pdf(file = paste0(plotDir, "/", Type, "_", sName, "_vs_FT__pcalaneNormGC.pdf"), height = 15, width = 20)
    plotPCA(nSet, cex=0.7)
    dev.off()
  }
  
  countsN <- as(nSet, "CountDataSet")
  sizeFactors(countsN) <- rep(1, length(sampleNames(phenoData(countsN))))
  countsN <- estimateDispersions(countsN)
  res <- nbinomTest(countsN, "FT", colnames(HGSOC)[i])
  
  
  
  counts <- as(dataNorm, "CountDataSet")
  sizeFactors(counts) <- rep(1,4)
  counts <- estimateDispersions(counts)
  res <- nbinomTest(counts, "wt", "mut")
  
} else {
  
  # make Counts a matrix for DESeq2 and order so HRD ctl columns are first
  cMatrix <- as.matrix(Counts)
  
  # input matrix with coldata to DESeq data set:
  dds <- DESeqDataSetFromMatrix(
    countData = cMatrix,
    colData = coldata,
    design = ~condition
  )
  
  # relevel condition so HRD will be used as control:
  dds$condition <- factor(dds$condition, levels = c("HRD", "CCNEamp"))
  dds <- DESeq(dds)
  res <- results(dds)
  
  res <- as.data.frame(res)

  saveRDS(res, paste0(newRobjectDir, "/DESeq_res.rds"))
}


### 5. Calculate differential expression values of repeats ###

comp <- "CCNEamp_vs_HRD"	

if ("repeats" %in% resultTypes) {
  # define repeat and sig DE repeat dfs:
	repGenes <- as.data.frame(res[grep("ENS",  rownames(res), invert = T),])
	print(repGenes)
	  
	if ( FCthresh == 0 ) {
	  sigGenes <- dplyr::filter(repGenes, padj < padjthresh)
	  repGenes$threshold <- as.factor(repGenes$padj < padjthresh)
	} else {
	  sigGenes <- dplyr::filter(repGenes, (padj < padjthresh & log2FoldChange < -(FCthresh))|(padj < DRthresh & log2FoldChange > FCthresh))
	  repGenes$threshold <- as.factor((repGenes$padj < padjthresh & repGenes$log2FoldChange < -(FCthresh))|(repGenes$padj <  padjthresh & repGenes$log2FoldChange > FCthresh))
	}
	
	sig <- subset(repGenes, threshold == T)
	
	# include the control genes for labelling:
	for (j in 1:length(posGeneIDs)) {
	  if (j==1) {
	    posGenes <- res[ posGeneIDs[j],]
	  } else {
	    posGenes <- rbind(posGenes,   res[posGeneIDs[j],])
	  }
	}
	rownames(posGenes) <- posGeneNames
	
	for (j in 1:length(negGeneIDs)) {
	  if (j==1) {
	    negGenes <- res[ negGeneIDs[j],]
	  } else {
	    negGenes <- rbind(negGenes,   res[negGeneIDs[j],])
	  }
	}
	rownames(negGenes) <- negGeneNames
	
	# set default threshold statuses  for control genes:
	posGenes$threshold <- "POSITIVE"
	if (nrow(posGenes[posGenes$padj< padjthresh,])>0) {
	  posGenes[posGenes$padj<  padjthresh,]$threshold <- "POSSIG"
	}
	
	negGenes$threshold = "NEGATIVE"
	if (nrow(negGenes[negGenes$padj< padjthresh,])>0) {
	  negGenes[negGenes$padj<  padjthresh,]$threshold <-  "NEGSIG"
	}
	
	lab <- rbind(rbind(sig,   posGenes), negGenes)
	repGenes <- rbind(rbind(repGenes,   posGenes), negGenes)
	lab$genes <- rownames(lab)
	
	# plot on volcano plot:
	p <- ggplot(data=repGenes, aes( x=log2FoldChange, y=-log10(padj), color=threshold))
	p <- p + geom_point(data=repGenes)
	p <- p + geom_text_repel(data=lab, aes(label=genes))
	p <- p + theme(legend.position =  "none")
	p <- p + labs(x="log2 fold change   vs FT control", y="-log10   padj")
	p <- p +  xlim(c(-4, 4))
	if (FCthresh == 0) {
	  if (file.exists(paste0(plotDir,   "/", Type,  "_volcano_padj_",   padjthresh, "_", comp, ".pdf"))) {
	    print(paste0(plotDir, "/",  Type,  "_volcano_padj_",   padjthresh, "_", comp, ".pdf"))
      p
	  } else {
	   print(paste0("Creating  ", plotDir, "/", Type,  "_volcano_padj_", padjthresh, "_", comp, ".pdf"))
	   pdf(file = paste0(plotDir, "/",   Type,  "_volcano_padj_",  padjthresh, "_", comp, ".pdf"))
	   print(p)
	   dev.off()
    }
  } else {
    if (file.exists(paste0(plotDir, "/",  Type,  "_volcano_padj",   padjthresh, "_FC", FCthresh, "_", comp, ".pdf"))) {
	    print(paste0(plotDir, "/",  Type,  "_volcano_padj",   padjthresh, "_FC", FCthresh, "_", comp, ".pdf already exists"))
	    p
	  } else {
	    print(paste0("Creating  ", plotDir, "/",  Type,  "_volcano_padj", padjthresh, "_FC", FCthresh, "_", comp, ".pdf"))
	    pdf(file = paste0(plotDir, "/",  Type,  "_volcano_padj",   padjthresh, "_FC", Cthresh, "_", comp, ".pdf"))
	   print(p)
	   dev.off()
	  }
  }
}

if ("all" %in% resultTypes) {
  
  if (FCthresh == 0) {
    sigGenes <- filter(res, padj < padjthresh)
    res$threshold <- as.factor(res$padj < padjthresh)
  } else {
    sigGenes <- filter(res, (padj < padjthresh & log2FoldChange < -(FCthresh))|(padj < padjthresh & log2FoldChange > FCthresh))
    res$threshold <- as.factor((res$adj < padjthresh & res$log2FoldChange < -(FCthresh))|(res$padj <  padjthresh & res$log2FoldChange > FCthresh))
  }
  
  sig <- subset(res, threshold == T)
  
  # include the control genes for labelling:
  for (j in 1:length(posGeneIDs)) {
    if (j==1) {
      posGenes <- res[ posGeneIDs[j],]
    } else {
      posGenes <- rbind(posGenes,   res[posGeneIDs[j],])
    }
  }
  rownames(posGenes) <- posGeneNames
  for (j in 1:length(negGeneIDs)) {
    if (j==1) {
      negGenes <- res[ negGeneIDs[j],]
    } else {
      negGenes <- rbind(negGenes,   res[negGeneIDs[j],])
    }
  }
  rownames(negGenes) <- negGeneNames
  
  # set default threshold statuses  for control genes:
  posGenes$threshold <- "POSITIVE"
  if (nrow(posGenes[posGenes$padj< padjthresh,])>0) {
    posGenes[posGenes$padj<  padjthresh,]$threshold <- "POSSIG"
  }
  
  negGenes$threshold = "NEGATIVE"
  if (nrow(negGenes[negGenes$padj< padjthresh,])>0) {
    negGenes[negGenes$padj<  padjthresh,]$threshold <-  "NEGSIG"
  }
  
  lab <- rbind(rbind(sig, posGenes), negGenes)
  lab <- lab[( lab$log2FoldChange > 10 | lab$log2FoldChange < -10 | lab$padj < 5e-4 ),]
  res <- rbind(rbind(res,   posGenes), negGenes)
  lab$genes <- rownames(lab)
  lab <- as.data.frame(lab)
  
  # add gene symbol annotations where relevant:
  egENSEMBL <- toTable(org.Hs.egENSEMBL)
  egSYMBOL <- toTable(org.Hs.egSYMBOL)
  
  # for rows with ensembl ids, annotate entrez ids and symbols in separate columns
  # for lab and sig:
  lab$gene_id <- egENSEMBL$gene_id[match(rownames(lab), egENSEMBL$ensembl_id)]
  lab$symbol <- egSYMBOL$symbol[match(lab$gene_id, egSYMBOL$gene_id)]
  
  sig$gene_id <- egENSEMBL$gene_id[match(rownames(sig), egENSEMBL$ensembl_id)]
  sig$symbol <- egSYMBOL$symbol[match(sig$gene_id, egSYMBOL$gene_id)]
  
  res$gene_id <- egENSEMBL$gene_id[match(rownames(res), egENSEMBL$ensembl_id)]
  res$symbol <- egSYMBOL$symbol[match(res$gene_id, egSYMBOL$gene_id)]
  # save CCNEamp_vs_HRD for comparison with Patch 2015 data:
  saveRDS(res, paste0(newRobjectDir, "/CCNEamp_vs_HRD_res.rds"))
  saveRDS(sig, paste0(newRobjectDir, "/CCNEamp_vs_HRD_sig.rds"))
  
  # plot on volcano plot:
  p <- ggplot(data=res, aes( x=log2FoldChange, y=-log10(padj), color=threshold) )
  p <- p + geom_point(data=res)
  p <- p + geom_text_repel(data=lab, aes(label=symbol))
  p <- p + theme(legend.position =  "none")
  p <- p + labs(x="log2 fold change   vs FT control", y="-log10   padj")
  if (length(FCthresh) == 0) {
    if (file.exists(paste0(plotDir,   "/", Type,  "_volcano_padj_5e_neg4_", comp, "_all.pdf"))) {
      print(paste0(plotDir, "/",  Type,  "_volcano_padj_5e_neg4_", comp, "_all.pdf"))
      p
    } else {
      print(paste0("Creating  ",plotDir, "/", Type,    "_volcano_padj_5e_neg4_", comp, "_all.pdf"))
      pdf(file = paste0(plotDir, "/",   Type,  "_volcano_padj_5e_neg4_", comp, "_all.pdf"))
      print(p)
      dev.off()
    }
  } else {
    if (file.exists(paste0(plotDir, "/",  Type,  "_volcano_padj_5e_neg4_", comp, "_all.pdf"))) {
      print(paste0(plotDir, "/",  Type,  "_volcano_padj_5e_neg4_", comp, "_all.pdf already exists"))
      p
    } else {
      print(paste0("Creating  ", plotDir, "/",  Type,  "_volcano_padj_10e_neg15_", comp, "_res.pdf"))
      pdf(file = paste0(plotDir, "/",  Type,  "_volcano_padj_10e_neg15_", comp, "_res.pdf"))
      print(p)
      dev.off()
    }
  }
  
  # save CCNEamp_vs_HRD for comparison with Patch 2015 data:
  if ( comp == "CCNEamp_vs_HRD") {
    saveRDS(res, paste0(newRobjectDir, "/CCNEamp_vs_HRD_allGenes.rds"))
    saveRDS(sig, paste0(newRobjectDir, "/CCNEamp_vs_HRD_sig.rds"))
    write.table(DEs, file=paste0(plotDir, "/DE_CCNEamp_vs_HRD.txt"), row.names = F, col.names = F, quote = F, sep = "\t")
  }
  
}




