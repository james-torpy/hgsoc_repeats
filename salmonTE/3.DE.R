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
expName <- "salmonTE"

Type <- "custom3"
sTypes <- c("FT", "HGSOC")
descript <- "SalmonTE_HGSOC_vs_FT"
#sTypes <- c("FT", "recurrent_ascites", "primary_resistant", "refractory", 
#            "acquired_resistant", "multiple_responder", "extreme_responder", 
#            "metastatic")
#descrip <- "htseq_hgsoc_split_more_by_drug_response_vs_FT"

# define sample groups to compare:
sGroups <- list("FT", "HGSOC")
names(sGroups) <- sTypes
cat_by_driver <- FALSE

# define sample group to use as control:
ctl <- "FT"

# specify what combination of repeat genes (repeats), epigenetic modulators (epiMods),
# RNAi genes (RNAi) and protein-coding genes (pCoding) should contribute to the results:
#resultTypes <- c("repeats", "epiMods")
resultTypes <- c("repeats")

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

# if necessary, divide samples by driver mutation:
# load in sample key for categories homologous repair deficient (HRD) and cyclin E gain/amplification (CCNE):
HRDkey <- read.table(paste0(rawDir, "/hrd_samples.txt"), header=F, sep="\t")
HRDnos <- gsub("AOCS_", "", HRDkey[,1])
CCNEkey <- read.table(paste0(rawDir, "/ccne_gain_or_amp_samples.txt"), header=F, sep="\t")
CCNEnos <- gsub("AOCS_", "", CCNEkey[,1])

# re-categorize samples as HRD, CCNE_amp, both_drivers or unknown_drivers:
if (cat_by_driver) {
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

# change the order of columns of Counts to alphabetical order:
Counts <- Counts[,order(
  gsub(
    "AOCS.*_[0-9][0-9][0-9]_", "", colnames(Counts)
  )
)]

# define sample groups:
splt <- unlist(
  lapply(
    split(
      colnames(Counts), gsub(
        "AOCS.*_[0-9][0-9][0-9]_", "", colnames(Counts)
      )
    ), length
  )
)

for (i in 1:length(splt)) {
  if (i==1) {
    typeF <- c(rep(names(splt)[i], splt[i]))
  } else {
    typeF <- c(typeF, rep(names(splt)[i], splt[i]))
  }
}
levels(typeF) <- sTypes

# save number of samples in each group:
saveRDS(splt, file = paste0(newRobjectDir, "/sample_no_per_cat.rds"))

# convert Counts into SeqExpressionSet - elements need to be delisted and changed to integers first:
Counts <- apply(Counts, 2, unlist)
storage.mode(Counts) <- "integer"
set <- newSeqExpressionSet(Counts, phenoData = data.frame(typeF, row.names=colnames(Counts)))

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
  plotPCA(set, cex=0.7)
  dev.off()
}


### 3. perform normalisation on counts using RUVseq:

# perform between lane full normalisation:
nSet <- betweenLaneNormalization(set, which="full")

# create post-norm RLE plot:
if (file.exists(paste0(plotDir, "/", Type, "_RLElaneNormGC.pdf"))) {
  print(paste0(plotDir, "/", Type, "_RLElaneNormGC.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", plotDir, "/", Type, "_RLElaneNormGC.pdf"))
  pdf(file = paste0(plotDir, "/", Type, "_RLElaneNormGC.pdf"))
  plotRLE(nSet, outline=FALSE, ylim=c(-4, 4))
  dev.off()
}

# create RUVseq post-norm PCA:
if (file.exists(paste0(plotDir, "/", Type, "_pcalaneNormGC.pdf"))) {
  print(paste0(plotDir, "/", Type, "_pcalaneNormGC.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", plotDir, "/", Type, "_pcalaneNormGC.pdf"))
  pdf(file = paste0(plotDir, "/", Type, "_pcalaneNormGC.pdf"), height = 15, width = 20)
  plotPCA(nSet, cex=0.7)
  dev.off()
}


### 4. Perform differential expression comparing normalised FT controls to cancer samples ###

# design matrix labelling all sample types:
design <- model.matrix(~0+typeF, data=pData(nSet))

# convert set into a DGElist format which specifies library size and normalisation factors:
y <- DGEList(counts=counts(nSet), group=typeF)

# create MDS plot:
if (file.exists(paste0(plotDir, Type, "_mdslaneNormGC.pdf"))) {
  print(paste0(plotDir, Type, "_mdslaneNormGC.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", plotDir, Type, "_mdslaneNormGC.pdf"))
  pdf(file = paste0(plotDir, Type, "_mdslaneNormGC.pdf"), height = 15, width = 20)
  plotMDS(y)
  dev.off()
}

# calculate the deviance residuals from a first-pass GLM regression of the counts on the co-variates of interest (p8 RUVseq manual):
# estimate dispersion:
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
# adjust values using dispersion:
fit <- glmFit(y, design)

saveRDS(fit, file=paste0(newRobjectDir, "/", Type, "DEfit.rds"))
save.image(paste0(newRobjectDir, "/", Type, "DEdone.rds"))

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
        if (file.exists(paste0(plotDir,   "/", Type,  "_volcano_FDR_",   FDRthresh, "_", comp, ".pdf"))) {
          print(paste0(plotDir, "/",  Type,  "_volcano_FDR_",   FDRthresh, "_", comp, ".pdf"))
          p
        } else {
          print(paste0("Creating  ",plotDir, "/", Type,    "_volcano_FDR_", FDRthresh, "_", comp, ".pdf"))
          pdf(file = paste0(plotDir, "/",   Type,  "_volcano_FDR_",  FDRthresh, "_", comp, ".pdf"))
          print(p)
          dev.off()
        }
      } else {
        if (file.exists(paste0(plotDir, "/",  Type,  "_volcano_FDR",   FDRthresh, "_FC", FCthresh, "_", comp, ".pdf"))) {
          print(paste0(plotDir, "/",  Type,  "_volcano_FDR",   FDRthresh, "_FC", FCthresh, "_", comp, ".pdf already exists"))
          p
        } else {
          print(paste0("Creating  ", plotDir, "/",  Type,  "_volcano_FDR", FDRthresh, "_FC", FCthresh, "_", comp, ".pdf"))
          pdf(file = paste0(plotDir, "/",  Type,  "_volcano_FDR",   FDRthresh, "_FC", FCthresh, "_", comp, ".pdf"))
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
          if (file.exists(paste0(plotDir,   "/", Type,  "_volcano_FDR_",   FDRthresh, "_", comp, "_epi.pdf"))) {
            print(paste0(plotDir, "/",  Type,  "_volcano_FDR_",   FDRthresh, "_", comp, "_epi.pdf"))
            p
          } else {
            print(paste0("Creating  ",plotDir, "/", Type,    "_volcano_FDR_", FDRthresh, "_", comp, "_epi.pdf"))
            pdf(file = paste0(plotDir, "/",   Type,  "_volcano_FDR_",  FDRthresh, "_", comp, "_epi.pdf"))
            print(p)
            dev.off()
          }
        } else {
          if (file.exists(paste0(plotDir, "/",  Type,  "_volcano_FDR",   FDRthresh, "_FC", FCthresh, "_", comp, "_epi.pdf"))) {
            print(paste0(plotDir, "/",  Type,  "_volcano_FDR",   FDRthresh, "_FC", FCthresh, "_", comp, "_epi.pdf already exists"))
            p
          } else {
            print(paste0("Creating  ", plotDir, "/",  Type,  "_volcano_FDR", FDRthresh, "_FC", FCthresh, "_", comp, "_epi.pdf"))
            pdf(file = paste0(plotDir, "/",  Type,  "_volcano_FDR",   FDRthresh, "_FC", FCthresh, "_", comp, "_epi.pdf"))
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
allReps <- allReps[-ctlInd]
sigReps <- sigReps[-ctlInd]

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
