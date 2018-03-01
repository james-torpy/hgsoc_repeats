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

Type <- "custom3"

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
                     expName, "/")
newRobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
                     expName, "/HGSOC_vs_HGSOC_bowtell_cats/")
plotDir <- paste0(resultsDir, "/R/", expName,
                  "/plots/DEplots/HGSOC_vs_HGSOC_bowtell_cats/")

system(paste0("mkdir -p ", plotDir))
system(paste0("mkdir -p ", RobjectDir))
system(paste0("mkdir -p ", newRobjectDir))

### 1. Load in all counts ###

custom3Counts <- readRDS(paste0(RobjectDir, "/", Type,
                                "_allcounts.htseq.rds"))
gcCounts <- readRDS(paste0(RobjectDir,
                           "/gc_allcounts.htseq.rds"))

sTypes <- unique(
  grep(
    "id", gsub(
      "^.*\\_", "", colnames(custom3Counts)
    ), value=T, invert = T
  )
)

# append gcCounts to custom3Counts:
Counts <- rbind(custom3Counts, gcCounts)

# remove FT and pAF samples:
Counts <- Counts[,-grep("FT", colnames(Counts))]
Counts <- Counts[,-grep("pAF", colnames(Counts))]
sTypes <- sTypes[-grep("FT", sTypes)]
sTypes <- sTypes[-grep("pAF", sTypes)]

# make rownames gene_id, get rid of latter column and change
# storage mode from factor to integer:
rownames(Counts) <- Counts$gene_id
Counts <- subset(Counts, select=-gene_id)
#storage.mode(Counts <- "integer")

# save Counts:
saveRDS(Counts, file=paste0(RobjectDir, "/", Type, "_counts.RData"))


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

if (file.exists(paste0(plotDir, "htseq_", Type, "_pcaCompsPrenormGC.pdf"))) {
  print(paste0(plotDir, "htseq_", Type, "_pcaCompsPrenormGC.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", plotDir, "htseq_", Type, "_pcaCompsPrenormGC.pdf"))
  pdf(file = paste0(plotDir, "htseq_", Type, "_pcaCompsPrenormGC.pdf"))
  plot(pca)
  dev.off()
}

# change the order of columns of Counts to alphabetical order:
colOrder <- unlist(split(colnames(Counts), gsub("^.*_*_", "", colnames(Counts))))
Counts <- Counts[,colOrder]

splt <- unlist(lapply(split(colnames(Counts), gsub("^.*_*_", "", colnames(Counts))), length))
for (i in 1:length(splt)) {
  if (i==1) {
    typeF <- c(rep(names(splt)[i], splt[i]))
  } else {
    typeF <- c(typeF, rep(names(splt)[i], splt[i]))
  }
}
levels(typeF) <- sTypes

# convert Counts into SeqExpressionSet - elements need to be delisted and changed to integers first:
Counts <- apply(Counts, 2, unlist)
storage.mode(Counts) <- "integer"
set <- newSeqExpressionSet(Counts, phenoData = data.frame(typeF, row.names=colnames(Counts)))

# create pre-norm RLE plot:
if (file.exists(paste0(plotDir, "htseq_", Type, "_RLEPrenormGC.pdf"))) {
  print(paste0(plotDir, "htseq_", Type, "_RLEPrenormGC.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", plotDir, "htseq_", Type, "_RLEPrenormGC.pdf"))
  par(mar=c(1,1,1,1))
  pdf(file = paste0(plotDir, "htseq_", Type, "_RLEPrenormGC.pdf"))
  plotRLE(set)
  dev.off()
}

# create RUVseq pre-norm PCA:
if (file.exists(paste0(plotDir, "htseq_", Type, "_pcaPrenormGC.pdf"))) {
  print(paste0(plotDir, "htseq_", Type, "_pcaPrenormGC.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", plotDir, "htseq_", Type, "_pcaPrenormGC.pdf"))
  pdf(file = paste0(plotDir, "htseq_", Type, "_pcaPrenormGC.pdf"), height = 10, width = 12)
  plotPCA(set, cex=0.7)
  dev.off()
}

### 3. perform normalisation on counts using RUVseq:

# perform between lane full normalisation:
nSet <- betweenLaneNormalization(set, which="full")

# create post-norm RLE plot:
if (file.exists(paste0(plotDir, "htseq_", Type, "_RLElaneNormGC.pdf"))) {
  print(paste0(plotDir, "htseq_", Type, "_RLElaneNormGC.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", plotDir, "htseq_", Type, "_RLElaneNormGC.pdf"))
  pdf(file = paste0(plotDir, "htseq_", Type, "_RLElaneNormGC.pdf"))
  plotRLE(nSet, outline=FALSE, ylim=c(-4, 4))
  dev.off()
}

# create RUVseq post-norm PCA:
if (file.exists(paste0(plotDir, "htseq_", Type, "_pcalaneNormGC.pdf"))) {
  print(paste0(plotDir, "htseq_", Type, "_pcalaneNormGC.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", plotDir, "htseq_", Type, "_pcalaneNormGC.pdf"))
  pdf(file = paste0(plotDir, "htseq_", Type, "_pcalaneNormGC.pdf"), height = 15, width = 20)
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

saveRDS(fit, file=paste0(RobjectDir, "/", Type, "DEfit.rds"))

# order sTypes:
sTypes <- sTypes[order(sTypes)]

# determine which column has FT control:
for (ctlInd in 1:length(sTypes)) {
  con <- c(rep(0, (ctlInd - 1) ), -1, rep(0, (ncol(design) - ctlInd)))
  
  for (i in 1:ncol(design)) {
    if (i!=ctlInd) {
      comp <- paste0(sTypes[ctlInd], "_vs_", sTypes[i])
      
      # perform likelihood ratio test:
      con[i] <- 1
      lrt <- glmLRT(fit, contrast = con)
      
      # determine the top DE genes:
      topTags(lrt)
      
      if (file.exists(paste0(RobjectDir, sTypes[ctlInd], "_vs_", sTypes[i],   "_lrt.rds"))) {
        print(paste0(RobjectDir, sTypes[ctlInd], "_vs_", sTypes[i],   "_lrt.rds already exists, no need to create"))
      } else {
        print(paste0("Creating ", RobjectDir, sTypes[ctlInd], "_vs_", sTypes[i],   "_lrt.rds"))
        saveRDS(lrt, file = paste0(RobjectDir, sTypes[ctlInd], "_vs_", sTypes[i],   "_lrt.rds"))
      }
      
      
      ### 5. Calculate differential expression values with FDR threshold only ###
      
      # fetch summary of differentially expressed genes (those  with FDR =< 0.05:
      DEs <- summary(result <- decideTestsDGE((lrt)))
      
      # fetch all gene DE info, 
      thresh <- 0.05
      allGenes <- as.data.frame(topTags(lrt, n=Inf))
      repGenes <- allGenes[grep("ENS", rownames(allGenes), invert   = T),]
      sigGenes <- filter(repGenes, FDR<thresh)
      print(repGenes)
      
      # plot on volcano plot:
      repGenes$threshold <- as.factor(repGenes$FDR < thresh)
      sig <- subset(repGenes, threshold == T)
      
      # include the control genes for labelling:
      # add positive and negative control genes CD47, CCNE1,  GAPDH, b-actin:
      posGenes <- rbind(allGenes["ENSG00000196776",],   allGenes["ENSG00000105173",])
      negGenes <- rbind(allGenes["ENSG00000075624",],   allGenes["ENSG00000204574",],
                        allGenes["ENSG00000104904",], allGenes["ENSG00000089157",])
      rownames(posGenes) <- c("CD47", "CCNE1")
      rownames(negGenes) <- c("beta-actin", "ABCF1",   "OAZ1", "RPLP0")
      
      posGenes$threshold <- "POSITIVE"
      if (nrow(posGenes[posGenes$FDR<thresh,])>0) {
        posGenes[posGenes$FDR<thresh,]$threshold <- "POSSIG"
      }
      
      negGenes$threshold = "NEGATIVE"
      if (nrow(negGenes[negGenes$FDR<thresh,])>0) {
        negGenes[negGenes$FDR<thresh,]$threshold <- "NEGSIG"
      }
      
      lab <- rbind(rbind(sig, posGenes), negGenes)
      repGenes <- rbind(rbind(repGenes, posGenes), negGenes)
      lab$genes <- rownames(lab)
      
            p <- ggplot(data=repGenes, aes(x=logFC, y=-log10(FDR),  color=threshold))
      p <- p + geom_point(data=repGenes)
      p <- p + geom_text_repel(data=lab, aes(label=genes))
      p <- p + theme(legend.position = "none")
      p <- p + labs(x="log2 fold change vs FT control", y="-log10   FDR")
      p <- p +  xlim(c(-7, 7))
      if (file.exists(paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, comp, sTypes[ctlInd], "_control.pdf"))) {
        print(paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, comp, sTypes[ctlInd], "_control.pdf"))
        p
      } else {
        print(paste0("Creating ",plotDir, "/", Type,  "_volcano_FDR", thresh, comp, sTypes[ctlInd], "_control.pdf"))
        pdf(file = paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, comp, sTypes[ctlInd], "_control.pdf"))
        print(p)
        dev.off()
      }
      
      ### 6. Calculate differential expression values with FDR and FC thresholds ###
      
      # PCA and RLE plots of RUVseq RUVr-normalised data looked   best with etween lane normalisation = 	'full',  will go  with this #
      
      # fetch summary of differentially expressed genes (those  with FDR =< 0.05:
      DEs <- summary(result <- decideTestsDGE((lrt)))
      
      # fetch all gene DE info, 
      thresh <- 0.05
      FCthresh <- 1
      allGenes <- as.data.frame(topTags(lrt, n=Inf))
      repGenes <- allGenes[grep("ENS", rownames(allGenes), invert = T),]
      sigGenes <- filter(repGenes, (FDR < thresh & logFC < -(FCthresh))|(FDR < thresh & logFC > FCthresh))
      print(repGenes)
      
      # plot on volcano plot:
      repGenes$threshold <- as.factor((repGenes$FDR < thresh & repGenes$logFC < -(FCthresh))|(repGenes$FDR <  thresh & repGenes$logFC > FCthresh))
      sig <- subset(repGenes, threshold == T)
      
      # include the control genes for labelling:
      # add positive and negative control genes CD47, CCNE1,  GAPDH, b-actin:
      
      posGenes <- rbind(allGenes["ENSG00000091831",],   allGenes["ENSG00000196776",],   allGenes["ENSG00000105173",])
      negGenes <- rbind(allGenes["ENSG00000111640",],   allGenes["ENSG00000075624",],   allGenes["ENSG00000204574",],
                        allGenes["ENSG00000104904",],   allGenes["ENSG00000089157",])
      rownames(posGenes) <- c("ESR1", "CD47", "CCNE1")
      rownames(negGenes) <- c("GAPDH", "beta-actin", "ABCF1",   "OAZ1", "RPLP0")
      
      posGenes$threshold <- "POSITIVE"
      if (nrow(posGenes[posGenes$FDR<thresh,])>0) {
        posGenes[posGenes$FDR<thresh,]$threshold <- "POSSIG"
      }
      
      negGenes$threshold = "NEGATIVE"
      if (nrow(negGenes[negGenes$FDR<thresh,])>0) {
        negGenes[negGenes$FDR<thresh,]$threshold <- "NEGSIG"
      }
      
      lab <- rbind(rbind(sig, posGenes), negGenes)
      repGenes <- rbind(rbind(repGenes, posGenes), negGenes)
      lab$genes <- rownames(lab)
      
      p <- ggplot(data=repGenes, aes(x=logFC, y=-log10(FDR),  color=threshold))
      p <- p + geom_point(data=repGenes)
      p <- p + geom_text_repel(data=lab, aes(label=genes))
      p <- p + theme(legend.position = "none")
      p <- p + labs(x="log2 fold change vs FT control", y="-log10   FDR")
      p <- p +  xlim(c(-7, 7))
      if (file.exists(paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, "_logFC", FCthresh, "_", comp,  ".pdf"))) {
        print(paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, "_logFC", FCthresh, "_", comp, sTypes[ctlInd], "_control.pdf"))
        p
      } else {
        print(paste0("Creating ",plotDir, "/", Type,  "_volcano_FDR", thresh, "_logFC", FCthresh, "_", comp,  ".pdf"))
        pdf(file = paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, "_logFC", FCthresh, "_", comp, sTypes[ctlInd], "_control.pdf"))
        print(p)
        dev.off()
      }
      
      ### 7. Identify epigenetic and RNAi regulators and assess expression: ###
      
      epiIDs <- c("ENSG00000130816", "ENSG00000119772", "ENSG00000088305",  
                  "ENSG00000276043", "ENSG00000138336", "ENSG00000168769", "ENSG00000187605", 
                  "ENSG00000100697", "ENSG00000092847", "ENSG00000123908", "ENSG00000126070", 
                  "ENSG00000134698", "ENSG00000101945")
      
      epiSym <- c("DNMT1", "DNMT3A", "DNMT3B", "UHRF1", "TET1", "TET2", "TET3", 
                  "Dicer", "AGO1", "AGO2", "AGO3", "AGO4", "SUV39H1")
      
      epiGenesDF <- allGenes[epiIDs,]
      rownames(epiGenesDF) <- epiSym
      
      epiThresh <- 0.1
      epiGenesDF$threshold <- as.factor((epiGenesDF$FDR < epiThresh & epiGenesDF$logFC < -(FCthresh))|(epiGenesDF$FDR <  epiThresh & epiGenesDF$logFC > FCthresh))
      epiSig <- subset(epiGenesDF, threshold == T)
      
      epiGenesDF$genes <- rownames(epiGenesDF)
      
      p <- ggplot(data=epiGenesDF, aes(x=logFC, y=-log10(FDR),color=threshold))
      p <- p + geom_point(data=epiGenesDF)
      p <- p + geom_text_repel(data=epiGenesDF, aes(label=genes))
      p <- p + theme(legend.position = "none")
      p <- p + labs(x="log2 fold change vs FT control", y="-log10   FDR")
      if (file.exists(paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, "_logFC", FCthresh, "_", comp,  "_epi.pdf"))) {
        print(paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, "_logFC", FCthresh, "_", comp, "_control_epi.pdf"))
        p
      } else {
        print(paste0("Creating ",plotDir, "/", Type,  "_volcano_FDR", thresh, "_logFC", FCthresh, "_", comp, "_control_epi.pdf"))
        pdf(file = paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, "_logFC", FCthresh, "_", comp, "_control_epi.pdf"))
        print(p)
        dev.off()
      }
      
      saveRDS(repGenes, file=paste0(RobjectDir, "/", Type, "_allDErep", thresh, "_logFC", FCthresh, "_", comp, "_control.rds"))
      
      if (!(ctlInd==1)) {
        if (i==1) {
          allReps <- list(repGenes)
          sigReps <- list(sig)
          epiGenes <- list(epiGenesDF)
        } else {
          allReps[[i]] <- repGenes
          sigReps[[i]] <- sig
          epiGenes[[i]] <- epiGenesDF
        }
      } else {
        if (i==2) {
          allReps <- list(NULL)
          allReps[[i]] <- repGenes
          sigReps <- list(NULL)
          sigReps[[i]] <- sig
          epiGenes <- list(NULL)
          epiGenes[[i]] <- epiGenesDF
        } else {
          allReps[[i]] <- repGenes
          sigReps[[i]] <- sig
          epiGenes[[i]] <- epiGenesDF
        }
      }
      
      con <- c(rep(0, (ctlInd - 1) ), -1, rep(0, (ncol(design) - ctlInd)))
    }
  }
  
  ### 8. Save all repeat DE, and all repeat DE that is significant and log2 FC > 1 for at least one comparison:
  
  # remove element 3 of lists as corresponds to FT_vs_FT:
  allReps <- c(allReps[-ctlInd])
  names(allReps) <- paste0(sTypes[ctlInd], "_vs_", sTypes[-ctlInd])
  
  sigReps <- c(sigReps[-ctlInd])
  names(sigReps) <- paste0(sTypes[ctlInd], "_vs_", sTypes[-ctlInd])

  epiGenes <- c(epiGenes[-ctlInd])
  names(epiGenes) <- paste0(sTypes[ctlInd], "_vs_", sTypes[-ctlInd])
  
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
  
  if (ctlInd==1) {
    allRepsTotal <- list(allReps)
    allSigTotal <- list(allSig)
    epiGenesTotal <- list(epiGenes)
  } else {
    allRepsTotal[[ctlInd]] <- allReps
    allSigTotal[[ctlInd]] <- allSig
    epiGenesTotal[[ctlInd]] <- epiGenes
  }
}

# fetch names of all significant repeats across all elements of allSigsTotal:
j=1
for (e in allSigTotal) {
  print(j)
  for (i in 1:length(e)) {
    print(i)
    if (j==1 & i==1) {
      totalNames <- c(rownames(e[[i]]))
    } else {
      totalNames <- append(totalNames, rownames(e[[i]]))
    }
  }
  j <<- j+1
}

totalNames <- unique(totalNames)

# standardise all elements of allSigsTotal across all elements in terms of repeats included:
allSigTotalNew <- lapply(allRepsTotal, function(x) {
  return(lapply(x, function(y) {
    return(y[totalNames,])
  }))
})

allSigTotalNewReduced <- allSigTotalNew


# delete duplicate comparisons:
######

# for each element of parent list, delete occurances of the same comparisons in subsequent elements:
#i=1
#for (x in allGene) {
#  for (i in 1:length(x)) {
#    t1 <- strsplit(names(x)[i], "_")[[1]][1]
#    t2 <- strsplit(names(x)[i], "_")[[1]][3]
#    x <- lapply(allGene[-i], function(y) {
#      index1 <- grep(t1, names(y))
#      index2 <- grep(t2, names(y))
#      delete_index <- grep(index1, index2)
#      return(y[-delete_index])
#    })
#    allGene[[1]] <- x
#  }
#i <<- i+1
#}

######

# delete duplicate comparisons:
allSigTotalNewReduced[[2]] <- allSigTotalNewReduced[[2]][2:6]
allSigTotalNewReduced[[3]] <- allSigTotalNewReduced[[3]][3:6]
allSigTotalNewReduced[[4]] <- allSigTotalNewReduced[[4]][4:6]
allSigTotalNewReduced[[5]] <- allSigTotalNewReduced[[5]][5:6]
allSigTotalNewReduced[[6]] <- allSigTotalNewReduced[[6]][6]
allSigTotalNewReduced <- allSigTotalNewReduced[1:6]

allSigTotalNewReduced <- unlist(allSigTotalNewReduced, recursive=F)

if (!file.exists(paste0(newRobjectDir, "/", Type, "_DEallReps.rds"))) {
  saveRDS(allRepsTotal, file=paste0(newRobjectDir, "/", Type, "_DEallReps.rds"))
}

if (!file.exists(paste0(newRobjectDir, "/", Type, "_DEsigReps.rds"))) {
  saveRDS(allSigTotalNew, file=paste0(newRobjectDir, "/", Type, "_DEsigReps.rds"))
}

if (!file.exists(paste0(newRobjectDir, "/", Type, "_DEsigRepsReduced.rds"))) {
  saveRDS(allSigTotalNewReduced, file=paste0(newRobjectDir, "/", Type, "_DEsigRepsReduced.rds"))
}

if (!file.exists(paste0(newRobjectDir, "/", Type, "_DEepiGenes.rds"))) {
  saveRDS(epiGenesTotal, file=paste0(newRobjectDir, "/", Type, "_DEepiGenes.rds"))
}

if (!file.exists(paste0(newRobjectDir, "DEImg_", expName, "_HGSOC_vs_HGSOC_bowtell_cats.RData"))) {
  save.image(file = paste0(newRobjectDir, "DEImg_", expName, "_HGSOC_vs_HGSOC_bowtell_cats.RData"))
}

# save number of samples in each category:
readRDS(paste0(newRobjectDir, "/sample_no_per_cat.rds"))

if (!file.exists(paste0(newRobjectDir, "/sample_no_per_cat.rds"))) {
  saveRDS(splt, file=paste0(newRobjectDir, "/sample_no_per_cat.rds"))
}