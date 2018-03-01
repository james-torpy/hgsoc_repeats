### 13c.DE_hrd_vs_ccne_unknown_and_ccne_vs_hrd_unknown.R ###

# This script takes a list of dfs of different classes and types of repeat counts for
# HGSOC and HRD control RNA-seq data sets and performs DE analysis:

#module load briglo/R/3.4.2

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
ctl <- "HRD"

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
                     expName, "/")
plotDir <- paste0(resultsDir, "/R/", expName,
                  "/plots/DEplots/hrd_vs_ccne_unknown/")

homeDir2 <- "/Users/jamestorpy/clusterHome2/"
projectDir2 <- paste0(homeDir2, "/projects/", project)
rawDir <- paste0(projectDir2, "/RNA-seq/raw_files/")
newRobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
                        expName, "/hrd_vs_ccne_unknown_and_ccne_vs_hrd_unknown/")

system(paste0("mkdir -p ", newRobjectDir))
system(paste0("mkdir -p ", plotDir))


### 1. Load in all counts ###

custom3Counts <- readRDS(paste0(RobjectDir, "/", Type,
                                "_allcounts.htseq.rds"))
gcCounts <- readRDS(paste0(RobjectDir,
                           "/gc_allcounts.htseq.rds"))

sTypes <- c("unknownDrivers", "HRD", "CCNEAmp", "bothDrivers")

# append gcCounts to custom3Counts:
Counts <- rbind(custom3Counts, gcCounts)

# remove FT samples from Counts:
Counts <- Counts[,-grep("FT", colnames(Counts))]

# make rownames gene_id, get rid of latter column and change
# storage mode from factor to integer:
rownames(Counts) <- Counts$gene_id
Counts <- subset(Counts, select=-gene_id)

# load in sample key for categories homologous repair deficient (HRD) and cyclin E gain/amplification (CCNE):
HRDkey <- read.table(paste0(rawDir, "/hrd_samples.txt"), header=F, sep="\t")
HRDnos <- gsub("AOCS_", "", HRDkey[,1])
CCNEkey <- read.table(paste0(rawDir, "/ccne_gain_or_amp_samples.txt"), header=F, sep="\t")
CCNEnos <- gsub("AOCS_", "", CCNEkey[,1])

# re-categorize samples as HRD, CCNE_amp, both_drivers or unknown_drivers:
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
      colnames(Counts)[i] <- paste0(colnames(Counts)[i], "_CCNEAmp")
    } else if (!(no %in% HRDnos | no %in% CCNEnos)) {
      colnames(Counts)[i] <- paste0(colnames(Counts)[i], "_unknownDrivers")
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

# record number of samples in each category:
saveRDS(splt, file=(paste0(newRobjectDir, "/sample_no_per_cat.rds")))

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


### 4. Perform differential expression comparing normalised HRD controls to cancer samples ###

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

# assign control group:
ctl <- "HRD"
ctlInd <- grep(ctl, colnames(design))
con <- c(rep(0, (ctlInd - 1) ), -1, rep(0, (ncol(design) - ctlInd)))

sTypes <- sTypes[order(sTypes)]

for (i in 1:length(sTypes)) {
  if (i!=ctlInd) {
    comp <- paste0("unknownDrivers_vs_", sTypes[i])
    
    # perform likelihood ratio test:
    con[i] <- 1
    lrt <- glmLRT(fit, contrast = con)
    
    # determine the top DE genes:
    topTags(lrt)
    
    if (file.exists(paste0(newRobjectDir, "HRD_vs_", sTypes[i],   "_lrt.rds"))) {
      print(paste0(newRobjectDir, "HRD_vs_", sTypes[i],   "_lrt.rds already exists, no need to create"))
    } else {
      print(paste0("Creating ", newRobjectDir, "HRD_vs_", sTypes[i],   "_lrt.rds"))
      saveRDS(lrt, file = paste0(newRobjectDir, "HRD_vs_", sTypes[i],   "_lrt.rds"))
    }
    
    
    ### 5. Calculate differential expression values with FDR threshold only ###
    
    # fetch summary of differentially expressed genes (those  with FDR =< 0.05:
    DEs <- summary(result <- decideTestsDGE((lrt)))
    
    # fetch all gene DE info, 
    thresh <- 0.1
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
    
    p <- ggplot(data=repGenes, aes(x=logFC, y=-log10(FDR),  color=threshold))
    p <- p + geom_point(data=repGenes)
    p <- p + geom_text_repel(data=lab, aes(label=genes))
    p <- p + theme(legend.position = "none")
    p <- p + labs(x="log2 fold change vs HRD control", y="-log10   FDR")
    p <- p +  xlim(c(-7, 7))
    if (file.exists(paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, comp, ".pdf"))) {
      print(paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, comp, ".pdf"))
      p
    } else {
      print(paste0("Creating ",plotDir, "/", Type,  "_volcano_FDR", thresh, comp, ".pdf"))
      pdf(file = paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, comp, ".pdf"))
      print(p)
      dev.off()
    }
    
    
    ### 6. Identify epigenetic and RNAi regulators and assess expression: ###
    
    epiIDs <- c("ENSG00000130816", "ENSG00000119772", "ENSG00000088305",  
                "ENSG00000276043", "ENSG00000138336", "ENSG00000168769", "ENSG00000187605", 
                "ENSG00000100697", "ENSG00000092847", "ENSG00000123908", "ENSG00000126070", 
                "ENSG00000134698", "ENSG00000101945")
    
    epiSym <- c("DNMT1", "DNMT3A", "DNMT3B", "UHRF1", "TET1", "TET2", "TET3", 
                "Dicer", "AGO1", "AGO2", "AGO3", "AGO4", "SUV39H1")
    
    epiGenes <- allGenes[epiIDs,]
    rownames(epiGenes) <- epiSym
    
    epiThresh <- 0.1
    epiGenes$threshold <- as.factor((epiGenes$FDR < epiThresh & epiGenes$logFC < -(thresh))|(epiGenes$FDR <  epiThresh & epiGenes$logFC > thresh))
    epiSig <- subset(epiGenes, threshold == T)
    
    epiGenes$genes <- rownames(epiGenes)
    
    if (i==1) {
      epiGenesHrd <- list(epiGenes)
    } else {
      epiGenesHrd[[i]] <- epiGenes
    }

    p <- ggplot(data=epiGenes, aes(x=logFC, y=-log10(FDR),color=threshold))
    p <- p + geom_point(data=epiGenes)
    p <- p + geom_text_repel(data=epiGenes, aes(label=genes))
    p <- p + theme(legend.position = "none")
    p <- p + labs(x="log2 fold change vs HRD control", y="-log10   FDR")
    if (file.exists(paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, "_logFC", thresh, "_", comp,  "_epi.pdf"))) {
      print(paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, "_logFC", thresh, "_", comp, "_epi.pdf"))
      p
    } else {
      print(paste0("Creating ",plotDir, "/", Type,  "_volcano_FDR", thresh, "_logFC", thresh, "_", comp,  "_epi.pdf"))
      pdf(file = paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, "_logFC", thresh, "_", comp, "_epi.pdf"))
      print(p)
      dev.off()
    }
    
    saveRDS(repGenes, file=paste0(newRobjectDir, "/", Type, "_allDErep", thresh, "_logFC", thresh, "_", comp,  ".rds"))
    
    con <- c(rep(0, (ctlInd - 1) ), -1, rep(0, (ncol(design) - ctlInd)))
  }
}


### 7. Save all repeat DE, and all repeat DE that is significant and log2 FC > 1 for at least one comparison:

allRepsHrd <- append(allReps[1:2], allReps[4])
names(allRepsHrd) <- c(paste0("HRD_vs_", sTypes[1:2]), paste0("HRD_vs_", sTypes[4]))

sigRepsHrd <- append(sigReps[1:2], sigReps[4])
names(sigRepsHrd) <- c(paste0("HRD_vs_", sTypes[1:2]), paste0("HRD_vs_", sTypes[4]))

epiGenesHrd <- epiGenesHrd[-3]
names(epiGenesHrd) <- c(paste0("HRD_vs_", sTypes[1:2]), paste0("HRD_vs_", sTypes[4]))

# don't just want significant genes for each individual group, so find names of
# all the sig genes from all groups by taking the rownames and unlisting into one vector:
totalSigHrd <- unique(unlist(lapply(sigReps, function(x) {
  return(rownames(x))
})))
names(totalSigHrd) <- NULL


### 8. Perform differential expression comparing normalised CCNE controls to cancer samples ###

# assign CCNE as the new control:
ctl <- "CCNEAmp"

# assign control group:
ctlInd <- grep(ctl, colnames(design))
con <- c(rep(0, (ctlInd - 1) ), -1, rep(0, (ncol(design) - ctlInd)))

sTypes <- sTypes[order(sTypes)]

for (i in 1:length(sTypes)) {
  if (i!=ctlInd) {
    comp <- paste0(ctl, "_vs_", sTypes[i])
    
    # perform likelihood ratio test:
    con[i] <- 1
    lrt <- glmLRT(fit, contrast = con)
    
    # determine the top DE genes:
    topTags(lrt)
    
    if (file.exists(paste0(newRobjectDir, "CCNE_vs_", sTypes[i],   "_lrt.rds"))) {
      print(paste0(newRobjectDir, "CCNE_vs_", sTypes[i],   "_lrt.rds already exists, no need to create"))
    } else {
      print(paste0("Creating ", newRobjectDir, "CCNE_vs_", sTypes[i],   "_lrt.rds"))
      saveRDS(lrt, file = paste0(newRobjectDir, "CCNE_vs_", sTypes[i],   "_lrt.rds"))
    }
    
    
    ### 5. Calculate differential expression values with FDR threshold only ###
    
    # fetch summary of differentially expressed genes (those  with FDR =< 0.05:
    DEs <- summary(result <- decideTestsDGE((lrt)))
    
    # fetch all gene DE info, 
    thresh <- 0.1
    allGenes <- as.data.frame(topTags(lrt, n=Inf))
    repGenes <- allGenes[grep("ENS", rownames(allGenes), invert   = T),]
    sigGenes <- filter(repGenes, FDR<thresh)
    print(repGenes)
    
    # plot on volcano plot:
    repGenes$threshold <- as.factor(repGenes$FDR < thresh)
    sig <- subset(repGenes, threshold == T)
    
    
    lab <- sig
    lab$genes <- rownames(lab)
    
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
    
    p <- ggplot(data=repGenes, aes(x=logFC, y=-log10(FDR),  color=threshold))
    p <- p + geom_point(data=repGenes)
    p <- p + geom_text_repel(data=lab, aes(label=genes))
    p <- p + theme(legend.position = "none")
    p <- p + labs(x="log2 fold change vs CCNE control", y="-log10   FDR")
    p <- p +  xlim(c(-7, 7))
    if (file.exists(paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, comp, ".pdf"))) {
      print(paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, comp, ".pdf"))
      p
    } else {
      print(paste0("Creating ",plotDir, "/", Type,  "_volcano_FDR", thresh, comp, ".pdf"))
      pdf(file = paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, comp, ".pdf"))
      print(p)
      dev.off()
    }
    
    
    ### 6. Identify epigenetic and RNAi regulators and assess expression: ###
    
    epiIDs <- c("ENSG00000130816", "ENSG00000119772", "ENSG00000088305",  
                "ENSG00000276043", "ENSG00000138336", "ENSG00000168769", "ENSG00000187605", 
                "ENSG00000100697", "ENSG00000092847", "ENSG00000123908", "ENSG00000126070", 
                "ENSG00000134698", "ENSG00000101945")
    
    epiSym <- c("DNMT1", "DNMT3A", "DNMT3B", "UHRF1", "TET1", "TET2", "TET3", 
                "Dicer", "AGO1", "AGO2", "AGO3", "AGO4", "SUV39H1")
    
    epiGenes <- allGenes[epiIDs,]
    rownames(epiGenes) <- epiSym
    
    epiThresh <- 0.1
    epiGenes$threshold <- as.factor((epiGenes$FDR < epiThresh & epiGenes$logFC < -(thresh))|(epiGenes$FDR <  epiThresh & epiGenes$logFC > thresh))
    epiSig <- subset(epiGenes, threshold == T)
    
    epiGenes$genes <- rownames(epiGenes)
    
    if (i==1) {
      epiGenesCcne <- list(epiGenes)
    } else {
      epiGenesCcne[[i]] <- epiGenes
    }

    p <- ggplot(data=epiGenes, aes(x=logFC, y=-log10(FDR),color=threshold))
    p <- p + geom_point(data=epiGenes)
    p <- p + geom_text_repel(data=epiGenes, aes(label=genes))
    p <- p + theme(legend.position = "none")
    p <- p + labs(x="log2 fold change vs CCNE control", y="-log10   FDR")
    if (file.exists(paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, "_logFC", thresh, "_", comp,  "_epi.pdf"))) {
      print(paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, "_logFC", thresh, "_", comp, "_epi.pdf"))
      p
    } else {
      print(paste0("Creating ",plotDir, "/", Type,  "_volcano_FDR", thresh, "_logFC", thresh, "_", comp,  "_epi.pdf"))
      pdf(file = paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, "_logFC", thresh, "_", comp, "_epi.pdf"))
      print(p)
      dev.off()
    }
    
    saveRDS(repGenes, file=paste0(newRobjectDir, "/", Type, "_allDErep", thresh, "_logFC", thresh, "_", comp,  ".rds"))
    
    con <- c(rep(0, (ctlInd - 1) ), -1, rep(0, (ncol(design) - ctlInd)))
  }
}


### 9. Save all repeat DE, and all repeat DE that is significant and log2 FC > 1 for at least one comparison:

allRepsCcne <- append(allReps[1], allReps[3:4])
names(allRepsCcne) <- c(paste0("CCNEAmp_vs_", sTypes[1]), paste0("CCNEAmp_vs_", sTypes[3:4]))

sigRepsCcne <- append(sigReps[1], sigReps[3:4])
names(sigRepsCcne) <- c(paste0("CCNEAmp_vs_", sTypes[1]), paste0("CCNEAmp_vs_", sTypes[3:4]))

epiGenesCcne <- epiGenesCcne[-2]
names(epiGenesCcne) <- c(paste0("CCNEAmp_vs_", sTypes[1]), paste0("CCNEAmp_vs_", sTypes[3:4]))

# don't just want significant genes for each individual group, so find names of
# all the sig genes from all groups by taking the rownames and unlisting into one vector:
totalSigCcne <- unique(unlist(lapply(sigRepsCcne, function(x) {
  return(rownames(x))
})))
names(totalSigCcne) <- NULL

# append totalSigHrd onto totalSigCcne to include in allSig:
totalSig <- unique(append(totalSigHrd, totalSigCcne))

# fetch all the entries for the above significant repeat names and save:
allSigHrd <- lapply(allRepsHrd, function(x) {
  return(x[totalSig,])
})

allSigCcne <- lapply(allRepsCcne, function(x) {
  return(x[totalSig,])
})

allSig <- c(allSigHrd, list(allSigCcne$CCNEAmp_vs_bothDrivers), list(allSigCcne$CCNEAmp_vs_unknownDrivers))
names(allSig)[4:5] <- c(names(allSigCcne)[1], names(allSigCcne[3]))

if (!file.exists(paste0(newRobjectDir, "/", Type, "_DEsigReps.rds"))) {
	saveRDS(allSig, file=paste0(newRobjectDir, "/", Type, "_DEsigReps.rds"))
}

epiGenesAll <- c(epiGenesHrd, epiGenesCcne)
epiGenesAll <- epiGenesAll[-5]

if (!file.exists(paste0(newRobjectDir, "/", Type, "_DEepiGenes.rds"))) {
  saveRDS(allSig, file=paste0(newRobjectDir, "/", Type, "_DEepiGenes.rds"))
}

if (!file.exists(paste0(newRobjectDir, "/", Type, "_DEepiGenes.rds"))) {
	saveRDS(epiGenesAll, file=paste0(newRobjectDir, "/", Type, "_DEepiGenes.rds"))
}

if (!file.exists(paste0(newRobjectDir, "DEImg_", expName, ".RData"))) {
save.image(file = paste0(newRobjectDir, "DEImg_", expName, ".RData"))
}


