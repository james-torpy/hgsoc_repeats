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
sTypes <- c("FT", "all_primary", "acquired resistance", "ascites")
descrip <- "hgsoc_split_by_biopsy_time_vs_FT"

# define sample groups to compare:
sGroups <- list("FT", c("erPT", "mrPT", "prPT", "rfPT"), "arPT", "rcAF")
names(sGroups) <- sTypes

# define sample group to use as control:
ctl <- "FT"

# specify what combination of repeat genes (repeats), epigenetic modulators (epiMods),
# RNAi genes (RNAi) and protein-coding genes (pCoding) should contribute to the results:
resultTypes <- c("repeats", "epiMods")

# specify what FDR and log2 fold change thresholds to use:
FDRthresh <- "0.05"
FCthresh <- "1"

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
                     expName, "/")
newRobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
                        expName, "/allHGSOC/")
plotDir <- paste0(resultsDir, "/R/", expName,
                  "/plots/DEplots/allHGSOC/")

system(paste0("mkdir -p ", plotDir))
system(paste0("mkdir -p ", RobjectDir))
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

# change all types except for FT to 'HGSOC':
newNames <- lapply(strsplit(colnames(Counts), "_"), function(x) {
  if (x[3] == "FT") {
    print("Sample is FT, not changing sample name")
  } else {
    print(paste0("Replacing sample type ", x[3], " with 'HGSOC' in sample name"))
    x[3] <- "HGSOC"
  }
  new <- paste0(x[1], "_", x[2], "_", x[3])
  return(new)
})

colnames(Counts) <- unlist(newNames)

# change the order of columns of Counts to alphabetical order:
colOrder <- unlist(split(colnames(Counts), gsub("^.*_*_", "", colnames(Counts))))
Counts <- Counts[,colOrder]

splt <- unlist(
  lapply(
    split(
      colnames(Counts), gsub(
        "HGSOC.[0-9]", "HGSOC", gsub("^.*_*_", "", colnames(Counts))
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

saveRDS(fit, file=paste0(newRobjectDir, "/", Type, "DEfit.rds"))
save.image(paste0(newRobjectDir, "/", Type, "DEdone.rds"))

# determine which column has FT control:
ctlInd <- grep("FT", colnames(design))
con <- c(rep(0, (ctlInd - 1) ), -1, rep(0, (ncol(design) - ctlInd)))

sTypes <- sTypes[order(sTypes)]

for (i in 1:ncol(design)) {
  print(i)
  if (i!=ctlInd) {
    comp <- paste0("FT_vs_", sTypes[i])
    
    # perform likelihood ratio test:
    con[i] <- 1
    lrt <- glmLRT(fit, contrast = con)
    
    # determine the top DE genes:
    topTags(lrt)
    
    #if (file.exists(paste0(newRobjectDir, "FT_vs_", sTypes[i],   "_lrt.rds"))) {
    #  print(paste0(newRobjectDir, "FT_vs_", sTypes[i],   "_lrt.rds already exists, no need to create"))
    #} else {
    #  print(paste0("Creating ", newRobjectDir, "FT_vs_", sTypes[i],   "_lrt.rds"))
    #  saveRDS(lrt, file = paste0(newRobjectDir, "FT_vs_", sTypes[i],   "_lrt.rds"))
    #}
    
    
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
    # add positive and negative control genes GAPDH, CD47, b-actin, GUSB:
    posGenes <- rbind(allGenes["ENSG00000111640",],   allGenes["ENSG00000196776",])
    negGenes <- rbind(allGenes["ENSG00000075624",],   allGenes["ENSG00000169919",])
    rownames(posGenes) <- c("GAPDH", "CD47")
    rownames(negGenes) <- c("beta-actin", "GUSB")
    
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
    
    if (!(ctlInd==1)) {
      if (i==1) {
        allReps <- list(repGenes)
      } else {
        allReps[[(i-1)]] <- repGenes
      }
      
      if (i==1) {
        sigReps <- list(sig)
      } else {
        sigReps[[(i-1)]] <- sig
      }
    } else {
      if (i==2) {
        allReps <- list(repGenes)
      } else {
        allReps[[(i-1)]] <- repGenes
      }
      
      if (i==2) {
        sigReps <- list(sig)
      } else {
        sigReps[[(i-1)]] <- sig
      }
    }
    
    
    p <- ggplot(data=repGenes, aes(x=logFC, y=-log10(FDR),  color=threshold))
    p <- p + geom_point(data=repGenes)
    p <- p + geom_text_repel(data=lab, aes(label=genes))
    p <- p + theme(legend.position = "none")
    p <- p + labs(x="log2 fold change vs FT control", y="-log10   FDR")
    p <- p +  xlim(c(-4, 4))
    if (file.exists(paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, comp, ".pdf"))) {
      print(paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, comp, ".pdf"))
      p
    } else {
      print(paste0("Creating ",plotDir, "/", Type,  "_volcano_FDR", thresh, comp, ".pdf"))
      pdf(file = paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, comp, ".pdf"))
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
    # add positive and negative control genes GAPDH, CD47, b-actin, GUSB:
    posGenes <- rbind(allGenes["ENSG00000111640",],   allGenes["ENSG00000196776",],
                      allGenes["FAM",], allGenes["REP522",])
    negGenes <- rbind(allGenes["ENSG00000075624",],   allGenes["ENSG00000169919",])
    rownames(posGenes) <- c("GAPDH", "CD47")
    rownames(negGenes) <- c("beta-actin", "GUSB")
    
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
      print(paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, "_logFC", FCthresh, "_", comp, ".pdf"))
      p
    } else {
      print(paste0("Creating ",plotDir, "/", Type,  "_volcano_FDR", thresh, "_logFC", FCthresh, "_", comp,  ".pdf"))
      pdf(file = paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, "_logFC", FCthresh, "_", comp, ".pdf"))
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
    
    epiGenes <- allGenes[epiIDs,]
    rownames(epiGenes) <- epiSym
    
    epiThresh <- 0.1
    epiGenes$threshold <- as.factor((epiGenes$FDR < epiThresh & epiGenes$logFC < -(FCthresh))|(epiGenes$FDR <  epiThresh & epiGenes$logFC > FCthresh))
    epiSig <- subset(epiGenes, threshold == T)
    
    epiGenes$genes <- rownames(epiGenes)
    
    
    p <- ggplot(data=epiGenes, aes(x=logFC, y=-log10(FDR),color=threshold))
    p <- p + geom_point(data=epiGenes)
    p <- p + geom_text_repel(data=epiGenes, aes(label=genes))
    p <- p + theme(legend.position = "none")
    p <- p + labs(x="log2 fold change vs FT control", y="-log10   FDR")
    p <- p +  xlim(c(-2, 2))
    #if (file.exists(paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, "_logFC", FCthresh, "_", comp,  "_epi.pdf"))) {
     # print(paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, "_logFC", FCthresh, "_", comp, "_epi.pdf"))
      #p
    #} else {
      print(paste0("Creating ",plotDir, "/", Type,  "_volcano_FDR", thresh, "_logFC", FCthresh, "_", comp,  "_epi.pdf"))
      pdf(file = paste0(plotDir, "/", Type,  "_volcano_FDR", thresh, "_logFC", FCthresh, "_", comp, "_epi.pdf"))
      print(p)
      dev.off()
    #}
    
    saveRDS(repGenes, file=paste0(newRobjectDir, "/", Type, "_allDErep", thresh, "_logFC", FCthresh, "_", comp,  ".rds"))
    
    con <- c(rep(0, (ctlInd - 1) ), -1, rep(0, (ncol(design) - ctlInd)))
  }
}

### 8. Save all repeat DE, and all repeat DE that is significant and log2 FC > 1 for at least one comparison:

names(allReps) <- paste0("FT_vs_", sTypes[-1])
names(sigReps) <- paste0("FT_vs_", sTypes[-1])
names(epiGenes) <- paste0("FT_vs_", sTypes[-1])

saveRDS(allReps, file=paste0(newRobjectDir, "/", Type, "_DEreps.rds"))
saveRDS(epiGenes, file=paste0(newRobjectDir, "/", Type, "_DEepiGenes.rds"))


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
