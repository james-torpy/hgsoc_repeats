### normaliseCounts.R ###

# This script takes a list of dfs of different classes and types of repeat counts and
# normalises using RUVseq, EdgeR and DEseq, then compares methods

### 0. Define variables/paths ###

# load packages needed:
library(tibble)
library(dplyr)
library(RColorBrewer)
library(RUVSeq)
library(ggplot2)
library(ggrepel)
library(preprocessCore)
library(reshape2)
library(Rmisc)

# define starting variables:
project <- "hgsoc_repeats"
expName <- "exp5"
STypes <- c("arPT5", "rcAF6")
Type <- "c"

# define directories:
#homeDir <- "/Users/jamestorpy/clusterHome/"
homeDir <- "/share/ScratchGeneral/jamtor/"
#homeDir <- "/Users/jamestorpy/Documents/Garvan/phd"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/", expName, "/")
plotDir <- paste0(resultsDir, "/R/", expName, "/plots/DEplots/arPT5_vs_rcAF6/")

system(paste0("mkdir -p ", plotDir))


### 1. Load in all counts ###

cCounts <- readRDS(paste0(RobjectDir, "/", Type, "_allcounts.htseq.rds"))
cCounts$gene_id <- rownames(cCounts)
gcCounts <- readRDS(paste0(RobjectDir, "/gc_allcounts.htseq.rds"))

# append GCcountsDF to each GRanges object of custom1Counts:
Counts <- rbind(cCounts, gcCounts)

# remove additional info rows:
Counts <- Counts[grep("__", Counts$gene_id, invert=T),]

# make rownames gene_id, get rid of latter column and change
# storage mode from factor to integer:
Rnames <- Counts$gene_id
Counts <- sapply(subset(Counts, select=-gene_id), unlist)
storage.mode(Counts) <- "integer"
Counts <- as.data.frame(Counts)
rownames(Counts) <- Rnames

# save Counts:
#saveRDS(Counts, file=paste0(RobjectDir, "/", Type, "_counts.RData"))
saveRDS(Counts, file=paste0(RobjectDir, "/", Type, "_counts.RData"))


### 3. Perform pre-normalisation PCA and RLE plots ###

# eliminate lowly expressed genes (rows where there are less than 3 counts where df > 4):
print(paste0("No. rows before filtering is: ", nrow(Counts)))
Counts <- Counts %>%
  rownames_to_column('gene_id') %>%
  dplyr::filter(rowSums(Counts > 5) >= 1) %>%
  column_to_rownames('gene_id')
print(paste0("No. rows after  filtering: ", nrow(Counts)))

# create pre-normalised PCA plot from counts and plot:
Counts <- apply(Counts, 2, unlist)
if (ncol(Counts) > nrow(Counts)) {
  pca <- prcomp(Counts)
} else {
  pca <- princomp(Counts)	  
}

splt <- unlist(lapply(split(colnames(Counts), gsub("^.*_*_", "", colnames(Counts))), length))
for (i in 1:length(splt)) {
  if (i==1) {
    typeF <- c(rep(names(splt)[i], splt[i]))
  } else {
    typeF <- c(typeF, rep(names(splt)[i], splt[i]))
  }
}
levels(typeF) <- c("arPT5", "rcAF6")

# convert matrix into SeqExpressionSet:
set <- newSeqExpressionSet(Counts, phenoData = data.frame(typeF, row.names=colnames(Counts)))


### 4. perform normalisation on counts using RUVseq:

# perform between lane full normalisation:
nSet <- betweenLaneNormalization(set, which="full")


### 5. Perform differential expression comparing normalised FT controls to cancer samples ###

genes <- rownames(Counts)
# normalise using upper-quartile normalisation (http://inaykmittal.blogspot.com.au/2013/10/	pkmrpkm-normalization-caveat-and-upper.html)
# design matrix specifying all samples as the thing to be compared to:
design <- model.matrix(~0+typeF, data=pData(nSet))
# convert set into a DGElist format which specifies library size and normalisation factors:
y <- DGEList(counts=counts(nSet), group=typeF)

lSizes <- y$samples$lib.size
saveRDS(lSizes, file=paste0(RobjectDir, "/", Type, "_libsizes.RData"))


# calculate CPMs:
CPM <- as.data.frame(t(t(Counts)/lSizes)*1000000)

# subset to include only control genes and repeats:
ctls <- c("ENSG00000111640", "ENSG00000012048", "ENSG00000156414", 
  "ENSG00000187605", "ENSG00000077800", "ENSG00000173809", 
  "ENSG00000105173")
names(ctls) <- c("GAPDH", "BRCA1", "TDRD9", "TET3", "FKBP6", "TDRD12", 
  "CCNE1")
#ctls <- c("ENSG00000111640", "ENSG00000204574",
#	"ENSG00000012048", "ENSG00000173809")
#names(ctls) <- c("GAPDH", "CCNE1", "BRCA1", "TDRD12")
#ctls <- c("ENSG00000111640",
#  "ENSG00000012048", "ENSG00000173809")
#names(ctls) <- c("GAPDH", "BRCA1", "TDRD12")
#CPM$id <- rownames(CPM)
ctlCPM <- CPM[rownames(CPM) %in% ctls,]
rnam <- names(ctls)[match(rownames(ctlCPM), ctls)]
repeatCPM <- CPM[grep("ENSG", rownames(CPM), invert=T),]
allCPM <- rbind(ctlCPM, repeatCPM)
allnam <- c(rnam, rownames(allCPM)[(length(ctls)+1):nrow(allCPM)])

# calculate total repeat count size:
#rSizes <- apply(Counts, 2, sum)
# calculate CPMRs
#CPMR <- as.data.frame(t(t(Counts)/rSizes)*1000000)

# prepare each df for plotting:
df <- allCPM
sampleNames <- colnames(df)
df <- sapply(df, as.numeric)
rownames(df) <- allnam

for (i in 1:4) {
  if (i==1) {
      sp <- c(rep(i, 5))
  } else {
    sp <- c(sp, rep(i, 5))
  }
}

dfs <- as.data.frame(split(df[,1], sp))
i=1
dfs <- lapply(dfs, function(x) {
  result <- cbind(x, split(df[,2], sp)[[i]])
  colnames(result) <- colnames(allCPM)
  i <<- i+1
return(result)
})

pDFs <- lapply(dfs, function(x) {
  melt(x, varnames = c("id", "sample"), value.name = "CPM")
})

# sort scpCounts data frames according to repeat_id:
sort_rID <- function(x) {
  if (colnames(x)[1] == "repeat_id") {
    return(x[with(x, order(repeat_id)),])
  } else {
    return(x)
  }
}
pDFs <- lapply(pDFs, sort_rID)


# relevel factors to put negative ctls first:
pDFs <- lapply(pDFs, function(x) {
  x$id <- factor(x$id,
  levels(x$id)[c(2:3, 1, 4:nrow(x))])
  return(x)
})

# add 1 to all values so they are loggable:
add2zero <- function(x) {
  if (ncol(x) == 3) {
    x[,3] <- x[,3]+1
    return(x)
  } else {
    x[,2] <- x[,2]+1
    return(x)
  }
}
pDF4Logs <- lapply(pDFs, add2zero)


# rename samples:
pDFs <- lapply(pDFs, function(x) {
  x$sample <- gsub(
    "rcAF", "recurrent_chemo_resistant_tumour",gsub(
    "arPT", "acquired_resistant_primary_tumour", x$sample
    )
  )
  return(x)
})



# create stats summary of the data:
#statDF4Log <- summarySE(pDF, measurevar = "CPMR", groupvar = c("repeat_id", "sample"))
#p <- ggplot(pDF, aes(x=sample, y=CPM))
#p <- p + geom_bar(stat = "identity", aes(fill = id), position = "dodge")
#p <- p + scale_y_log10()
#p <- p + xlab("samples") + ylab("log10_CPM")
#pdf(file = paste0(plotDir, "/ctlCPM_barplot.pdf"), width = 10, height=10)
##print(p)
#dev.off()

#saveRDS(pDF, file=paste0(RobjectDir, "/", Type, "_arPT5_vs_rcAF6_log10_plotDF.rds"))


for (i in 1:length(pDFs)) {
  p <- ggplot(pDFs[[i]], aes(x=sample, y=CPM, group=id, colour=id))
  p <- p + geom_line()
  p <- p + theme(axis.text.x = element_text(angle = 90), legend.position = "bottom")
  p <- p + xlab("samples") + ylab("CPM")
  pdf(file = paste0(plotDir, "/", Type, "_CPMlineplot", i, ".pdf"), width = 10, height=10)
  print(p)
  dev.off()
}



# rename samples:
pDF4Logs <- lapply(pDF4Logs, function(x) {
  x$sample <- gsub(
    "rcAF", "recurrent_chemo_resistant_tumour",gsub(
    "arPT", "acquired_resistant_primary_tumour", x$sample
    )
  )
  return(x)
})


# relevel factors to put negative ctls first:
pDF4Logs <- lapply(pDF4Logs, function(x) {
  x$id <- factor(x$id,
  levels(x$id)[c(2:3, 1, 4:nrow(x))])
  return(x)
})

# create stats summary of the data:
#statDF4Log <- summarySE(pDF, measurevar = "CPMR", groupvar = c("repeat_id", "sample"))
#p <- ggplot(pDF, aes(x=sample, y=CPM))
#p <- p + geom_bar(stat = "identity", aes(fill = id), position = "dodge")
#p <- p + scale_y_log10()
#p <- p + xlab("samples") + ylab("log10_CPM")
#pdf(file = paste0(plotDir, "/ctlCPM_barplot.pdf"), width = 10, height=10)
##print(p)
#dev.off()

#saveRDS(pDF, file=paste0(RobjectDir, "/", Type, "_arPT5_vs_rcAF6_log10_plotDF.rds"))


for (i in 1:length(pDF4Logs)) {
  p <- ggplot(pDF4Logs[[i]], aes(x=sample, y=CPM, group=id, colour=id))
  p <- p + geom_line()
  p <- p + theme(axis.text.x = element_text(angle = 90), legend.position = "bottom")
  p <- p + scale_y_log10()
  p <- p + xlab("samples") + ylab("CPM")
  pdf(file = paste0(plotDir, "/", Type, "_CPM_log10_lineplot", i, ".pdf"), width = 10, height=10)
  print(p)
  dev.off()
}





######


ctlDF <- df[1:length(ctls),]

pCtlDF <- melt(ctlDF, varnames = c("id", "sample"), value.name = "CPM")

# sort scpCounts data frames according to repeat_id:
sort_rID <- function(x) {
  if (colnames(x)[1] == "repeat_id") {
    return(x[with(x, order(repeat_id)),])
  } else {
    return(x)
  }
}
pCtlDF <- sort_rID(pCtlDF)


# relevel factors to put negative ctls first:
pCtlDF$id <- factor(pCtlDF$id,
levels(pCtlDF$id)[c(4, 1, 5, 7, 2, 6, 3)])


# add 1 to all values so they are loggable:
add2zero <- function(x) {
  if (ncol(x) == 3) {
    x[,3] <- x[,3]+1
    return(x)
  } else {
    x[,2] <- x[,2]+1
    return(x)
  }
}
pCtlDF4Log <- add2zero(pCtlDF)

# rename samples:

pCtlDF4Log$sample <- gsub(
  "rcAF", "recurrent_chemo_resistant_tumour",gsub(
  "arPT", "acquired_resistant_primary_tumour", pCtlDF4Log$sample
  )
)


p <- ggplot(pCtlDF4Log, aes(x=sample, y=CPM))
p <- p + geom_bar(stat = "identity", aes(fill = id), position = "dodge")
p <- p + scale_y_log10()
p <- p + xlab("samples") + ylab("log10_CPM")
pdf(file = paste0(plotDir, "/ctlCPM_barplot4.pdf"), width = 10, height=10)
print(p)
dev.off()


######




plot_itdf <- function(x, logNo="nope") {
  if (ncol(x) == 7) {
    p <- ggplot(x, aes(x=sample, y=x[,3], group=repeat_id, colour=repeat_id))
    p <- p + geom_errorbar(aes(ymin=eval(parse(text="CPMR"))-se, ymax=eval(parse(text="CPMR"))+se), width=0.1)
    p <- p + geom_line()
    p <- p + theme(axis.text.x = element_text(angle = 90), legend.position = "bottom")
    p <- p + ylab("CPMR")
    if (logNo=="log10") {
      p <- p + scale_y_log10()
      return(p)
    } else {
      return(p)
    }
  } else {
    p <- ggplot(x, aes(x=sample, y=x[,4], group = 1))
    #p <- p + geom_errorbar(aes(ymin=eval(parse(text="CPMR"))-se-se, ymax=eval(parse(text="CPMR"))-se+se), width=0.1)
    p <- p + geom_line()
    p <- p + theme(axis.text.x = element_text(angle = 90), legend.position = "bottom")
    p <- p + ylab("CPMR")
    if (logNo=="log10") {
      p <- p + scale_y_log10()
      return(p)
    } else {
      return(p)
    }
  }
}

pDFPlot <- plot_itdf(statDF)
pDFPlotLog10 <- plot_itdf(statDF4Log, logNo="log10")

if (file.exists(paste0(plotDir, "/", Type, "_", names(countDFs)[j], "_", "CPMR", ".pdf"))) {
  print(paste0(plotDir, "/", Type, "_", names(countDFs)[j], "_", "CPMR", ".pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", plotDir, "/", Type, "_", names(countDFs)[j], "_", "CPMR", ".pdf"))
  pdf(file = paste0(plotDir, "/", Type, "_", "CPMR", ".pdf"), width = 10, height=10)
  print(pDFPlot)
  dev.off()
}

if (file.exists(paste0(plotDir, "/", Type, "_", names(countDFs)[j], "_", "CPMR", "Log10.pdf"))) {
  print(paste0(plotDir, "/", Type, "_", names(countDFs)[j], "_", "CPMR", "Log10.pdf already exists, no need to create"))
} else {
  print(paste0("Creating ", plotDir, "/", Type, "_", names(countDFs)[j], "_", "CPMR", "Log10.pdf"))
  pdf(file = paste0(plotDir, "/", Type, "_", names(countDFs)[j], "_", "CPMR", "Log10.pdf"), width = 10, height=10)
  print(pDFPlotLog10)
  dev.off()
}
n=n+1
}
j=j+1
}

names(CPMRcount) <- names(countDFs)
saveRDS(CPMRcount, file=paste0(RobjectDir, "/", expName, "/", Type, "CPMR.rds"))

names(CPMcount) <- names(countDFs)
saveRDS(CPMcount, file=paste0(RobjectDir, "/", expName, "/", Type, "CPM.rds"))
