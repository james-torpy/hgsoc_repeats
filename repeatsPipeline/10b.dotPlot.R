# load packages needed:
library(ggplot2)
library(reshape2)
library(dplyr)
library(Rmisc)
library(RColorBrewer)

# define starting variables:
project <- "hgsoc_repeats"
expName <- "exp2"

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
#homeDir <- "/Users/jamestorpy/Documents/Garvan/phd/"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/")
inDir <- paste0(resultsDir, "/R/", expName)
plotDir <- paste0(inDir, "/plots/candidatePlots")

system(paste0("mkdir -p ", plotDir))


### 1. Load in inputs ###
custom3CPMR <- readRDS(file=paste0(RobjectDir, "/", expName, "/custom3CPMR.rds"))
cCPMR <- readRDS(file=paste0(RobjectDir, "/", expName, "/cCPMR.rds"))

# remove any samples present in c but not custom3, and vice versa:
#cCPMR <- lapply(cCPMR, function(x) {
#  return(x[-grep("prPT8", x$sample),])
#})

# grep specific repeats from custom3CMPR:
spec1 <- c("HSATII", "ACRO1", "L1MD3")
for (i in 1:length(spec1)) {
  rep <- lapply(custom3CPMR, function(x) {
    return(x[grep(spec1[i], x$repeat_id),])
  })
  
  # remove empty elements from list:
  rep <- as.data.frame(rep[unlist(lapply(rep, function(x) {
    return(nrow(x)>0)
  }))
  ])
  
  # remove repeat type ids from colnames:
  colnames(rep) <- gsub("^.*\\.", "", colnames(rep))
  
  # bind result dfs together:
  if (i==1) {
    rp1 <- rep
  } else {
    rp1 <- rbind(rp1, rep)
  }
}

# grep specific repeats from cCMPR:
spec2 <- c("RC/Helitron", "RC\\?/Helitron\\?")
for (i in 1:length(spec2)) {
  rep <- lapply(cCPMR, function(x) {
    return(x[grep(spec2[i], x$repeat_id),])
  })
  
  # remove empty elements from list:
  rep <- as.data.frame(rep[unlist(lapply(rep, function(x) {
    return(nrow(x)>0)
  }))
  ])
  
  # remove repeat type ids from colnames:
  colnames(rep) <- gsub("^.*\\.", "", colnames(rep))
  
  # bind result dfs together:
  if (i==1) {
    rp2 <- rep
  } else {
    rp2 <- rbind(rp2, rep)
  }
}      
         
# bind rp1 and 2 together:
rp <- rbind(rp1, rp2)
# create sample type column:
rp$sType <- gsub("AOCS_[0-9][0-9][0-9]_", "", rp$sample)

# convert rp$sType to factor:
#rp$sType <- factor(rp$sType, levels = c("FT", "erPT", "mrPT", "arPT", "prPT", "rfPT", "lrcT", "msST", "pAF", "rcAF"))
rp$sType <- factor(rp$sType, levels = c("FT", "prPT"))

# add 1 to each CPMR value so they are loggable:
rp4Log <- rp
rp4Log$CPMR <- rp4Log$CPMR + 1

# perform Kruskal-Wallis Test (non-parametric non-normal ANOVA):
kwTest <- c("HSATII", "ACRO1", "L1MD3", "RC/Helitron", "RC?/Helitron?")
for (i in 1:length(kwTest)) {
  df <- rp[rp$repeat_id==kwTest[i],]
  splt <- split(df, df$sType)
  dati <- lapply(splt, function(x) { return(x[,3]) })
  kw <- kruskal.test(dati)
  if (i==1) {
    kwRes <- list(kw)
    p1 <- c(kw$p.value)
  } else {
    kwRes[[i]] <- kw
    p1[i] <- kw$p.value
  }
}
names(kwRes) <- kwTest
names(p1) <- kwTest

for (i in 1:length(p1)) {
  if (i==1) {
    if (p1[i]>0.05) {
      pVals <- c("p > 0.05")
    } else if (p1[i]==0.05) {
      pVals <- c("p = 0.05")
    } else {
      pVals <- c(paste0("p = ", round(p1[i], 7)))
    }
  } else {
    if (p1[i]>0.05) {
      pVals[i] <- "p > 0.05"
    } else if (p1[i]==0.05) {
      pVals[i] <- "p = 0.05"
    } else {
      pVals[i] <- paste0("p = ", round(p1[i], 7))
    }
  }
}


# Combine RC/Helitron and RC?/Helitron? values together:
rpC <- rp
rpC$repeat_id <- gsub("\\?", "", rpC$repeat_id)


rpC4Log <- rp4Log
rpC4Log$repeat_id <- gsub("\\?", "", rpC4Log$repeat_id)

dfC <- rpC[rpC$repeat_id=="RC/Helitron",]
spltC <- split(dfC, dfC$sType)
datiC <- lapply(spltC, function(x) { return(x[,3]) })
kwResC <- kruskal.test(datiC)
p2 <- kwResC$p.value

cols <- brewer.pal(10, "Paired")
# create scatter plot with 2 Helitron values:
p <- ggplot(rp4Log, aes(x=repeat_id, y=CPMR, group=sType, colour=sType))
p <- p + geom_point(aes(fill=sType), position = position_dodge(0.5))
p <- p + scale_y_log10()
p <- p + annotate("text", x=1, y=90000, label=pVals[1])
p <- p + annotate("text", x=2, y=900000, label=pVals[2])
p <- p + annotate("text", x=3, y=900000, label=pVals[3])
p <- p + annotate("text", x=4, y=9000, label=pVals[4])
p <- p + annotate("text", x=5, y=9000, label=pVals[5])
p <- p + scale_colour_manual(values = cols)
p
if (file.exists(paste0(plotDir, "/", "candidate_log10CPMRscatter.pdf"))) {
  print(paste0(plotDir, "/", "candidate_log10CPMRscatter.pdf"))
} else {
  print(paste0("Creating ", plotDir, "/", "candidate_log10CPMRscatter.pdf"))
  pdf(file = paste0(plotDir, "/", "candidate_log10CPMRscatter.pdf"), width = 10, height=10)
  print(p)
  dev.off()
}

# create boxplot with 2 Helitron values:
p <- ggplot(rp4Log, aes(x=repeat_id, y=CPMR, fill=sType))
p <- p + geom_boxplot()
p <- p + scale_y_log10()
p <- p + annotate("text", x=1, y=90000, label=pVals[1])
p <- p + annotate("text", x=2, y=900000, label=pVals[2])
p <- p + annotate("text", x=3, y=900000, label=pVals[3])
p <- p + annotate("text", x=4, y=9000, label=pVals[4])
p <- p + annotate("text", x=5, y=9000, label=pVals[5])
p
if (file.exists(paste0(plotDir, "/", "candidate_log10CPMRboxplot.pdf"))) {
  print(paste0(plotDir, "/", "candidate_log10CPMRboxplot.pdf"))
} else {
  print(paste0("Creating ", plotDir, "/", "candidate_log10CPMRboxplot.pdf"))
  pdf(file = paste0(plotDir, "/", "candidate_log10CPMRboxplot.pdf"), width = 10, height=10)
  print(p)
  dev.off()
}

# create scatter plot with 1 Helitron value:
p <- ggplot(rpC4Log, aes(x=repeat_id, y=CPMR, group=sType, colour=sType))
p <- p + geom_point(aes(fill=sType), position = position_dodge(0.5))
p <- p + scale_y_log10()
p <- p + annotate("text", x=1, y=900000, label=pVals[1])
p <- p + annotate("text", x=2, y=900000, label=pVals[2])
p <- p + annotate("text", x=3, y=900000, label=pVals[3])
p <- p + annotate("text", x=4, y=9000, label=paste0("p = ", round(p2, 7)))
p <- p + scale_colour_manual(values = cols)
p
if (file.exists(paste0(plotDir, "/", "candidate_log10CPMRscatter_helitron_combined.pdf"))) {
  print(paste0(plotDir, "/", "candidate_log10CPMRscatter_helitron_combined.pdf"))
} else {
  print(paste0("Creating ", plotDir, "/", "candidate_log10CPMRscatter_helitron_combined.pdf"))
  pdf(file = paste0(plotDir, "/", "candidate_log10CPMRscatter_helitron_combined.pdf"), width = 10, height=10)
  print(p)
  dev.off()
}

# create boxplot with 1 Helitron value:
p <- ggplot(rpC4Log, aes(x=repeat_id, y=CPMR, fill=sType))
p <- p + geom_boxplot()
p <- p + scale_y_log10()
p <- p + annotate("text", x=1, y=900000, label=pVals[1])
p <- p + annotate("text", x=2, y=900000, label=pVals[2])
p <- p + annotate("text", x=3, y=900000, label=pVals[3])
p <- p + annotate("text", x=4, y=9000, label=paste0("p = ", round(p2, 7)))
p
if (file.exists(paste0(plotDir, "/", "candidate_log10CPMRboxplot_helitron_combined.pdf"))) {
  print(paste0(plotDir, "/", "candidate_log10CPMRboxplot_helitron_combined.pdf"))
} else {
  print(paste0("Creating ", plotDir, "/", "candidate_log10CPMRboxplot_helitron_combined.pdf"))
  pdf(file = paste0(plotDir, "/", "candidate_log10CPMRboxplot_helitron_combined.pdf"), width = 10, height=10)
  print(p)
  dev.off()
}
