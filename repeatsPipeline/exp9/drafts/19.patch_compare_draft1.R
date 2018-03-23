### 19.patch_compare.R ###

# This script takes a CCNEamp vs HRD comparison and compares it with DE genes from 
# Patch et. al. 2015:

# Run on cluster with:
#briR
#qsub -N patchcomp -b y -wd \
#/share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/logs/exp9/DE \
#-j y -R y -pe smp 2 -V 
#"Rscript /share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/exp9/19.patch_compare.R"

### 0. Define variables/paths ###

# load packages needed:
library(tibble)
library(dplyr)
library(org.Hs.eg.db)
library(pheatmap)

# define starting variables:
project <- "hgsoc_repeats"
expName <- "exp9"

Type <- "all"
descrip <- "htseq_EdgeR_sig_vs_patch"

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project)
rawDir <- paste0(projectDir, "/RNA-seq/raw_files/fullsamples/bowtell_primary/")
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
                     expName, "/DEcompare/htseq_EDA_EdgeR_HGSOC_CCNEamp_unknown_vs_HRD/")
newRobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
                        expName, "/DEcompare/", descrip, "/")
plotDir <- paste0(resultsDir, "/R/", expName,
                  "/plots/DEplots/DEcompare/", descrip, "/")

system(paste0("mkdir -p ", newRobjectDir))
system(paste0("mkdir -p ", plotDir))


### 1. Load in data ###

# load in my DE data:
sig <- readRDS(paste0(RobjectDir, "/CCNEamp_vs_HRD_sig.rds"))

# load in Patch data:
patch <- read.csv(file=paste0(rawDir, "/CCNEamp_vs_HRD_Patch_DE.csv"), header=T, fill=T)
colnames(patch)[1] <- "symbol"


### 2. Create volcano plot comparing DE genes common to my DE and Patch DE analysis ###

patch_comp <- dplyr::select(merge(sig, patch, by="symbol"), symbol, logFC.x, FDR, logFC.y, Adjusted.P.value)
colnames(patch_comp) <- c("symbol", "my_logFC", "my_FDR", "Patch_logFC", "Patch_adj_Pvalue")

temp1 <- patch_comp[,c(1, 4, 5)]
colnames(temp1) <- c("symbol", "logFC", "adjPvalue_FDR")
temp1$type <- "patch"
temp1$threshold <- FALSE
temp1$threshold[temp1$adjPvalue_FDR < 0.005] <- TRUE
temp1$typethresh <- paste0(temp1$type, temp1$threshold)


temp2 <- patch_comp[,c(1, 2, 3)]
colnames(temp2) <- c("symbol", "logFC", "adjPvalue_FDR")
temp2$type <- "sig"
temp2$threshold <- FALSE
temp2$threshold[temp2$adjPvalue_FDR < 0.005] <- TRUE
temp2$typethresh <- paste0(temp2$type, temp2$threshold)

patch_comp_sep <- rbind(temp1, temp2)

lab3 <- patch_comp_sep[patch_comp_sep$threshold == TRUE,]

# plot on volcano plot:
p <- ggplot(data=patch_comp_sep, aes( x=logFC, y=-log10(adjPvalue_FDR), color=type) )
p <- p + geom_point(data=patch_comp_sep)
p <- p + geom_text_repel(data=lab3, aes(label=symbol))
p <- p + theme(legend.position =  "none")
p <- p + labs(x="log2 fold change   vs FT control", y="-log10   FDR")
#p <- p +  ylim(c(0, 45))
if (file.exists(paste0(plotDir, "/",  Type, "_CCNEamp_vs_HRD_common_sig_and_patch_volcano.pdf"))) {
  print(paste0(plotDir, "/",  Type, "_CCNEamp_vs_HRD_common_sig_and_patch_volcano.pdf already exists"))
  p
} else {
  print(paste0("Creating  ", plotDir, "/",  Type, "_CCNEamp_vs_HRD_common_sig_and_patch_volcano.pdf"))
  pdf(file = paste0(plotDir, "/",  Type, "_CCNEamp_vs_HRD_common_sig_and_patch_volcano.pdf"))
  print(p)
  dev.off()
}


### 3. Create volcano plot with all significant DE genes from my analysis and Patch 2015 analysis:

temp1 <- patch[,c(1, 2, 6)]
colnames(temp1) <- c("symbol", "logFC", "adjPvalue_FDR")
temp1$type <- "patch"
temp1$threshold <- FALSE
temp1$threshold[temp1$adjPvalue_FDR < 0.005] <- TRUE
temp1$typethresh <- paste0(temp1$type, temp1$threshold)

temp2 <- sig[,c(8, 1, 5)]
colnames(temp2) <- c("symbol", "logFC", "adjPvalue_FDR")
temp2$type <- "sig"
temp2$threshold <- FALSE
temp2$threshold[temp2$adjPvalue_FDR < 10e-20] <- TRUE
temp2$typethresh <- paste0(temp2$type, temp2$threshold)

patch_and_sig <- rbind(temp1, temp2)
lab2 <- patch_and_sig[patch_and_sig$threshold == TRUE,]

# plot on volcano plot:
p <- ggplot(data=patch_and_sig, aes( x=logFC, y=-log10(adjPvalue_FDR), color=typethresh) )
p <- p + geom_point(data=patch_and_sig)
p <- p + geom_text_repel(data=lab2, aes(label=symbol))
p <- p + theme(legend.position =  "none")
p <- p + labs(x="log2 fold change   vs FT control", y="-log10   FDR")
p <- p +  ylim(c(0, 45))
if (file.exists(paste0(plotDir, "/",  Type, "_CCNEamp vs HRD_sig_vs_patch.pdf"))) {
  print(paste0(plotDir, "/",  Type, "_CCNEamp vs HRD_sig_vs_patch.pdf already exists"))
  p
} else {
  print(paste0("Creating  ", plotDir, "/",  Type, "_CCNEamp vs HRD_sig_vs_patch.pdf"))
  pdf(file = paste0(plotDir, "/",  Type, "_CCNEamp vs HRD_sig_vs_patch.pdf"))
  print(p)
  dev.off()
}


### 4. Create heatmap with top 20 significant DE genes from my analysis and Patch 2015 analysis:

# take top 20 DE genes for patch and sig and save for heatmap:
vec1 <- as.character(head(temp1$symbol[order(temp1$adjPvalue_FDR)], 20))
vec2 <- head(na.omit(temp2$symbol[order(temp2$adjPvalue_FDR)]), 20)
vec_total <- c(vec1, vec2)

top_20s <- patch_comp[patch_comp$symbol %in% vec_total,]

# split into logFC and FDR dfs:
rownames(top_20s) <- top_20s$symbol

FC <- dplyr::select(top_20s, my_logFC, Patch_logFC)

FDR <- dplyr::select(top_20s, my_FDR, Patch_adj_Pvalue)
  
# create logFC and FDR heatmaps:
# resize margins for plots:
par(mar=c(4,4,4,4))
  
library(grid)
  
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}
  
## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                    ns=asNamespace("pheatmap"))

pdf(file = paste0(plotDir, "/FC_common_to_patch.pdf"))
pheatmap(FC, color = colorRampPalette(c("white", "firebrick3"))(50),
         display_numbers = as.matrix(ifelse(FDR < 0.1, "*", "")), fontsize = 8)
dev.off()