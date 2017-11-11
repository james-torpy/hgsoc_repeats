library(GenomicRanges)
library(ShortRead)
library(rtracklayer)
library("BSgenome.Hsapiens.UCSC.hg38")
library(ggplot2)

# define starting variables:
project <- "hgsoc_repeats"
methodName <- "method1"
annot <- "t"

# define directories:
#homeDir <- "/Users/jamestorpy/clusterHome/"
homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project)
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/")
bamDir <- paste0(projectDir, "/RNA-seq/results/star/", methodName)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
plotDir <- paste0(resultsDir, "/R/plots")


### 1. Create custom 3 annotation and save ###

Genes <- readRDS(file = paste0(RobjectDir, "/", annot, "_RepeatGenes.rds"))

# split Genes into custom3 annotation:
# key: RC/Helitron + RC?/Helitron?, DNA/TcMar-Tc2, all satellites, LINE/L1, SINE/Alu
Genes <- list(c(Genes[[43]], Genes[[44]]), Genes[[21]], c(Genes[[1]], Genes[[48]], Genes[[49]], Genes[[50]], Genes[[51]]), Genes[[27]], Genes[[57]])
names(Genes) <- c("RC/Helitron_all", "DNA/TcMar-Tc2", "Satellites_all", "LINE1/L1", "SINE/Alu")

# split L1s and Alus:
L1Simp <- gsub("[0-9]$", "", substr(names(Genes[[4]]), 1, 4))
L1Names <- split(names(Genes[[4]]), L1Simp)

for (i in 1:length(L1Names)) {
  Genes[[(5+i)]] <- Genes[[4]][L1Names[[i]]]
  print(i)
}
names(Genes)[6:(5+length(L1Names))] <- names(L1Names)

AluSimp <- substr(names(Genes[[5]]), 1, 4)
AluNames <- split(names(Genes[[5]]), AluSimp)

for (i in 1:length(AluNames)) {
  Genes[[(18+i)]] <- Genes[[5]][AluNames[[i]]]
  print(i)
}
names(Genes)[19:(18+length(AluNames))] <- names(AluNames)

# remove old L1 and Alu lists:
Genes <- c(Genes[1:3], Genes[6:length(Genes)])

saveRDS(Genes, file=paste0(RobjectDir, "/custom3_RepeatGenes.rds"))
#Genes <- readRDS(file=paste0(RobjectDir, "/custom3_RepeatGenes.rds"))
