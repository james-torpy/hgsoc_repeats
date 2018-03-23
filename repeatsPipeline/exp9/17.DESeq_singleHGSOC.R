
# Run on cluster with:
#briR
#qsub -N DESeq -b y -wd \
#/share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/logs/exp9/DE \
#-j y -R y -pe smp 4 -V "Rscript /share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/exp9/17a.DESeq_indivHGSOC_vs_FT.R"


### 0. Define variables/paths ###

# load packages needed:
library(DESeq)
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
descrip <- "htseq_DESeq_CCNEamp_vs_HRD"

sTypes <- c("bothDrivers", "FT", "HRD", "CCNEamp", "unknown_driver")
descrip <- "htseq_DESeq_HGSOC_CCNEamp_unknown_vs_HRD"

# define sample group to use as control:
ctl <- "CCNEamp"

# define sample groups to compare:
allHGSOC_vs_FT <- FALSE
cat_by_driver <- TRUE
EDAnormalise <- TRUE

# specify what combination of repeat genes (repeats), epigenetic modulators (epiMods),
# RNAi genes (RNAi) and protein-coding genes (pCoding) should contribute to the results:
#resultTypes <- c("repeats", "epiMods")
resultTypes <- c("repeats", "all")

# specify what FDR and log2 fold change thresholds to use:
FDRthresh <- 0.05
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
                  "/plots/DEplots/", descrip, "/")

system(paste0("mkdir -p ", plotDir))
system(paste0("mkdir -p ", newRobjectDir))


if (file.exists(paste0(RobjectDir, "/post_DESeq.RData"))) {
	load(paste0(RobjectDir, "/post_DESeq.RData"))
} else {
  
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
	
	# remove pAF and remove categorisations of HGSOC samples:
	Counts <- Counts[,grep("pAF", colnames(Counts), invert=T)]
	#colnames(Counts) <- gsub("_[a-z][a-z][A-Z][A-Z]$", "", colnames(Counts))
	# append '.2' onto duplicate patient IDs:
	#colnames(Counts)[duplicated(colnames(Counts))] <- paste0(colnames(Counts)[duplicated(colnames(Counts))], ".2")
	
	# change the order of columns of Counts to 0.1betical order of subtypes:
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
	
	# save number of samples in each group:
	saveRDS(splt, file = paste0(newRobjectDir, "/sample_no_per_cat.rds"))
	
	# delist elements need to be delisted and change to integers:
	Counts <- apply(Counts, 2, unlist)
	storage.mode(Counts) <- "integer"
	
	if ( allHGSOC_vs_FT ) {
	  # add '.2' to IDs of samples with duplicate names:
	  colnames(Counts)[duplicated(colnames(Counts))] <- gsub("_HGSOC", ".2_HGSOC", colnames(Counts)[duplicated(colnames(Counts))])
	  
	  #HGSOC <- colnames(Counts)[grep("FT", colnames(Counts), invert=T)]
	  HGSOC <- Counts[,grep("FT", colnames(Counts), invert=T)]
	}

	
	### 3. Split counts into FT vs each HGSOC data set and process separately ###
	
	# define upTo:
	upTo <- 1

	# load image if loop below has previously been partially completed:
	if (file.exists(paste0(RobjectDir, "/DElist.rds"))) {
		DElist <- readRDS(paste0(RobjectDir, "/DElist.rds"))
		upTo <- length(DElist)
	}

	# define conditions:
	for (i in upTo:length(colnames(HGSOC))) {
		writeLines("\n")
		
		# define samplename and remove any brackets:
		sName <- gsub("\\(|\\)", "", colnames(HGSOC)[i])
		print(paste0("Processing sample number ", i, ", i.e. ", sName))
		# build counts df:
		countDF <- cbind(Counts[,grep("FT", colnames(Counts))], HGSOC[,i])
		colnames(countDF)[length(grep("FT", colnames(Counts)))+1] <- colnames(HGSOC)[i]
		
		conditions <- factor(gsub("^.*FT", "FT", colnames(countDF)), levels = c("FT", colnames(countDF)[(grep("FT", colnames(countDF), invert=T))]))
		
		# convert countDF into SeqExpressionSet
		set <- newSeqExpressionSet(countDF, phenoData = data.frame(conditions, row.names=colnames(countDF)))
		
		# create pre-norm RLE plot:
		if (file.exists(paste0(plotDir, "/", Type, "_", sName, "_vs_FT_RLEPrenormGC.pdf"))) {
		  print(paste0(plotDir, "/", Type, "_", sName, "_vs_FT_RLEPrenormGC.pdf already exists, no need to create"))
		} else {
		  print(paste0(plotDir, "/", Type, "_", sName, "_vs_FT_RLEPrenormGC.pdf"))
		  par(mar=c(1,1,1,1))
		  pdf(file = paste0(plotDir, "/", Type, "_", sName, "_vs_FT_RLEPrenormGC.pdf"))
		  plotRLE(set)
		  dev.off()
		}
		
		# create RUVseq pre-norm PCA:
		if (file.exists(paste0(plotDir, "/", Type, "_", sName, "_vs_FT__pcaPrenormGC.pdf"))) {
		  print(paste0(plotDir, "/", Type, "_", sName, "_vs_FT__pcaPrenormGC.pdf already exists, no need to create"))
		} else {
		  print(paste0("Creating ", plotDir, "/", Type, "_", sName, "_vs_FT__pcaPrenormGC.pdf"))
		  pdf(file = paste0(plotDir, "/", Type, "_", sName, "_vs_FT__pcaPrenormGC.pdf"), height = 10, width = 12)
		  plotPCA(set, cex=0.7)
		  dev.off()
		}
		
		
		### 4. perform normalisation on countDF using EDAseq:
		
		# perform between lane full normalisation:
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
		
		
		### 5. Perform DESeq on countDF ###
		# convert to DESeq data set (from EDAseq manual page 12):
		# rename typeF to 'conditions' in nSet, for DESeq:
		countsN <- as(nSet, "CountDataSet")
		sizeFactors(countsN) <- rep(1, length(sampleNames(phenoData(countsN))))
		#varLabels(phenoData(countsN)) <- c("sizeFactor", "condition")
		#varMetadata(phenoData(countsN))[2,1] <- "experimental condition, treatment or phenotype"
		
		# design matrix labelling all sample types:
		#design <- model.matrix(~0+conditions, data=pData(nSet))
		
		countsN <- estimateDispersions(countsN)
		res <- nbinomTest(countsN, "FT", colnames(HGSOC)[i])
		
		# save to results list:
		if (i==1) {
		  DElist <- list(res)
		} else {
		  DElist[[i]] <- res
		}
		names(DElist) <- paste0(colnames(HGSOC), "_vs_FT")
		saveRDS(DElist, file=paste0(RobjectDir, "/", descrip, "_DElist.rds"))
		# save progress to image:
		
		save.image(file=paste0(RobjectDir, "/post_DESeq.RData"))
	}
}	




