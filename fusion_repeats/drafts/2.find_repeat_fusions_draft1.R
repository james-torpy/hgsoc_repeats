### 3.find_repeat_fusions.R ###

# This script takes filtered STAR-chimera output and creates plots
# of repeat fusions that have occurred


### 0. Define variables/paths ###

library(rtracklayer)
library(GenomicRanges)
library("BiocParallel")
library(org.Hs.eg.db)

project <- "hgsoc_repeats"
expName <- "exp9"

Type <- "fusion_repeats"
descrip <- "fusion_repeats"


# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project, "/")
resultsDir <- paste0(projectDir, "/RNA-seq/results")
starDir <- paste0(resultsDir, "/star/", descrip, "/")
rpAnnotDir <- paste0(projectDir, "/RNA-seq/refs/repeats/")
gcAnnotDir <- paste0(projectDir, "/RNA-seq/refs/")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
                     expName, "/", descrip, "/")
plotDir <- paste0(resultsDir, "/R/", expName,
                  "/plots/", descrip)

system(paste0("mkdir -p ", plotDir))
system(paste0("mkdir -p ", RobjectDir))


### 1. load in filtered chimera file and repeat annotation ###

chimeras <- read.table(file=paste0(starDir, "/prPT4/chimeras_filtered.txt"))
colnames(chimeras) <- c("readNo", "donorChr", "donorBP1", "donorStrand", 
                        "acceptorChr", "acceptorBP1", "acceptorStrand", 
                        "junction_type", "leftRepeatLength", "rightRepeatLength")
chimericRanges <- GRanges(
  seqnames = Rle(chimeras$donorChr),
  ranges = IRanges(chimeras$donorBP1, end = chimeras$donorBP1),
  strand = Rle(chimeras$donorStrand),
  readNo = chimeras$readNo,
  donorStrand = chimeras$donorStrand,
  chimeras$acceptorChr,
  chimeras$acceptorBP1,
  chimeras$acceptorStrand,
  chimeras$junction_type
)

repeatAnnot <- import(paste0(rpAnnotDir, "/human-89.repeats.gtf"))
gcAnnot <- import(paste0(gcAnnotDir, "/gencode_v24_hg38_annotation.gtf"))


### 2. Find rOverlaps of chimeras with repeats annotation ranges ###

rOverlaps <- findOverlaps(chimericRanges, repeatAnnot)


# find maximum number of hits for the same chimera:
#q <- rle(queryHits(rOverlaps))
#maxHits <- max(q$lengths)

# add more rows onto chimericRanges to input repeat annotations:
#df <- data.frame(matrix(NA, nrow = maxHits, ncol=length(chimericRanges)))

#for (n in 1:maxHits) {
#  chimericRanges <- (chimericRanges, rep(NA, length(chimericRanges)))
#}

chimericRanges$repeatHits <- rep(NA, length(chimericRanges))

for (i in 1:length(queryHits(rOverlaps))) {
  qHit <- queryHits(rOverlaps)[i]
  sHit <- subjectHits(rOverlaps)[i]
  print(qHit)
  
  if (i==1) {
    record <- c(qHit)
  } else {
    record[i] <- qHit
  }
  
  dup <- table(record)
  dupNo <- dup[names(dup)==qHit]
  
  if (dupNo==1) {
    chimericRanges$repeatHits[qHit] <- as.character(repeatAnnot[sHit]$type)
  } else {
    chimericRanges[qHit]$repeatHits <- paste0(chimericRanges[qHit]$repeatHits, ",", as.character(repeatAnnot[sHit]$type))
  }
}

repeat_chimeras <- chimericRanges[!is.na(chimericRanges$repeatHits)]

# remove 'chimera$' from acceptor colnames so these columns can later be referred
# to without problems:
colnames(values(repeat_chimeras))[3:6] <- gsub(
  "chimeras\\$", "", colnames(values(repeat_chimeras))[3:6]
)

# find acceptor sequence gene names from both gencode and repeat annotation GR:
# merge gc and repeats annotations:
allAnnot <- c(
  GRanges(
    seqnames=seqnames(repeatAnnot),
    ranges=ranges(repeatAnnot),
    strand=strand(repeatAnnot),
    source=repeatAnnot$source,
    gene_id=repeatAnnot$type,
    gene_type=NA,
    gene_name=NA,
    transcript_id=NA,
    transcript_type=NA,
    transcript_name=NA,
    exon_number=NA,
    exon_id=NA
  ),

  GRanges(
    seqnames=seqnames(gcAnnot),
    ranges=ranges(gcAnnot),
    strand=strand(gcAnnot),
    source=gcAnnot$source,
    gene_id=gcAnnot$gene_id,
    gene_type=gcAnnot$gene_type,
    gene_name=gcAnnot$gene_name,
    transcript_id=gcAnnot$transcript_id,
    transcript_type=gcAnnot$transcript_type,
    transcript_name=gcAnnot$transcript_name,
    exon_number=gcAnnot$exon_number,
    exon_id=gcAnnot$exon_id
  )
)

# make GR object out of acceptor co-ordinates from repeat_chimeras:
acceptorRanges <- GRanges(
  seqnames=repeat_chimeras$acceptorChr,
  ranges=IRanges(repeat_chimeras$acceptorBP1, repeat_chimeras$acceptorBP1),
  strand=repeat_chimeras$acceptorStrand
)

# find overlaps between acceptor co-ordinates and allAnnotation:
acceptorOverlaps <- findOverlaps(acceptorRanges, allAnnot)

# add acceptor location information to repeat_chimeras:
repeat_chimeras$acceptorHits <- rep(NA, length(repeat_chimeras))

for (j in 1:length(queryHits(acceptorOverlaps))) {
  qHit <- queryHits(acceptorOverlaps)[j]
  sHit <- subjectHits(acceptorOverlaps)[j]
  print(qHit)
  
  if (j==1) {
    record <- c(qHit)
  } else {
    record[j] <- qHit
  }
  
  dup <- table(record)
  dupNo <- dup[names(dup)==qHit]
  
  if (dupNo==1) {
    repeat_chimeras$acceptorHits[qHit] <- as.character(allAnnot[sHit]$gene_id)
  } else {
    repeat_chimeras[qHit]$acceptorHits <- paste0(repeat_chimeras[qHit]$acceptorHits, 
    	",", as.character(allAnnot[sHit]$gene_id))
  }
}

save.image(paste0(RobjectDir, "/acceptors_annotated.RData"))

# delete duplicates within repeatHits or acceptorHits entries:
for (n in 1:length(repeat_chimeras)) {
  print(n)
  newRhits <- unique(strsplit(repeat_chimeras$repeatHits[n], ",")[[1]])
  if (length(newRhits) > 1) {
    repeat_chimeras$repeatHits[n] <- paste(newRhits, collapse=",")
  } else {
    repeat_chimeras$repeatHits[n] <- newRhits
  }
  
  newAcHits <- unique(strsplit(repeat_chimeras$acceptorHits[n], ",")[[1]])
  if (length(newAcHits) > 1) {
    repeat_chimeras$acceptorHits[n] <- paste(newAcHits, collapse=",")
  } else {
    repeat_chimeras$acceptorHits[n] <- newAcHits
  }
}

# isolate repeat fusions with more than 6 supporting reads:
strongFusions <- repeat_chimeras[repeat_chimeras$readNo > 5]
strongFusions <- strongFusions[!(strongFusions$repeatHits=="trf" | strongFusions$repeatHits=="dust")]
strongDF <- data.frame(strongFusions$repeatHits, strongFusions$readNo)
colnames(strongDF) <- c("repeat_id", "no_of_reads")

# create vector of donor and acceptor hit names:
both_names <- paste0(strongDF$repeat_id, " -> ", strongFusions$symbol)
both_names <- gsub("NA", "unknown", both_names)

par(mar=c(4,4,4,4))
pdf(file = paste0(plotDir, "/strong_repeats_fusion_read_nos.pdf"))
barplot(strongDF$no_of_reads, names.arg = both_names, cex.names = 0.5, 
        ylab = "no. reads supporting fusion", las=2)
dev.off()

save.image(paste0(RobjectDir, "/strong_fusions_created.RData"))

# place ensembl_ids in separate columns in repeat_chimeras:
# find out the biggest number of hits within acceptorHits values:
for (m in 1:length(repeat_chimeras)) {
  if (m==1) {
    hit_nos <- c(length(strsplit(repeat_chimeras$acceptorHits[m], ",")[[1]]))
  } else {
    hit_nos[m] <- length(strsplit(repeat_chimeras$acceptorHits[m], ",")[[1]])
  }
}
max_no <- max(hit_nos)

ensemblID_cols <- data.frame(matrix(ncol = max_no, nrow = length(repeat_chimeras)))






for (m in 1:length(repeat_chimeras)) {
  spl <- strsplit(repeat_chimeras$acceptorHits[m], ",")[[1]]
  for (s in 1:length(spl)) {
    print(spl[s])
    if (length(grep("ENSG", spl[s]))>0) {
      values(repeat_chimeras)[,(ncol(values(repeat_chimeras))+1)] <- NA
      
    }
  }
  
  
  
  if (length(grep("ENSG", repeat_chimeras$acceptorHits[m])) > 1) {
    print(m)
    repeat_chimeras$ensembl_id[m] <- repeat_chimeras$acceptorHits[m]
  } else if (length(grep("ENSG", repeat_chimeras$acceptorHits[m])) > 0) {
    print(m)
    id <- gsub(
      "^.*,", "", gsub(
        ",*,$", "",repeat_chimeras$acceptorHits[m]
      )
    )
    repeat_chimeras$ensembl_id[m] <- id
  }
}

# convert ensembl_ids to gene_ids and add to additional column in repeat_chimeras:
egENSEMBL <- toTable(org.Hs.egENSEMBL)
egSYMBOL <- toTable(org.Hs.egSYMBOL)

repeat_chimeras$symbol <- NA
for (o in 1:length(repeat_chimeras)) {
  if (!is.na(repeat_chimeras$ensembl_id[o])) {
    print(o)
    ensembl_id <- gsub("\\..*$", "", repeat_chimeras$ensembl_id[o])
    gene_id <- egENSEMBL$gene_id[match(ensembl_id, egENSEMBL$ensembl_id)]
    symbol <- egSYMBOL$symbol[match(gene_id, egSYMBOL$gene_id)]
    repeat_chimeras$symbol[o] <- symbol
  }
  if(length(grep("ENSG", repeat_chimeras$acceptorHits[o])) < 1) {
    repeat_chimeras$symbol[o] <- repeat_chimeras$acceptorHits[o]
  }
}


customList <- rownames(readRDS(file=paste0(projectDir, 
	"/RNA-seq/Robjects/exp9/htseq_allHGSOC_vs_FT/custom3_DEsigReps.rds"))[[1]])

for (k in 1:length(customList)) {
	print(k)
	rp <- grep(customList[k], repeat_chimeras$repeatHits)
	if (length(rp) > 0) {
		print(paste0("Adding ", customList[k], " to intFusions"))
		if (k==1) {
			intFusions <- repeat_chimeras[rp]
		} else {
			intFusions <- c(intFusions, repeat_chimeras[rp])
		}
	}
}

intFusions <- intFusions[intFusions$readNo > 1]


intDF <- data.frame(intFusions$repeatHits, intFusions$readNo, intFusions$acceptorHits)
colnames(intDF) <- c("repeat_id", "no_of_reads")
pdf(file = paste0(plotDir, "/interesting_repeats_fusion_read_nos.pdf"))
barplot(intDF$no_of_reads, names.arg = intDF$repeat_id, cex.names = 0.4, ylab = "no. reads supporting fusion", xlab = "repeat id of fusion")
dev.off()




# isolate repeat fusions with more than 2 supporting reads:
trustFusions <- repeat_chimeras[repeat_chimeras$readNo > 2]

## isolate L1s and Alus:
#L1s <- repeat_chimeras[grep("L1", repeat_chimeras$repeatHits)]
#L1s <- GRanges(
#  seqnames = seqnames(L1s),
#  ranges = ranges(L1s),
#  strand = strand(L1s),
#  repeatHits = L1s$repeatHits,
#  readNo = L1s$readNo
#)
#
#for (n in 1:length(L1s)) {
#  print(L1s[n]$repeatHits)
#  rp <- strsplit(L1s[n]$repeatHits, ",")[[1]]
#  L1s[n]$repeatHits <- paste(rp[!duplicated(rp)], sep=",", collapse=",")
#}
#
#Alus <- repeat_chimeras[grep("Alu", repeat_chimeras$repeatHits)]
#Alus <- GRanges(
#  seqnames = seqnames(Alus),
#  ranges = ranges(Alus),
#  strand = strand(Alus),
#  repeatHits = Alus$repeatHits,
#  readNo = Alus$readNo
#)
#
#for (n in 1:length(Alus)) {
#  print(Alus[n]$repeatHits)
#  rp <- strsplit(Alus[n]$repeatHits, ",")[[1]]
#  Alus[n]$repeatHits <- paste(rp[!duplicated(rp)], sep=",", collapse=",")
#}
#
#
#L1sDF <- data.frame(L1s$repeatHits, L1s$readNo)
#colnames(L1sDF) <- c("repeat_id", "no_of_reads")
#L1sDF <- L1sDF[L1sDF$no_of_reads > 2,]
#pdf(file = paste0(plotDir, "/L1_fusion_read_nos.pdf"))
#barplot(L1sDF$no_of_reads, names.arg = L1sDF$repeat_id, cex.names = 0.2, ylab = "no. reads supporting fusion", xlab = "repeat id of fusion")
#dev.off()
#
#AlusDF <- data.frame(Alus$repeatHits, Alus$readNo)
#colnames(AlusDF) <- c("repeat_id", "no_of_reads")
#AlusDF <- AlusDF[AlusDF$no_of_reads > 2,]
#pdf(file = paste0(plotDir, "/Alu_fusion_read_nos.pdf"))
#barplot(AlusDF$no_of_reads, names.arg = AlusDF$repeat_id, cex.names = 0.2, ylab = "no. reads supporting fusion", xlab = "repeat id of fusion")
#dev.off()
