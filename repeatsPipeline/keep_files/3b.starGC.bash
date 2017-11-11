#!/bin/bash

echo -e
echo -e
echo "### 3b.starGC.bash ###"
echo -e

# make directory hierachy:
projectname="hgsoc_repeats"
sampleType="fullsamples/bowtell_primary"

homeDir="/share/ClusterShare/thingamajigs/jamtor/"
projectDir="$homeDir/projects/$projectname/RNA-seq/"
rawDir="$projectDir/raw_files/"

# scripts/logs directories:
scriptsPath="/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/"

# set in/outPaths:
fastqPath="$rawDir/fullsamples/bowtell_primary/fastq"

genomeName="hg38_ercc"
genomeDir="/share/ScratchGeneral/jamtor/genomes/$genomeName/starRef/"

echo "This is the genomeDir:"
echo $genomeDir
echo -e

# fetch input variables and define outDir:
uID=$1
fastqPath=$2
projectDir=$3
numcores=$4

#input/output:
outPath="$projectDir/results/star/GC"
  
echo "This is the uID:"
echo $uID
echo -e
echo "This is the fastqPath:"
echo $fastqPath
echo -e
echo "This is the outPath:"
echo $outPath
echo -e

# Fetch inFiles:
inFile1="$fastqPath/$uID.1.fastq.gz"
inFile2="$fastqPath/$uID.2.fastq.gz"
outDir=$outPath/$uID/
  
mkdir -p $outDir
echo -e
echo This is the uID:
echo $uID
echo -e
echo This is the outDir:
echo $outDir
echo -e
echo This is the inFile1:
echo $inFile1
echo -e
echo This is the inFile2:
echo $inFile2
echo -e

# align reads of input files with STAR, output into .bam files:
starJobName="star."$uID
     bamJobName="bam."$uID
     sortJobName="sort."$uID
     indexJobName="index."$uID
     indexStatsJobName="indexstats."$uID
     outSam=$outDir"Aligned.out.sam"
     outBam=$outDir"$uID.bam"
     outSortedBam=$outDir"$uID.sorted.bam"
/home/jamtor/local/bin/STAR --runMode alignReads \
    	--readFilesCommand zcat \
    	--genomeDir $genomeDir \
 		--outFilterType BySJout \
   	--outSAMattributes NH HI AS NM MD\
   	--outFilterMultimapNmax 999 \
    --outMultimapperOrder Random \
    --runRNGseed 666 \
    --outSAMmultNmax 1 \
   	--outFilterMismatchNmax 999 \
   	--outFilterMismatchNoverReadLmax 0.04 \
   	--alignIntronMin 20 \
   	--alignIntronMax 1500000 \
   	--alignMatesGapMax 1500000 \
   	--alignSJoverhangMin 6 \
   	--alignSJDBoverhangMin 1 \
   	--readFilesIn $inFile1 $inFile2 \
   	--outFileNamePrefix $outDir \
   	--runThreadN $numcores \
   	--outFilterMatchNmin 76 \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM 95927520246

