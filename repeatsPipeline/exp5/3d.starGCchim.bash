#!/bin/bash

module load phuluu/samtools/1.4

echo -e
echo -e
echo "### 3b.starGC.bash ###"
echo -e

# make directory hierachy:
projectname="hgsoc_repeats"
sampleType="fullsamples/bowtell_primary"
expName="exp5"

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

#input/output:
outPath="$projectDir/results/star/GC/$expName"
uID=$1
numcores=$2

echo "This is the uID:"
echo $uID
echo -e
echo "This is the numcores:"
echo $numcores
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

#align reads of input files with STAR, output into .bam files:
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
  --chimSegmentMin 25 \
    --chimJunctionOverhangMin 25 \
    --chimScoreMin 0 \
    --chimScoreDropMax 20 \
    --chimScoreSeparation 10 \
    --chimScoreJunctionNonGTAG -1 \
    --outSAMtype BAM Unsorted