#!/bin/bash

echo -e
echo -e
echo "### 3b.starGC.bash ###"
echo -e

date

# fetch input variables:
cores=$1
fqFile1=$2
fqFile2=$3
outPath=$4

# define refDir:
refDir="/home/jamtor/genomes/hg38_ercc/starRef/"

# Fetch uID and create outDir:
uID=`basename fqFile1 | sed "s/\\..*$//"`
outDir="$outPath/$uID/"
    
mkdir -p $outDir

echo "This is the fqFile1:"
echo $fqFile1
echo -e
echo "This is the fqFile2:"
echo $fqFile2
echo -e
echo "This is the outPath:"
echo $outPath
echo -e
echo "This is the refDir:"
echo -e
echo "This is the uID:"
echo $uID
echo -e
echo "This is the outDir:"
echo $outDir

#align reads of input files with STAR, output into .bam files:
starJobName="star."$uID
     bamJobName="bam."$uID
     sortJobName="sort."$uID
     indexJobName="index."$uID
     indexStatsJobName="indexstats."$uID
     outSam=$outDir"Aligned.out.sam"
     outBam=$outDir"$uID.bam"
     outSortedBam=$outDir"$uID.sorted.bam"

starline="/home/jamtor/local/bin/STAR --runMode alignReads \
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
    --outSAMtype BAM Unsorted"

echo -e
echo "This is the starline:"
echo $starline

# run command:
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

date