#!/bin/bash

echo -e
echo -e
echo "### 3a.starRibo.bash ###"
echo -e

date

module load gi/zlib/1.2.8
module load phuluu/samtools/1.4

# fetch input variables:
cores=$1
fqFile1=$2
fqFile2=$3
outFile=$4

# define refDir:
refDir="/home/jamtor/genomes/ribosome"

# Fetch uID and create outDir:
uID=`basename $fqFile1 | sed "s/\\..*$//"`

outDir=`echo $outFile | sed "s/Log.final.out//"`
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

# align reads of input files with STAR, output into .bam files:
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
     	--genomeDir $refDir \
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
    	--readFilesIn $fqFile1 $fqFile2 \
    	--outFileNamePrefix $outDir \
    	--runThreadN $cores \
    	--outFilterMatchNmin 76 \
      --outSAMtype BAM Unsorted"

echo -e
echo "This is the starline:"
echo $starline
echo -e

# run command:
/home/jamtor/local/bin/STAR --runMode alignReads \
      --readFilesCommand zcat \
      --genomeDir $refDir \
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
      --readFilesIn $fqFile1 $fqFile2 \
      --outFileNamePrefix $outDir \
      --runThreadN $cores \
      --outFilterMatchNmin 76 \
      --outSAMtype BAM Unsorted

 dl=( $(ls $outDir | grep -v Log.final.out) )
 for file in ${dl[@]}; do
   rm -r $outDir/$file
 done;

date