#!/bin/bash

echo -e
echo "### 4.htseq_gc.bash ###"
echo -e
echo -e

source /home/jamtor/.bashrc

module load gi/boost/1.53.0
module load gi/zlib/1.2.8
module load phuluu/samtools/1.4
module load borgue/HTSeq/0.6.1

date

cores=$1
bamFile=$2
outFile=$3

# define genome reference file:
gffFile="/home/jamtor/genomes/hg38_ercc/gencode_v24_hg38_annotation.gff"

# create outDir:
outDir=`echo ${outFile%/*}`/
mkdir -p $outDir

echo -e
echo This is the bamFile:
echo $bamFile
echo -e
echo This is the outFile:
echo $outFile
echo -e
echo This is the gffFile:
echo $gffFile
echo -e
	
#perform analysis on bam files using the reference genome
htseq_line="htseq-count -f bam -i ID -t exon $bamFile $gffFile >> $outFile"
echo This is the htseq_line:
echo $htseq_line
echo -e

#/usr/bin/python -m HTSeq.scripts.count -f bam -i ID -t exon $inFile $gffFile >> $outFile
htseq-count -f bam -i ID -t exon -a 0 -m intersection-strict $bamFile $gffFile >> $outFile
#python -m htseq-count -f bam -i ID -t exon $inFile $gffFile >> $outFile


date