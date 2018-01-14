#!/bin/bash

echo -e
echo "### 4.htseq_gc.bash ###"
echo -e
echo -e

module load gi/boost/1.53.0
module load gi/zlib/1.2.8
module load phuluu/samtools/1.4

date

cores=6
bamFile="/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/results/star/GC/exp9/test/FT1/Aligned.novosortedByName.out.bam"
outFile="/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/results/htseq/exp9/test/FT1/FT1.custom3.htseq.txt"

# define gff file name:
gff="/home/jamtor/genomes/repeats/custom3rep.final.gff"

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
echo $gff
echo -e

#perform analysis on bam files using the reference genome:
htseq_line="htseq-count -f bam -i ID -t exon --stranded=no $bamFile $gff >> $outFile"
echo This is the htseq_line:
echo $htseq_line
echo -e

#/usr/bin/python -m HTSeq.scripts.count -f bam -i ID -t exon $inFile $gffFile >> $outFile
#python -m HTSeq.scripts.count -f bam -i ID -t exon $inFile $gffFile >> $outFile
htseq-count -f bam -i ID -t exon --stranded=no -a 0 -m intersection-strict $bamFile $gff >> $outFile


date
