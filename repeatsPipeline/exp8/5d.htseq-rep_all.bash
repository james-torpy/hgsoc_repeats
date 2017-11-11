#!/bin/bash

module load gi/boost/1.53.0

numcores=4

# make directory hierachy
projectname="hgsoc_repeats"
samplename="GC"
expName="exp8"
Type="all"
draft=$2

date

homeDir="/share/ClusterShare/thingamajigs/jamtor/"
projectDir="$homeDir/projects/$projectname"
resultsDir="$projectDir/RNA-seq/results"
refDir="/share/ScratchGeneral/jamtor/projects/$projectname/RNA-seq/refs/$Type"

# input/output directories
inType="star"
outType="htseq"

inPath="$resultsDir/$inType/$samplename/$expName/$draft/"
outPath="$resultsDir/$outType/$expName/$draft/"

# load in variable from calling script:
uID=$1

# scripts/log directory
scriptsDir="$projectDir/RNA-seq/scripts/repeatsPipeline/$expName"
logDir="$scriptsDir/logs"
mkdir -p $logDir

echo -e
echo "The logDir is:"
echo $logDir
echo -e
echo "The inPath is:"
echo $inPath

inFile=$inPath/$uID/Aligned.novosortedByName.out.bam
outDir=$outPath/$uID
out_filePrefix=$outDir/$uID

mkdir -p $outDir

echo -e
echo This is the uID:
echo $uID
echo -e
echo This is the inFile:
echo $inFile
echo -e
echo This is the outDir:
echo $outDir
echo -e

# define gtf file name:
gtf="$refDir/all_rep.final.gtf"
gtfID=`basename $gtf | sed 's/.gtf//' | sed 's/.final//'`

echo This is the gtfFile:
echo $gtf
echo -e

#perform analysis on bam files using the reference genome
htseq_line="htseq-count -f bam -i ID -t exon --stranded=no $inFile $gtf >> $outDir/$uID.$Type.htseq.txt"
echo This is the htseq_line:
echo $htseq_line
echo -e

htseq-count -f bam -i ID -t exon --stranded=no $inFile $gtf >> $outDir/$uID.$Type.htseq.txt

echo "Removing $inFile"
echo -e
rm $inFile

date
