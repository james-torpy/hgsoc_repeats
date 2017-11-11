#!/bin/bash

module load gi/boost/1.53.0

numcores=4

# make directory hierachy
projectname="hgsoc_repeats"
samplename="GC"
expName="exp5"
Type="gc"

homeDir="/share/ClusterShare/thingamajigs/jamtor/"
projectDir="$homeDir/projects/$projectname"
resultsDir="$projectDir/RNA-seq/results"

# genome reference file
genomeName="hg38_ercc"
gffFile="/share/ScratchGeneral/jamtor/genomes/$genomeName/gencode_v24_hg38_annotation.gff"

# input/output directories
inType="star"
outType="htseq"

inPath="$resultsDir/$inType/$samplename/$expName"
outPath="$resultsDir/$outType/$expName"

# load in variable from calling script:
uID=$1

# scripts/log directory
scriptsDir="$projectDir/RNA-seq/scripts/repeatsPipeline/$expName"
logDir="$scriptsDir/logs"
mkdir -p $logDir

echo -e
echo "The logDir is:"
echo $logDir

inFile=$inPath/$uID/Aligned.sortedByName.out.bam
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
echo This is the gffFile:
echo $gffFile
echo -e
echo This is the outDir:
echo $outDir
echo -e
	
#perform analysis on bam files using the reference genome
htseq_line="htseq-count -f bam -i ID -t exon -o $outDir/$uID.$Type.out.sam $inFile $gffFile >> $outDir/$uID.$Type.htseq.txt"
echo This is the htseq_line:
echo $htseq_line
echo -e

htseq-count -f bam -i ID -t exon -o $outDir/$uID.$Type.out.sam $inFile $gffFile >> $outDir/$uID.$Type.htseq.txt