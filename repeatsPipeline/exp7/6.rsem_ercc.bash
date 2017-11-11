#!/bin/bash

# make directory hierachy
projectname="hgsoc_repeats"
expName="exp7"

homeDir="/share/ClusterShare/thingamajigs/jamtor/"
projectDir="$homeDir/projects/$projectname"
resultsDir="$projectDir/RNA-seq/results"

#genome reference file
genomeName="hg38_ercc"
genome_refFile="/share/ScratchGeneral/jamtor/genomes/$genomeName/rsem_ref"

# input/output directories
inType="star/GC"
outType="rsem"

inPath="$resultsDir/$inType/$expName"
outPath="$resultsDir/$outType/$expName"

# load in variable from calling script:
uID=$1
numcores=$2

# scripts/log directory
scriptsDir="$projectDir/RNA-seq/scripts/repeatsPipeline/$expName"
logDir="$scriptsDir/logs"
mkdir -p $logDir

echo -e
echo "The logDir is:"
echo $logDir

inFile="$inPath/$uID/Aligned.sortedByName.out.bam"
outDir="$outPath/$uID/"
out_filePrefix="$outDir/$uID"

mkdir -p $outDir

echo -e
echo This is the uID:
echo $uID
echo -e
echo This is the numcores:
echo $numcores
echo -e
echo This is the inFile:
echo $inFile
echo -e
echo This is the out_filePrefix:
echo $out_filePrefix
echo -e

mkdir -p $outDir

	
#perform analysis on bam files using the reference genome
rsem_line="rsem-calculate-expression -p $numcores --paired-end --bam $inFile $genome_refFile $out_filePrefix"
echo This is the rsem_line:
echo $rsem_line
echo -e
	
rsem-calculate-expression -p $numcores --paired-end --bam $inFile $genome_refFile $out_filePrefix

