#!/bin/bash

module load gi/boost/1.53.0

# make directory hierachy
projectname="hgsoc_repeats"
samplename="GC"
expName="exp8"
Type="gc"

homeDir="/share/ClusterShare/thingamajigs/jamtor/"
projectDir="$homeDir/projects/$projectname/"
resultsDir="$projectDir/RNA-seq/results/"
rawDir="$projectDir/RNA-seq/raw_files/fullsamples/bowtell_primary/"

# input/output directories
inType="star"

inPath="$resultsDir/$inType/$samplename/$expName"

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
inFile2=$inPath/$uID/Aligned.sortedByCoord.out.bam
inFile3=$inPath/$uID/Aligned.sortedByCoord.out.bam.bai
inFile4=$inPath/$uID/Aligned.out.bam
inFile5=$rawDir/fastq/$uID.1.fastq.gz
inFile6=$rawDir/fastq/$uID.2.fastq.gz
inFile7=$rawDir/fastq/$uID.unpaired.1.fastq.gz
inFile8=$rawDir/fastq/$uID.unpaired.2.fastq.gz
inFile9=$rawDir/$uID.bam

#echo -e
#echo "Removing $inFile..."
#if [ -a $inFile ]; then
#	rm $inFile
#fi;

echo -e
echo "Removing $inFile2..."
if [ -a $inFile2 ]; then
	rm $inFile2
fi;
echo -e
echo "Removing $inFile3..."
if [ -a $inFile3 ]; then
	rm $inFile3
fi;

echo -e
echo "Removing $inFile4..."
if [ -a $inFile4 ]; then
	rm $inFile4
fi;

echo -e
echo "Removing $inFile5..."
if [ -a $inFile5 ]; then
	rm $inFile5
fi;

echo -e
echo "Removing $inFile6..."
if [ -a $inFile6 ]; then
	rm $inFile6
fi;

echo -e
echo "Removing $inFile7..."
if [ -a $inFile7 ]; then
	rm $inFile7
fi;

echo -e
echo "Removing $inFile8..."
if [ -a $inFile8 ]; then
	rm $inFile8
fi;

echo -e
echo "Removing $inFile9..."
if [ -a $inFile9 ]; then
	rm $inFile9
fi;