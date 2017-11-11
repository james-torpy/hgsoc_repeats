#!/bin/bash

module load gi/boost/1.53.0

# make directory hierachy
projectname="hgsoc_repeats"
samplename="GC"
expName="exp7"
Type="gc"

homeDir="/share/ClusterShare/thingamajigs/jamtor/"
projectDir="$homeDir/projects/$projectname"
resultsDir="$projectDir/RNA-seq/results"

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

#echo -e
#echo "Removing $inFile..."
#if [ -a $inFile ]; then
#	rm $inFile
#fi;

#echo -e
#echo "Removing $inFile2..."
#if [ -a $inFile2 ]; then
#	rm $inFile2
#fi;

#echo -e
#echo "Removing $inFile3..."
#if [ -a $inFile3 ]; then
#	rm $inFile3
#fi;

echo -e
echo "Removing $inFile4..."
if [ -a $inFile4 ]; then
	rm $inFile4
fi;