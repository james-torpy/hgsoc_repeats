#!/bin/bash

echo -e
echo -e
echo "### 3c.novosort.bash ###"
echo -e

cores=$1
bamFile=$2
outBamFile=$3
mem=$4

date

module load gi/novosort/precompiled/1.03.08

# make directory hierachy:
projectname="hgsoc_repeats"
expName="exp8"

homeDir="/share/ClusterShare/thingamajigs/jamtor/"
projectDir="$homeDir/projects/$projectname/RNA-seq/"

# scripts/logs directories:
scriptsPath="/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/"

uID=$1

#input/output:
starDir="$projectDir/results/star/GC/$expName/$draft/$uID"
tempDir="$starDir/$uID/temp"
mkdir -p $tempDir

echo "This is the uID:"
echo $uID
echo -e
echo "This is the starDir:"
echo $starDir
echo -e

/home/jamtor//local/bin/novocraft/novosort -t $tempDir -n -c 6 -m 22G \
	$starDir/Aligned.out.bam > $starDir/Aligned.novosortedByName.out.bam

date