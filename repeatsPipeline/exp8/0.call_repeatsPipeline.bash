#!/bin/bash

### 0.call_repeatsPipeline.bash ###

# Takes paired read bam files and converts them to fastq files, one for 1st
# paired reads, one for 2nd paired reads
# Calls rest of repeatsPipeline

echo -e
echo "### 0.call_repeatsPipeline.bash ###"
echo -e

module load phuluu/samtools/1.4

# specify cores needed:
fqcores=2
ribocores=10
gccores=6
sortcores=7
countcores=12
mem=24G
subType="long"
expName="exp8"
draft=""

# make directory hierachy:
projectname="hgsoc_repeats"
sampleType="fullsamples/bowtell_primary"
RDir="/home/jamtor/local/lib/r/R-3.2.2/bin/R"

homeDir="/share/ClusterShare/thingamajigs/jamtor/"
projectDir="$homeDir/projects/$projectname/RNA-seq/"
rawDir="$projectDir/raw_files/$sampleType"

echo -e
echo This is the rawDir:
echo $rawDir

# input/output types:
inType="bam"
outType="fastq"

# extension of files to be used:
inExt=".bam"

# scripts/logs directories:
scriptsPath="/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/$expName"
logDir="$scriptsPath/logs/"
mkdir -p $logDir

echo -e
echo "The logDir is:"
echo $logDir

# set in/outPaths:
inPath="$rawDir"

# fetch input bamfiles:
files=( $(ls $inPath/*.bam | grep -v sorted) )

echo -e
echo "The inFiles are: "
echo ${files[@]}
echo -e
echo "The inFiles text file is: "
echo "$rawDir/files.txt"
echo ${files[@]} > "$rawDir/files.txt"

export fqcores=$fqcores
export ribocores=$ribocores
export gccores=$gccores
export sortcores=$sortcores
export countcores=$countcores
export draft=$draft
export subType=$subType
export expName=$expName
export mem=$mem

source $scriptsPath/1.bamtofastq.bash



