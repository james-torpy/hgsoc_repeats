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
smallcores=1
largecores=12
countcores=12
subType="long"
expName="exp7"

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
outPath="$rawDir/fastq"

mkdir -p $outPath

# fetch input bamfiles:
files=( $(ls $inPath/*.bam | grep -v sorted) )

echo -e
echo "The inFiles are: "
echo ${files[@]}
echo -e
echo "The inFiles text file is: "
echo "$rawDir/files.txt"
echo ${files[@]} > "$rawDir/files.txt"

export largecores=$largecores
export countcores=$countcores

source $scriptsPath/1.bam_to_fastq.bash

# remove files.txt
#rm "$rawDir/files.txt"


