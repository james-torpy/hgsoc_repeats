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
countcores=15
subType="short"

# make directory hierachy:
projectname="hgsoc_repeats"
sampleType="fullsamples/bowtell_primary"
RDir="/home/jamtor/local/lib/r/R-3.2.2/bin/R"

homeDir="/share/ClusterShare/thingamajigs/jamtor/"
projectDir="$homeDir/projects/$projectname/RNA-seq/"
rawDir="$projectDir/raw_files/"

echo -e
echo This is the rawDir:
echo $rawDir

# input/output types:
inType="bam"
outType="fastq"

# extension of files to be used:
inExt=".bam"

# scripts/logs directories:
scriptsPath="/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/"
logDir="$scriptsPath/logs/bam_to_fastq/"
mkdir -p $logDir

echo -e
echo "The logDir is:"
echo $logDir

# set in/outPaths:
inPath="$rawDir/$sampleType/"
outPath="$rawDir/fullsamples/bowtell_primary/fastq"

mkdir -p $outPath

# fetch input bamfiles:
files=( $(ls $inPath/*.bam | grep -v sorted | grep -e 'FT[1-3]\|prPT[1-3]') )

echo -e
echo "The inFiles are: "
echo ${files[@]}
echo -e
echo "The inFiles text file is: "
echo "$rawDir/fullsamples/bowtell_primary/files.txt"
echo ${files[@]} > "$rawDir/fullsamples/bowtell_primary/files.txt"

export largecores=12
export countcores=15

source $scriptsPath/repeatsPipeline/1.bam_to_fastq.bash

# remove files.txt
rm "$rawDir/fullsamples/bowtell_primary/files.txt"


