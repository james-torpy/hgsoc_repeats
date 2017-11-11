#!/bin/bash

### 1a.call_count_overlaps.bash ###
# This script calls count_overlaps.R and submits it with appropriate to the
# variables to the cluster

# define starting variables:
project="hgsoc_repeats"
methodName="method1"
gtfName="human-89.repeats4.gtf"
numcores="10"

# define directories:
homeDir="/share/ScratchGeneral/jamtor/"
projectDir="$homeDir/projects/$project"
bamDir="$projectDir/RNA-seq/results/star/$methodName"
refDir"$projectDir/RNA-seq/refs/"
resultsDir="$projectDir/RNA-seq/results")
RobjectDir="$projectDir/RNA-seq/Robjects/")
outDir="$resultsDir/R/$methodName"
scriptDir="$projectDir/RNA-seq/scripts/R"
logDir="$scriptDir/logs"

# create outDir, logDir:
mkdir -p outDir
mkdir -p logDir

echo -e
echo "The bamDir is:"
echo $bamDir
echo -e
echo "The refDir is:"
echo $refDir
echo -e
echo "The RobjectDir is:"
echo $RobjectDir
echo "The outDir is:"
echo $outDir

# define R script to be called:
callScript="$scriptDir/1b.count_overlaps_draft10.R"

# fetch GTF file:
gtf="$refDir/$gtfName"

echo -e
echo "The gtf file is:"
echo $gtf

# fetch the inFiles (bams):
inFiles=( $(ls $bamDir/**/Aligned.sortedByCoord.out.bam | grep -V bai | grep -V subset ))

echo -e
echo "The inFiles are:"
echo $inFiles

# make R script specifying variables for count_overlaps.R:
echo "methodName='$methodName'" >> "$scriptDir/$methodName_vars.R"
echo "bamDir='$bamDir'" >> "$scriptDir/$methodName_vars.R"
echo "RobjectDir='$RobjectDir'" >> "$scriptDir/$methodName_vars.R"
echo "outDir='$outDir'" >> "$scriptDir/$methodName_vars.R"
echo "gtf='$gtf'" >> "$scriptDir/$methodName_vars.R"
echo "inFiles='$inFiles'" >> "$scriptDir/$methodName_vars.R"

#Submit the job to the processing script, passing on arguments:
  #qsub -b y -j y -N count$methodName -wd $logDir -pe smp $numcores -V $RDir "--vanilla < $callScript $scriptDir/$methodName"_vars.R








