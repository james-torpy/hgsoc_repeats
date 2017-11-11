#!/bin/bash

numcores=6

module load gi/boost/1.53.0

#genome directories:
genomeName="hg38_ercc"
annotationName="gencode_and_repeat"

projectname="hgsoc_repeats"

homeDir="/share/ScratchGeneral/jamtor"
projectDir="$homeDir/projects/$projectname"
genomeDir="$homeDir/genomes/gencode_and_repeat"
genomeFile="$genomeDir/$genomeName.fa"
annotationFile="$genomeDir/$annotationName.gtf"
outDir="$genomeDir/rsem_ref"

mkdir -p $outDir

#log directory:
logDir="/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/mapping/logs"
mkdir -p $logDir

echo This is the genomeDir:
echo $genomeDir
echo -e
echo This is the genomeFile:
echo $genomeFile
echo -e
echo This is the annotationFile:
echo $annotationFile
echo -e
echo This is the outDir
echo $outDir
echo -e

#generate the rsem reference files:

#Set up conditions
rsem_ref_line="rsem-prepare-reference --gtf $annotationFile $genomeFile $outDir"

echo This is the rsem_ref_line:
echo $rsem_ref_line

#qsub -N RSEM_ref_$genome -wd $logDir -b y -j y -R y -pe smp $numcores -V $rsem_ref_line
