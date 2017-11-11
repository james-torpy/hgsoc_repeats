#!/bin/bash

#make modules to load here

module load gi/samtools/1.2

#number of cores
numcores=6

#genome directories
genomeName="hg38_ercc_L1"
annotName="gencode_L1"
projectname="hgsoc_repeats"

homeDir="/share/ScratchGeneral/jamtor"
projectDir="$homeDir/projects/$projectname"
genomeDir="$homeDir/genomes/hg38_ercc/$genomeName"
genomeFile="$genomeDir/$genomeName.fa"
annotFile="$genomeDir/$annotName.gtf"
outDir="$genomeDir/star_ref_L1_"

mkdir -p $outDir

#log directory

logDir="$projectDir/RNA-seq/scripts/mapping/logs"
mkdir -p $logDir

echo This is the genomeDir:
echo $genomeDir
echo -e
echo This is the genomeFile:
echo $genomeFile
echo -e
echo This is the annotFile:
echo $annotFile
echo -e
echo This is the outDir
echo $outDir
echo -e
echo "This is the logDir:"
echo $logDir

#generate the star reference files:
star_ref_line="STAR --runMode genomeGenerate \
	--sjdbGTFfile $annotFile --genomeDir $genomeDir \
	--genomeFastaFiles $genomeFile --genomeSAindexNbases 5 --runThreadN $numcores --outFileNamePrefix \
	$outDir"

echo -e
echo This is the star_ref_line:
echo $star_ref_line

#submit job with name 'RSEM_count_$sample' to 15 cluster cores:
qsub -N STAR_ref_$genomeName -wd $logDir -b y -j y -R y -pe smp $numcores -V $star_ref_line
