#!/bin/bash

# make directory hierachy:
projectname="hgsoc_repeats"
sampleType=""
expName="exp4"

homeDir="/share/ScratchGeneral/jamtor/"
projectDir="$homeDir/projects/$projectname/RNA-seq/"
rawDir="$projectDir/raw_files/"

# scripts/logs directories:
scriptsPath="/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/$expName"

# set in/outPaths:
fastqPath="$rawDir/$expName/fastq"

genomeName="hg38_ercc"
genomeDir="/share/ScratchGeneral/jamtor/genomes/$genomeName/starRef/"

logDir="$scriptsPath/logs"
mkdir -p $logDir

echo "This is the logDir:"
echo $logDir
echo -e
echo "This is the genomeDir:"
echo $genomeDir
echo -e

# fetch input variables and define outDir:
uIDs=( FT1_subset )
numcores=8

#input/output:
outPath="$projectDir/results/star/GC/$expName"
mkdir -p $outPath  

for uID in ${uIDs[@]}; do
    echo "This is the uID:"
    echo $uID
    echo -e
    echo "This is the fastqPath:"
    echo $fastqPath
    echo -e
    echo "This is the outPath:"
    echo $outPath
    echo -e

    # Fetch inFiles:
    inFile1="$fastqPath/$uID.1.fastq.gz"
    inFile2="$fastqPath/$uID.2.fastq.gz"
    outDir=$outPath/$uID/
  
	mkdir -p $outDir
    echo -e
    echo This is the uID:
    echo $uID
    echo -e
    echo This is the outDir:
    echo $outDir
    echo -e
    echo This is the inFile1:
    echo $inFile1
    echo -e
    echo This is the inFile2:
    echo $inFile2
    echo -e
    
    # align reads of input files with STAR, output into .bam files:
    starJobName="star."$uID
         bamJobName="bam."$uID
         sortJobName="sort."$uID
         indexJobName="index."$uID
         indexStatsJobName="indexstats."$uID
         outSam=$outDir"Aligned.out.sam"
         outBam=$outDir"$uID.bam"
         outSortedBam=$outDir"$uID.sorted.bam"
    star_line="/home/jamtor/local/bin/STAR --runMode alignReads \
            --readFilesCommand zcat \
            --genomeDir $genomeDir \
                    --outFilterType BySJout \
            --outSAMattributes NH HI AS NM MD\
            --outFilterMultimapNmax 999 \
        --outMultimapperOrder Random \
        --runRNGseed 666 \
        --outSAMmultNmax 1 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverReadLmax 0.04 \
            --alignIntronMin 20 \
            --alignIntronMax 1500000 \
            --alignMatesGapMax 1500000 \
            --alignSJoverhangMin 6 \
            --alignSJDBoverhangMin 1 \
            --readFilesIn $inFile1 $inFile2 \
            --outFileNamePrefix $outDir \
            --runThreadN $numcores \
            --outFilterMatchNmin 76 \
        --quantMode TranscriptomeSAM
    --outSAMtype BAM SortedByCoordinate \
    --outWigType wiggle \
    --outWigNorm None \
    --limitBAMsortRAM 80000000000"
    
    echo -e
    echo "The star_line is: "
    echo $star_line
    
    qsub -N STAR_$uID -b y -wd $logDir -j y -R y -pe smp $numcores -V $star_line   

done; 
