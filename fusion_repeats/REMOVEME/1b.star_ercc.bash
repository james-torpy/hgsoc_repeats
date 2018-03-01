#!/bin/bash

module load gi/samtools/1.2

numcores=6

# make directory hierachy:
projectname="hgsoc_repeats"

homeDir="/share/ScratchGeneral/jamtor"
projectDir="$homeDir/projects/$projectname"
resultsDir="$projectDir/RNA-seq/results"
genomeDir="/home/jamtor/genomes/hg38_ercc/starRef"

echo -e
echo "This is the genomeDir:"
echo $genomeDir
echo -e

# input/output types:
inExt="fastq.gz"
outType="star"

#scripts/logs directory
scriptsPath="$projectDir/RNA-seq/scripts/fusion_repeats"
logDir="/share/ClusterShare/thingamajigs/jamtor/projects/$projectname/RNA-seq/logs/fusion_repeats"
mkdir -p $logDir

#input/output:
inPath="$projectDir/RNA-seq/raw_files//fullsamples/bowtell_primary/fastq"
outPath="$resultsDir/$outType/fusion_repeats/"
	
echo This is the inPath:
echo $inPath
echo -e
echo This is the logDir:
echo $logDir
echo e

#fetch file names of all projects and put into an array:
i=0
#files=( $(ls $inPath/*$inExt | grep -v unpaired) )
files=( $inPath/prPT4_1.fastq.gz $inPath/prPT4_2.fastq.gz)
for file in ${files[@]}; do
	echo The file used is: $file
	echo -e
	filesTotal[i]=$file
	let i++;
done;

#fetch the inFiles and create an outDir based on their uniqueID:
j=0
echo Total files = ${#filesTotal[@]}
echo -e
while [ $j -lt ${#filesTotal[@]} ]; do
	inFile1="${filesTotal[$j]}"
	inFile2="${filesTotal[$(($j+1))]}"
	uniqueID=`basename $inFile1 | sed s/\\_1.$inExt//`
	outDir=$outPath/$uniqueID/
		
	mkdir -p $outDir
	echo -e
	echo This is the uniqueID:
	echo $uniqueID
	echo -e
	echo This is the outDir:
	echo $outDir
	echo -e

#align reads of input files with STAR, output into .bam files:
  starJobName="star."$uniqueID
       bamJobName="bam."$uniqueID
       sortJobName="sort."$uniqueID
       indexJobName="index."$uniqueID
       indexStatsJobName="indexstats."$uniqueID
       outSam=$outDir"Aligned.out.sam"
       outBam=$outDir"$uniqueID.bam"
       outSortedBam=$outDir"$uniqueID.sorted.bam"
      star_line="STAR \
        --readFilesCommand zcat \
        --genomeDir $genomeDir \
        --twopassMode Basic \
      --outFilterMultimapNmax 999 \
      --outMultimapperOrder Random \
      --runRNGseed 666 \
      --outSAMmultNmax 1 \
      --outReadsUnmapped None \
      --alignIntronMax 100000 \
      --alignMatesGapMax 100000 \
      --chimSegmentMin 12 \                                                                                                    
      --chimJunctionOverhangMin 12 \
      --alignSJDBoverhangMin 10 \
      --chimSegmentReadGapMax parameter 3 \                                                                                    
      --alignSJstitchMismatchNmax 5 -1 5 5 \
      --readFilesIn $inFile1 $inFile2 \
      --outFileNamePrefix $outDir \
      --runThreadN $numcores \
      --limitBAMsortRAM 31532137230
      --outSAMtype BAM SortedByCoordinate"

       echo $star_line

#submit jobs to the cluster, creating a log in $logDir which includes reported errors:                	
  qsub -N STAR_$uniqueID -b y -wd $logDir -j y -R y -pe smp $numcores -V $star_line
	j=$(($j+2))

done;
