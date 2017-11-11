#!/bin/bash

module load gi/samtools/1.2

numcores=6

# make directory hierachy:
projectname="hgsoc_repeats"
samplename="subset"

homeDir="/share/ScratchGeneral/jamtor"
projectDir="$homeDir/projects/$projectname"
resultsDir="$projectDir/RNA-seq/results"
genomeDir="$homeDir/genomes/hg38_ercc/repeats"

echo -e
echo "This is the genomeDir:"
echo $genomeDir
echo -e

# input/output types:
inExt=".fq.gz"
outType="star"

#scripts/logs directory
scriptsPath="$projectDir/RNA-seq/scripts/mapping"
logDir="$scriptsPath/logs"
mkdir -p $logDir

echo This is the logDir:
echo $logDir
echo -e

#input/output:
inPath="$projectDir/RNA-seq/raw_files/$samplename"
outPath="$resultsDir/$outType/L1subset7"
	
echo This is the inPath:
echo $inPath
echo -e

#fetch file names of all projects and put into an array:
i=0
#files=( $(ls $inPath/**/*$inExt | grep -v unpaired | grep -v subset) )
files=( bowtell_FT1_R1_subset.fq.gz bowtell_FT1_R2_subset.fq.gz kur_v1_R1_subset.fq.gz kur_v1_R2_subset.fq.gz )
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
	inFile1="$inPath/${filesTotal[$j]}"
	inFile2="$inPath/${filesTotal[$(($j+1))]}"
	#uniqueID=`basename $inFile1 | sed s/_R1$inExt//`
	uniqueID=`basename $inFile1 | sed s/_R1_subset$inExt//`
  outDir=$outPath/$uniqueID/
		
	mkdir -p $outDir
	echo -e
	echo This is the uniqueID:
	echo $uniqueID
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

#align reads of input files with STAR, output into .bam files:
	starJobName="star."$uniqueID
       bamJobName="bam."$uniqueID
       sortJobName="sort."$uniqueID
       indexJobName="index."$uniqueID
       indexStatsJobName="indexstats."$uniqueID
       outSam=$outDir"Aligned.out.sam"
       outBam=$outDir"$uniqueID.bam"
       outSortedBam=$outDir"$uniqueID.sorted.bam"
     	star_line="STAR --runMode alignReads \
      	--readFilesCommand zcat \
      	--genomeDir $genomeDir \
   		--outFilterType BySJout \
     	--outSAMattributes NH HI AS NM MD\
     	--outFilterMultimapNmax 999 \
      --outMultimapperOrder Random
      --runRNGseed 666 \
     	--outFilterMismatchNmax 999 \
     	--outFilterMismatchNoverReadLmax 0.1 \
     	--alignIntronMin 20 \
     	--alignIntronMax 1500000 \
     	--alignMatesGapMax 1500000 \
     	--alignSJoverhangMin 6 \
     	--alignSJDBoverhangMin 1 \
     	--readFilesIn $inFile1 $inFile2 \
     	--outFileNamePrefix $outDir \
     	--runThreadN $numcores \
     	--outFilterMatchNmin 76 \
		--outSAMtype BAM SortedByCoordinate \
		--limitBAMsortRAM 80000000000"
       echo $star_line

#submit jobs to the cluster, creating a log in $logDir which includes reported errors:                	
  qsub -N STAR_$uniqueID -hold_jid TRIMGALORE_$UniqueID -b y -wd $logDir -j y -R y -pe smp $numcores -V $star_line
	j=$(($j+2))

done;
