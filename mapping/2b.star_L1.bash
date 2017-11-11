#!/bin/bash

module load gi/samtools/1.2

numcores=6

# make directory hierachy:
projectName="hgsoc_repeats"
sampleName="fulltest"
genomeName="L1"

homeDir="/share/ScratchGeneral/jamtor"
projectDir="$homeDir/projects/$projectName"
rawDir="$projectDir/RNA-seq/raw_files"
genomeDir="/home/jamtor/genomes/$genomeName/star_ref_$genomeName"
resultsDir="$projectDir/RNA-seq/results"

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

#input/output:
inPath="$rawDir/$sampleName"
outPath="$resultsDir/$outType/L1/$sampleName"
	
echo This is the inPath:
echo $inPath
echo -e

#fetch file names of all projects and put into an array:
i=0
files=( $(ls $inPath/*L1s$inExt | grep -v unpaired) )
#files=( bowtell_FT1_R1_subset.fq.gz bowtell_FT1_R2_subset.fq.gz kur_v1_R1_subset.fq.gz kur_v1_R2_subset.fq.gz)
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
	uniqueID=`basename $inFile1 | sed s/_R1// | sed s/_L1s$inExt//`
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
     	star_line="STAR --runMode alignReads \
      	--readFilesCommand zcat \
      	--genomeDir $genomeDir \
   		--outFilterType BySJout \
     	--outSAMattributes NH HI AS NM MD\
     	--outFilterMultimapNmax 999 \
      --outMultimapperOrder Random
      --runRNGseed 666
      --outSAMmultNmax 1
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
		--quantMode TranscriptomeSAM \
     	--outFilterMatchNmin 76 \
		--outSAMtype BAM SortedByCoordinate \
		--outWigType wiggle \
		--outWigNorm None \
		--limitBAMsortRAM 80000000000"
       echo $star_line

# run job on current qsub:
$star_line
#submit jobs to the cluster, creating a log in $logDir which includes reported errors:                	
  #qsub -N STAR_$uniqueID -hold_jid TRIMGALORE_$UniqueID -b y -wd $logDir -j y -R y -pe smp $numcores -V $star_line
	j=$(($j+2))

done;
