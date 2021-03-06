#!/bin/bash

module load gi/samtools/1.2

numcores=6

# make directory hierachy:
projectname="hgsoc_repeats"
sampleName="ribosome"
genomeName="ribosome_12"
methodName="method1"
mapL=0.5

homeDir="/share/ScratchGeneral/jamtor"
projectDir="$homeDir/projects/$projectname"
rawDir="$projectDir/RNA-seq/raw_files"
genomeDir="$homeDir/genomes/$genomeName"
resultsDir="$projectDir/RNA-seq/results"

echo -e
echo "This is the genomeDir:"
echo $genomeDir
echo -e

# input/output types:
inExt=".fq.gz"
#inType="trimgalore"
outType="star"

#scripts/logs directory
scriptsPath="$projectDir/RNA-seq/scripts/mapping"
logDir="$scriptsPath/logs"
mkdir -p $logDir

echo This is the logDir:
echo $logDir
echo -e

#input/output:
inPath="$rawDir/$sampleName"
outPath="$resultsDir/$outType/$methodName/$genomeName"
	
echo This is the inPath:
echo $inPath
echo -e

#fetch file names of all projects and put into an array:
i=0
#files=( $(ls $inPath/**/*$inExt | grep -v unpaired) )
files=( $(ls $inPath/*$inExt | grep -v unpaired) )
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
	#uniqueID=`basename $inFile1 | sed s/_R1$inExt//`
  uniqueID=`basename $inFile1 | sed 's/_R1//' | sed "s/$inExt//"`
  #| sed 's/_val_1$inExt//'`
  outDir="$outPath/basic2/"
		
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
     	--readFilesIn $inFile1 $inFile2 \
     	--outFileNamePrefix $outDir \
     	--runThreadN $numcores \
      --outReadsUnmapped Fastx  \
		--outSAMtype BAM SortedByCoordinate \
		--limitBAMsortRAM 70000000000"
       echo $star_line

#submit jobs to the cluster, creating a log in $logDir which includes reported errors:                	
  qsub -N basic_rSTAR -b y -wd $logDir -j y -R y -pe smp $numcores -V $star_line

	j=$(($j+2))

done;
