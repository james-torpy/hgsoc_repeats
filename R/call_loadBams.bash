#!/bin/bash

### call_loadBams.bash ###

# Call loadBams.R to load bams and save as RDS files

# define variables:
projectName="hgsoc_repeats"
methodName="method1"
bamExt=".bam"
numcores=3

# define directories:
homeDir="/share/ScratchGeneral/jamtor"
projectDir="$homeDir/projects/$projectName"
resultsDir="$projectDir/RNA-seq/results"
bamPath="$resultsDir/star/$methodName"
scriptDir="$projectDir/RNA-seq/scripts/R"
logDir="$scriptDir/logs"

mkdir -p $logDir

RDir="/home/jamtor/local/lib/r/R-3.2.2/bin/R"
outScript="$scriptDir/loadBams.R"

echo -e
echo "This is the bamPath:"
echo $bamPath

# fetch directory names for files:
i=0
directory_names=`ls $bamPath | grep subset | grep -v ribosome | grep -v bowtell_FT1`
for directory in ${directory_names[@]}; do
	echo -e
	echo The directory used is: $directory
	directoriesTotal[i]=$directory
	let i++
done;

echo -e
echo The total number of directories is: ${#directoriesTotal[@]}

# fetch inFiles:
j=0
while [ $j -lt ${#directoriesTotal[@]} ]; do
	uniqueDir="$bamPath/${directoriesTotal[$j]}"
	uID=`basename $uniqueDir`
	inFile=$uniqueDir/Aligned.sortedByCoord.out.bam

	echo -e
	echo "The file used is:"
	echo $inFile

	# create command to send bam file to loadBams.R for loading:
	submitLine="--vanilla --args $inFile $uID < $outScript"

	echo -e
	echo "The submitLine is:"
	echo $submitLine

	# submit method_reports.R to the cluster for each sample, with associated variables:
	qsub -b y -j y -N load$uID -wd $logDir -pe smp $numcores -V $RDir $submitLine

	j=$(($j+1))
done;





