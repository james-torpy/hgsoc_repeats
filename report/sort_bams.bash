#!/bin/bash

### bam_report_gen.bash ###

# For each sample subset counts alignments in mapped bam,primary alignments in bam,
# primary alignments reported by STAR and multimappers reported by STAR and
# generates a report to be inputted into method_report_gen.R
# Also sorts each bam for R script.

# define variables:
projectName="hgsoc_repeats"
methodName="method1"
bamExt=".bam"

# define directories:
homeDir="/share/ScratchGeneral/jamtor"
projectDir="$homeDir/projects/$projectName"
resultsDir="$projectDir/RNA-seq/results"

bamPath="$resultsDir/star/$methodName"

echo -e
echo "This is the bamPath:"
echo $bamPath

# fetch directory names for files:
i=0
directory_names=`ls $bamPath`
for directory in ${directory_names[@]}; do
	echo -e
	echo The directory used is: $directory;
	echo -e
	directoriesTotal[i]=$directory
	let i++
done;

echo The total number of directories is: ${#directoriesTotal[@]}
echo -e

# fetch inFiles:
j=0
while [ $j -lt ${#directoriesTotal[@]} ]; do
	uniqueDir="$bamPath/${directoriesTotal[$j]}"
	uID=`basename $uniqueDir`
	inFile=$uniqueDir/Aligned.sortedByCoord.out.bam
	outPrefix=`echo $inFile | sed 's/bam/sorted/'`

	echo -e
	echo "The inFile is:"
	echo $inFile
	echo -e
	echo "The outPrefix is:"
	echo $outPrefix

	# sort inFiles:
	samtools sort $inFile $outPrefix

	j=$(($j+1))
done;





