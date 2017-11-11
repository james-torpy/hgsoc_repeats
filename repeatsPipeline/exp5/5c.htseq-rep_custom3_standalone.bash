#!/bin/bash

module load gi/boost/1.53.0

numcores=4

# make directory hierachy
projectname="hgsoc_repeats"
samplename="GC"
expName="exp5"

homeDir="/share/ScratchGeneral/jamtor"
projectDir="$homeDir/projects/$projectname"
resultsDir="$projectDir/RNA-seq/results"
refDir="$projectDir/RNA-seq/refs/$Type"

# input/output directories
inType="star"
outType="htseq"

inPath="$resultsDir/$inType/$samplename/$expName"
outPath="$resultsDir/$outType/$expName"

# scripts/log directory
scriptsDir="$projectDir/RNA-seq/scripts/repeatsPipeline/$expName"
logDir="$scriptsDir/logs"
mkdir -p $logDir

echo -e
echo "The logDir is:"
echo $logDir

# fetch directory names for files:
i=0
directory_names=`ls $inPath | grep -v subset | grep -v genome | grep prPT10`
for directory in ${directory_names[@]};do
	echo The directory used is: $directory;
	echo -e
	directoriesTotal[i]=$directory
	let i++
done;


# set up conditions to perform analysis on bam files from STAR alignment output:
j=0
echo The total number of directories is: ${#directoriesTotal[@]}
echo -e

#set up directories specific to each file being analysed
while [ $j -lt ${#directoriesTotal[@]} ]; do

	unique_inDir=$inPath/${directoriesTotal[$j]}
	inFile=$unique_inDir/Aligned.sortedByName.out.bam
	outDir=$outPath/${directoriesTotal[$j]}
	out_filePrefix=$outDir/${directoriesTotal[$j]}
	uID=`basename $unique_inDir`

	mkdir -p $outDir

	echo -e
	echo This is the uID:
	echo $uID
	echo -e
	echo This is the inFile:
	echo $inFile
	echo -e
	echo This is the outDir:
	echo $outDir
	echo -e
	# define gff file name:
	gff="$refDir/custom3/custom3rep.final.gff"
	gffID=`basename $gff | sed 's/.gff//' | sed 's/.final//'`

	#perform analysis on bam files using the reference genome
	htseq_line="htseq-count -f bam -i ID -t exon --stranded=no -o $outDir/$Type.$gffID.out.sam $inFile $gff >> $outDir/$uID.$gffID.htseq.txt"
	echo This is the htseq_line:
	echo $htseq_line
	echo -e
	
	qsub -q short.q -N cu.$uID.htseq -wd $logDir -b y -j y -R y -pe smp $numcores -V $htseq_line
	
	j=$(($j+1))

done;
