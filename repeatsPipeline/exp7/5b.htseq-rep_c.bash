#!/bin/bash

module load gi/boost/1.53.0

numcores=4

# make directory hierachy
projectname="hgsoc_repeats"
samplename="GC"
expName="exp7"
Type="c"

homeDir="/share/ClusterShare/thingamajigs/jamtor/"
projectDir="$homeDir/projects/$projectname"
resultsDir="$projectDir/RNA-seq/results"
refDir="$projectDir/RNA-seq/refs/$Type"

# genome reference file
genomeName="hg38_ercc"
gtfFile="$refDir/human-89.repeats.htseq.gtf"

# input/output directories
inType="star"
outType="htseq"

inPath="$resultsDir/$inType/$samplename/$expName"
outPath="$resultsDir/$outType/$expName"

# load in variable from calling script:
uID=$1

# scripts/log directory
scriptsDir="$projectDir/RNA-seq/scripts/repeatsPipeline/$expName"
logDir="$scriptsDir/logs"
mkdir -p $logDir

echo -e
echo "The logDir is:"
echo $logDir
echo -e
echo "The inPath is:"
echo $inPath

# fetch directory names for files:
i=0
directory_names=`ls $inPath | grep -v subset | grep -v genome | grep FT[2-3]`
for directory in ${directory_names[@]};do
	echo The directory used is: $directory;
	echo -e
	directoriesTotal[i]=$directory
	let i++
done;

# fetch gff file names:
GFFs=( $(ls $refDir/*.gff | grep final) )

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
	uniqueID=`basename $unique_inDir`

	 mkdir -p $outDir

	echo -e
	echo This is the uniqueID:
	echo $uniqueID
	echo -e
	echo This is the inFile:
	echo $inFile
	echo -e
	echo This is the outDir:
	echo $outDir
	echo -e
	
	# loop through each gff file:
	for gff in ${GFFs[@]}; do
		gffID=`basename $gff | sed 's/.gff//' | sed 's/.final//'`

		#perform analysis on bam files using the reference genome
		htseq_line="htseq-count -f bam -i ID -t exon --stranded=no $inFile $gff >> $outDir/$uniqueID.$Type.$gffID.htseq.txt"
		echo This is the htseq_line:
		echo $htseq_line
		echo -e
	
		#submit job with name 'RSEM_count_$samplename' to 10 cluster cores:
		qsub -q long.q -N $uniqueID.$Type.$gffID.htseq -wd $logDir -b y -j y -R y -P GenomeInformatics -pe smp $numcores -V $htseq_line
	done;
    j=$(($j+1))

done;


