#!/bin/bash

module load gi/boost/1.53.0

numcores=6

# make directory hierachy
projectname="hgsoc_repeats"
samplename="GC"
expName="exp4"

homeDir="/share/ScratchGeneral/jamtor"
projectDir="$homeDir/projects/$projectname"
resultsDir="$projectDir/RNA-seq/results"

# genome reference file
genomeName="hg38_ercc"
genome_refFile="$homeDir/genomes/$genomeName/rsem_ref"

# input/output directories
inType="star"
outType="rsem"

inPath="$resultsDir/$inType/$samplename/$expName"
outPath="$resultsDir/$outType/$expName"

# scripts/log directory
scriptsDir="$projectDir/scripts/mapping"
logDir="$scriptsDir/logs"
mkdir -p $logDir

echo -e
echo "The logDir is:"
echo $logDir

# fetch directory names for files:
i=0
directory_names=`ls $inPath | grep -v genome | grep -v subset`
for directory in ${directory_names[@]};do
	echo The directory used is: $directory;
	echo -e
	directoriesTotal[i]=$directory
	let i++
done;

# set up conditions to perform analysis on bam files from STAR alignment output
j=0
echo The total number of directories is: ${#directoriesTotal[@]}
echo -e

#set up directories specific to each file being analysed
while [ $j -lt ${#directoriesTotal[@]} ]; do

	unique_inDir=$inPath/${directoriesTotal[$j]}
	inFile=$unique_inDir/Aligned.toTranscriptome.out.bam
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
	echo This is the genome_refFile:
	echo $genome_refFile
	echo -e
	echo This is the out_filePrefix:
	echo $out_filePrefix
	echo -e
	
	#perform analysis on bam files using the reference genome
	rsem_line="rsem-plot-transcript-wiggles -p $numcores $uniqueID $out_filePrefix"
	echo This is the rsem_line:
	echo $rsem_line
	echo -e
	
	#submit job with name 'RSEM_count_$samplename' to 10 cluster cores:
	qsub -N RSEMwig_FT1_subset -wd $logDir -b y -j y -R y -pe smp $numcores -V $rsem_line
qsub -N RSEMwig_FT1_subset -wd /share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/exp4/logs -b y -j y -R y -pe smp 6 -V \
"rsem-bam2wig /share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/results/rsem/exp4/FT1_subset/FT1_subset.genes.results FT1_subset.wig"
plot-transcript-wiggles -p 6 FT1_subset_wiggle /share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/results/rsem/exp4/FT1_subset/FT1_subset_tIDs.txt
        j=$(($j+1))

done;


