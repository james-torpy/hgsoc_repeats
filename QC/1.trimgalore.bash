#!/bin/bash

numcores=1

# define variables:
projectname="hgsoc_repeats"
samplename="fullsamples"
inExt=".fq.gz"
outType="trimgalore"

# define directories:
homeDir="/share/ScratchGeneral/jamtor"
projectDir="$homeDir/projects/$projectname"
resultsDir="$projectDir/RNA-seq/results"
scriptsPath="$projectDir/RNA-seq/scripts/QC"

inDir="$projectDir/RNA-seq/raw_files/$samplename"
outPath="$resultsDir/$outType/$samplename"
logDir="$scriptsPath/logs"

# check directories:
echo -e
echo The inDir is:
echo $inDir
echo -e
echo The outPath is:
echo $outPath
echo -e
echo The logDir is:
echo $logDir

# create directories:
mkdir -p outPath
mkdir -p $logDir


# fetch file names:

i=0
files=( $(ls $inDir/*$inExt | grep bowtell_FT1) )
for file in ${files[@]} ;do
	echo -e
	echo The file used is: $file
	echo -e
	filesTotal[i]=$file;
	let i++;
done;

# print how many files there are:
j=0
echo -e
echo The total number of files is: ${#filesTotal[@]}
echo -e

#fetch directory names for files in pairs and set up directories specific to each file pair being analysed:
while [ $j -lt ${#filesTotal[@]} ]; do
	inFile1=${files[$j]}
	inFile2=`echo $inFile1 | sed s/_R1$inExt/_R2$inExt/g`
	uniqueID=`basename $inFile1 | sed s/_R1$inExt//g`
	outDir=$outPath/$uniqueID/

	# create outDir:
	mkdir -p $outDir
	
	# check variables/directories:
	echo -e
	echo This is the uniqueID:
	echo $uniqueID
	echo -e
	echo This is inFile1:
	echo $inFile1
	echo -e
	echo This is inFile2:
	echo $inFile2
	echo -e
	echo This is the outDir: $outDir

	# check trimgalore_line:
	trimgalore_line="trim_galore $inFile1 $inFile2 --fastqc --paired --retain_unpaired --length 16 -o $outDir"
	echo -e
	echo The trimgalore_line is:
	echo $trimgalore_line
	echo -e

#submit the job to the cluster as a binary file with name trimgalore_$samplename, creating a log in $logDir
#which includes reported errors:
#	qsub -N TRIMGALORE_$uniqueID -wd $logDir -b y -j y -R y -P GenomeInformatics -pe smp $numcores -V $trimgalore_line

	j=$(($j+2))

done;

