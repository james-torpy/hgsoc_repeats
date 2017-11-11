#!/bin/bash

### 2.fastq_report_gen.bash ###

module load phuluu/samtools/1.4

# For each sample subset, counts reads in each of pairs of fastqs and
# generates a report to be inputted into method_report_gen.R

echo -e
echo "### 2.fastq_report_gen.bash ###"
echo -e

# fetch input variables and define outDir:
uID=$1
fastqPath=$2
projectDir=$3
bamFile=$4

outDir="$projectDir/reports/exp1/"

echo "The uID is:"
echo $uID
echo -e
echo "The bamFile is:"
echo $bamFile
echo -e
echo "The fastqPath is:"
echo $fastqPath
echo -e
echo "The projectDir is:"
echo $projectDir
echo -e
echo "The outDir is:"
echo $outDir
echo -e

mkdir -p $outDir

echo "Removing $bamFile..."
echo -e
rm $bamFile

# add the number of reads from each fastq file into report:
# fetch input files, counting the number of files:
i=1

fastq=( $(ls $fastqPath/*.fastq.gz | grep $uID | grep -v unpaired) )

for file in ${fastq[@]}; do
	echo -e
	echo "The file used is:"
	echo $file
	echo -e
	filesTotal[i]=$file
	let i++;
done;

echo -e
echo "The total number of files is:"
echo ${#filesTotal[@]}

# define each inFile:
j=1
while [ $j -lt ${#filesTotal[@]} ]; do
	inFile1="${filesTotal[$j]}"
	inFile2="${filesTotal[$(($j+1))]}"
	# count the number of lines of each inFile and put into report:
	count1=`zcat $inFile1 | wc -l`
	echo -e
	echo "Number of lines in $inFile1:"
	echo $count1

	echo -e Number of lines in fastq1:\\t$count1 > $outDir/$uID.report.txt
	
	count2=`zcat $inFile2 | wc -l`
	echo -e
	echo "Number of lines in $inFile2:"
	echo $count2

	echo -e Number of lines in fastq2:\\t$count2 >> $outDir/$uID.report.txt

	j=$(($j+2))
done;
