#!/bin/bash

### fastq_report_gen.bash ###

# For each sample subset, counts reads in each of pairs of fastqs and
# generates a report to be inputted into method_report_gen.R

# define variables:
projectName="hgsoc_repeats"
sampleName="subset"
methodName="method1"
fastqExt=".fq.gz"

# define directories:
homeDir="/share/ScratchGeneral/jamtor"
projectDir="$homeDir/projects/$projectName"
resultsDir="$projectDir/RNA-seq/results"

fqPath="$resultsDir/trimgalore/$sampleName"
outDir="$resultsDir/method_reports"

mkdir -p $outDir

echo -e
echo "This is the fqPath:"
echo $fqPath
echo -e
echo "This is the outDir:"
echo $outDir

# add the number of reads from each fastq file into report:
# fetch input files, counting the number of files:
i=0
#fastqs=( $(ls $fqPath/**/$fastqExt | grep -v unpaired))
fastqs=( $(ls $fqPath/**/*$fastqExt | grep -v unpaired | grep 1 | grep -v grant | grep -v gtx | grep -v PR | grep -v FT[2:3] | grep -v v[2:3]))
for file in ${fastqs[@]}; do
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

# define each inFile and unique IDs:
j=0
while [ $j -lt ${#filesTotal[@]} ]; do
	inFile1="${filesTotal[$j]}"
	inFile2="${filesTotal[$(($j+1))]}"
	uID=`basename $inFile1 | sed 's/R1_//' | sed "s/_val_1\$fastqExt//"`
	# count the number of lines of each inFile and put into report:
	count1=$(wc -l < $inFile1)
	echo -e "Number of reads in fastq1:\t$count1" > "$outDir/$methodName"_$uID\_report.txt
	count2=$(wc -l < $inFile2)
	echo -e "Number of reads in fastq2:\t$count2" >> "$outDir/$methodName"_$uID\_report.txt

	j=$(($j+2))
done;


