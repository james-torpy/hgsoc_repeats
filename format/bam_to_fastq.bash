#!/bin/bash

### bam_to_fastq.bash ###

# Takes paired read bam files and converts them to fastq files, one for 1st
# paired reads, one for 2nd paired reads

module load phuluu/samtools/1.4

# specify cores needed:
numcores=1

# make directory hierachy:
projectname="hgsoc_repeats"
sampleType="fullsamples/bowtell_primary"

homeDir="/share/ScratchGeneral/jamtor/"
projectDir="$homeDir/projects/$projectname/"
rawDir="$projectDir/RNA-seq/raw_files/"

echo -e
echo This is the rawDir:
echo $rawDir

# input/output types:
inType="bam"
outType="fastq"

# extension of files to be used:
inExt=".bam"

# scripts/logs directories:
scriptsPath="$projectDir/RNA-seq/scripts/"
logDir="$scriptsPath/logs"
mkdir -p $logDir

echo -e
echo "The logDir is:"
echo $logDir

# set in/outPaths:
inPath="$rawDir/$sampleType/"
outPath="$rawDir/$sampleType/$outType"

mkdir -p $outPath

# fetch input bamfiles:
files=( $(ls $inPath/*.bam | grep -v sorted | grep FT1) )

echo -e
echo "The bam files used are:"
echo ${files[@]}

# convert files:
for inFile in $files; do

	# create unique IDs for each bamfile:
	uID=`basename $inFile | sed 's/.bam//g'`
	sortedFile="$inPath/$uID.sorted"
	fastqFile="$outPath/$uID.fastq"

	echo -e
	echo "The inFile is:"
	echo $inFile
	echo -e
	echo "The sorted bam file is:"
	echo $sortedFile
	echo -e
	echo "The first outFile is:"
	echo $outPath/$uID.1.fastq
	echo -e
	echo "The second outFile is:"
	echo $outPath/$uID.2.fastq

	# sort bamfiles by read name to put pairs together for bam2fastx:
	#sort_line="samtools sort -n $inFile -o $sortedFile\.bam"
	# index bamfiles:
	#index_line="samtools index $sortedFile\.bam"
	# convert bamfiles to fastq file containing both 1st and 2nd paired reads
	fastq_line="java -jar /home/jamtor/local/lib/picard-2.7.1/picard/build/libs/picard.jar SamToFastq I=$inPath/$uID\.bam F=$outPath/$uID.1.fastq F2=$outPath/$uID.2.fastq FU=$outPath/$uID\_unpaired.fastq VALIDATION_STRINGENCY=LENIENT"
	# gzip resulting fastqs:
	gzip_line1="gzip $outPath/$uID_R1.fastq"
	gzip_line2="gzip $outPath/$uID_R2.fastq"
	gzip_line3="gzip FU=$outPath/$uID\_unpaired.fastq"

	echo -e
	echo "The sort_line is: "
	echo $sort_line
	echo -e
	echo "The index_line is:"
	echo $index_line
	echo -e
	echo "The fastq_line is: "
	echo $fastq_line
	echo -e
	echo "The gzip lines are:"
	echo $gzip_line1
	echo $gzip_line2
	echo $gzip_line3

	# submit all jobs to the cluster:
	#qsub -N srt_$uID -b y -wd $logDir -j y -R y -pe smp $numcores -V $sort_line
	#qsub -N ind_$uID -hold_jid srt_$uID -b y -wd $logDir -j y -R y -pe smp $numcores -V $index_line
#	qsub -N fq_$uID -hold_jid ind_$uID -b y -wd $logDir -j y -R y -pe smp $numcores -V $fastq_line
#	qsub -N gz1_$uID -hold_jid fq_$uID -b y -wd $logDir -j y -R y -pe smp $numcores -V $gzip1_line
#	qsub -N gz2_$uID -hold_jid gz1_$uID -b y -wd $logDir -j y -R y -pe smp $numcores -V $gzip2_line
#	qsub -N gz3_$uID -hold_jid gz2_$uID -b y -wd $logDir -j y -R y -pe smp $numcores -V $gzip3_line

done;
