#!/bin/bash

### 1.bam_to_fastq.bash ###

# Takes paired read bam files and converts them to fastq files, one for 1st
# paired reads, one for 2nd paired reads
# Calls rest of repeatsPipeline

echo -e
echo "### 1.bam_to_fastq.bash ###"
echo -e

module load phuluu/samtools/1.4

# specify cores needed:
smallcores=1
#largecores=12
#countcores=12
subType="long"
expName="exp7"

# make directory hierachy:
projectname="hgsoc_repeats"
sampleType="fullsamples/bowtell_primary"
RDir="/home/jamtor/local/lib/r/R-3.2.2/bin/R"

homeDir="/share/ClusterShare/thingamajigs/jamtor/"
projectDir="$homeDir/projects/$projectname/RNA-seq/"
rawDir="$projectDir/raw_files/$sampleType"
origRawDir="/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/raw_files/$sampleType/"

echo -e
echo This is the rawDir:
echo $rawDir

# input/output types:
inType="bam"
outType="fastq"

# extension of files to be used:
inExt=".bam"

# scripts/logs directories:
scriptsPath="/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline"
logDir="$scriptsPath/$expName/logs"
mkdir -p $logDir

echo -e
echo "The logDir is:"
echo $logDir
echo -e
echo "The scriptsPath is:"
echo $scriptsPath

# set in/outPaths:
inPath="$rawDir/"
outPath="$rawDir/fastq/$expName"

mkdir -p $outPath

# fetch input bamfiles:
files=( $(cat "$rawDir/files.txt") )
#files=( $rawDir/FT4.bam $rawDir/FT5.bam $rawDir/FT6.bam \
#	$rawDir/prPT4.bam $rawDir/prPT5.bam $rawDir/prPT6.bam \
#	$rawDir/rfPT4.bam $rawDir/rfPT5.bam $rawDir/rfPT6.bam  )
#files=$rawDir/bowtell_FT3_subset.bam
# remove files.txt
rm "$rawDir/files.txt"

echo ${files[@]}

# convert files:
for inFile in ${files[@]}; do

	# create unique IDs for each bamfile:
	uID=`basename $inFile | sed 's/.bam//g'`
	# copy the bam file across if it's not in the rawDir:
#	if [ ! -e "$inPath/$uID.bam" ]; then
#		qsub -q $subType.q -N cp_$uID -b y -wd $logDir -j y -R y -P DSGClinicalGenomics -pe smp $smallcores -V cp "$origRawDir/$uID.bam" $inPath
#	fi;
	fastqFile="$outPath/$expName/$uID.fastq"

	echo -e
	echo "The inFile is:"
	echo $inFile
	echo -e
	echo "The first outFile is:"
	echo $outPath/$uID.1.fastq
	echo -e
	echo "The second outFile is:"
	echo $outPath/$uID.2.fastq
	echo -e

	fastqcores=$(($largecores-1))
	# convert bamfiles to fastq file containing both 1st and 2nd paired reads
	fastq_line="java -XX:ParallelGCThreads=$fastqcores -jar \
	/home/jamtor/local/lib/picard-2.7.1/picard/build/libs/picard.jar SamToFastq \
	I=$inPath/$uID\.bam F=$outPath/$uID.1.fastq F2=$outPath/$uID.2.fastq 
	FU=$outPath/$uID.unpaired.fastq VALIDATION_STRINGENCY=LENIENT"
	# gzip resulting fastqs:
	gzip_line1="gzip $outPath/$uID.1.fastq"
	gzip_line2="gzip $outPath/$uID.2.fastq"
	gzip_line3="gzip $outPath/$uID.unpaired.fastq"

	echo "The fastq_line is: "
	echo $fastq_line
	echo -e
	echo "The gzip lines are:"
	echo $gzip_line1
	echo $gzip_line2
	echo $gzip_line3

	# submit fastq job to the cluster:
	qsub -q $subType.q -N fq_$uID -hold_jid cp_$uID -b y -wd $logDir -j y -R y \
	-P DSGClinicalGenomics -pe smp $largecores -V java -XX:ParallelGCThreads=$fastqcores -jar \
	/home/jamtor/local/lib/picard-2.7.1/picard/build/libs/picard.jar SamToFastq \
	I=$inPath/$uID\.bam F=$outPath/$uID.1.fastq F2=$outPath/$uID.2.fastq \
	FU=$outPath/$uID.unpaired.fastq VALIDATION_STRINGENCY=LENIENT

	# gzip the fastq files:
	qsub -q $subType.q -N gz1_$uID -hold_jid fq_$uID -b y -wd $logDir -j y -R y \
	-P DSGClinicalGenomics -pe smp $smallcores -V gzip $outPath/$uID.1.fastq

	qsub -q $subType.q -N gz2_$uID -hold_jid gz1_$uID -b y -wd $logDir -j y -R y \
	-P DSGClinicalGenomics -pe smp $smallcores -V gzip $outPath/$uID.2.fastq

	qsub -q $subType.q -N gz3_$uID -hold_jid gz2_$uID -b y -wd $logDir -j y -R y \
	-P DSGClinicalGenomics -pe smp $smallcores -V gzip $outPath/$uID.unpaired.fastq

	# pass variables on to and submit rest of scripts:
	# 2.fastq_report_gen.bash: 
	qsub -q $subType.q -N j2_$uID -hold_jid gz1_$uID,gz2_$uID,gz3_$uID -b y -wd $logDir -j y -R y \
	-P DSGClinicalGenomics -pe smp $smallcores \
	-V $scriptsPath/2.fastq_report_gen.bash $uID $outPath $projectDir $inFile

	# 3a.starRibo.bash:
	qsub -q $subType.q -N j3a_$uID -hold_jid gz1_$uID,gz2_$uID,gz3_$uID -b y -wd $logDir -j y -R y \
	-P DSGClinicalGenomics -pe smp $largecores \
	-V $scriptsPath/$expName/3a.starRibo.bash $uID $largecores

	#3b.starGCchim.bash:
	qsub -q $subType.q -N j3b_$uID -hold_jid gz1_$uID,gz2_$uID,gz3_$uID -b y -wd $logDir -j y -R y \
	-P DSGClinicalGenomics -pe smp $largecores \
	-V $scriptsPath/$expName/3b.starGCchim.bash $uID $largecores
	
	#3c sort bam:
	qsub -q $subType.q -N j3c_$uID -hold_jid j3b_$uID -b y -wd $logDir -j y -R y \
	-P DSGClinicalGenomics -pe smp $largecores -V samtools sort -n $projectDir/results/star/GC/$expName/$uID/Aligned.out.bam -o $projectDir/results/star/GC/$expName/$uID/Aligned.sortedByName.out.bam

	# 4.bam_report_gen.bash:
	qsub -q $subType.q -N j4_$uID -hold_jid j3b_$uID,j3c_$uID -b y -wd $logDir -j y -R y \
	-P DSGClinicalGenomics -pe smp $smallcores \
	-V $scriptsPath/$expName/4.bam_report_gen.bash $uID $outPath
	
	# 5a.htseq.bash:
	qsub -q $subType.q -N j5a_$uID -hold_jid j3c_$uID -b y -wd $logDir -j y -R y \
	-P DSGClinicalGenomics -pe smp $countcores -v inFile=$inFile,uID=$uID,projectDir=$projectDir \
	-V $scriptsPath/$expName/5a.htseq.bash $uID
 
	# 5c.htseq-rep_custom3.bash:
	qsub -q $subType.q -N j5c_$uID -hold_jid j3c_$uID -b y -wd $logDir -j y -R y \
	-P DSGClinicalGenomics -pe smp $countcores \
	-V $scriptsPath/$expName/5c.htseq-rep_custom3.bash $uID

 	# 5d.htseq-rep_all.bash:
	qsub -q $subType.q -N j5d_$uID -hold_jid j3c_$uID -b y -wd $logDir -j y -R y \
	-P DSGClinicalGenomics -pe smp $countcores \
	-V $scriptsPath/$expName/5d.htseq-rep_all.bash $uID
#
#	# 6.rsem_ercc.bash:
#	qsub -q $subType.q -N j6_$uID -hold_jid j3c_$uID -b y -wd $logDir -j y -R y \
#	-P DSGClinicalGenomics -pe smp $largecores \
#	-V $scriptsPath/$expName/6.rsem_ercc.bash $uID $largecores
#
#	# 7.remove_files.bash:
#	qsub -q $subType.q -N j7_$uID -hold_jid j5a_$uID,j5c_$uID,j5d_$uID,j6_$uID -b y -wd $logDir -j y -R y \
#	-P DSGClinicalGenomics -pe smp 1 \
#	-V $scriptsPath/$expName/7.remove_files.bash $uID $largecores

done;