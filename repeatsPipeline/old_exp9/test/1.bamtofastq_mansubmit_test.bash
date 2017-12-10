#!/bin/bash

source /home/jamtor/.bashrc
echo $SGE_ROOT

### 1.bam_to_fastq.bash ###

# Takes paired read bam files and converts them to fastq files, one for 1st
# paired reads, one for 2nd paired reads
# Calls rest of repeatsPipeline

module load gi/zlib/1.2.8

echo -e
echo "### 1.bam_to_fastq.bash ###"
echo -e

module load phuluu/samtools/1.4

date

# specify cores needed:
fqcores=2
ribocores=10
gccores=6
sortcores=7
countcores=12
mem=45G
subType="long"
expName="exp9"
draft=""

# make directory hierachy:
projectname="hgsoc_repeats"
sampleType="fullsamples/bowtell_primary"
RDir="/home/jamtor/local/lib/r/R-3.2.2/bin/R"

homeDir="/share/ScratchGeneral/jamtor/"
projectDir="$homeDir/projects/$projectname/RNA-seq/"
rawDir="$projectDir/raw_files/$sampleType"

echo -e
echo This is the rawDir:
echo $rawDir

# input/output types:
inType="bam"
outType="fastq"

# extension of files to be used:
inExt=".bam"

# scripts/logs directories:
scriptsPath="/share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/$expName"
logDir="/share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/$expName/logs"
mkdir -p $logDir

echo -e
echo "The logDir is:"
echo $logDir
echo -e
echo "The scriptsPath is:"
echo $scriptsPath

# set in/outPaths:
inPath="$rawDir/"
outPath="$rawDir/fastq/test/"

mkdir -p $outPath

# fetch input bamfiles:
files=( $(cat "$rawDir/files.txt") )
#files=( $rawDir/FT4.bam $rawDir/FT5.bam $rawDir/FT6.bam \
#	$rawDir/prPT4.bam $rawDir/prPT5.bam $rawDir/prPT6.bam \
#	$rawDir/rfPT4.bam $rawDir/rfPT5.bam $rawDir/rfPT6.bam  )
#files=$rawDir/bowtell_FT3_subset.bam
# remove files.txt
rm "$rawDir/files.txt"

# convert files:
for inFile in ${files[@]}; do

	# create unique IDs for each bamfile:
	uID=`basename $inFile | sed 's/.bam//g'`
	# copy the bam file across if it's not in the rawDir:
#	if [ ! -e "$inPath/$uID.bam" ]; then
#		/opt/gridengine/bin/linux-x64/qsub -q $subType.q -N cp_$uID -b y -wd $logDir -j y -R y -pe smp 1 -V cp "$origRawDir/$uID.bam" $inPath
#	fi;
	fastqFile="$outPath/$uID.fastq"

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

	# convert bamfiles to fastq file containing both 1st and 2nd paired reads
	fastq_line="/home/jamtor/local/bin/bamtofastq F=$outPath/$uID.1.fastq.gz \
	F2=$outPath/$uID.2.fastq.gz O=$outPath/$uID.unmatched.1.fastq.gz \
	O2=$outPath/$uID.unmatched.2.fastq.gz filename=$inPath/$uID\.bam gz=1"

	echo -e
	echo "The fastq_line is: "
	echo $fastq_line

	# submit fastq job to the cluster:
	/opt/gridengine/bin/linux-x64/qsub -q $subType.q -N fq_$uID -hold_jid cp_$uID -b y -wd $logDir -j y -R y \
	-pe smp $fqcores -V bamtofastq \
	F=$outPath/$uID.1.fastq.gz F2=$outPath/$uID.2.fastq.gz \
	O=$outPath/$uID.unpaired.1.fastq.gz O2=$outPath/$uID.unpaired.2.fastq.gz \
	filename=$inPath/$uID.bam gz=1

#	# pass variables on to and submit rest of scripts:
#	# 2.fastq_report_gen.bash: 
#	/opt/gridengine/bin/linux-x64/qsub -q $subType.q -N j2_$uID -hold_jid fq_$uID -b y -wd $logDir -j y -R y \
#	-pe smp 1 -V $scriptsPath/2.fastq_report_gen.bash $uID $outPath \
#	$projectDir $inFile $draft
#
#	# 3a.starRibo.bash:
#	/opt/gridengine/bin/linux-x64/qsub -q $subType.q -N j3a_$uID -hold_jid fq_$uID -b y -wd \
#	$logDir -j y -R y -pe smp $ribocores -V $scriptsPath/3a.starRibo.bash \
#	$uID $ribocores $draft
#
#	#3b.starGCchim.bash:
#	/opt/gridengine/bin/linux-x64/qsub -q $subType.q -N j3b_$uID -hold_jid fq_$uID -b y -wd \
#	$logDir -j y -R y -pe smp $gccores -V $scriptsPath/3b.starGCchim.bash \
#	$uID $gccores $draft $outPath
#	
#	#3c sort bam:
#	/opt/gridengine/bin/linux-x64/qsub -q $subType.q -N j3c_$uID -hold_jid j3b_$uID -b y -wd $logDir -j y -R y \
#	-pe smp $sortcores -V $scriptsPath/3c.novosort.bash $uID $draft
#	
#	# 4.bam_report_gen.bash:
#	/opt/gridengine/bin/linux-x64/qsub -q $subType.q -N j4_$uID -hold_jid j3b_$uID,j3c_$uID -b y -wd $logDir -j y -R y \
#	-pe smp 1 \
#	-V $scriptsPath/4.bam_report_gen.bash $uID $outPath $draft
	
#	# 5a.htseq.bash:
#	/opt/gridengine/bin/linux-x64/qsub -q $subType.q -N j5a_$uID -hold_jid j3c_$uID -b y -wd $logDir -j y -R y \
#	-pe smp $countcores -l mem_free=$mem,h_vmem=$mem -v inFile=$inFile,uID=$uID,projectDir=$projectDir \
#	-V $scriptsPath/5a.htseq.bash $uID $draft
# 
#	# 5c.htseq-rep_custom3.bash:
#	/opt/gridengine/bin/linux-x64/qsub -q $subType.q -N j5c_$uID -hold_jid j3c_$uID -b y -wd $logDir -j y -R y \
#	-pe smp $countcores -l mem_free=$mem,h_vmem=$mem \
#	-V $scriptsPath/5c.htseq-rep_custom3.bash $uID $draft
#
# 	# 5d.htseq-rep_all.bash:
#	/opt/gridengine/bin/linux-x64/qsub -q $subType.q -N j5d_$uID -hold_jid j3c_$uID -b y -wd $logDir -j y -R y \
#	-pe smp $countcores -l mem_free=$mem,h_vmem=$mem \
#	-V $scriptsPath/5d.htseq-rep_all.bash $uID $draft

	# 6.rsem_ercc.bash:
#	/opt/gridengine/bin/linux-x64/qsub -q $subType.q -N j6_$uID -hold_jid j3c_$uID -b y -wd $logDir -j y -R y \
#	-pe smp $gccores \
#	-V $scriptsPath/6.rsem_ercc.bash $uID $gccores

	# 7.remove_files.bash:
#	/opt/gridengine/bin/linux-x64/qsub -q $subType.q -N j7_$uID -hold_jid j3a_$uID,j3b_$uID,j3c_$uID,j5a_$uID,j5c_$uID,j5d_$uID,j6_$uID -b y -wd $logDir -j y -R y \
#	-pe smp 1 \
#	-V $scriptsPath/7.remove_files.bash $uID $gccores $draft

done;
