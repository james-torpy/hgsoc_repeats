#!/usr/bin/env python

### 1.bam_to_fastq.bash ###

# Takes paired read bam files and converts them to fastq files, one for 1st
# paired reads, one for 2nd paired reads
# Calls rest of repeatsPipeline

import os
import fileinput

print '\n### 1.bam_to_fastq.bash ###\n'

smallcores=1
largecores=12
countcores=15
subType='long'

projectname='hgsoc_repeats'
sampleType=''
RDir='/home/jamtor/local/lib/r/R-3.2.2/bin/R'

homeDir='/share/ClusterShare/thingamajigs/jamtor/'
projectDir=homeDir + '/projects/' + projectname + '/RNA-seq/python/'
rawDir=projectDir + '/raw_files/'
origRawDir='/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/raw_files/fullsamples/bowtell_primary/'

print '\n### 1.bam_to_fastq.bash ###\n'

print 'This is the rawDir:'
print rawDir

# input/output types:
inType='bam'
outType='fastq'

# extension of files to be used:
inExt='.bam'

# scripts/logs directories:
scriptsPath='/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/'
logDir=scriptsPath + '/logs/bam_to_fastq/'

if not os.path.exists(logDir):
	os.makedirs(logDir)

print '\nThe logDir is:\n'
print logDir

# set in/outPaths:
inPath=rawDir + '/' + sampleType + '/'
outPath=rawDir + '/fullsamples/bowtell_primary/fastq'

if not os.path.exists(outPath):
	os.makedirs(outPath)


#f= open("guru99.txt","w+")
#f.write("This is line %d\r\n" % (i+1))
#f.close()


# fetch input bamfiles:
fls = fileinput.input(rawDir + '/files.txt')

for line in fls:
	print 'The file being processed is: '
	print line,

fls.close()


files=[]

# fetch input bamfiles:
files=( $(cat '$rawDir/fullsamples/bowtell_primary/files.txt') )

# remove files.txt
rm '$rawDir/fullsamples/bowtell_primary/files.txt'

# convert files:
for inFile in ${files[@]}; do

	# create unique IDs for each bamfile:
	uID=`basename $inFile | sed 's/.bam//g'`
	# copy the bam file across if it's not in the rawDir:
#	if [ ! -e '$inPath/$uID.bam' ]; then
#		qsub -q $subType.q -N cp_$uID -b y -wd $logDir -j y -R y -P \
#	GenomeInformatics -pe smp $smallcores -V cp '$origRawDir/$uID.bam' $inPath
#	fi;

	fastqFile='$outPath/$uID.fastq'

	echo -e
	echo 'The inFile is:'
	echo $inFile
	echo -e
	echo 'The first outFile is:'
	echo $outPath/$uID.1.fastq
	echo -e
	echo 'The second outFile is:'
	echo $outPath/$uID.2.fastq
	echo -e

	fastqcores=$(($largecores-1))
	# convert bamfiles to fastq file containing both 1st and 2nd paired reads
	fastq_line='java -XX:ParallelGCThreads=$fastqcores -jar \
	/home/jamtor/local/lib/picard-2.7.1/picard/build/libs/picard.jar SamToFastq \
	I=$inPath/$uID\.bam F=$outPath/$uID.1.fastq F2=$outPath/$uID.2.fastq 
	FU=$outPath/$uID.unpaired.fastq VALIDATION_STRINGENCY=LENIENT'
	# gzip resulting fastqs:
	gzip_line1='gzip $outPath/$uID.1.fastq'
	gzip_line2='gzip $outPath/$uID.2.fastq'
	gzip_line3='gzip $outPath/$uID.unpaired.fastq'

	echo 'The fastq_line is: '
	echo $fastq_line
	echo -e
	echo 'The gzip lines are:'
	echo $gzip_line1
	echo $gzip_line2
	echo $gzip_line3

	# submit fastq job to the cluster:
#	qsub -q $subType.q -N fastq_$uID -hold_jid cp_$uID -b y -wd $logDir -j y -R y -P GenomeInformatics \
#	-pe smp $largecores -V java -jar \
#	/home/jamtor/local/lib/picard-2.7.1/picard/build/libs/picard.jar \
#	SamToFastq I=$inPath/$uID\.bam \
#	F=$outPath/$uID.1.fastq F2=$outPath/$uID.2.fastq \
#	FU=$outPath/$uID.unpaired.fastq VALIDATION_STRINGENCY=LENIENT
#	# gzip the fastq files:
#	qsub -q $subType.q -N gz1_$uID -hold_jid fastq_$uID -b y -wd $logDir -j y -R y -P \
#	GenomeInformatics -pe smp $smallcores -V gzip $outPath/$uID.1.fastq
#
#	qsub -q $subType.q -N gz2_$uID -hold_jid fastq_$uID -b y -wd $logDir -j y -R y -P \
#	GenomeInformatics -pe smp $smallcores -V gzip $outPath/$uID.2.fastq
#
#	qsub -q $subType.q -N gz3_$uID -hold_jid gz2_$uID -b y -wd $logDir -j y -R y -P \
#	GenomeInformatics -pe smp $smallcores -V gzip $outPath/$uID.unpaired.fastq

	# pass variables on to and submit rest of scripts:
	# 2.fastq_report_gen.bash: 
	qsub -q $subType.q -N sc2_$uID -hold_jid gz1_$uID,gz2_$uID,gz3_$uID -b y -wd $logDir -j y -R y -P \
	GenomeInformatics -pe smp $smallcores \
	$scriptsPath/repeatsPipeline/2.fastq_report_gen.bash $uID $outPath $projectDir \
	$inFile
	# 3a.starRibo.bash:
	qsub -q $subType.q -N sc3a_$uID -hold_jid sc2_$uID -b y -wd $logDir -j y -R y -P \
	GenomeInformatics -pe smp $largecores \
	$scriptsPath/repeatsPipeline/3a.starRibo.bash $uID $outPath $projectDir \
	$largecores

	# 3b.starGC.bash:
	qsub -q $subType.q -N sc3b_$uID -hold_jid sc2_$uID -b y -wd $logDir -j y -R y -P \
	GenomeInformatics -pe smp $largecores \
	$scriptsPath/repeatsPipeline/3b.starGC.bash $uID $outPath $projectDir $largecores
	

	# 4.bam_report_gen.bash:
	qsub -q $subType.q -N sc4_$uID -hold_jid sc3a_$uID,sc3b_$uID -b y -wd $logDir -j y -R y -P \
	GenomeInformatics -pe smp $smallcores \
	$scriptsPath/repeatsPipeline/4.bam_report_gen.bash $uID $outPath
	
	annot='c'
	# 5a.countOverlaps_c1.R:
	qsub -q $subType.q -N isc5a_$uID -hold_jid sc4_$uID -b y -wd $logDir -j y -R y -P \
	GenomeInformatics -pe smp $countcores -V '$RDir --vanilla --args $inFile $uID \
	$projectDir $annot < $scriptsPath/repeatsPipeline/5a.countOverlaps_c1.R'
	# 5b.countOverlaps_c2.R:
	qsub -q $subType.q -N isc5b_$uID -hold_jid sc4_$uID -b y -wd $logDir -j y -R y -P \
	GenomeInformatics -pe smp $countcores -V '$RDir --vanilla --args $inFile $uID \
	$projectDir $annot < $scriptsPath/repeatsPipeline/5b.countOverlaps_c2.R'
	
	annot='custom3'
	# 5a.countOverlaps_c1.R:
	qsub -q $subType.q -N iisc5a_$uID -hold_jid sc4_$uID -b y -wd $logDir -j y -R y -P \
	GenomeInformatics -pe smp $countcores -V '$RDir --vanilla --args $inFile $uID \
	$projectDir $annot < $scriptsPath/repeatsPipeline/5a.countOverlaps_c1.R'
	# 5b.countOverlaps_c2.R:
	qsub -q $subType.q -N iisc5b_$uID -hold_jid sc4_$uID -b y -wd $logDir -j y -R y -P \
	GenomeInformatics -pe smp $countcores -V '$RDir --vanilla --args $inFile $uID \
	$projectDir $annot < $scriptsPath/repeatsPipeline/5b.countOverlaps_c2.R'
	# 6.compGen.R:
	qsub -q $subType.q -N sc6_$uID -hold_jid isc5a_$uID,isc5b_$uID,iisc5a_$uID,iisc5b_$uID -b y -wd $logDir -j y -R y -P \
	GenomeInformatics -pe smp $countcores -V '$RDir --vanilla --args $uID $projectDir \
	< $scriptsPath/repeatsPipeline/6.compGen.R'

done;
