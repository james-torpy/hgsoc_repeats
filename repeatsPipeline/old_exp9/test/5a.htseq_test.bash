#!/bin/bash

module load gi/boost/1.53.0

numcores=12

# make directory hierachy
projectname="hgsoc_repeats"
samplename="GC"
expName="exp9"
Type="gc"

date

homeDir="/share/ScratchGeneral/jamtor/"
projectDir="$homeDir/projects/$projectname"
resultsDir="$projectDir/RNA-seq/results"

# genome reference file
genomeName="hg38_ercc"
gffFile="/home/jamtor/genomes/$genomeName/gencode_v24_hg38_annotation.gff"

# input/output directories
inType="star"
outType="htseq"

inPath="$resultsDir/$inType/$samplename/$expName/"
outPath="$resultsDir/$outType/$expName/"

# load in variable from calling script:
uID="prPT9"

# scripts/log directory
scriptsDir="$projectDir/RNA-seq/scripts/repeatsPipeline/$expName"
logDir="$projectDir/logs"
mkdir -p $logDir

echo -e
echo "The logDir is:"
echo $logDir

inFile=$inPath/$uID/Aligned.novosortedByName.out.bam
outDir=$outPath/$uID
out_filePrefix=$outDir/$uID

mkdir -p $outDir

echo -e
echo This is the uID:
echo $uID
echo -e
echo This is the inFile:
echo $inFile
echo -e
echo This is the gffFile:
echo $gffFile
echo -e
echo This is the outDir:
echo $outDir
echo -e
	
#perform analysis on bam files using the reference genome
# -o $outDir/$uID.$Type.out.sam
htseq_line="htseq-count -f bam -i ID -t exon $inFile $gffFile >> $outDir/$uID.$Type.htseq.txt"
echo This is the htseq_line:
echo $htseq_line
echo -e

/opt/gridengine/bin/linux-x64/qsub -q short.q -N j5a_$uID -hold_jid 3c_$uID -b y -wd \
/share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/logs/exp9 -j y -R y -pe smp $numcores -V \
$htseq_line


date