#!/bin/bash

module load gi/boost/1.53.0

numcores=6

# make directory hierachy
projectname="hgsoc_repeats"
samplename="GC"
expName="exp9"
Type="all"
draft="intersection_test"

date

homeDir="/share/ScratchGeneral/jamtor/"
projectDir="$homeDir/projects/$projectname"
resultsDir="$projectDir/RNA-seq/results"
refDir="/home/jamtor/genomes/repeats/"

# input/output directories
inType="star"
outType="htseq"

inPath="$resultsDir/$inType/$samplename/$expName/"
outPath="$resultsDir/$outType/$expName/$draft/"

# load in variable from calling script:
uID="prPT9"

# scripts/log directory
scriptsDir="$projectDir/RNA-seq/scripts/repeatsPipeline/$expName"
logDir="$projectDir/logs"
mkdir -p $logDir

echo -e
echo "The logDir is:"
echo $logDir
echo -e
echo "The inPath is:"
echo $inPath

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
echo This is the outDir:
echo $outDir
echo -e

# define gff file name:
gff="$refDir/all_rep.final.gff"
gffID=`basename $gff | sed 's/.gff//' | sed 's/.final//'`

echo This is the gffFile:
echo $gff
echo -e

#perform analysis on bam files using the reference genome
htseq_line="htseq-count -f bam -i ID -t exon --stranded=no -o $outDir/$uID.align_out.sam  -m intersection-strict $inFile $gff >> $outDir/$uID.$Type.htseq.txt"
echo This is the htseq_line:
echo $htseq_line
echo -e

/opt/gridengine/bin/linux-x64/qsub -q short.q -N intj5d_$uID -hold_jid j3c_$uID -b y -wd \
/share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/logs/exp9 -j y -R y -pe smp $numcores -V \
"htseq-count -f bam -i ID -t exon --stranded=no -o $outDir/$uID.align_out.sam  -m intersection-strict $inFile $gff >> $outDir/$uID.$Type.htseq.txt"


date
