#!/bin/bash

date

#numcores=6

# make directory hierachy
projectname="hgsoc_repeats"
samplename="GC"
expName="exp9"
Type="custom3"
draft=$1
numcores=$2
uID=$3


homeDir="/share/ScratchGeneral/jamtor/"
projectDir="$homeDir/projects/$projectname"
resultsDir="$projectDir/RNA-seq/results"
refDir="/home/jamtor/genomes/repeats"

# genome reference file
genomeName="hg38_ercc"

# input/output directories
inType="star"
outType="htseq"

inPath="$resultsDir/$inType/$samplename/$expName/"
outPath="$resultsDir/$outType/$expName/"

# load in variable from calling script:
#uID="prPT9_subset"

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
outDir=$outPath/$uID/$draft/
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
gff="$refDir/custom3rep.final.gff"
gffID=`basename $gff | sed 's/.gff//' | sed 's/.final//'`

echo This is the gffFile:
echo $gff
echo -e

#perform analysis on bam files using the reference genome

htseq_line="python -m HTSeq.scripts.count -f bam -i ID -t exon --stranded=no -o $outDir/$uID.align_out.sam $inFile $gff >> $outDir/$uID.$Type.htseq.txt"
echo This is the htseq_line:
echo $htseq_line
echo -e

python -m HTSeq.scripts.count -f bam -i ID -t exon --stranded=no -o $outDir/$uID.align_out.sam $inFile $gff >> $outDir/$uID.$Type.htseq.txt

#/opt/gridengine/bin/linux-x64/qsub -q long.q -N $draft.j5c_$uID -hold_jid 3c_$uID -b y -wd \
#/share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/logs/exp9 -j y -R y -P \
#DSGClinicalGenomics -pe smp $numcores -V "python -m HTSeq.scripts.count -f bam -i ID -t exon --stranded=no $inFile $gff >> $outDir/$uID.$Type.htseq.txt"

date
