#!/bin/bash

echo -e
echo -e
echo "### 3c.novosort.bash ###"
echo -e

date

numcores=7

module load gi/novosort/precompiled/1.03.08

# make directory hierachy:
projectname="hgsoc_repeats"
expName="exp9"

homeDir="/share/ScratchGeneral/jamtor/"
projectDir="$homeDir/projects/$projectname/RNA-seq/"

# scripts/logs directories:
scriptsPath="/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/"

uID="rcAF16"

#input/output:
starDir="$projectDir/results/star/GC/$expName/$draft/$uID"
tempDir="$starDir/$uID/temp"
mkdir -p $tempDir

echo "This is the uID:"
echo $uID
echo -e
echo "This is the starDir:"
echo $starDir
echo -e

sortline="/home/jamtor//local/bin/novocraft/novosort -t $tempDir -n -c 6 -m 22G \
	$starDir/Aligned.out.bam > $starDir/Aligned.novosortedByName.out.bam"

echo $sortline

/opt/gridengine/bin/linux-x64/qsub -q long.q -N j3c_$uID -hold_jid j3b_$uID -b y -wd \
/share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/logs/exp9 -j y \
-R y -pe smp $numcores -V $sortline

date