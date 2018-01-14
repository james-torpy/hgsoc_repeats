#!/bin/bash

echo -e
echo -e
echo "### 3c.novosortCoord.bash ###"
echo -e

source /home/jamtor/.bashrc

module load gi/zlib/1.2.8
module load phuluu/samtools/1.4
module load gi/novosort/precompiled/1.03.08

date

# fetch input variables:
cores=$1
bamFile=$2
outBamFile=$3
mem=$4

echo "This is the bamFile:"
echo $bamFile
echo -e
echo "This is the outBamFile:"
echo $outBamFile
echo -e

# define tempDir:
tempDir=`echo $bamFile | sed "s/Aligned.out.bam/temp/"`
mkdir -p $tempDir

echo "This is the tempDir:"
echo $tempDir
echo -e

# sort bamFile by co-ordinate:
/home/jamtor//local/bin/novocraft/novosort -t $tempDir -n -c 6 -m $mem \
	$bamFile > $outBamFile

date