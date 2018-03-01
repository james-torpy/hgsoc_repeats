#!/bin/bash

echo -e
echo -e
echo "### 2.salmon.bash ###"
echo -e

date

# fetch input variables:
cores=$1
fqFile1=$2
fqFile2=$3
outFile=$4
annot=$5

# define and create outPath and genomeDir:
outPath=`echo $outFile | sed "s/EXPR.csv//"`
mkdir -p $outPath

# define the genomeDir:
genomeDir="/home/jamtor/genomes/hg38_ercc/salmon/"

# define salmonTE line:
salmon_line="/home/jamtor/local/bin/Salmon-latest_linux_x86_64/bin/salmon \
quant -i $genomeDir/hg38_transcriptome -l A -1 $fqFile1 -2 $fqFile2 -o \
$outPath"

# echo variables:
echo "This is the outPath:"
echo $outPath
echo -e

echo "This is the genomeDir:"
echo $genomeDir
echo -e

echo -e
echo "This is the salmon_line:"
echo $salmon_line

# run command:
/home/jamtor/local/bin/Salmon-latest_linux_x86_64/bin/salmon \
quant -i $genomeDir/hg38_transcriptome -l A -1 $fqFile1 -2 $fqFile2 -o \
$outPath

date