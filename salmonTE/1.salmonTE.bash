#!/bin/bash

echo -e
echo -e
echo "### 1.salmonTE.bash ###"
echo -e

date

module load briglo/R/3.4.2

# fetch input variables:
cores=$1
fqFile1=$2
fqFile2=$3
outFile=$4
annot=$5

# define and create outPath:
outPath=`echo $outFile | sed "s/EXPR.csv//"`
mkdir -p $outPath

# define salmonTE line:
salmonTEline="SalmonTE.py quant --reference=$annot --outpath=$outPath --num_threads=$cores \
--exprtype=count $fqFile1 $fqFile2"

# echo variables:
echo "This is the outPath:"
echo $outPath
echo -e

echo -e
echo "This is the salmonTEline:"
echo $salmonTEline

# run command:
SalmonTE.py quant --reference=$annot --outpath=$outPath --num_threads=$cores \
--exprtype=count $fqFile1 $fqFile2

date