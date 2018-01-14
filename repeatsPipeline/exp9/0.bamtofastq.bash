#!/bin/bash

### 0.bamtofastq.bash ###

# Takes paired read bam files and converts them to fastq files, one for 1st
# paired reads, one for 2nd paired reads, and two for unpaired reads

source /home/jamtor/.bashrc

module load gi/zlib/1.2.8
module load phuluu/samtools/1.4

echo -e
echo "### 1.bamtofastq.bash ###"
echo -e

date

# read in inputs/outputs:
cores=$1
bamFile=$2
fqFile1=$3
fqFile2=$4
fqUnpaired1=$5
fqUnpaired2=$6

echo "This is the bamFile:"
echo $bamFile
echo -e
echo "This is the fqFile1:"
echo $fqFile1
echo -e
echo "This is the fqFile2:"
echo $fqFile2
echo -e
echo "This is the fqUnpaired1:"
echo $fqUnpaired1
echo -e
echo "This is the fqUnpaired2:"
echo $fqUnpaired2
echo -e

# define tool path:
f2bPath="/home/jamtor/local/bin/bamtofastq"

# convert bamfiles to fastq file containing both 1st and 2nd paired reads
fastq_line="$f2bPath F=fqFile1 F2=fqFile2 O=fqUnpaired1 O2=fqUnpaired2 \
filename=$bamFile gz=1"

echo "the fastq line is:"
echo $fastq_line
echo -e

$f2bPath F=$fqFile1 F2=$fqFile2 O=$fqUnpaired1 O2=$fqUnpaired2 filename=$bamFile gz=1

date