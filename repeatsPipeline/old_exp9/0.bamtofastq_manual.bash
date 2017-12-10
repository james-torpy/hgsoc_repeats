#!/bin/bash

### 0.bamtofastq.bash ###

# Takes paired read bam files and converts them to fastq files, one for 1st
# paired reads, one for 2nd paired reads, and two for unpaired reads

module load gi/zlib/1.2.8

echo -e
echo "### 0.bamtofastq.bash ###"
echo -e

date

# read in inputs/outputs:
#cores=$1
#bamFile=$2
#fqFile1=$3
#fqFile2=$4
#fqUnpaired1=$5
#fqUnpaired2=$6

cores=2
bamFile="/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/raw_files/fullsamples/bowtell_primary/FT1.bam"
fqFile1="/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/raw_files/fullsamples/bowtell_primary/fastq/FT1.1.fastq.gz"
fqFile2="/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/raw_files/fullsamples/bowtell_primary/fastq/FT1.2.fastq.gz"
fqUnpaired1="/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/raw_files/fullsamples/bowtell_primary/fastq/FT1.unpaired.1.fastq.gz"
fqUnpaired2="/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/raw_files/fullsamples/bowtell_primary/fastq/FT1.unpaired.2.fastq.gz"

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

q="short"
logDir="/share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/logs/exp9"
uID="FT1"

# define tool path:
f2bPath="/home/jamtor/local/bin/bamtofastq"

# convert bamfiles to fastq file containing both 1st and 2nd paired reads
fastq_line="$f2bPath F=$fqFile1 F2=$fqFile2 O=$fqUnpaired1 O2=$fqUnpaired2 \
filename=$bamFile gz=1"

echo "the fastq line is:"
echo $fastq_line
echo -e

/opt/gridengine/bin/linux-x64/qsub -q $q.q -N fq_$uID -hold_jid cp_$uID -b y -wd $logDir -j y -R y \
	-pe smp $cores -V "$f2bPath F=$fqFile1 F2=$fqFile2 O=$fqUnpaired1 O2=$fqUnpaired2 filename=$bamFile gz=1"

date