#!/bin/bash

echo -e
echo "### 2.rsem_ercc.bash ###"
echo -e
echo -e

source /home/jamtor/.bashrc

module load /share/ClusterShare/Modules/modulefiles/contrib/gi/zlib/1.2.8
module load /share/ClusterShare/Modules/modulefiles/contrib/phuluu/samtools/1.4
module load /share/ClusterShare/Modules/modulefiles/contrib/gi/boost/1.53.0

date

cores=$1
bamFile=$2
outFile=$3

# define genome reference file:
genomeName="hg38_ercc"
genome_refPrefix="/home/jamtor/genomes/$genomeName/rsem_ref"

# define outPrefix:
uID=`echo $bamFile | sed "s/^.*mappers//" | sed "s/Aligned.toTranscriptome.out.bam//" | sed "s/\/\///" | sed "s/\///"`

# define outDir:
outDir=`echo $outFile | sed "s/.genes.results//"`

mkdir -p $outDir

echo -e
echo This is the bamFile:
echo $bamFile
echo -e
echo This is the outDir:
echo $outDir
echo -e
echo This is the genome_refPrefix:
echo $genome_refPrefix
echo -e

#perform analysis on bam files using the reference genome
rsem_line="/home/jamtor/local/lib/RSEM-1.3.0/rsem-calculate-expression -p $cores --paired-end --bam $bamFile $genome_refPrefix $outDir"
echo This is the rsem_line:
echo $rsem_line
echo -e

/home/jamtor/local/lib/RSEM-1.3.0/rsem-calculate-expression -p $cores --paired-end --bam $bamFile $genome_refPrefix $outDir
