#!/bin/bash

echo -e
echo -e
echo "### 1.star_ercc.bash ###"
echo -e

date

module load /share/ClusterShare/Modules/modulefiles/contrib/gi/zlib/1.2.8
module load /share/ClusterShare/Modules/modulefiles/contrib/phuluu/samtools/1.4
module load /share/ClusterShare/Modules/modulefiles/contrib/gi/boost/1.53.0

# fetch input variables:
cores=$1
fqFile1=$2
fqFile2=$3
outFile=$4

# define refDir:
refDir="/home/jamtor/genomes/hg38_ercc/starRef/"

# fetch uID and create outDir:
uID=`basename $fqFile1 | sed "s/\\..*$//"`

outDir=`echo ${outFile%/*}`/
mkdir -p $outDir

echo "This is the fqFile1:"
echo $fqFile1
echo -e
echo "This is the fqFile2:"
echo $fqFile2
echo -e
echo "This is the outFile:"
echo $outFile
echo -e
echo "This is the refDir:"
echo $refDir
echo -e
echo "This is the uID:"
echo $uID
echo -e
echo "This is the outDir:"
echo $outDir

#align reads of input files with STAR, output into .bam files:
		starJobName="star."$uID
        bamJobName="bam."$uID
        sortJobName="sort."$uID

        indexJobName="index."$uID
        indexStatsJobName="indexstats."$uID
        outSam=$outDir"Aligned.out.sam"
        outBam=$outDir"$uID.bam"
        outSortedBam=$outDir"$uID.sorted.bam"

      	star_line="/home/jamtor/local/bin/STAR --runMode alignReads \
       	--readFilesCommand zcat \
       	--genomeDir $refDir \
    	--outFilterType BySJout \
      	--outSAMattributes NH HI AS NM MD\
      	--outFilterMultimapNmax 20 \
      	--outFilterMismatchNmax 999 \
      	--outFilterMismatchNoverReadLmax 0.04 \
      	--alignIntronMin 20 \
      	--alignIntronMax 1500000 \
      	--alignMatesGapMax 1500000 \
      	--alignSJoverhangMin 6 \
      	--alignSJDBoverhangMin 1 \
      	--readFilesIn $fqFile1 $fqFile2 \
      	--outFileNamePrefix $outDir \
      	--runThreadN $cores \
		--quantMode TranscriptomeSAM \
      	--outFilterMatchNmin 76 \
    	--outSAMtype BAM Unsorted"

echo -e
echo "This is the star_line:"
echo $star_line

# run command:
/home/jamtor/local/bin/STAR --runMode alignReads \
       	--readFilesCommand zcat \
       	--genomeDir $refDir \
    	--outFilterType BySJout \
      	--outSAMattributes NH HI AS NM MD\
      	--outFilterMultimapNmax 20 \
      	--outFilterMismatchNmax 999 \
      	--outFilterMismatchNoverReadLmax 0.04 \
      	--alignIntronMin 20 \
      	--alignIntronMax 1500000 \
      	--alignMatesGapMax 1500000 \
      	--alignSJoverhangMin 6 \
      	--alignSJDBoverhangMin 1 \
      	--readFilesIn $fqFile1 $fqFile2 \
      	--outFileNamePrefix $outDir \
      	--runThreadN $cores \
		--quantMode TranscriptomeSAM \
      	--outFilterMatchNmin 45 \
    	--outSAMtype BAM SortedByCoordinate \
        --limitBAMsortRAM 75000000000

date
