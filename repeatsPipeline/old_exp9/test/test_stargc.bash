#!/bin/bash

mkdir -p /share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/results/star/GC/exp9/rfPT1

	starJobName="star."rfPT1
      bamJobName="bam."rfPT1
      sortJobName="sort."rfPT1
      indexJobName="index."rfPT1
      indexStatsJobName="indexstats."rfPT1
      outSam=$outDir"Aligned.out.sam"
      outBam=$outDir"rfPT1.bam"
      outSortedBam=$outDir"rfPT1.sorted.bam"
starline="/home/jamtor/local/bin/STAR --runMode alignReads \
      --readFilesCommand zcat \
      --genomeDir /home/jamtor/genomes/hg38_ercc/starRef/ \
    --outFilterType BySJout \
    --outSAMattributes NH HI AS NM MD\
    --outFilterMultimapNmax 999 \
    --outMultimapperOrder Random \
    --runRNGseed 666 \
    --outSAMmultNmax 1 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1500000 \
    --alignMatesGapMax 1500000 \
    --alignSJoverhangMin 6 \
    --alignSJDBoverhangMin 1 \
    --readFilesIn /share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/raw_files/fullsamples/bowtell_primary/fastq/rfPT1.1.fastq.gz /share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/raw_files/fullsamples/bowtell_primary/fastq/rfPT1.1.fastq.gz \
    --outFileNamePrefix /share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/results/star/GC/exp9/rfPT1/ \
    --runThreadN 6 \
    --outFilterMatchNmin 76 \
  --chimSegmentMin 25 \
    --chimJunctionOverhangMin 25 \
    --chimScoreMin 0 \
    --chimScoreDropMax 20 \
    --chimScoreSeparation 10 \
    --chimScoreJunctionNonGTAG -1 \
    --outSAMtype BAM Unsorted"


/opt/gridengine/bin/linux-x64/qsub -q short.q -N j3b_rfPT1 -b y -wd /share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/logs/exp9 -j y -R y -pe smp 6 -V $starline
