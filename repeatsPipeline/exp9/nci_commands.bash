#!/bin/bash

ssh jt3341@raijin.nci.org.au "touch test.txt; \
	mdss put /short/ku3/jt3341/projects/hgsoc_repeats/RNA-seq/raw_files/fastq/*.fastq.gz \
	jt3341/projects/hgsoc_repeats/RNA-seq/raw_files/fastq; rm /short/ku3/jt3341/projects/hgsoc_repeats/RNA-seq/raw_files/fastq/prPT20.unpaired.1.fastq.gz; \
	rm /short/ku3/jt3341/projects/hgsoc_repeats/RNA-seq/raw_files/fastq/*.fastq.gz"
