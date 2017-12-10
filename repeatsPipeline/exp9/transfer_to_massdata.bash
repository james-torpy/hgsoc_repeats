#!/bin/bash

### transfer_fq_to_massdata.bash ###

shortDir="/short/ku3/jt3341/projects/hgsoc_repeats/RNA-seq/"
massDir="jt3341/projects/hgsoc_repeats/RNA-seq/"

ssh jt3341@raijin.nci.org.au "mdss put $shortDir/raw_files/*bam.gz $massDir/raw_files/; \
	rm $shortDir/raw_files/*bam.gz"