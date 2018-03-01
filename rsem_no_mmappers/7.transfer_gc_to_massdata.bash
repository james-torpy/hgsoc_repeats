#!/bin/bash

### transfer_fq_to_massdata.bash ###

ssh jt3341@raijin.nci.org.au "mdss put \
	/short/ku3/jt3341/projects/hgsoc_repeats/RNA-seq/results/star/GC/rsem_no_mmappers/*.tar.gz \
	jt3341/projects/hgsoc_repeats/RNA-seq/results/star/GC/rsem_no_mmappers; \
	rm /short/ku3/jt3341/projects/hgsoc_repeats/RNA-seq/results/star/GC/rsem_no_mmappers/*.tar.gz"
