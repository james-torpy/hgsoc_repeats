# import libraries and set dir:
import HTSeq
import os
os.chdir('/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/results/star/GC/exp5')
os.getcwd()

# read in bam:
bam_reader = HTSeq.BAM_Reader( "bowtell_FT3_subset/Aligned.sortedByCoord.out.bam" )
# check first 5 lines of bam:
import itertools
for a in itertools.islice( bam_reader, 5 ):
     print a

# read in gencode annotation:
homeDir='/share/ScratchGeneral/jamtor/'
gc=homeDir + '/genomes/hg38_ercc/gencode_v24_hg38_annotation.gtf'
gtf_file=HTSeq.GFF_Reader( gc, end_included=True )
# check first 10 lines of annotation:
for feature in itertools.islice( gtf_file, 10 ):
     print feature

# initiate a GenomicArrayOfSets object and fill with exons only from annotation:
exons = HTSeq.GenomicArrayOfSets( "auto", stranded=True )

for feature in gtf_file:
    if feature.type == "exon":
       exons[ feature.iv ] += feature.name

for e in itertools.islice( gtf_file, 10 ):
	print e