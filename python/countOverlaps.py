 #!/usr/bin/env python

### countOverlaps.py ###

# This script takes a bam file of mapped reads and counts overlaps with
# annotations for 17 repeats classes and a gencode annotation, outputting a
# repeat counts and a gencode counts .tab file


### 0. Set up variables and directories ###

import HTSeq

# define starting variables:
project = 'hgsoc_repeats'
dataSource = 'dx'

# define directories:
homeDir = '/share/ScratchGeneral/jamtor/'
projectDir = homeDir + '/projects/' + project
annotDir = projectDir + '/RNA-seq/refs/c_GenesTab'
bamDir <- projectDir + '/RNA-seq/results/star/' + methodName
resultsDir = projectDir, '/RNA-seq/results'