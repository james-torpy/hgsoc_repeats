#!/bin/bash

source /home/jamtor/.bashrc
echo $SGE_ROOT

numcores=6
draft="stest"
uID="prPT11"

/opt/gridengine/bin/linux-x64/qsub -q long.q -N $draft.j5c_$uID -hold_jid 3c_$uID -b y -wd \
/share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/logs/exp9 -j y -R y -P \
DSGClinicalGenomics -pe smp $numcores -V /share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/exp9/test/5c.htseq-rep_custom3_test.bash $draft $numcores $uID