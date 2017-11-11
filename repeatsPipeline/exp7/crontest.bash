#!/bin/bash
source /home/jamtor/.bashrc
echo $SGE_ROOT
exec >> /tmp/jamtor.log
exec 2>&1
echo "works" >> /share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/exp7/works.txt
/opt/gridengine/bin/linux-x64/qsub -V /share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/exp7/test.bash
