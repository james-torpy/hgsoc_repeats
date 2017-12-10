module load gi/novosort/precompiled/1.03.08

/opt/gridengine/bin/linux-x64/qsub -q short.q -N j3c_prPT11 -b y -wd /share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/logs/exp9 \
-j y -R y -pe smp 7 -V /home/jamtor//local/bin/novocraft/novosort -t /share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/results/star/GC/exp9/prPT11/temp/ \
-n -c 6 -m 22G /share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/results/star/GC/exp9/prPT11/Aligned.out.bam \
> /share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/results/star/GC/exp9/prPT11/Aligned.novosortedByName.out.bam