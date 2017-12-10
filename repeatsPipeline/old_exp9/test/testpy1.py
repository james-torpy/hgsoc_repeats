#!/share/ClusterShare/software/contrib/marsmi/python/3.5.1/python

import os

var = 'pass'

os.system('export var=' + var)

os.system('/share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/exp9/test/testpy2.py')

#os.system('/opt/gridengine/bin/linux-x64/qsub -q short.q -N testpy -b y -wd /share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/logs/exp9/ -j y -R y -V /share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/exp9/test/testpy2.py')
