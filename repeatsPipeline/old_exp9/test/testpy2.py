#!/share/ClusterShare/software/contrib/marsmi/python/3.5.1/python

import os

os.environ

#print(os.environ['var'])

file = open('/share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/exp9/test/testpyout.txt', 'w')
file.write('%s\n' % 'pass')
