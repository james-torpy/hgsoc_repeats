#!/home/jamtor/local/bin/Python-3.6.3/bin/python3

import os

os.system('source /home/jamtor/.bashrc')

os.system('echo pass > /share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/logs//exp9_test2/crontab_test1.txt')

testline = '/home/jamtor/local/bin/Python-3.6.3/bin/python3.6 /share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/exp9_test2/crontab_test2.py'

#os.system('echo /opt/gridengine/bin/linux-x64/qsub -b y -V "/home/jamtor/local/bin/Python-3.6.3/bin/python3.6 /share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/exp9_test2/crontab_test2.py" > /share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/exp9_test2/crontab_test3.txt')

#os.system('echo /opt/gridengine/bin/linux-x64/qsub -b y -V /home/jamtor/local/bin/Python-3.6.3/bin/python3.6 /share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/exp9_test2/crontab_test2.py > /share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/exp9_test2/crontab_test3.txt')

#os.system('echo /opt/gridengine/bin/linux-x64/qsub -b y -V ' + testline ' > /share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/exp9_test2/crontab_test3.txt')

#os.system('echo /opt/gridengine/bin/linux-x64/qsub -b y -V ' + testline + ' > /share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/logs/exp9_test2/crontab_test3.txt')

os.system('/opt/gridengine/bin/linux-x64/qsub -b y -V ' + testline)

os.system('/opt/gridengine/bin/linux-x64/qsub -b y -V "/home/jamtor/local/bin/Python-3.6.3/bin/python3.6 /share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/exp9_test2/crontab_test2.py"')

os.system('/opt/gridengine/bin/linux-x64/qsub -b y -V /home/jamtor/local/bin/Python-3.6.3/bin/python3.6 /share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/exp9_test2/crontab_test2.py')

os.system('"/opt/gridengine/bin/linux-x64/qsub -b y -V /home/jamtor/local/bin/Python-3.6.3/bin/python3.6 /share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/exp9_test2/crontab_test2.py"')
