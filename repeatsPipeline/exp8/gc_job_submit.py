#!/share/ClusterShare/software/contrib/fabbus/python/2.7.3/bin/python

### add_gc.py ###

# This script detects whether the number of 
# cores I'm using drops below 80, and if so
# submits the next GC job #

# run marsmi/python/3.5.1 before executing!

### 0. Define starting variables ###

import os
#import numpy as np
#import pandas as pd
import re

home_dir = '/share/ScratchGeneral/jamtor/'
project_dir = home_dir + '/projects/hgsoc_repeats/RNA-seq/'
script_dir = project_dir + '/scripts/repeatsPipeline/exp8/'
raw_dir = '/share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/raw_files/fullsamples/bowtell_primary'



### 1. Fetch number of cores being used ###

# save qstat output as file:
os.system("qstat -u '*' | grep jamtor | grep -e 'q@' | grep -v QRLOGIN | awk '{print $9}' > " + script_dir + "/qstat.txt")

# fetch core number:
slots = 0
with open(script_dir + '/qstat.txt') as f:
    for line in f:
    	slots = slots + int(re.sub('\\n', '', line))


### 2. Initiate counter for gc jobs ###

if slots < 61:
	print 'Less than 61 cores being used!'
	if os.path.exists(raw_dir + '/htseq_count.txt'):
		print('Count file exists, adding to this')
		print ''
		with open(raw_dir + '/htseq_count.txt') as f:
			for line in f:
				num = int(line) + 1
				file = open((raw_dir + '/htseq_count.txt'), 'w')
				file.write("%s\n" % num)
				file.close()
	else:
		print('Count file does not exist, creating')
		print ''
		num=4
		file = open((raw_dir + '/htseq_count.txt'), 'w')
		file.write('%s\n' % '4')
		file.close()


	### 3. Submit next job in gc_redo.txt

	# define next file:
	with open(raw_dir + '/htseq_redo.txt') as f:
	   	for i, subfile in enumerate(f, 1):
	  			if i == num:
					break
	print 'File to be submitted is: ' + subfile
	print ''
	
	# write this file into files.txt:
	file = open(raw_dir + '/files.txt', 'w')
	file.write('%s\n' % subfile)
	file.close()
	# call 1.bamtofastq_mansubmit.bash:
	os.system(script_dir + '/1.bamtofastq_mansubmit.bash')
