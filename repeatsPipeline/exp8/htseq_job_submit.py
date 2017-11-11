#!/share/ClusterShare/software/contrib/fabbus/python/2.7.3/bin/python

### add_gc.py ###

# This script detects whether the number of 
# cores I'm using drops below 80, and if so
# submits the next GC job #

# run marsmi/python/3.5.1 before executing!

### 0. Define starting variables ###

import os
import re
import itertools
import pandas as pd

home_dir = '/share/ScratchGeneral/jamtor/'
project_dir = home_dir + '/projects/hgsoc_repeats/RNA-seq/'
script_dir = project_dir + '/scripts/repeatsPipeline/exp8/'
raw_dir = '/share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/raw_files/fullsamples/bowtell_primary'



### 1. Fetch number of cores being used ###

# save qstat cores output as file:
os.system("qstat -u '*' | grep jamtor | grep -e 'q@' | grep -v QRLOGIN | awk '{print $9}' > " + script_dir + "/qstat.txt")

# fetch core number:
slots = 0
with open(script_dir + '/qstat.txt') as f:
	for line in f:
		slots = slots + int(re.sub('\\n', '', line))

# save entire qstat output for section 3:
os.system("qstat -u '*' | grep jamtor > " + script_dir + "/qstat2.txt")


### 2. Initiate counter for htseq jobs ###

if slots < 72:
	print 'Less than 72 cores being used!'
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
		num=1
		file = open((raw_dir + '/htseq_count.txt'), 'w')
		file.write('%s\n' % '1')
		file.close()


	### 3. Submit next job in htseq_redo.txt

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


### 4. Limit jobs to never take up more than 100 cores ###

elif slots > 100:
	print('Too many cores being used, removing all except first 6 jobs')
	qstat = pd.read_csv(script_dir + '/qstat2.txt', header = None)
	for line in qstat.iterrows():
		j_id = str.split(str(line[1]))[1]
	   	print(j_id)
	   	os.system('qdel ' + j_id)

