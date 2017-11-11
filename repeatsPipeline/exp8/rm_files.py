#!/share/ClusterShare/software/contrib/fabbus/python/2.7.3/bin/python

### rm_files.py ###

# This script detects whether intermediate files are 
# still necessary, and if not deletes them: 

### 0. Define starting variables ###

import os
from path import path
import re

home_dir = '/share/ScratchGeneral/jamtor/'
project_dir = home_dir + '/projects/hgsoc_repeats/RNA-seq/'
script_dir = project_dir + '/scripts/repeatsPipeline/exp8/'
raw_dir = '/share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/raw_files/fullsamples/bowtell_primary'


### 1. Fetch bam files in raw_dir ###

bam_ids = []

d = path(raw_dir)
for f in d.walkfiles('*.bam'):
	bam_ids.append(re.sub('.bam', '', os.path.basename(f)))


### 2. Grep bamfiles against qstat jobs ###

os.system("qstat -u '*' | grep jamtor > " + \
	script_dir + "/qstat2.txt")

skiplist = []
removelist = []
for f in bam_ids:
	print ''
	print ''
	print 'Bam id is: ' + f
	with open(script_dir + '/qstat2.txt') as qstat:
		for line in qstat:
			if line.split(" ")[3].endswith(f) and f not in skiplist:
				print ''
				print f + ' is in ' + line
				if 'fq' in line or 'cp' in line:
					print 'Not touching ' + f + \
					', dependent job still running'
					skiplist.append(f)
				else:
					print 'Removing ' + raw_dir + '/' + f + '.bam'
    				os.system('rm ' + raw_dir + '/' + f + '.bam')
    				removelist.append(f)
			elif line.split(" ")[3].endswith(f) and f in skiplist:
				print ''
				print f + ' already taken care of'
    		
	if f not in skiplist and f not in removelist:
		print "Removing " + raw_dir + "/" + f + ".bam as it's not in qstat list"
		os.system('rm ' + raw_dir + '/' + f + '.bam')


### 3. Fetch fastq files ###

#fq_ids = []
#
#d = path(raw_dir + '/fastq/')
#for f in d.walkfiles('*[0-9].[1-2].fastq.gz'):
#	fq_ids.append(re.sub('.[1-2].fastq.gz', '', os.path.basename(f)))
#
#skiplist2 = []
#for f in fq_ids:
#	print ''
#	print ''
#	print 'Fastq id is: ' + f
#	with open(script_dir + '/qstat2.txt') as qstat:
#		for line in qstat:
#			if line.split(" ")[3].endswith(f) and f not in skiplist2:
#				print ''
#				print f + ' is in ' + line
#				if 'j2' in line or 'j3' in line:
#					print 'Not touching ' + f + \
#					', dependent job still running'
#					skiplist2.append(f)
#				else:
#					print 'Removing ' + raw_dir + '/' + f + '.1.fastq.gz'
#    				#os.system('rm ' + raw_dir + '/' + f + '.1.fastq.gz')
#    				print 'Removing ' + raw_dir + '/' + f + '.2.fastq.gz'
#    				#os.system('rm ' + raw_dir + '/' + f + '.2.fastq.gz')
#    				print 'Removing ' + raw_dir + '/' + f + '.unpaired.1.fastq.gz'
#    				#os.system('rm ' + raw_dir + '/' + f + '.unpaired.1.fastq.gz')
#    				print 'Removing ' + raw_dir + '/' + f + '.unpaired.2.fastq.gz'
#    				#os.system('rm ' + raw_dir + '/' + f + '.unpaired.2.fastq.gz')
#    		else:
#    			print ''
#    			print f + ' not in qstat list or already taken care of'
#
#
#
#