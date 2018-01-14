#!/home/jamtor/local/bin/Python-3.6.3/bin/python3

### transfer_fq_to_massdata.py ###

# This script moves any fastq files to mass data on NCI and
# removes the original files

import os
import glob
import getpass
import pexpect
import sys

# qsub line:
# qsub -q short.q -N gctr -b y -wd \
# /share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/logs/exp9 \
# -j y -R y -pe smp 1 -V /share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/exp9/transfer_fq_to_massdata.py


os.system('source /home/jamtor/.bashrc')

# make directory hierachy:
projectname = 'hgsoc_repeats'
sample_type = 'fullsamples/bowtell_primary'
exp_name = 'exp9'
nci_username = 'jt3341@raijin.nci.org.au'

# home_dir is where input/intermediate files are located:
home_dir = '/share/ScratchGeneral/jamtor/'
# home_dir2 is where scripts are located and final outputs will be stored:
home_dir2 = '/share/ClusterShare/thingamajigs/jamtor'

project_dir = home_dir + '/projects/' + projectname + '/RNA-seq/'
raw_dir = project_dir + '/raw_files/' + sample_type + '/'

fq_dir = raw_dir + '/fastq/'

print('')
print('This is the fq_dir:')
print(fq_dir)

nci_home = '/short/ku3/jt3341/'
nci_project = nci_home + '/projects/hgsoc_repeats/RNA-seq/'

nci_fq = nci_project + '/raw_files/fastq/'

print('This is the nci_fq:')
print(nci_fq)

# scripts/logs directories:
script_dir = home_dir2 + \
	'/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/' + exp_name +\
	 '/'

print('')
print('The script_dir is:')
print(script_dir)

for f in os.listdir(fq_dir):
	print(f)
	if 'fastq.gz' in f:
		print('')
		print('Copying ' + f + ' to NCI short')
		os.system('rsync -avPS ' + fq_dir + '/' + f + ' ' + nci_username + \
			':' + nci_fq)
	
		os.system(script_dir + '/transfer_fq_to_massdata.bash')
		os.remove(fq_dir + '/' + f)


