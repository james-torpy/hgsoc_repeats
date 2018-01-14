#!/home/jamtor/local/bin/Python-3.6.3/bin/python3

### transfer_gc_to_massdata.py ###

# This script tars and moves any star gc file directories to
# mass data on NCI and removes the original files

# qsub line:
# qsub -q short.q -N gctr -b y -wd \
# /share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/logs/exp9 \
# -j y -R y -pe smp 1 -V /share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/exp9/transfer_gc_to_massdata.py

import os
import glob
import getpass
import pexpect
import sys
import shutil
import re

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
results_dir = project_dir + '/results/'

gc_dir = results_dir + '/star/GC/' + exp_name + '/' 
print('')
print('This is the gc_dir:')
print(gc_dir)

nci_home = '/short/ku3/jt3341/'
nci_project = nci_home + '/projects/hgsoc_repeats/RNA-seq/'

nci_gc = nci_project + '/results/star/GC/' + exp_name + "/"

print('')
print('This is the nci_gc:')
print(nci_gc)

# scripts/logs directories:
script_dir = home_dir2 + \
	'/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/' + exp_name + \
	 '/'

print('')
print('The script_dir is:')
print(script_dir)

# define files as those that match the regex pattern specified:
files = os.listdir(gc_dir)

reg = re.compile('.*[A-Z][A-Z].*[0-9]$')
selected_files = filter(reg.search, files)

files = []
for f in selected_files:
	files.append(f)

print('Files are: ' + str(files))

# copy logs of each sample to separate log directory:
for f in files:
	print('')
	print('Copying ' + f + ' log file to log directory')
	os.system('mkdir -p ' + gc_dir + '/logs/' + f)
	os.system('cp ' + gc_dir + '/' + f + '/Log.final.out ' + gc_dir + \
		'/logs/' + f + '/')
	# if bam file exists, tar its directory, submit it to
	# transfer_gc_to_massdata.bash and remove the original directory:
	if os.path.isfile(gc_dir + '/' + f + '/Aligned.novosortedByName.out.bam'):
		print('')
		print('Tarring ' + f)
		os.system('tar -zcvf ' + gc_dir + '/' + f + '.tar.gz ' + gc_dir + '/' + f)
		if os.path.isfile(gc_dir + '/' + f + '.tar.gz'):
			print('')
			print('Removing original directory ' + f)
			shutil.rmtree(gc_dir + '/' + f)
			print('')
			print('Copying ' + gc_dir + '/' + f + '.tar.gz to NCI short')
			os.system('rsync -avPS ' + gc_dir + '/' + f + '.tar.gz ' + nci_username + \
				':' + nci_gc)
			os.system(script_dir + '/transfer_gc_to_massdata.bash')
			os.remove(gc_dir + '/' + f + '.tar.gz')
