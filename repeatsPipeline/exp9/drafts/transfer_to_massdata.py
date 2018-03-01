#!/home/jamtor/local/bin/Python-3.6.3/bin/python3

### transfer_to_massdata.py ###

# This script moves files to mass data on NCI and
# removes the original files

import os
import glob
import getpass
import pexpect
import sys

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

project_dir = home_dir2 + '/projects/' + projectname + '/RNA-seq/'
in_dir = project_dir + '/raw_files/'

print('')
print('This is the in_dir:')
print(in_dir)

nci_home = '/short/ku3/jt3341/'
nci_project = nci_home + '/projects/hgsoc_repeats/RNA-seq/'

nci_out = nci_project + '/raw_files/'

print('This is the nci_out:')
print(nci_out)

# scripts/logs directories:
script_dir = home_dir2 + \
	'/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/' + exp_name +\
	 '/'

print('')
print('The script_dir is:')
print(script_dir)

for f in os.listdir(in_dir):
	print('')
	print('Copying ' + f + ' to NCI short')
	os.system('rsync -avPS ' + in_dir + '/' + f + ' ' + nci_username + \
		':' + nci_out)

	os.system(script_dir + '/transfer_to_massdata.bash')
	os.remove(out_dir + '/' + f)


