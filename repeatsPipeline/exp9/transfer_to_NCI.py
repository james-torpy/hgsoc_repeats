#!/home/jamtor/local/bin/Python-3.6.3/bin/python3

### transfer_GC_to_NCI.py ###

# This script copies log files to new directories,
# tars and gzips each star gc directory, and sends
# it to NCI short with fastq files


### 1. Define directories and import packages ###

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
q = 'short'
nci_username = 'jt3341@raijin.nci.org.au'
#nci_pass = getpass.getpass()

# home_dir is where input/intermediate files are located:
home_dir = '/share/ScratchGeneral/jamtor/'
# home_dir2 is where scripts are located and final outputs will be stored:
home_dir2 = '/share/ClusterShare/thingamajigs/jamtor'

project_dir = home_dir + '/projects/' + projectname + '/RNA-seq/'
raw_dir = project_dir + '/raw_files/' + sample_type + '/'
results_dir = project_dir + '/results/'

fq_dir = raw_dir + '/fastq/'
gc_dir = results_dir + '/star/GC/' + exp_name

print('')
print('This is the raw_dir:')
print(raw_dir)
print('')
print('This is the fq_dir:')
print(fq_dir)
print('')
print('This is the gc_dir:')
print(gc_dir)

nci_home = '/short/ku3/jt3341/'
nci_project = nci_home + '/projects/hgsoc_repeats/RNA-seq/'
nci_raw = nci_project + '/raw_files/'
nci_results = nci_project + '/results/'

nci_fq = nci_raw + '/fastq/'
nci_gc = nci_results + '/star/GC/' + exp_name + "/"

print('This is the nci_fq:')
print(nci_fq)
print('')
print('This is the nci_gc:')
print(nci_gc)

mass_project = 'jt3341//projects/hgsoc_repeats/RNA-seq/'
mass_raw = mass_project + '/raw_files/'
mass_results = mass_project + '/results/'
mass_script = mass_project + '/scripts/'

mass_fq = mass_raw + '/fastq/'
mass_gc = mass_results + '/star/GC/' + exp_name + "/"

print('This is the mass_fq:')
print(mass_fq)
print('')
print('This is the mass_gc:')
print(mass_gc)

# scripts/logs directories:
script_dir = home_dir2 + \
	'/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/' + exp_name +\
	 '/'
log_dir = project_dir + '/logs/' + exp_name + '/'
os.makedirs(log_dir, exist_ok=True)

print('')
print('The log_dir is:')
print(log_dir)
print('')
print('The script_dir is:')
print(script_dir)


#f = os.listdir(fq_dir)[17]
#print('')
#print('Copying ' + f + ' to NCI short')
#os.system('rsync -avPS ' + fq_dir + '/' + f + ' ' + nci_username + \
#	':' + nci_fq)

child = pexpect.spawn('ssh jt3341@raijin.nci.org.au')
#child.expect("jt3341@raijin.nci.org.au's password: ")
#child.sendline(nci_pass)
child.expect('[jt3341@raijin1 ~]$')
child.sendline('touch test.txt')
#child.expect('[jt3341@raijin1 ~]$')
#child.sendline('for f in `ls /short/ku3/jt3341//projects/hgsoc_repeats/RNA-seq//raw_files//fastq/`; do')
#child.expect('>')
#child.sendline('mdss put /short/ku3/jt3341//projects/hgsoc_repeats/RNA-seq//raw_files//fastq/$f jt3341/projects/hgsoc_repeats/RNA-seq/raw_files/fastq')
#child.expect('>')
#child.sendline('rm -f /short/ku3/jt3341//projects/hgsoc_repeats/RNA-seq//raw_files//fastq/$f')
#child.expect('>')
#child.sendline('done;')

### 1. Send fastq files to NCI short ###

#for f in os.listdir(fq_dir):
#	print('')
#	print('Copying ' + f + ' to NCI short')
##	os.system('rsync -avPS ' + fq_dir + '/' + f + ' ' + nci_username + \
##		':' + nci_fq)
#
#	# login to nci and run script to transfer file to mass data:
#	child = pexpect.spawn('ssh jt3341@raijin.nci.org.au')
#	child.expect("jt3341@raijin.nci.org.au's password:")
#	child.sendline(nci_pass)
#	child.expect('[jt3341@raijin1 ~]$')
#	child.sendline('rm test.txt')
#	child.expect('[jt3341@raijin1 ~]$')
#	child.sendline('touch blarg.txt')


	#os.system("/short/ku3/jt3341/projects/hgsoc_repeats/RNA-seq/scripts/mdss_transfer.bash")
	#logout


### 2. Tar star GC dirs and send to NCI short

#dirs = os.listdir(gc_dir)
#for d in dirs:
#	if os.path.isdir(gc_dir + '/' + d) and '.tar.gz' not in d and \
#		'Store' not in d and 'log' not in d:
#
#		# tar directory:
#		print('')
#		print('Tarring ' + gc_dir + '/' + d)
#		os.system('tar -zvcf ' + gc_dir + '/' + d + '.tar.gz ' + gc_dir + \
#			'/' + d)
#
#		# send tar file to NCI:
#		print('')
#		print('Copying ' + gc_dir + '/' + d + '.tar.gz' + ' to NCI short')
#		os.system('rsync -avPS ' + gc_dir + '/' + d + '.tar.gz' + ' ' + \
#			nci_username + ':' + nci_gc)
#
#		# login to nci and run script to transfer file to mass data:
#		os.system('ssh jt3341@raijin.nci.org.au')
#		os.system("/short/ku3/jt3341/projects/hgsoc_repeats/RNA-seq/scripts/mdss_transfer.bash")
#		logout

