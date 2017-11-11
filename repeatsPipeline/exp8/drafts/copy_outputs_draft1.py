#!/share/ClusterShare/software/contrib/fabbus/python/2.7.3/bin/python

### copy_outputs.py ###

# This script copies output files from repeatsPipeline which are in ClusterShare
# directories but not ScratchGeneral directories

### 0. Define starting variables ###

import os
from path import path

exp_name = "exp8"

home_dir = '/share/ClusterShare/thingamajigs/jamtor/'
project_dir = home_dir + '/projects/hgsoc_repeats/RNA-seq/'
results_dir = project_dir + '/results/'
gc_dir = results_dir + '/star/GC/' + exp_name + '/'
htseq_dir = results_dir + '/htseq/' + exp_name + '/'
ribo_dir = results_dir + '/star/ribo/' + exp_name + '/'

dest_home = '/share/ScratchGeneral/jamtor/'
dest_project = dest_home + '/projects/hgsoc_repeats/RNA-seq/'
dest_results = dest_project + '/results/'
dest_gc = dest_results + '/star/GC/' + exp_name + '/'
dest_htseq = dest_results + '/htseq/' + exp_name + '/'
dest_ribo = dest_results + '/star/ribo/' + exp_name + '/'


### 1. Fetch names of gc STAR dirs in ClusterShare and check ScratchGeneral for
# these ###

d = path(gc_dir)
for f in d.walkfiles('Log.final.out'):
	print('Log file is ' + f)

	temp = str.split(f, '/')
	uID = temp[len(temp)-2]
	print('uID is ' + uID)

	if os.path.exists(dest_gc + '/' + uID):
		if os.path.exists(dest_gc + '/' + uID + '/Log.final.out'):
			print('Log.final.out exists for '+ uID + ', assuming all files are already copied')
		else:		
			print('Copying ' + f + ' to ' + dest_gc + '/' + uID)
			print('cp ' + f + ' ' + dest_gc + '/' + uID)
			#os.system('cp ' + f + ' ' + dest_gc)
			print('Copying ' + gc_dir + '/' + uID + '/Chimeric* to ' + dest_gc + '/' + uID)
			print('cp ' + gc_dir + '/' + uID + '/Chimeric* ' + dest_gc + '/' + uID)
			#os.system('cp ' + gc_dir + '/' + uID + '/Chimeric* ' + dest_gc + '/' + uID)
			print('Copying ' + gc_dir + '/' + uID + '/SJ.out.tab to ' + dest_gc + '/' + uID)
			print('cp ' + gc_dir + '/' + uID + '/SJ.out.tab ' + dest_gc + '/' + uID)
			#os.system('cp ' + gc_dir + '/' + uID + '/SJ.out.tab ' + dest_gc + '/' + uID)
	else:
		print('Creating directory for ' + uID)
		os.system('mkdir ' + dest_gc + '/' + uID)
		print('Copying ' + f + ' to ' + dest_gc + '/' + uID)
		print('cp ' + f + ' ' + dest_gc + '/' + uID)
		#os.system('cp ' + f + ' ' + dest_gc)
		print('Copying ' + gc_dir + '/' + uID + '/Chimeric* to ' + dest_gc + '/' + uID)
		print('cp ' + gc_dir + '/' + uID + '/Chimeric* ' + dest_gc + '/' + uID)
		#os.system('cp ' + gc_dir + '/' + uID + '/Chimeric* ' + dest_gc + '/' + uID)
		print('Copying ' + gc_dir + '/' + uID + '/SJ.out.tab to ' + dest_gc + '/' + uID)
		print('cp ' + gc_dir + '/' + uID + '/SJ.out.tab ' + dest_gc + '/' + uID)
		#os.system('cp ' + gc_dir + '/' + uID + '/SJ.out.tab ' + dest_gc + '/' + uID)

### 1. Fetch names of htseq files in ClusterShare and check ScratchGeneral for
# these ###

d = path(htseq_dir)
for f in d.walkfiles('*.gc.htseq.txt'):
	statinfo1 = os.stat(f)
	if statinfo1.stsize > 0:
		print('htseq file is ' + f)
	
		temp = str.split(f, '/')
		uID = temp[len(temp)-2]
		print('uID is ' + uID)
		if os.path.exists(dest_gc + '/' + uID):
			if os.path.exists(dest_gc + '/' + uID + '/Log.final.out'):
				print('Log.final.out exists for uID, assuming all files are already copied')
			else:		
				print('Copying ' + f + ' to ' + dest_gc + '/' + uID)
				print('cp ' + f + ' ' + dest_gc + '/' + uID)
				#os.system('cp ' + f + ' ' + dest_gc + '/' + uID)
				print('Copying ' + gc_dir + '/' + uID + '/Chimeric* to ' + dest_gc + \
				'/' + uID)
				print('cp ' + gc_dir + '/' + uID + '/Chimeric* ' + dest_gc + '/' + uID)
				#os.system('cp ' + gc_dir + '/' + uID + '/Chimeric* ' + dest_gc + '/' + uID)
				print('Copying ' + gc_dir + '/' + uID + '/SJ.out.tab to ' + dest_gc + \
				'/' + uID)
				print('cp ' + gc_dir + '/' + uID + '/SJ.out.tab ' + dest_gc + '/' + uID)
				#os.system('cp ' + gc_dir + '/' + uID + '/SJ.out.tab ' + dest_gc + '/' + uID)
		else:
			print('Creating directory for ' + uID)
			os.system('mkdir ' + dest_gc + '/' + uID)
			print('Copying ' + f + ' to ' + dest_gc + '/' + uID)
			print('cp ' + f + ' ' + dest_gc + '/' + uID)
			#os.system('cp ' + f + ' ' + dest_gc + '/' + uID)
			print('Copying ' + gc_dir + '/' + uID + '/Chimeric* to ' + dest_gc + '/' \
			+ uID)
			print('cp ' + gc_dir + '/' + uID + '/Chimeric* ' + dest_gc + '/' + uID)
			#os.system('cp ' + gc_dir + '/' + uID + '/Chimeric* ' + dest_gc + '/' + uID)
			print('Copying ' + gc_dir + '/' + uID + '/SJ.out.tab to ' + dest_gc + \
			'/' + uID)
			print('cp ' + gc_dir + '/' + uID + '/SJ.out.tab ' + dest_gc + '/' + uID)
			#os.system('cp ' + gc_dir + '/' + uID + '/SJ.out.tab ' + dest_gc + '/' + uID)
	else:
		print(f 'is 0 bytes in size')





temp = str.split(f, '/')
uID = temp[len(temp)-2]
print('uID is ' + uID)
for t in ['gc', 'custom3', 'all']:
	if os.path.exists(dest_htseq + '/' + uID):
		statinfo2 = os.stat(dest_htseq + '/' + uID):
		if statinfo2.stsize > 0:
			if os.path.exists(dest_htseq + '/' + uID + '/' + t + '.' + t + \
				'.htseq.txt'):
				print(t + '.htseq.txt exists for ' + uID + ', assuming all files are already copied')
			else:		
				print('Copying ' + htseq_dir + '/' + uID + "/" + uID + '.' + t + \
				'.htseq.txt' + ' to ' + dest_htseq + '/' + uID)
				print('cp ' + htseq_dir + '/' + uID + "/" + uID + '.' + t + '.htseq.txt ' + \
				dest_htseq + '/' + uID)
				#os.system('cp ' htseq_dir + '/' + uID + "/" + uID + '.' + t + '.htseq.txt ' + \
				#dest_htseq + '/' + uID)
	else:
		os.system('mkdir ' + dest_htseq + '/' + uID)
		print('Copying ' + htseq_dir + '/' + uID + "/" + uID + '.' + t + \
		'.htseq.txt' + ' to ' + dest_htseq + '/' + uID)
		print('cp ' + htseq_dir + '/' + uID + "/" + uID + '.' + t + '.htseq.txt ' + \
		dest_htseq + '/' + uID)
		#os.system('cp ' htseq_dir + '/' + uID + "/" + uID + '.' + t + '.htseq.txt ' + \
		#dest_htseq + '/' + uID)

