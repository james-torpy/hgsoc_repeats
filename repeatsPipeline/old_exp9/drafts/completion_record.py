#!/share/ClusterShare/software/contrib/phuluu/python/2.7.8/bin/python

### completion_record ###

# This python script returns which samples have been processed at what stage of the
# repeatsPipeline #

### 0. Define starting variables ###

from os import listdir
from os.path import isfile, join
import re
import os
import pandas as pd
from path import path

exp_name = 'exp9'

# define directories:
home_dir = '/share/ScratchGeneral/jamtor/'
project_dir = home_dir + '/projects/hgsoc_repeats/RNA-seq/'
results_dir = project_dir + '/results/'

report_dir = project_dir + '/reports/' + exp_name + '/'
ribo_dir = results_dir + '/star/ribo/' + exp_name + '/'
gc_dir = results_dir + '/star/GC/' + exp_name + '/'
htseq_dir = path(results_dir + '/htseq/' + exp_name + '/')
rkey_dir = project_dir + '/raw_files/fullsamples/bowtell_primary/'
record_dir = project_dir + '/record/'

if not os.path.exists(record_dir):
	os.makedirs(record_dir)

print('The report_dir is ' + report_dir)
print('The ribo_dir is ' + ribo_dir)
print('The gc_dir is ' + gc_dir)
print('The htseq_dir is ' + htseq_dir)
print('The record_dir is ' + record_dir)


### 1. Fetch directory names in reports directory ###

rep = []
for f in listdir(report_dir):
	if 'subset' not in f:
		rep.append(re.sub('\.report.txt$', '', f))


### 2. Fetch directory names in ribo directory ###

ribo = []
for d in listdir(ribo_dir):
	if 'subset' not in d:
		for file in listdir(ribo_dir + d):
			if file == 'Log.final.out':
				ribo.append(d)


### 3. Fetch directory names in star GC directory ###

gc = []
for d in listdir(gc_dir):
	if 'subset' not in d:
		for file in listdir(gc_dir + d):
			print(file)
			if file == 'Log.final.out':
				gc.append(d)


### 4. Fetch names of the htseq counts all, custom3 and gc:

all_htseq=[]
custom3_htseq=[]
gc_htseq=[]

for f in htseq_dir.walk():
	if '.htseq' in str(f) and 'subset' not in str(f):
			print(f)
			if 'all' in str(f) and f.size>1:
				all_htseq.append(re.sub('\.all.htseq.txt$', '', os.path.basename(str(f))))
			elif 'custom3' in str(f) and f.size>1:
				custom3_htseq.append(re.sub('\.custom3.htseq.txt$', '', os.path.basename(str(f))))
			elif 'gc' in str(f) and f.size>1:
				gc_htseq.append(re.sub('\.gc.htseq.txt$', '', os.path.basename(str(f))))


### 5. Create empty data frame for each sample ###

# create sorted list with all unique samples from above:
all_samples_list = list(set(rep + ribo + gc + all_htseq + custom3_htseq +
	gc_htseq))

all_samples_list.sort()

# for each list item, iterate through each above lists to record whether it is
# present:
record = {'rep':'', 'ribo':'', 'gc':'', 'all_htseq':'', 'custom3_htseq':'',
	'gc_htseq':''}


for i in range((len(all_samples_list)-1)+1):
	print('for ' + all_samples_list[i])
	for t in ['rep', 'ribo', 'gc', 'all_htseq', 'custom3_htseq',
	'gc_htseq']:
		if i==0:
			if all_samples_list[i] in eval(t):
				print('sample is in ' + t + ', appending to ' + t)
				record[t] = ['True']
			else:
				print('sample not in ' + t)
				record[t] = ['False']
		else:
			if all_samples_list[i] in eval(t):
				print('sample is in ' + t + ', appending to ' + t)
				record[t] = record[t] + ['True']
			else:
				print('sample not in ' + t)
				record[t] = record[t] + ['False']

# add sample names to dictionary:
record['samples'] = all_samples_list

# convert record dictionary into data frame:
record_df = pd.DataFrame(record)

# add sample names as row:
record_df.set_index('samples', inplace=True)

# rearrange the column names of df:
colnames=["rep","ribo", "gc", "all_htseq", "custom3_htseq",
"gc_htseq"]
record_df=record_df.reindex(columns=colnames)

# write dataframe to record directory:
print('')
print('writing ' + record_dir + '/completion_record.tab')
record_df.to_csv(record_dir + '/completion_record.tab', sep='\t')

# drop all rows containing 'False' extremely inefficiently:
record_df = record_df[record_df.rep != 'False']
record_df = record_df[record_df.ribo != 'False']
record_df = record_df[record_df.gc != 'False']
record_df = record_df[record_df.all_htseq != 'False']
record_df = record_df[record_df.custom3_htseq != 'False']
record_df = record_df[record_df.gc_htseq != 'False']

# get the sample names of completely processed files:
complete = list(record_df.index)

# write list to file:
print('')
print('writing ' + record_dir + '/completed_files.txt')

compfile = open(record_dir + '/completed_files.txt', 'w')
for item in complete:
  compfile.write('%s\n' % item)

