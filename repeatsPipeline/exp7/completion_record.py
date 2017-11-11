#!/share/ClusterShare/software/contrib/fabbus/python/2.7.3/bin/python

### completion_record ###

# This python script returns which samples have been processed at what stage of the
# repeatsPipeline #

### 0. Define starting variables ###

from os import listdir
from os.path import isfile, join
import re
import os
import pandas as pd

exp_name = 'exp7_and_8'

# define directories:
home_dir = '/share/ScratchGeneral/jamtor/'
project_dir = home_dir + '/projects/hgsoc_repeats/RNA-seq/'
results_dir = project_dir + '/results/'

report_dir = results_dir + '/reports/' + exp_name
ribo_dir = results_dir + '/star/ribo/' + exp_name
gc_dir = results_dir + '/star/GC/' + exp_name
htseq_dir = results_dir + '/htseq/' + exp_name
rkey_dir = project_dir + '/raw_files/fullsamples/bowtell_primary/'

print 'The report_dir is ' + report_dir
print 'The ribo_dir is ' + ribo_dir
print 'The gc_dir is ' + gc_dir
print 'The htseq_dir is ' + htseq_dir


### 1. Fetch directory names in reports directory ###

rep = []
for f in listdir(report_dir):
	if 'subset' not in f and '.2' not in f:
		rep.append(re.sub('\.report.txt$', '', f))


### 2. Fetch directory names in ribo directory ###

ribo = []
for d in listdir(ribo_dir):
	if 'subset' not in d and '.2' not in d:
		ribo.append(d)


### 3. Fetch directory names in star GC directory ###

gc = []
for d in listdir(gc_dir):
	if 'subset' not in d and 'DS_Store' not in d and '.2' not in d:
		gc.append(d)

### 4. Fetch names of the htseq counts all, custom3 and gc:

all_htseq=[]
custom3_htseq=[]
gc_htseq=[]

for f in os.walk(htseq_dir):
	for x in f[2]:
		if 'htseq' in x and 'subset' not in x and '.2' not in x:
			print x
			if 'all' in x:
				all_htseq.append(re.sub('\.all.htseq.txt$', '', x))
			elif 'custom3' in x:
				custom3_htseq.append(re.sub('\.custom3.htseq.txt$', '', x))
			elif 'gc' in x:
				gc_htseq.append(re.sub('\.gc.htseq.txt$', '', x))


### 5. Create empty data frame for each sample ###

# create sorted list with all unique samples from above:
all_samples_list = list(set(rep + ribo + gc + all_htseq + custom3_htseq +
	gc_htseq))

test = set(rep + ribo + gc + all_htseq + custom3_htseq + gc_htseq)

all_samples_list.sort()

# for each list item, iterate through each above lists to record whether it is
# present:
record = {'rep':'', 'ribo':'', 'gc':'', 'all_htseq':'', 'custom3_htseq':'',
	'gc_htseq':''}

for i in range((len(all_samples_list)-1)+1):
	print 'for ' + all_samples_list[i]
	for t in ['rep', 'ribo', 'gc', 'all_htseq', 'custom3_htseq',
	'gc_htseq']:
		if i==0:
			if all_samples_list[i] in eval(t):
				print 'sample is in ' + t + ', appending to ' + t
				record[t] = ['True']
			else:
				print 'sample not in ' + t
				record[t] = ['False']
		else:
			if all_samples_list[i] in eval(t):
				print 'sample is in ' + t + ', appending to ' + t
				record[t] = record[t] + ['True']
			else:
				print 'sample not in ' + t
				record[t] = record[t] + ['False']

# add sample names to dictionary:
record['samples'] = all_samples_list

# convert record dictionary into data frame:
record_df = pd.DataFrame(record)

# add sample names as row:
record_df.set_index('samples', inplace=True)

# write dataframe to htseq directory:
record_df.to_csv(htseq_dir + '/completion_record.tab', sep='\t')

# create list of sample names with all output files present and save as
# text file:
completed_files = list(record_df[(record_df != 'False').all(1)].index)

comp_file = open(htseq_dir + '/completed_files.txt', 'w')
for item in completed_files:
	print item
	comp_file.write("%s\n" % item)



