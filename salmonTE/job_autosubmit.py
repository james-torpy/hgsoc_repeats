#!/home/jamtor/local/bin/Python-3.6.3/bin/python3


### job_autosubmit.py ###

# This script detects which output files still must be created from the repeatsPipeline
# and submits jobs for them #

# run once and then input the following command into crontab for automation:
# */10 * * * * source /etc/environment; source /usr/bin; source /home/jamtor/.bashrc; \
# source activate salmonTE; \
# python /share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/salmonTE/job_autosubmit.py

### 0. Define starting variables and environment###

import os

os.system('source /etc/profile')
os.system('source /etc/environment')
os.system('source /home/jamtor/.bashrc')

import pandas as pd
import glob
import re
import time

# make directory hierachy:
projectname = 'hgsoc_repeats'
sample_type = 'fullsamples/bowtell_primary'
exp_name = 'salmonTE'
total_no_jobs = 3
q = 'short'
upper_core_lim = 56
lower_core_lim = 46

# home_dir is where input/intermediate files are located:
home_dir = '/share/ScratchGeneral/jamtor/'
#home_dir = '/Users/jamestorpy/clusterHome/'
# home_dir2 is where scripts are located and final outputs will be stored:
home_dir2 = '/share/ClusterShare/thingamajigs/jamtor/'
#home_dir2 = '/Users/jamestorpy/clusterHome2/'

project_dir = home_dir + 'projects/' + projectname + '/RNA-seq/'
raw_dir = project_dir + 'raw_files/' + sample_type + '/fastq/'
results_dir = project_dir + 'results/'

print('')
print('This is the raw_dir:')
print(raw_dir)

os.makedirs(raw_dir, exist_ok=True)

# scripts/logs directories:
script_dir = home_dir2 + 'projects/hgsoc_repeats/RNA-seq/scripts/' + exp_name + '/'
log_dir = project_dir + 'logs/' + exp_name + '/'
old_log_dir = project_dir + 'logs/' + exp_name + '/old/'
os.makedirs(log_dir, exist_ok=True)
os.makedirs(old_log_dir, exist_ok=True)

print('')
print('The log_dir is:')
print(log_dir)
print('')
print('The script_dir is:')
print(script_dir)

# add indicator that crontab is running:
ind = open(script_dir + '/indicator.txt', 'w')
ind.write(str(time.strftime('%H:%M:%S')))

# define other output directories:
salmonTE_dir = results_dir + exp_name + '/'
os.makedirs(salmonTE_dir, exist_ok=True)

print('')
print('The salmonTE_dir is:')
print(salmonTE_dir)


### 1. Set up core number checkpoints ###

# save entire qstat output for 1st core number checkpoint:
os.system("qstat -u '*' | grep jamtor > " + script_dir + "/qstat.txt")

# save qstat job output as file:
os.system("qstat -u '*' | grep jamtor | grep -v QRLOGIN | awk '{print $3}' > " + script_dir \
    + "/qstat_jobs.txt")

# fetch job ids currently running or queued:
qstat_jobs = []
with open(script_dir + '/qstat_jobs.txt') as f:
    for line in f:
        qstat_jobs.append(
            re.sub(
                '\n', '', line
            )
        )
print('')
print('Jobs already running or queued:')
print(qstat_jobs)

# save qstat cores output as file:
os.system("qstat -u '*' | grep jamtor | grep -e 'q@' | grep -v QRLOGIN | awk '{print $9}' > " 
    + script_dir + "/qstat_cores.txt")

# fetch core number:
slots = 0
with open(script_dir + '/qstat_cores.txt') as f:
    for line in f:
        slots = slots + int(re.sub('\\n', '', line))

# if too many cores are being used, remove all but running jobs:
if slots > upper_core_lim:
    print('Too many cores being used, removing all except first 14 jobs')
    qstat = pd.read_csv(script_dir + '/qstat.txt', header = None)
    for line in qstat[14:].iterrows():
        j_id = str.split(str(line[1]))[1]
        print(j_id)
        os.system('qdel ' + j_id)

# create list of jobs currently running or queued in qstat:
running_jobs = []
if not os.stat(script_dir + '/qstat.txt').st_size == 0:
    index = 0
    qstat = pd.read_csv(script_dir + '/qstat.txt', header = None)
    for line in qstat.iterrows():
        if index == 0:
            running_jobs = [str.split(str(line[1]))[1]]
        else:
            running_jobs[index] = str.split(str(line[1]))[1]


if slots < lower_core_lim:
    print('')
    print('Less than ' + str(lower_core_lim) + 
    	' cores being used! Submitting more jobs...')
    
    ### 2. Check for any errors in the log files and delete output files
    # from erroneous jobs, rrecording the sample involved in redo_log.txt: ###
    
    sterms = ['error', 'ERROR', 'Error', 'halted', 'unexpected', 'Aborted', 'bad_alloc', 'Please check', 'failed']
    redo = []
    deleting_jobs = []
    for log in os.listdir(log_dir):
        if 'core' in log:
        	print('Removing ' + log)
        	os.system('rm ' + log_dir + '/' + log)
        elif '.o' in log:
            for term in sterms:
            	if os.path.isfile(log_dir + '/' + log):
            		try:
            			if term in open(log_dir + '/' + log).read():
            				redo_id = re.sub('.o[0-9].*', '', os.path.basename(log))
            				redo.append(redo_id)
            				os.system('mv ' + log_dir + '/' + log + ' ' + old_log_dir)
            				with open(script_dir + '/qstat.txt') as f:
            					for line in f:
            						if redo_id in line:
            							job_id = line.split()[0]
            							print("Removing " + redo_id + " with job ID " + job_id 
            								+ " from job queue as it's logfile contains an error")
            							os.system('qdel ' + job_id)
            							deleting_jobs.append(redo_id)
            		except ValueError:
            			print(log + " couldn't be read")

    redo = list(set(redo))

    print('')
    print('The following jobs will be redone:')

    if os.path.isfile(script_dir + '/redo_log.txt'):
        redo_log = open(script_dir + '/redo_log.txt', 'a')
    else:
        redo_log = open(script_dir + '/redo_log.txt', 'w')

    for r in redo:
        print(r)
        redo_log.write("%s\n" % r)


    ### 3. Check whether job has been previously completed ###

    if os.path.isfile(script_dir + '/done_samples.txt'):
        with open(script_dir + '/done_samples.txt') as f:
            done_ids = f.readlines()
        # remove whitespace characters like `\n` at the end of each line
        done_ids = [x.strip() for x in done_ids]
        print('')
        print(done_ids)
    else:
        done_ids = []


    ### 4. Define specific job parameters for each sample ###
    
    # read in raw files and fetch ids:
    #in_files = glob.glob(raw_dir + '*.bam')
    #print(in_files)
    #in_files.extend(glob.glob(raw_dir + 'FT[1-7].bam'))
	#in_files.extend(glob.glob(raw_dir + 'prPT1[2-7].bam'))
	#in_files.extend(glob.glob(raw_dir + 'rfPT[1-9].bam'))
	#in_files.extend(glob.glob(raw_dir + 'rcAF[1-9].bam'))
	#in_files.extend(glob.glob(raw_dir + 'arPT[1-9].bam'))
	#in_files.extend(glob.glob(raw_dir + 'msST[1-5].bam'))
	#in_files.extend(glob.glob(raw_dir + 'mrPT[1-8].bam'))
	#in_files.extend(glob.glob(raw_dir + 'erPT[1-8].bam'))
    in_files = [raw_dir + 'FT3_1.fastq.gz', raw_dir + 'FT4_1.fastq.gz', 
        raw_dir + 'FT5_1.fastq.gz', raw_dir + 'prPT1_1.fastq.gz',  
        raw_dir + 'prPT2_1.fastq.gz', raw_dir + 'prPT3_1.fastq.gz']
    #in_files = [raw_dir + '/FT1.bam', raw_dir + '/FT2.bam']
    #in_files = [raw_dir + '/erPT6.sub.bam']
    
for infile in in_files:
    if slots < lower_core_lim:
        print('')
        print('')
        print('The infile is ' + infile)
        u_id = re.sub(
            '.1.fastq.gz', '', os.path.basename(infile)
        )
        print('')
        print('The u_id is: ' + u_id)

        if u_id not in done_ids:


            ### Job 0 inputs ###
            # SalmonTE on all repeats #
            
            cores = ['4']
        
            inputs1 = [raw_dir + u_id + '_1.fastq.gz']
            inputs2 = [raw_dir + u_id + '_2.fastq.gz']
            inputs3 = ['none']
            inputs4 = ['none']
        
            outputs1 = [salmonTE_dir + u_id + '/all/EXPR.csv']
            outputs2 = ['none']
            outputs3 = ['none']
            outputs4 = ['none']
    
            extra_params = ['all']
            
            scripts = [script_dir + '1.salmonTE.bash']
            
            qsub_params = ['none']
    
            dependent_job = ['none']
    
            # define the minimum size the output files should be before
            # the job is accepted as complete
            min_output_size1 = [1000]
            #min_output_size1 = ['none']
            min_output_size2 = ['none']
            min_output_size3 = ['none']
            min_output_size4 = ['none']
            

            ### Job 1 inputs ###
            # SalmonTE on custom3 repeats #
                
            cores.append('4')
                
            inputs1.append(raw_dir + u_id + '_1.fastq.gz')
            inputs2.append(raw_dir + u_id + '_2.fastq.gz')
            inputs3.append('none')
            inputs4.append('none')
                
            outputs1.append(salmonTE_dir + u_id + '/custom3/EXPR.csv')
            outputs2.append('none')
            outputs3.append('none')
            outputs4.append('none')
        
            extra_params.append('custom3')
        
            scripts.append(script_dir + '1.salmonTE.bash')
                
            qsub_params.append('none')
        
            dependent_job.append('none')
        
            # define the minimum size the output files should be before
            # the job is accepted as complete
            min_output_size1.append(1000)
            min_output_size2.append('none')
            min_output_size3.append('none')
            min_output_size4.append('none')


            ### Job 2 inputs ###
            # SalmonTE on hg38 transcriptome #
                
            cores.append('4')
                
            inputs1.append(raw_dir + u_id + '_1.fastq.gz')
            inputs2.append(raw_dir + u_id + '_2.fastq.gz')
            inputs3.append('none')
            inputs4.append('none')
                
            outputs1.append(salmonTE_dir + u_id + '/hg38_transcriptome//EXPR.csv')
            outputs2.append('none')
            outputs3.append('none')
            outputs4.append('none')
        
            extra_params.append('hg38_transcriptome')
        
            scripts.append(script_dir + '2.salmon.bash')
                
            qsub_params.append('none')
        
            dependent_job.append('none')
        
            # define the minimum size the output files should be before
            # the job is accepted as complete
            min_output_size1.append('none')
            min_output_size2.append('none')
            min_output_size3.append('none')
            min_output_size4.append('none')
           

            for i in range(total_no_jobs):
                print('')
                print('')
                print('For job number ' + str(i) + ':')
    
                # set up job holding if necessary:
                if dependent_job[i] == 'none':
                    holdj = ' '
                else:
                    holdj = ' -hold_jid ' + dependent_job[i] + '_' + u_id
    
                print('')
                print('This is the hold string:')
                print(holdj)
                    
                # define extra_params:
                if extra_params:
                    if extra_params[i] == 'none':
                        extra_params[i] = ''
    
                # define qsub_params:
                if qsub_params:
                    if qsub_params[i] == 'none':
                        qsub_params[i] = ''
    
                # define dependent_job:
                if dependent_job:
                    if dependent_job[i] == 'none':
                        dependent_job[i] = ''
    
                
                inputs = [inputs1[i]]
                for input in [inputs2[i], inputs3[i], inputs4[i]]:
                    if input !='none':
                        inputs.append(input)
                
                outputs = [outputs1[i]]
                for output in [outputs2[i], outputs3[i], outputs4[i]]:
                    if output !='none':
                        outputs.append(output)
    
                min_output_sizes = [min_output_size1[i]]
                for size in [min_output_size2[i], min_output_size3[i], min_output_size4[i]]:
                        min_output_sizes.append(size)
    
                print('')
                print('min_output_sizes are:')
                print(min_output_sizes)
    
                # determine whether output files are min file size required:
                min_size_pass = []
                for j, out in enumerate(outputs):
                    if os.path.isfile(out):
                        if min_output_sizes[j] != 'none':
                            min_size_pass.append(os.path.getsize(out) > min_output_sizes[j])
    
                # if all outputs are minimum size, record size_pass as True, if not, False:
                size_pass = all(min_size_pass)
                print('')
                print('sizepass is ' + str(size_pass))
                
                # add u_id to log_dir and create:
                u_log_dir = log_dir + u_id
                os.makedirs(u_log_dir, exist_ok=True)

    
                # define qsub command for jobs:
                q_command = '/opt/gridengine/bin/linux-x64/qsub ' \
                          + '-q ' + q + '.q ' \
                          + '-N j' + str(i) + '_' + u_id \
                          + holdj \
                          + ' -b y -wd ' \
                          + u_log_dir \
                          + ' -j y -R y ' \
                          + '-P DSGClinicalGenomics -pe smp ' + cores[i] + ' ' \
                          + qsub_params[i] \
                          + ' -V ' + scripts[i] + ' ' \
                          + cores[i] + ' ' \
                          + ' '.join(inputs) + ' ' \
                          + ' '.join(outputs) + ' ' \
                          + extra_params[i]
    

                ### 5. Run jobs if they meet the conditions required ###
    
                print('')
                if not all([os.path.isfile(inp) for inp in inputs]):
                    print('At least one job ' + str(i) + ' input file does not exist, holding job')
    
                elif 'j' + str(i) + '_' + u_id in redo and not 'j' + str(i) + '_' + u_id in qstat_jobs:

                    # count how many times this job has been redone before:
                    redo_count = 0

                    redo_file = open(script_dir + '/redo_log.txt', 'r')
                    for line in redo_file:
                        if 'j' + str(i) + '_' + u_id in line:
                            redo_count +=1

                    print('Job has been redone ' + str(redo_count) + ' times before')

                    cores[i] = str(int(cores[i]) + (2 * (redo_count + 1)))
                    print('Job previously encountered an error!')
                    print('Submitting job ' + str(i) + ' with input/s or in dir/s and ' \
                    	+ str(cores[i]) + ' cores:')
                    print(inputs)
                    print('')
                    print('and output/s or out dir/s:')
                    print(outputs)

                    # redefine qsub command for jobs:
                    q_command = '/opt/gridengine/bin/linux-x64/qsub ' \
                          + '-q ' + q + '.q ' \
                          + '-N j' + str(i) + '_' + u_id \
                          + holdj \
                          + ' -b y -wd ' \
                          + log_dir \
                          + ' -j y -R y ' \
                          + '-P DSGClinicalGenomics -pe smp ' + cores[i] + ' ' \
                          + qsub_params[i] \
                          + ' -V ' + scripts[i] + ' ' \
                          + cores[i] + ' ' \
                          + ' '.join(inputs) + ' ' \
                          + ' '.join(outputs) + ' ' \
                          + extra_params[i]
    
                    #print(q_command)
                    os.system(q_command)

                elif 'j' + str(i) + '_' + u_id in redo and 'j' + str(i) + '_' \
                	+ u_id in deleting_jobs:
                    
                    # count how many times this job has been redone before:
                    redo_count = 0

                    redo_file = open(script_dir + '/redo_log.txt', 'r')
                    for line in redo_file:
                        if 'j' + str(i) + '_' + u_id in line:
                            redo_count +=1

                    print('Job has been redone ' + str(redo_count) + ' times before')

                    cores[i] = str(int(cores[i]) + (2 * (redo_count + 1)))
                    print('Job previously encountered an error!')
                    print('Submitting job ' + str(i) + ' with input/s or in dir/s and ' \
                    	+ str(cores[i]) + ' cores:')
                    print(inputs)
                    print('')
                    print('and output/s or out dir/s:')
                    print(outputs)

                    # redefine qsub command for jobs:
                    q_command = '/opt/gridengine/bin/linux-x64/qsub ' \
                          + '-q ' + q + '.q ' \
                          + '-N j' + str(i) + '_' + u_id \
                          + holdj \
                          + ' -b y -wd ' \
                          + log_dir \
                          + ' -j y -R y ' \
                          + '-P DSGClinicalGenomics -pe smp ' + cores[i] + ' ' \
                          + qsub_params[i] \
                          + ' -V ' + scripts[i] + ' ' \
                          + cores[i] + ' ' \
                          + ' '.join(inputs) + ' ' \
                          + ' '.join(outputs) + ' ' \
                          + extra_params[i]
    
                    #print(q_command)
                    os.system(q_command)

                elif all([os.path.isfile(inp) for inp in inputs]) and not all([os.path.isfile(outp) \
                	for outp in outputs]) and not 'j' + str(i) + '_' + u_id in qstat_jobs:
                    print('Submitting job ' + str(i) + ' with input/s or in dir/s:')
                    print(inputs)
                    print('')
                    print('and output/s or out dir/s:')
                    print(outputs)

                    #print(q_command)
                    os.system(q_command)
     
                # run job if outputs exist but are not the minimum file size:
                elif all([os.path.isfile(outp) for outp in outputs]) and size_pass == False and not \
                    'j' + str(i) + '_' + u_id in qstat_jobs:
                    print('Job outputs are not the minimum size required!')
                    print('Submitting job ' + str(i) + ' with input/s or in dir/s and ' + str(cores[i]) + ' cores:')
                    print(inputs)
                    print('')
                    print('and output/s or out dir/s:')
                    print(outputs)
    
                    #print(q_command)
                    os.system(q_command)
    
                elif all([os.path.isfile(outp) for outp in outputs]) and size_pass == True:
                    print('Job ' + str(i) + ' has previously been completed, no need to run')
                
                elif 'j' + str(i) + '_' + u_id in qstat_jobs and 'j' + str(i) + '_' + u_id not in deleting_jobs:
                    os.system('echo Job ' + str(i) + ' is currently running or queued, no need to run again')
