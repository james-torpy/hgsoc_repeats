#!/home/jamtor/local/bin/Python-3.6.3/bin/python3

### 1.bamtofastq_autosubmit.py ###

# This script detects which output files still must be created from the repeatsPipeline
# and submits jobs for them #

# run once and then input the following command into crontab for automation:
# */10 * * * * source /etc/environment; source /usr/bin; source /home/jamtor/.bashrc; \
# /home/jamtor/local/bin/Python-3.6.3/bin/python3.6 \
# /share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/rsem_no_mmappers/job_autosubmit.py

### 0. Define starting variables ###

import os
import pandas as pd
import glob
import re
import time

os.system('source /etc/profile')
os.system('source /etc/environment')
os.system('source /home/jamtor/.bashrc')
os.system('module load /share/ClusterShare/Modules/modulefiles/contrib/gi/zlib/1.2.8')
os.system('module load /share/ClusterShare/Modules/modulefiles/contrib/phuluu/samtools/1.4')
os.system('module load /share/ClusterShare/Modules/modulefiles/contrib/gi/boost/1.53.0')

# make directory hierachy:
projectname = 'hgsoc_repeats'
sample_type = 'fullsamples/bowtell_primary'
exp_name = 'rsem_no_mmappers'
total_no_jobs = 2
q = 'long'
upper_core_lim = 200
lower_core_lim = 200

# home_dir is where input/intermediate files are located:
home_dir = '/share/ScratchGeneral/jamtor/'
#home_dir = '/Users/jamestorpy/clusterHome/'
# home_dir2 is where scripts are located and final outputs will be stored:
home_dir2 = '/share/ClusterShare/thingamajigs/jamtor'
#home_dir2 = '/Users/jamestorpy/clusterHome2/'

project_dir = home_dir + '/projects/' + projectname + '/RNA-seq/'
raw_dir = project_dir + '/raw_files/' + sample_type + '/'
results_dir = project_dir + '/results/'

print('')
print('This is the raw_dir:')
print(raw_dir)

# scripts/logs directories:
script_dir = home_dir2 + '/projects/hgsoc_repeats/RNA-seq/scripts/' + exp_name + '/'
log_dir = project_dir + '/logs/' + exp_name + '/'
old_log_dir = project_dir + '/logs/' + exp_name + '/old/'
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
fastq_dir=raw_dir + '/fastq/'
os.makedirs(fastq_dir, exist_ok=True)

ribo_dir = results_dir + '/star/ribo/' + exp_name + '/'
os.makedirs(ribo_dir, exist_ok=True)

gc_dir = results_dir + '/star/GC/' + exp_name + '/'
os.makedirs(gc_dir, exist_ok=True)

rsem_dir = results_dir + '/rsem/' + exp_name + '/'
os.makedirs(rsem_dir, exist_ok=True)

print('')
print('The fastq_dir is:')
print(fastq_dir)
print('')
print('The ribo_dir is:')
print(ribo_dir)
print('')
print('The gc_dir is:')
print(gc_dir)
print('')
print('The rsem_dir is:')
print(rsem_dir)


### 1. Set up core number checkpoints ###

# save entire qstat output for 1st core number checkpoint:
os.system("qstat -u '*' | grep jamtor > " + script_dir + "/qstat.txt")

# save qstat job output as file:
os.system("qstat -u '*' | grep jamtor | grep -v QRLOGIN | awk '{print $3}' > " + script_dir 
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
#if slots > upper_core_lim:
#    print('Too many cores being used, removing all except first 8 jobs')
#    qstat = pd.read_csv(script_dir + '/qstat.txt', header = None)
#    for line in qstat[8:].iterrows():
#        j_id = str.split(str(line[1]))[1]
#        print(j_id)
#        os.system('qdel ' + j_id)

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
    
    sterms = ['error', 'ERROR', 'Error', 'halted', 'unexpected', 'Aborted', 'bad_alloc', 'Please check']
    redo = []
    deleting_jobs = []
    for log in os.listdir(log_dir):
        if '.o' in log:
            for term in sterms:
                if os.path.isfile(log_dir + '/' + log):
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
    #in_files = [raw_dir + '/erPT6.bam', raw_dir + '/FT3.bam']
    #in_files = [raw_dir + '/rcAF14.bam', raw_dir + '/arPT8.bam']
    in_files = [raw_dir + '/erPT6.s.bam']
    
for infile in in_files:
    if slots < lower_core_lim:
        print('')
        print('')
        print('The infile is ' + infile)
        u_id = re.sub(
            '.bam', '', os.path.basename(infile)
        )
        print('')
        print('The u_id is: ' + u_id)

        if u_id not in done_ids:


            ### Job 0 inputs ###
            # alignment of fastqs to gc transcriptome reference #
            
            cores = ['7']
            
            inputs1 = [fastq_dir + '/' + u_id + '.1.fastq.gz']
            inputs2 = [fastq_dir + '/' + u_id + '.2.fastq.gz']
            inputs3 = ['none']
            inputs4 = ['none']
            
            outputs1 = [gc_dir + '/' + u_id + '/Log.final.out']
            outputs2 = [gc_dir + '/' + u_id + '/Aligned.toTranscriptome.out.bam']
            outputs3 = [gc_dir + '/' + u_id + '/SJ.out.tab']
            outputs4 = ['none']
    
            extra_params = ['none']
    
            scripts = [script_dir + '/1.star_ercc.bash']
            
            qsub_params = ['none']
    
            dependent_job = ['none']
    
            # define the minimum size the output files should be before
            # the job is accepted as complete
            min_output_size1 = ['none']
            min_output_size2 = ['none']
            min_output_size3 = ['none']
            min_output_size4 = ['none']
    
    
            ### Job 1 inputs ###
            # count gc transcripts from bam using rsem #
            
            cores.append('7')
            
            inputs1.append(gc_dir + '/' + u_id + '/Aligned.toTranscriptome.out.bam')
            inputs2.append('none')
            inputs3.append('none')
            inputs4.append('none')
            
            outputs1.append(rsem_dir + '/' + u_id + '/' + u_id + '.genes.results')
            outputs2.append('none')
            outputs3.append('none')
            outputs4.append('none')
    
            extra_params.append('none')
    
            scripts.append(script_dir + '/2.rsem_ercc.bash')
            
            qsub_params.append('none')
    
            dependent_job.append('j0')
    
            # define the minimum size the output files should be before
            # the job is accepted as complete
            min_output_size1.append(2000)
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
    
    
                # define qsub command for jobs:
                q_command = '/opt/gridengine/bin/linux-x64/qsub ' \
                          + '-q ' + q + '.q ' \
                          + '-N j' + str(i) + '_' + u_id \
                          + holdj \
                          + ' -b y -wd ' \
                          + log_dir \
                          + ' -j y -R y ' \
                          + '-pe smp ' + cores[i] + ' ' \
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
                    print('Submitting job ' + str(i) + ' with input/s or in dir/s and ' 
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
                          + '-pe smp ' + cores[i] + ' ' \
                          + qsub_params[i] \
                          + ' -V ' + scripts[i] + ' ' \
                          + cores[i] + ' ' \
                          + ' '.join(inputs) + ' ' \
                          + ' '.join(outputs) + ' ' \
                          + extra_params[i]
    
                    print(q_command)
                    #os.system(q_command)

                elif 'j' + str(i) + '_' + u_id in redo and 'j' + str(i) + '_' 
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
                    print('Submitting job ' + str(i) + ' with input/s or in dir/s and ' 
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
                          + '-pe smp ' + cores[i] + ' ' \
                          + qsub_params[i] \
                          + ' -V ' + scripts[i] + ' ' \
                          + cores[i] + ' ' \
                          + ' '.join(inputs) + ' ' \
                          + ' '.join(outputs) + ' ' \
                          + extra_params[i]
    
                    print(q_command)
                    #os.system(q_command)

                elif all([os.path.isfile(inp) for inp in inputs]) and not all([os.path.isfile(outp) \
                	for outp in outputs]) and not 'j' + str(i) + '_' + u_id in qstat_jobs:
                    print('Submitting job ' + str(i) + ' with input/s or in dir/s:')
                    print(inputs)
                    print('')
                    print('and output/s or out dir/s:')
                    print(outputs)
    
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
    
                    os.system(q_command)
    
                elif all([os.path.isfile(outp) for outp in outputs]) and size_pass == True:
                    print('Job ' + str(i) + ' has previously been completed, no need to run')
                
                elif 'j' + str(i) + '_' + u_id in qstat_jobs and 'j' + str(i) + '_' + u_id not in deleting_jobs:
                    os.system('echo Job ' + str(i) + ' is currently running or queued, no need to run again')
    
    
            ### 6. Remove files no longer needed ###
            rm_cores = '1'
            
            rm1 = raw_dir + '/' + u_id + '.bam'
            rm2 = gc_dir + '/' + u_id + '/Aligned.sortedByCoord.out.bam'
            rm3 = rsem_dir + '/' + u_id + '/' + u_id + '.transcript.bam'
    
            rm_script = (script_dir + '/3.remove_files.bash')
    
            rm_dependents = ['j0', 'j1', 'j2']
    
    
            rm_files = [rm1, rm2, rm3]
    
            holdj = '-hold_jid '
            for d in rm_dependents:
                holdj = holdj + d + '_' + u_id + ','
            
            # double check all outputs have been created:
            final_check = []
            for n, op in enumerate(outputs):
                print(op)
                final_check.append(os.path.isfile(op))
            print(final_check)

            final_size_check = []
            for n, op in enumerate(outputs):
                if os.path.isfile(op):
                    final_size_check.append(os.path.getsize(op) > 2000)
            print(final_size_check)

            if all(final_check):
                if all(final_size_check):
                    for k, rm in enumerate(rm_files):
                        if not 'r' + str(k) + '_' + u_id in qstat_jobs:
                            print('')
                            print('Removing ' + rm)
                            if os.path.isfile(rm):
                                os.system('/opt/gridengine/bin/linux-x64/qsub ' \
                                      + '-q ' + q + '.q ' \
                                      + '-N r' + str(k) + '_' + u_id + ' ' \
                                      + holdj \
                                      + ' -b y -wd ' \
                                      + log_dir \
                                      + ' -j y -R y ' \
                                      + '-pe smp ' + rm_cores
                                      + ' -V ' + rm_script + ' ' \
                                      + rm + ' ' + u_id + ' ' + script_dir
                                    )
