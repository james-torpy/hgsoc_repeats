#!/home/jamtor/local/bin/Python-3.6.3/bin/python3

### 1.bamtofastq_autosubmit.py ###

# This script detects which output files still must be created from the repeatsPipeline
# and submits jobs for them #

# run once and then input the following command into crontab for automation:
# */1 * * * * source /etc/environment; source /usr/bin; source /home/jamtor/.bashrc; \
# /home/jamtor/local/bin/Python-3.6.3/bin/python3.6 \
# /share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/exp9_test2/job_autosubmit.py


### 0. Define starting variables ###

import os
import pandas as pd
import glob
import re
import time

os.system('source /home/jamtor/.bashrc')
os.system('module load gi/zlib/1.2.8')
os.system('module load phuluu/samtools/1.4')

# make directory hierachy:
projectname = 'hgsoc_repeats'
sample_type = 'fullsamples/bowtell_primary'
exp_name = 'exp9_test2'
total_no_jobs = 4
q = 'long'

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
script_dir = home_dir2 + '/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/' + exp_name + '/'
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

print('')
print('The fastq_dir is:')
print(fastq_dir)
print('')
print('The ribo_dir is:')
print(ribo_dir)
print('')
print('The gc_dir is:')
print(gc_dir)


### 1. Set up core number checkpoints ###

# save entire qstat output for 1st core number checkpoint:
os.system("qstat -u '*' | grep jamtor > " + script_dir + "/qstat.txt")

# save qstat job output as file:
os.system("qstat -u '*' | grep jamtor | grep -v QRLOGIN | awk '{print $3}' > " + script_dir + "/qstat_jobs.txt")

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
os.system("qstat -u '*' | grep jamtor | grep -e 'q@' | grep -v QRLOGIN | awk '{print $9}' > " + script_dir + "/qstat_cores.txt")

# fetch core number:
slots = 0
with open(script_dir + '/qstat_cores.txt') as f:
    for line in f:
        slots = slots + int(re.sub('\\n', '', line))

# if too many cores are being used, remove all but running jobs:
if slots > 100:
    print('Too many cores being used, removing all except first 10 jobs')
    qstat = pd.read_csv(script_dir + '/qstat.txt', header = None)
    for line in qstat.iterrows():
        next(line)
        next(line)
        next(line)
        next(line)
        next(line)
        next(line)
        next(line)
        next(line)
        next(line)
        next(line)
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


if slots < 80 and slots <100:
    print('')
    print('Less than 80 cores being used! Submitting more jobs...')

    ### 2. Check for any errors in the log files and delete output files
    # from erroneous jobs: ###

    sterms = ['error', 'ERROR', 'Error', 'halted', 'unexpected']
    redo = []
    for log in os.listdir(log_dir):
        if '.o' in log:
            for term in sterms:
                if os.path.isfile(log_dir + '/' + log):
                    if term in open(log_dir + '/' + log).read():
                        redo.append(
                            re.sub(
                                '.o[0-9].*', '', os.path.basename(log)	
                            )
                        )
                        os.system('mv ' + log_dir + '/' + log + ' ' + old_log_dir)
    redo = list(set(redo))
    print('The following jobs will be redone:')
    for r in redo:
        print(r)


    ### 3. Define specific job parameters for each sample ###
    
    # read in raw files and fetch ids:
    #in_files = glob.glob(raw_dir + '*.bam')
    in_files = [raw_dir + '/prPT9.bam', raw_dir \
    + '/prPT11.bam', raw_dir + '/rfPT3.bam', raw_dir \
    + '/rfPT1.bam', raw_dir + '/prPT10.bam']
    #in_files = [raw_dir + '/prPT9.bam']
    
    for infile in in_files:
        print('')
        print('')
        print('The infile is ' + infile)
        u_id = re.sub(
            '.bam', '', os.path.basename(infile)
        )
        print('')
        print('The u_id is: ' + u_id)
    
        
        ### Job 0 inputs ###
        # fastq to bam conversion #
        
        cores = ['2']
    
        inputs1 = [raw_dir + '/' + u_id + '.bam']
        inputs2 = ['none']
        inputs3 = ['none']
        inputs4 = ['none']
    
        outputs1 = [fastq_dir + '/' + u_id + '.1.fastq.gz']
        outputs2 = [fastq_dir + '/' + u_id + '.2.fastq.gz']
        outputs3 = [fastq_dir + '/' + u_id + '.unpaired.1.fastq.gz']
        outputs4 = [fastq_dir + '/' + u_id + '.unpaired.2.fastq.gz']
        outputs5 = ['none']

        extra_params = ['none']
        
        scripts = [script_dir + '/0.bamtofastq.bash']
        
        qsub_params = ['none']

        dependent_job = ['none']
        
        
        ### Job 1 inputs ###
        # alignment of fastqs to ribosomal reference #
        
        cores.append('8')
        
        inputs1.append(fastq_dir + '/' + u_id + '.1.fastq.gz')
        inputs2.append(fastq_dir + '/' + u_id + '.2.fastq.gz')
        inputs3.append('none')
        inputs4.append('none')
        
        outputs1.append(ribo_dir + '/' + u_id + '/Log.final.out')
        outputs2.append('none')
        outputs3.append('none')
        outputs4.append('none')
        outputs5.append('none')

        extra_params.append('none')
        
        scripts.append(script_dir + '/1.starRibo.bash')
        
        qsub_params.append('none')

        dependent_job.append('j0')


        ### Job 2 inputs ###
        # alignment of fastqs to gc reference #
        
        cores.append('6')
        
        inputs1.append(fastq_dir + '/' + u_id + '.1.fastq.gz')
        inputs2.append(fastq_dir + '/' + u_id + '.2.fastq.gz')
        inputs3.append('none')
        inputs4.append('none')
        
        outputs1.append(gc_dir + '/' + u_id + '/Log.final.out')
        outputs2.append(gc_dir + '/' + u_id + '/Aligned.out.bam')
        outputs3.append(gc_dir + '/' + u_id + '/Chimeric.out.junction')
        outputs4.append(gc_dir + '/' + u_id + '/Chimeric.out.sam')
        outputs5.append(gc_dir + '/' + u_id + '/SJ.out.tab')

        extra_params.append('none')

        scripts.append(script_dir + '/2.starGC.bash')
        
        qsub_params.append('none')

        dependent_job.append('j0')


        ### Job 3 inputs ###
        # sort bams by name #
        
        cores.append('6')
        
        inputs1.append(gc_dir + '/' + u_id + '/Aligned.out.bam')
        inputs2.append('none')
        inputs3.append('none')
        inputs4.append('none')
        
        outputs1.append(gc_dir + '/' + u_id + '/Aligned.novosortedByName.out.bam')
        outputs2.append('none')
        outputs3.append('none')
        outputs4.append('none')
        outputs5.append('none')

        extra_params.append('22G')

        scripts.append(script_dir + '/3.novosortName.bash')
        
        qsub_params.append('none')

        dependent_job.append('j2')        

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
        
            print('')
            if not all([os.path.isfile(inp) for inp in inputs]):
                print('At least one job ' + str(i) + ' input file does not exist, holding job')

            elif all([os.path.isfile(inp) for inp in inputs]) and not all([os.path.isfile(outp) \
            	for outp in outputs]) and not 'j' + str(i) + '_' + u_id in qstat_jobs:
                print('Submitting job ' + str(i) + ' with input/s or in dir/s:')
                print(inputs)
                print('')
                print('and output/s or out dir/s:')
                print(outputs)

                os.system('/opt/gridengine/bin/linux-x64/qsub ' \
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
                    )

            elif 'j' + str(i) + '_' + u_id in redo and not 'j' + str(i) + '_' + u_id in qstat_jobs:
                cores[i] = int(cores[i]) + 2
                print('Submitting job ' + str(i) + ' with input/s or in dir/s and ' + str(cores[i]) + ' cores:')
                print(inputs)
                print('')
                print('and output/s or out dir/s:')
                print(outputs)

                os.system('/opt/gridengine/bin/linux-x64/qsub ' \
                      + '-q ' + q + '.q ' \
                      + '-N j' + str(i) + '_' + u_id \
                      + holdj \
                      + ' -b y -wd ' \
                      + log_dir \
                      + ' -j y -R y ' \
                      + '-P DSGClinicalGenomics -pe smp ' + str(cores[i]) + ' ' \
                      + qsub_params[i] \
                      + ' -V ' + scripts[i] + ' ' \
                      + str(cores[i]) + ' ' \
                      + ' '.join(inputs) + ' ' \
                      + ' '.join(outputs) + ' ' \
                      + extra_params[i]
                    )

            elif all([os.path.isfile(outp) for outp in outputs]):
                print('Job ' + str(i) + ' has previously been completed, no need to run')
            
            elif 'j' + str(i) + '_' + u_id in qstat_jobs:
                os.system('echo Job ' + str(i) + ' is currently running or queued, no need to run again')