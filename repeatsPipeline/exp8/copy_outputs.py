#!/share/ClusterShare/software/contrib/fabbus/python/2.7.3/bin/python

### copy_outputs.py ###

# This script copies output files from repeatsPipeline which are in ClusterShare
# directories but not ScratchGeneral directories

### 0. Define starting variables ###

import os
from path import path

exp_name = "exp8"

home_dir = '/share/ClusterShare/thingamajigs/jamtor/'
#home_dir = '/Users/jamestorpy/clusterHome2'
project_dir = home_dir + '/projects/hgsoc_repeats/RNA-seq/'
results_dir = project_dir + '/results/'
gc_dir = results_dir + '/star/GC/' + exp_name + '/'
htseq_dir = results_dir + '/htseq/' + exp_name + '/'
ribo_dir = results_dir + '/star/ribo/' + exp_name + '/'

dest_home = '/share/ScratchGeneral/jamtor/'
#dest_home = '/Users/jamestorpy/clusterHome'
dest_project = dest_home + '/projects/hgsoc_repeats/RNA-seq/'
dest_results = dest_project + '/results/'
dest_gc = dest_results + '/star/GC/' + exp_name + '/'
dest_htseq = dest_results + '/htseq/' + exp_name + '/'
dest_ribo = dest_results + '/star/ribo/' + exp_name + '/'


### 1. Fetch names of gc STAR dirs in ClusterShare and check ScratchGeneral for
# these ###

d = path(gc_dir)
for f in d.walkfiles('Log.final.out'):
    print('')
    print('Log file is ' + f)

    temp = str.split(f, '/')
    uID = temp[len(temp)-2]
    print('')
    print('uID is ' + uID)

    if os.path.exists(dest_gc + '/' + uID):
        if os.path.exists(dest_gc + '/' + uID + '/Log.final.out'):
            print('')
            print('Log.final.out exists for '+ uID + ', assuming all files are already copied')
        else:
            print('')
            print('Copying ' + f + ' to ' + dest_gc + '/' + uID)
            print('')
            print('cp ' + f + ' ' + dest_gc + '/' + uID)
            os.system('cp ' + f + ' ' + dest_gc + '/' + uID)
            print('')
            print('Copying ' + gc_dir + '/' + uID + '/Chimeric* to ' + dest_gc + '/' + uID)
            print('')
            print('cp ' + gc_dir + '/' + uID + '/Chimeric* ' + dest_gc + '/' + uID)
            os.system('cp ' + gc_dir + '/' + uID + '/Chimeric* ' + dest_gc + '/' + uID)
            print('')
            print('Copying ' + gc_dir + '/' + uID + '/SJ.out.tab to ' + dest_gc + '/' + uID)
            print('')
            print('cp ' + gc_dir + '/' + uID + '/SJ.out.tab ' + dest_gc + '/' + uID)
            os.system('cp ' + gc_dir + '/' + uID + '/SJ.out.tab ' + dest_gc + '/' + uID)
    else:
        print('')
        print('Creating directory for ' + uID)
        os.system('mkdir -p ' + dest_gc + '/' + uID)
        print('')
        print('Copying ' + f + ' to ' + dest_gc + '/' + uID)
        print('')
        print('cp ' + f + ' ' + dest_gc + '/' + uID)
        os.system('cp ' + f + ' ' + dest_gc + '/' + uID)
        print('')
        print('Copying ' + gc_dir + '/' + uID + '/Chimeric* to ' + dest_gc + '/' + uID)
        print('')
        print('cp ' + gc_dir + '/' + uID + '/Chimeric* ' + dest_gc + '/' + uID)
        os.system('cp ' + gc_dir + '/' + uID + '/Chimeric* ' + dest_gc + '/' + uID)
        print('')
        print('Copying ' + gc_dir + '/' + uID + '/SJ.out.tab to ' + dest_gc + '/' + uID)
        print('')
        print('cp ' + gc_dir + '/' + uID + '/SJ.out.tab ' + dest_gc + '/' + uID)
        os.system('cp ' + gc_dir + '/' + uID + '/SJ.out.tab ' + dest_gc + '/' + uID)

    print('')

    
### 2. Fetch names of htseq files in ClusterShare and check ScratchGeneral for
# these ###

d = path(htseq_dir)
for f in d.walkfiles('*.gc.htseq.txt'):
    statinfo1 = os.stat(f)
    if statinfo1.st_size > 0:
        print('')
        print('htseq file is ' + f)

        temp = str.split(f, '/')
        uID = temp[len(temp)-2]
        print('uID is ' + uID)
        for t in ['gc', 'custom3', 'all']:
            if os.path.exists(dest_htseq + '/' + uID):
                if os.path.exists(dest_htseq + '/' + uID + '/' + uID + '.' + t + '.htseq.txt'):
                    statinfo2 = os.stat(dest_htseq + '/' + uID + '/' + uID + '.' + t + '.htseq.txt')
                    if statinfo2.st_size > 0:
                        print('')
                        print(t + '.htseq.txt exists for ' + uID + ', file is already copied')
                    else:
                        print('')
                        print('Copying ' + htseq_dir + '/' + uID + '/' + uID + '.' + t + '.htseq.txt' + ' to ' + dest_htseq + '/' + uID)
                        print('')
                        print('cp ' + htseq_dir + '/' + uID + '/' + uID + '.' + t + '.htseq.txt ' + dest_htseq + '/' + uID)
                        os.system('cp ' + htseq_dir + '/' + uID + '/' + uID + '.' + t + '.htseq.txt ' + dest_htseq + '/' + uID)
                else:
                    print('')
                    print('Copying ' + htseq_dir + '/' + uID + '/' + uID + '.' + t + '.htseq.txt' + ' to ' + dest_htseq + '/' + uID)
                    print('')
                    print('cp ' + htseq_dir + '/' + uID + '/' + uID + '.' + t + '.htseq.txt ' + dest_htseq + '/' + uID)
                    os.system('cp ' + htseq_dir + '/' + uID + '/' + uID + '.' + t + '.htseq.txt ' + dest_htseq + '/' + uID)
            else:
                print('')
                print('Creating directory for ' + uID)
                os.system('mkdir -p -p ' + dest_htseq + '/' + uID)
                print('')
                print('Copying ' + htseq_dir + '/' + uID + '/' + uID + '.' + t + '.htseq.txt' + ' to ' + dest_htseq + '/' + uID)
                print('')
                print('cp ' + htseq_dir + '/' + uID + '/' + uID + '.' + t + '.htseq.txt ' + dest_htseq + '/' + uID)
                os.system('cp ' + htseq_dir + '/' + uID + '/' + uID + '.' + t + '.htseq.txt ' + dest_htseq + '/' + uID)
    print('')
    
    
### 3. Fetch names of ribo file in ClusterShare and check ScratchGeneral for
# this ###
 
d = path(ribo_dir)
for f in d.walkfiles('Log.final.out'):
    print('')
    print('Log file is ' + f)

    temp = str.split(f, '/')
    uID = temp[len(temp)-2]
    print('')
    print('uID is ' + uID)

    if os.path.exists(dest_ribo + '/' + uID):
        if os.path.exists(dest_ribo + '/' + uID + '/Log.final.out'):
            print('')
            print('Log.final.out exists for '+ uID + ', file is already copied')
        else:
            print('')
            print('Copying ' + f + ' to ' + dest_ribo + '/' + uID)
            print('')
            print('cp ' + f + ' ' + dest_ribo + '/' + uID)
            os.system('cp ' + f + ' ' + dest_ribo + '/' + uID)
    else:
        print('')
        print('Creating directory for ' + uID)
        os.system('mkdir -p ' + dest_ribo + '/' + uID)
        print('')
        print('Copying ' + f + ' to ' + dest_ribo + '/' + uID)
        print('')
        print('cp ' + f + ' ' + dest_ribo + '/' + uID)
        os.system('cp ' + f + ' ' + dest_ribo + '/' + uID)

    print('')

