#!/bin/bash

### jobSubmit.bash ###

# Takes qsub command from job_autosubmit.py and runs #

# set up bash environment:
source /home/jamtor/.bashrc
echo $SGE_ROOT

date

echo -e
echo "The submit_line is:"
echo $submit_line

$submit_line