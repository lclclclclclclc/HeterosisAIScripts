#!/bin/bash

# This is an example batch script for slurm on Oscar
# 
# The commands for slurm start with #SBATCH
# All slurm commands need to come before the program 
# you want to run.  In this example, 'echo "Hello World!"
# is the command we are running.
#
# This is a bash script, so any line that starts with # is
# a comment.  If you need to comment out an #SBATCH line 
# use ##SBATCH 
#
# To submit this script to slurm do:
#    sbatch batch.script
#
# Once the job starts you will see a file MySerialJob-****.out
# The **** will be the slurm JobID

# --- Start of slurm commands -----------

# Request an hour of runtime:
#SBATCH --time=2:00:00

# Default resources are 1 core with 2.8GB of memory.
# Use more memory (4GB):
##SBATCH --mem=4G
##SBATCH -c 10

# Specify a job name:
#SBATCH -J heterosis

# Specify an output file
# %j is a special variable that is replaced by the JobID when 
# job starts
#SBATCH -o het-%j.out
#SBATCH -e het-%j.out

#----- End of slurm commands ----

# Run a command
module load python/3.7.4 
module load slim/3.3

python3 
