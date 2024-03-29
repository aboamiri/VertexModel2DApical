#!/bin/bash

################
#
# Setting slurm options
#
################

# lines starting with "#SBATCH" define your jobs parameters

# requesting the type of node on which to run job
#SBATCH --partition medium

# telling slurm how many instances of this job to spawn (typically 1)
#SBATCH --ntasks 1

# setting number of CPUs per task (1 for serial jobs)
#SBATCH --cpus-per-task 1

# setting memory requirements
#SBATCH --mem-per-cpu 2048

# propagating max time for job to run
#SBATCH --time 30:00:00

# Setting the name for the job
#SBATCH --job-name edgecol

# setting notifications for job
# accepted values are ALL, BEGIN, END, FAIL, REQUEUE
#SBATCH --mail-type FAIL

# telling slurm where to write output and error
#SBATCH --output job-%j.out
#SBATCH --error job-%j.out

################
#
# copying your data to /scratch
#
################

# create local folder on ComputeNode
#scratch=/scratch/$USER/$SLURM_JOB_ID
#mkdir -p $scratch

# copy all your NEEDED data to ComputeNode
#cp /home/aamiri/Documents/VertexModel/y2019/overlap $scratch
#cd $scratch

# dont access /home after this line

# if needed load modules here
#module load <module_name>
#module load anaconda3
module load python3  

# if needed add export variables here

################
#
# run the program
#
################
srun python runcard.py edge_collapse 14 0.02 2000 0.1 0.02
#$scratch/ /home/aamiri/Documents/VertexModel/y2019/overlap

# copy results to data
#cp -r . /data/<my_dir>
#cd

# clean up scratch
#rm -rf $scratch
#unset scratch

exit 0
