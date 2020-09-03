#!/bin/bash
#!/bin/bash

################
# Setting slurm options
################
#SBATCH --partition medium

#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu 512
#SBATCH --time 10:00:00
#SBATCH --job-name "gamma"

#SBATCH --array=0-9
################
# parameter array
# arr=(p1 p2 ...), or arr[0]=p1 if set individual
# echo ${par[i]} returns i'th element
################
par_array=(0.018 0.019 0.020 0.021 0.022 0.023 0.024 0.025 0.026 0.027)

################
# run the program
################
module load anaconda3
module load python3
srun python runcard.py check_vertex_opening 18 0.02 10000 0.1 ${par_array[${SLURM_ARRAY_TASK_ID}]}
#$scratch/a.out ${SLURM_ARRAY_TASK_ID} ${par_array[${SLURM_ARRAY_TASK_ID}]}


exit 0
