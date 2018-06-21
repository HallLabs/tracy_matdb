#!/bin/bash

#SBATCH --time=4:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=4G   # memory per CPU core
#SBATCH -J "Pd DB"   # job name
#SBATCH --partition=physics
#SBATCH --array=1-1

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

cd /Users/trunks/codes/matdb-dist/tests/Pd/Pd.phonon-16/dynmatrix/W.$SLURM_ARRAY_TASK_ID
# Get the path to the executable; should be on user's path after the modules have been loaded.
vasp53s
