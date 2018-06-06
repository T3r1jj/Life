#! /bin/bash

#SBATCH --ntasks 2
#SBATCH -c 8

srun --mpi=pmi2 ./life 50 1000000 8 input50.txt
