#! /bin/bash

#SBATCH -n 16
#SBATCH -c 8

mpiexec -np 16 ./life 30000 10 8 input4000.txt
