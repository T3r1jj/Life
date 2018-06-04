#! /bin/bash

#SBATCH -n 2
#SBATCH -c 8

mpiexec -np 2 ./life 50 1000000 8 input50.txt
