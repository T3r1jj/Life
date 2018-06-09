#!/bin/bash

#SBATCH -n 16
#SBATCH -c 8

mpiexec -np 16 ./life 4000 100 8 input4000.txt
