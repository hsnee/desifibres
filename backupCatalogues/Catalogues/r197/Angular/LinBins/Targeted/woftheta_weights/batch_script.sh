#!/bin/bash -l
#
# Batch script for bash users
#
#BSUB -L /bin/bash
#BSUB -n 12
#BSUB -R "span[hosts=1]"
#BSUB -x
#BSUB -J CUTE
#BSUB -oo angA0r0noweight0.out
#BSUB -q cosma
#BSUB -P durham
#
# 1 hour wall clock time limit
#BSUB -W 10:00
#BSUB -M 60000
export OMP_NUM_THREADS=24

# Run the program
module load python/2.7.3
./CUTE param.txt
