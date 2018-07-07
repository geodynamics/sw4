#!/bin/bash
#SBATCH --account m2868
#SBATCH -N 1
#SBATCH -S 4
#SBATCH -t 0:30:00
#SBATCH -p debug 
#SBATCH -C knl,quad,cache
#SBATCH -L SCRATCH 

# Script for correctness testing of checkpoint restart

# Set number of threads per node
# Set number of OpenMP threads per node
export OMP_NUM_THREADS=2

# Always use these values:
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

# soft link to the local executable
SW4EXE=./sw4

all_cases='straight checkp straight-att checkp-att straight-mr checkp-mr straight-att-mr checkp-att-mr'
for CASE in $all_cases
do
    echo "------------------------- NEW RUN ------------------------------"
    echo CASE: $CASE
    srun -n 32 -c 8 --cpu_bind=cores ${SW4EXE} ${CASE}.in 
#    echo srun -N ${NODES} -n ${TASKS} -c ${NUMLC} --cpu_bind=cores ${SW4EXE} ${CASE}.in 
    echo "------------------------- END OF RUN ------------------------------"
done

    
