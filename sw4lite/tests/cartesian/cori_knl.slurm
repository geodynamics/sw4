#!/bin/bash
#SBATCH -N 1
#SBATCH -p debug
#SBATCH -A m2545
#SBATCH -S 4
#SBATCH -t 00:10:00 
#SBATCH -C knl,quad,flat
##SBATCH -C knl,quad,cache

export OMP_NUM_THREADS=32
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

#numactl -H

# quad,flat
# C-executable
srun -n 4 -c 64 --cpu_bind=cores numactl -m 1 ./sw4lite_c uni.in
# fortran executable
#srun -n 4 -c 64 --cpu_bind=cores numactl -m 1 ./sw4lite uni.in

#quad, cache
#srun -n 4 -c 64 --cpu_bind=cores ./sw4lite_c uni.in


