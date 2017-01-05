#!/bin/bash
#SBATCH -N 1
#SBATCH -p debug
#SBATCH -A m2545
# #SBATCH -S
#SBATCH -t 00:10:00 
#SBATCH -C knl,quad,flat

export OMP_NUM_THREADS=32
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
# export KMP_AFFINITY=verbose

numactl -H

srun -n 4 -c 64 --cpu_bind=cores numactl -m 1 ./sw4lite_c uni.in
# srun -n 4 -c 64 --cpu_bind=cores numactl -m 1 check-hybrid.intel.cori

