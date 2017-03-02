#!/bin/bash
#SBATCH -N 1
#SBATCH -p debug
#SBATCH -A m2545
#SBATCH -S 4
#SBATCH -t 00:10:00 
##SBATCH -C knl,quad,flat
#SBATCH -C knl,quad,cache

export OMP_PLACES=threads
export OMP_PROC_BIND=spread

# quad,flat
# C-executable
#srun -n 4 -c 64 --cpu_bind=cores numactl -m 1 ./sw4lite_c gaussianHill.in
# fortran executable
#srun -n 4 -c 64 --cpu_bind=cores numactl -m 1 ./sw4lite gaussianHill.in

#quad, cache, fortran exe

echo "64 MPI, 4 cores, 1 HT"
export OMP_NUM_THREADS=1
srun -n 64 -c 4 --cpu_bind=cores ./sw4 gauss-twi-2.in



