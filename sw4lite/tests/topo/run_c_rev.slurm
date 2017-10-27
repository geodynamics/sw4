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

#quad, cache
# echo "1 MPI, 64 cores, 1 HT"
# export OMP_NUM_THREADS=64
# srun -n 1 -c 256 --cpu_bind=cores ./sw4lite_c gaussianHill-rev.in

# echo "2 MPI, 32 cores, 1 HT"
# export OMP_NUM_THREADS=32
# srun -n 2 -c 128 --cpu_bind=cores ./sw4lite_c gaussianHill-rev.in

echo "4 MPI, 16 cores, 1 HT"
export OMP_NUM_THREADS=16
srun -n 4 -c 64 --cpu_bind=cores ./sw4lite_c gaussianHill-rev.in

# echo "4 MPI, 16 cores, 2 HT"
# export OMP_NUM_THREADS=32
# srun -n 4 -c 64 --cpu_bind=cores ./sw4lite_c gaussianHill-rev.in



