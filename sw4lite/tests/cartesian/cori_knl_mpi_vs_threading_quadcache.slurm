#!/bin/bash
#SBATCH -N 1
#SBATCH -p debug
#SBATCH -A m2545
# #SBATCH -S
#SBATCH -t 00:10:00 
#SBATCH -C knl,quad,cache

module load craype-hugepages2M

export KMP_AFFINITY=compact,granularity=thread,1
export KMP_BLOCKTIME=infinite

# do the loop
for numprocs in 1 2 4 8 16 32 64; do
# for numprocs in 1 2 4; do

# Number of ranks x threads is constant
numomp=$(echo ${numprocs} | awk '{print 64/$1}')
numht=2
numcores=$(( ${numomp}*${numht} ))

    echo "Running with ${numprocs} MPI ranks and ${numomp} cores, ${numht} threads per core."

    # Intel OpenMP runtime parameters
    export OMP_NUM_THREADS=${numcores}
    # export OMP_PLACES=threads
    # export OMP_PROC_BIND=spread
    # export KMP_PLACE_THREADS=1s${numomp}c${numht}t

# export OMP_NUM_THREADS=8
# srun --ntasks=64 --ntasks-per-node=8 --cpus-per-task=32  --cpu_bind=threads numactl -m 0 ./hpgmg-fv 8 1
# export OMP_NUM_THREADS=16
# srun --ntasks=32 --ntasks-per-node=4 --cpus-per-task=64  --cpu_bind=threads numactl -m 0 ./hpgmg-fv 8 2
# export OMP_NUM_THREADS=32
# srun --ntasks=16 --ntasks-per-node=2 --cpus-per-task=128 --cpu_bind=threads numactl -m 0 ./hpgmg-fv 8 4
# export OMP_NUM_THREADS=64
# srun --ntasks=8  --ntasks-per-node=1 --cpus-per-task=256 --cpu_bind=threads numactl -m 0 ./hpgmg-fv 8 8

    # Run the job with this MPI + OpenMP configuration
    MPI_COMMAND="srun --ntasks=${numprocs} --ntasks-per-node=${numprocs} --cpus-per-task=${numcores} --cpu_bind=cores" 
#    RUN_COMMAND="numactl -m 1 check-hybrid.intel.cori"
    RUN_COMMAND="numactl -m 0 ./sw4lite_c uni.in"
    COMMAND="${MPI_COMMAND} ${RUN_COMMAND}" 
    echo ${COMMAND}
#    echo ${COMMAND} > timing_mpi${numprocs}_omp${numomp}_ht${numht}.txt
    ${COMMAND} > mpivsomp_mpi${numprocs}_omp${numomp}_ht${numht}_qc.txt

done
