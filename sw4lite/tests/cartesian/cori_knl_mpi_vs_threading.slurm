#!/bin/bash
#SBATCH -N 1
#SBATCH -p debug
#SBATCH -A m2545
# #SBATCH -S
#SBATCH -t 00:10:00 
#SBATCH -C knl,quad,flat

# export OMP_NUM_THREADS=32
# export OMP_PLACES=threads
# export OMP_PROC_BIND=spread
# export KMP_AFFINITY=verbose

# ./sw4lite_c uni.in
# srun -n 4 -c 64 --cpu_bind=cores numactl -m 1 ./sw4lite_c uni.in
# srun -n 4 -c 64 --cpu_bind=cores numactl -m 1 check-hybrid.intel.cori


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
    export OMP_PLACES=threads
    export OMP_PROC_BIND=spread
    # export KMP_PLACE_THREADS=1s${numomp}c${numht}t

    # Run the job with this MPI + OpenMP configuration
    MPI_COMMAND="srun -n ${numprocs} -c ${numcores} --cpu_bind=cores" 
#    RUN_COMMAND="numactl -m 1 check-hybrid.intel.cori"
    RUN_COMMAND="numactl -m 1 ./sw4lite_c uni.in"
    COMMAND="${MPI_COMMAND} ${RUN_COMMAND}" 
    echo ${COMMAND}
#    echo ${COMMAND} > timing_mpi${numprocs}_omp${numomp}_ht${numht}.txt
    ${COMMAND} > mpivsomp_mpi${numprocs}_omp${numomp}_ht${numht}.txt

done
