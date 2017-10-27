#!/bin/bash
#SBATCH -N 12
#SBATCH -p debug
#SBATCH -A m2545
#SBATCH -S 4     # Special cores per nodes (4 are idle)
#SBATCH -t 00:30:00 
#SBATCH -C knl,quad,cache

numnodes=12

CASE=LOH.3-h50
RUN_COMMAND="./sw4 ${CASE}.in"
#    RUN_COMMAND="numactl -m 1 check-hybrid.intel.cori"

numprocspernode=64
numprocs=$(( ${numnodes}*${numprocspernode} ))
echo "numprocspernode=" $numprocspernode
# for numprocs in 1 2 4; do

# Number of ranks x threads is constant
# hc = hardware cores
numhc=$(echo ${numprocspernode} | awk '{print 64/$1}')
echo "numhc=" $numhc
# hyper-threads
numht=1
# logical cores = OMP_NUM_TASK
numcores=$(( ${numhc}*${numht} ))
# srun reservation logical cores per mpi-task  (-c arg)
numlc=$(( ${numhc}*4 ))

echo "Running on ${numnodes} nodes with ${numprocs} MPI ranks and OMP_NUM_THREADS=${numcores}."

# Intel OpenMP runtime parameters
export OMP_NUM_THREADS=${numcores}
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
# export KMP_PLACE_THREADS=1s${numhc}c${numht}t

# Run the job with this MPI + OpenMP configuration
MPI_COMMAND="srun -N ${numnodes} -n ${numprocs} -c ${numlc} --cpu_bind=cores" 

COMMAND="${MPI_COMMAND} ${RUN_COMMAND}" 
#echo "OMP_NUM_THREADS="${OMP_NUM_THREADS}
#echo ${COMMAND}

# actual sw4 command
${COMMAND} >& ${CASE}.out

