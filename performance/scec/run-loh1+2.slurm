#!/bin/bash
#SBATCH -N 4
#SBATCH -p debug
#SBATCH -A m2545
#SBATCH -S 4     # Special cores per nodes (4 are idle)
#SBATCH -t 00:30:00 
#SBATCH -C knl,quad,cache
##SBATCH -C knl,quad,flat

numnodes=4


#  icoe branch
#    RUN_COMMAND="numactl -m 1 check-hybrid.intel.cori"

#numprocspernode=16
for numprocspernode in 1 2 4 8; do

# developer branch
#  numprocspernode=64

numprocs=$(( ${numnodes}*${numprocspernode} ))
echo "numprocspernode=" $numprocspernode

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
# quad modes
MPI_COMMAND="srun -N ${numnodes} -n ${numprocs} -c ${numlc} --cpu_bind=cores" 
# flat mode needs the numactl command
#MPI_COMMAND="srun -N ${numnodes} -n ${numprocs} -c ${numlc} --cpu_bind=cores numactl -m 1" 

# 1st case
CASE=LOH.1-h50
RUN_COMMAND="./sw4 ${CASE}.in"
${MPI_COMMAND} ${RUN_COMMAND}  >& ${CASE}-sw4-${numprocs}mpi-${numcores}threads.out

# 2nd case
CASE=LOH.2-h50
RUN_COMMAND="./sw4 ${CASE}.in"
${MPI_COMMAND} ${RUN_COMMAND}  >& ${CASE}-sw4-${numprocs}mpi-${numcores}threads.out

done

