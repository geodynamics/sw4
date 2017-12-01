#!/bin/bash
#SBATCH -N 4
#SBATCH -p pdebug
#SBATCH -t 00:30:00 

numnodes=4

#numprocspernode = Number of MPI-tasks per node
for numprocspernode in 16 8 4; do

numprocs=$(( ${numnodes}*${numprocspernode} ))
echo "numprocspernode=" $numprocspernode

# Number of ranks x threads is constant
# hc = hardware cores
numhc=$(echo ${numprocspernode} | awk '{print 16/$1}')
echo "numhc=" $numhc
# hyper-threads
numht=1
# logical cores = OMP_NUM_TASK
numcores=$(( ${numhc}*${numht} ))
# srun reservation logical cores per mpi-task  (-c arg) (*1 for cab)
numlc=$(( ${numhc}*1 ))

echo "Running on ${numnodes} nodes with ${numprocs} MPI ranks and OMP_NUM_THREADS=${numcores}."

# Intel OpenMP runtime parameters
export OMP_NUM_THREADS=${numcores}
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

# Run the job with this MPI + OpenMP configuration
MPI_COMMAND="srun -N ${numnodes} -n ${numprocs} -c ${numlc} --cpu_bind=cores" 

# 1st case
CASE=berkeley-nomr-cart
RUN_COMMAND="./sw4_mp ${CASE}.in"
${MPI_COMMAND} ${RUN_COMMAND}  >& ${CASE}-sw4_mp-${numprocs}mpi-${numcores}threads.out

done

