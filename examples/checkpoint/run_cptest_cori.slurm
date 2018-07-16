#!/bin/csh
#SBATCH --account m2868
#SBATCH -N 4
#SBATCH -S 4
#SBATCH -t 0:30:00
#SBATCH -p debug 
#SBATCH -C knl,quad,cache
#SBATCH -L SCRATCH 
#DW persistentdw name=EQSIM

#NOTE: this script needs to be placed on the parallel file system before executing
# The files ${CASE}.in and ${RUPTURE} must be in the same directory on the parallel file system

# Run H=800m on 4 nodes
set RUPTURE = m7.0-50.0x20.0.s600.0.v5.1.curved.srf
set CASE = cptest4

# number of parallel writers = striping on the parallel file system
# this number should be the same as the fileio writers in the input file
set IOWRITERS=4

# Set total number of nodes request (must match -N above)
set NODES = 4 

# Set number of threads per node
# Set number of OpenMP threads per node
setenv OMP_NUM_THREADS 2
set PROCPERNODE = 32
#setenv OMP_NUM_THREADS 8
#set PROCPERNODE = 8 

# Always use these values:
setenv OMP_PLACES threads
setenv OMP_PROC_BIND spread

echo
echo NODES: $NODES
echo PROCPERNODE: $PROCPERNODE
echo OMP_NUM_THREADS: $OMP_NUM_THREADS

 set SW4BIN = /global/project/projectdirs/m2545/sw4/cori-knl/optimize
# This always links to the latest executable
set SW4FILE = sw4

set TASKS = ` echo $NODES $PROCPERNODE | awk '{ print $1 * $2 }' `
echo TASKS: $TASKS

# Number of Logical cores, for argument of -c in srun (below)
set NUMLC = ` echo $OMP_NUM_THREADS | awk '{ print $1 * 4 }' ` 
echo NUMLC: $NUMLC
echo
echo "Running on ${NODES} nodes with ${TASKS} MPI ranks and OMP_NUM_THREADS=${OMP_NUM_THREADS}"

echo CASE: $CASE

# insert path to burst buffer
echo DW_PERSISTENT_STRIPED_EQSIM: $DW_PERSISTENT_STRIPED_EQSIM
envsubst < $CASE.ENV >& $CASE.in

# Make and stripe output directories
if (! -d $CASE.sw4output ) then
  mkdir -p $CASE.sw4output
  lfs setstripe -c $IOWRITERS $CASE.sw4output
endif
#if (! -d $CASE.sw4restart ) then
#  mkdir -p $CASE.sw4restart
#  lfs setstripe -c $IOWRITERS  $CASE.sw4restart
#endif

# SBCAST files to assigned nodes
sbcast -f -F2 -t 300 --compress=lz4 $SW4BIN/$SW4FILE /tmp/$SW4FILE
sbcast -f -F2 -t 300 --compress=lz4 ./$CASE.in /tmp/$CASE.in
sbcast -f -F2 -t 300 --compress=lz4 ./$RUPTURE /tmp/$RUPTURE

date
srun -N $NODES -n $TASKS -c $NUMLC --cpu_bind=cores /tmp/$SW4FILE /tmp/$CASE.in >! $CASE.output
date
