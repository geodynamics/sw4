#!/usr/bin/env python3

import os, sys, argparse, subprocess, verify_hdf5

sw4_exe = '../optimize_mp/sw4'
sw4_input_file = './reference/hdf5/loh1-h100-mr-hdf5-1.in'
node_name = os.uname()[1]


if 'quartz' in node_name:
    omp_threads=2
    mpi_tasks = int(36/omp_threads)
    mpirun_cmd="srun -ppdebug -n " + str(mpi_tasks) + " -c " + str(omp_threads)
elif 'nid' in node_name: # the cori knl nodes are called nid
    omp_threads=4;
    mpi_tasks = int(32/omp_threads)# use 64 hardware cores per node
    sw_threads = 1*omp_threads # Cori uses hyperthreading by default
    mpirun_cmd="srun --cpu_bind=cores -n " + str(mpi_tasks) + " -c " + str(sw_threads)
else:
    #default mpi command
    omp_threads=4
    mpi_tasks = 4
    mpirun_cmd="mpirun -np " + str(mpi_tasks)

sw4_mpi_run = mpirun_cmd + ' ' + sw4_exe
if (omp_threads>0):
    os.putenv("OMP_NUM_THREADS", str(omp_threads))

run_cmd = mpirun_cmd.split() + [
    sw4_exe,
    sw4_input_file
]

print('MPI run command is: ', run_cmd)

status = subprocess.run(run_cmd)

success = verify_hdf5.verify('./', 1e-5)

if success:        
    print('HDF5 Test PASSED')
else:
    print('HDF5 Test FAILED')

