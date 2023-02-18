#!/bin/bash -l
#SBATCH -A GEO130
#SBATCH -J SW41KCOMPARE
#SBATCH -o %x-%j.out
#SBATCH -t 06:00:00
#SBATCH -p batch
#SBATCH -N 1000
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=pankajakshan1@llnl.gov

module load PrgEnv-amd cray-hdf5-parallel cray-fftw cray-python sqlite libtiff
module load rocm
module -t list

export FI_MR_CACHE_MAX_COUNT=0
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_OFI_NIC_POLICY=NUMA

date
stdbuf -o0 -e0 srun -N 1000 -n 8000 -c7 --gpus-per-task=1 --gpu-bind=closest -p batch ./sw4 Frontier-test-10Hz.sw4input
date
