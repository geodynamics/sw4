#!/bin/bash
#SBATCH -A GEO130_crusher
#SBATCH -J SW4PYTEST
#SBATCH -o %x-%j.out
#SBATCH -t 00:35:00
#SBATCH -p batch
#SBATCH -N 1

export FI_MR_CACHE_MAX_COUNT=0
export MPICH_GPU_SUPPORT_ENABLED=1
# Option below for multinode runs only
export MPICH_OFI_NIC_POLICY=NUMA
date
./test_sw4.py -d optimize_mp -l 2
date
