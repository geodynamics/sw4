#!/bin/bash
#SBATCH -N 1
#SBATCH -p debug
#SBATCH -A m2868
#SBATCH -S 4     # Special cores per nodes (4 are idle)
#SBATCH -t 00:30:00 
#SBATCH -C knl,quad,cache

# the srun commands are inside the python script
# 
#./test_sw4.py --verbose
# default level is 0; level 2 tests all cases
./test_sw4.py --level 0 -d optimize_mp

