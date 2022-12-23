#!/bin/bash
#SBATCH -N 1
#SBATCH -p pbatch
#SBATCH -t 120

date
stdbuf -e0 -o0 ./test_sw4.py -d optimize_mp -l 2

