#!/bin/bash
#BSUB -P lc
#BSUB -W 05:30
#BSUB -nnodes 1
#BSUB -q pdebug
#BSUB -J Test
#BSUB -o RunRajaMagma_%J.out
#BSUB -e RunRajaMagma_%J.err
date
lrun -T4 -M -gpu  ./sw4 hmr.in 
