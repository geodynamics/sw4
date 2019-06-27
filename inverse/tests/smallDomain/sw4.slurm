#!/bin/bash -l
#SBATCH --mail-user=kechen@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH -N 16                           # number of nodes
#SBATCH -q premium                        # submit to the regular QOS
#SBATCH -C knl,quad,cache               # note: non-quad,cache are via reservation only
#SBATCH -S 4                            # use core specialization
#SBATCH -t 00:30:00                     # time limit
#SBATCH -L SCRATCH,project              # specify license need for the file systems the job needs

export OMP_NUM_THREADS=2
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

srun -n 512 -c 8 --cpu_bind=cores sw4 alamo-sep-5-2008-forward.in
