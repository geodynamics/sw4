#!/bin/sh
#Submit using flux batch <filename>

#flux: --job-name=SW4
#flux: --output='SPX.{{id}}.out'
#flux: --error='SPX.{{id}}.err'
#flux: -N 1
#flux: -l # Label I/O
#flux: --setattr=thp=always # Enable Transparent Huge Pages (THP)
#flux: -t 10
#flux: -q pdebug # Other options are plarge and pdebug
#flux: --setattr=gpumode=SPX

module load PrgEnv-amd
module load rocm/6.4.2beta1
module list

export MPICH_GPU_SUPPORT_ENABLED=1
export HSA_XNACK=1
export LIBSCI_ARCH_OVERRIDE=genoa
export OMP_NUM_THREADS=1


# Check if THP are enabled
cat /sys/kernel/mm/transparent_hugepage/enabled
rocm-smi

for i in $(seq 1 1);
do
	flux run -N1 -n4 -x -l -g1 ./sw4 hmr3.in.lassen
done

