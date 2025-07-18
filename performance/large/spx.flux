#!/bin/sh
#Submit using flux batch <filename>

#flux: --job-name=SW4
#flux: --output='SPX.{{id}}.out'
#flux: --error='SPX.{{id}}.err'
#flux: -N 1
#flux: -l # Label I/O
#flux: --setattr=thp=always # Enable Transparent Huge Pages (THP)
#flux: -t 600
#flux: -q pbatch # Other options are plarge and pdebug
#flux: --requires=host:elcap3461
#flux: --setattr=gpumode=SPX

export MPICH_GPU_SUPPORT_ENABLED=1
export HSA_XNACK=1
export LIBSCI_ARCH_OVERRIDE=genoa
export OMP_NUM_THREADS=1


# Check if THP are enabled
cat /sys/kernel/mm/transparent_hugepage/enabled
rocm-smi

for i in $(seq 1 50);
do
	flux run -N1 -n4 -x -l -g1 ./sw4 hmr3.in.lassen
done

