#!/bin/sh
#Submit using flux batch <filename>

#flux: --job-name=SW4
#flux: --output='SW4EL.{{id}}.out'
#flux: --error='SW4EL.{{id}}.err'
#flux: -N 750
#flux: -l # Label I/O
#flux: --setattr=thp=always # Enable Transparent Huge Pages (THP)
#flux: -t 720
#flux: -q pbatch # Other options are plarge and pdebug

export MPICH_GPU_SUPPORT_ENABLED=1
export HSA_XNACK=1
export LIBSCI_ARCH_OVERRIDE=genoa
export OMP_NUM_THREADS=1
export FASTLOAD_VERBOSE=1

# Check if THP are enabled
cat /sys/kernel/mm/transparent_hugepage/enabled

#flux run -N750 -n3000 -x -l -g1 -c24 ./sw4 hmr3.in.750.large
flux run -N750 -n750 -x -l -g1 -c24 cat /sys/kernel/mm/transparent_hugepage/enabled
flux run -N750 -n3000 -x -l -g1 -c24 fastload ./sw4.org hmr3.in.750.large


