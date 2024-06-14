#!/bin/sh
#Submit using flux batch <filename>

#flux: --job-name=SW4
#flux: --output='SW4.{{id}}.out'
#flux: --error='SW4.{{id}}.err'
#flux: -N 1
#flux: -l # Label I/O
#flux: --setattr=thp=always # Enable Transparent Huge Pages (THP)
#flux: -t 20
#flux: -q pdev # Other options are plarge and pdebug

export MPICH_GPU_SUPPORT_ENABLED=1
export HSA_XNACK=1

# Check if THP are enabled
cat /sys/kernel/mm/transparent_hugepage/enabled

flux run -N1 -n4 -x -l -g1 -c24 ./sw4 hmr3.in
flux run -N1 -n4 -x -l -g1 -c24 ./sw4 h.in.rzansel


