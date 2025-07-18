#!/bin/sh
#Submit using flux batch <filename>

#flux: --job-name=SW4
#flux: --output='CPX.{{id}}.out'
#flux: --error='CPX.{{id}}.err'
#flux: -N 1
#flux: -l # Label I/O
#flux: --setattr=thp=always # Enable Transparent Huge Pages (THP)
#flux: -t 20
#flux: -q pdebug # Other options are plarge and pdebug
#flux: --setattr=gpumode=CPX

export MPICH_GPU_SUPPORT_ENABLED=1
export HSA_XNACK=1

# Check if THP are enabled
cat /sys/kernel/mm/transparent_hugepage/enabled
rocm-smi

flux run -N1 -n24 -x -l ./cpx_bind ./sw4 hmr3.in.lassen

