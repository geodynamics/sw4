#!/bin/sh
#Submit using flux batch <filename>

#flux: --job-name=SW4
#flux: --output='TPXNEW_4NODE.{{id}}.out'
#flux: --error='TPXNEW_4NODE.{{id}}.err'
#flux: -N 4
#flux: -l # Label I/O
#flux: --setattr=thp=always # Enable Transparent Huge Pages (THP)
#flux: -t 30
#flux: -q pdebug # Other options are plarge and pdebug
#flux: --setattr=gpumode=TPX

export MPICH_GPU_SUPPORT_ENABLED=1
export HSA_XNACK=1
export FLUX_MPIBIND_USE_TOPOFILE=1
export MPIBIND_TOPOFILE=/collab/usr/global/tools/mpi/mpibind/utils/tuo-fix-tpx.xml
export LIBSCI_ARCH_OVERRIDE=genoa
export OMP_NUM_THREADS=1


# Check if THP are enabled
cat /sys/kernel/mm/transparent_hugepage/enabled
flux run -N4 -n4 -x -l -c 48 rocm-smi
flux run -N 4 -n 48 -x -l -o initrc=/collab/usr/global/tools/mpi/mpibind/utils/mpibind-flux-250109.lua env | grep ROCR | sort -n | uniq
flux run -N4 -n48 -x -l -o initrc=/collab/usr/global/tools/mpi/mpibind/utils/mpibind-flux-250109.lua ./sw4 4node.in

