#!/bin/sh
#Submit using flux batch <filename>

#flux: --job-name=SW4
#flux: --output='TPXNEW.{{id}}.out'
#flux: --error='TPXNEW.{{id}}.err'
#flux: -N 1
#flux: -l # Label I/O
#flux: --setattr=thp=always # Enable Transparent Huge Pages (THP)
#flux: -t 240
#flux: -q pbatch  # Other options are plarge and pdebug
#flux: --setattr=gpumode=TPX

export MPICH_GPU_SUPPORT_ENABLED=1
export HSA_XNACK=1
export FLUX_MPIBIND_USE_TOPOFILE=1
export MPIBIND_TOPOFILE=/collab/usr/global/tools/mpi/mpibind/utils/tuo-fix-tpx.xml
export LIBSCI_ARCH_OVERRIDE=genoa
export OMP_NUM_THREADS=1


# Check if THP are enabled
cat /sys/kernel/mm/transparent_hugepage/enabled
rocm-smi
flux run -N 1 -n 12 -x -l -o initrc=/collab/usr/global/tools/mpi/mpibind/utils/mpibind-flux-250109.lua env | grep ROCR
for i in $(seq 1 50);
do
	flux run -N1 -n12 -x -l -o initrc=/collab/usr/global/tools/mpi/mpibind/utils/mpibind-flux-250109.lua ./sw4 hmr3.in.lassen
done

