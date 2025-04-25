#!/bin/sh
#Submit using flux batch <filename>

#flux: --job-name=SW4
#flux: --output='CPX_2NODE.{{id}}.out'
#flux: --error='CPX_2NODE.{{id}}.err'
#flux: -N 2
#flux: -l # Label I/O
#flux: --setattr=thp=always # Enable Transparent Huge Pages (THP)
#flux: -t 300
#flux: -q pbatch # Other options are plarge and pdebug
#flux: --requires=hostlist:elcap[11359-11361,11759]
#flux: --setattr=gpumode=CPX

export MPICH_GPU_SUPPORT_ENABLED=1
export HSA_XNACK=1
export FLUX_MPIBIND_USE_TOPOFILE=1
export MPIBIND_TOPOFILE=/collab/usr/global/tools/mpi/mpibind/utils/tuo-fix-cpx.xml
export LIBSCI_ARCH_OVERRIDE=genoa
export OMP_NUM_THREADS=1


# Check if THP are enabled
cat /sys/kernel/mm/transparent_hugepage/enabled
rocm-smi
flux run -N 2 -n 48 -x -l -o initrc=/collab/usr/global/tools/mpi/mpibind/utils/mpibind-flux-250109.lua env | grep ROCR
flux run -N 2 -n 48 -x -l -o initrc=/collab/usr/global/tools/mpi/mpibind/utils/mpibind-flux-250109.lua rocm-smi
for i in $(seq 1 50);
do
	flux run -N2 -n48 -x -l -o initrc=/collab/usr/global/tools/mpi/mpibind/utils/mpibind-flux-250109.lua ./sw4 hmr3.in.lassen
done
