#!/bin/sh
#Submit using flux batch <filename>

#flux: --job-name=SW4
#flux: --output='CPXNEW.{{id}}.out'
#flux: --error='CPXNEW.{{id}}.err'
#flux: -N 1
#flux: -l # Label I/O
#flux: --setattr=thp=always # Enable Transparent Huge Pages (THP)
#flux: -t 600
#flux: -q pbatch # Other options are plarge and pdebug
#flux: --requires=host:elcap3461
#flux: --setattr=gpumode=CPX

export MPICH_GPU_SUPPORT_ENABLED=1
export HSA_XNACK=1
# for proper binding, add next 2 lines + the -o initrc=/collab/usr/global/tools/mpi/mpibind/utils/mpibind-flux-250109.lua env option in the flux run command
export FLUX_MPIBIND_USE_TOPOFILE=1
export MPIBIND_TOPOFILE=/collab/usr/global/tools/mpi/mpibind/utils/tuo-fix-cpx.xml

export LIBSCI_ARCH_OVERRIDE=genoa
export OMP_NUM_THREADS=1


# Check if THP are enabled
cat /sys/kernel/mm/transparent_hugepage/enabled
rocm-smi
flux run -N 1 -n 24 -x -l -o initrc=/collab/usr/global/tools/mpi/mpibind/utils/mpibind-flux-250109.lua env | grep ROCR
for i in $(seq 1 50);
do
	flux run -N1 -n24 -x -l -o initrc=/collab/usr/global/tools/mpi/mpibind/utils/mpibind-flux-250109.lua ./sw4 hmr3.in.lassen
done
