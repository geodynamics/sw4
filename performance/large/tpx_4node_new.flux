#!/bin/sh
#Submit using flux batch <filename>

#flux: --job-name=SW4
#flux: --output='TPX_4NODE_NEW.{{id}}.out'
#flux: --error='TPX_4NODE_NEW.{{id}}.err'
#flux: -N 4
#flux: -l # Label I/O
#flux: --setattr=thp=always # Enable Transparent Huge Pages (THP)
#flux: -t 300
#flux: -q pbatch  # Other options are plarge and pdebug
#flux: --requires=hostlist:tuolumne[1848,1857-1858,1861]
#flux: --setattr=gpumode=TPX
#flux: --conf=resource.rediscover=true

export MPICH_GPU_SUPPORT_ENABLED=1
export HSA_XNACK=1
export LIBSCI_ARCH_OVERRIDE=genoa
export OMP_NUM_THREADS=1


# Check if THP are enabled
cat /sys/kernel/mm/transparent_hugepage/enabled
rocm-smi
flux run -N 4-n 48 -c 84 -x -l rocm-smi
flux run -N 4-n 48 -c 84 -x -l env | grep ROCR
for i in $(seq 1 50);
do
	flux run -N4 -n48 -x -l ./sw4 4node.in
done

