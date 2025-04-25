#!/bin/sh
#Submit using flux batch <filename>

#flux: --job-name=CORESPEC
#flux: --output='CORE.{{id}}.out'
#flux: --error='CORE.{{id}}.err'
#flux: -N 1
#flux: -l # Label I/O
#flux: --setattr=thp=always # Enable Transparent Huge Pages (THP)
#flux: -t 360
#flux: -q pbatch # Other options are plarge and pdebug

export MPICH_GPU_SUPPORT_ENABLED=1
export HSA_XNACK=1
export LIBSCI_ARCH_OVERRIDE=genoa

# Check if THP are enabled
cat /sys/kernel/mm/transparent_hugepage/enabled
# Warmup
flux run -N1 -n4 -x -g1 -c24 ./sw4.org hmr3.in.lassen

flux run -N1 -n4 -x -g1 -c24 ./sw4.org hmr3.in.lassen
flux run -N1 -n4 -x -g1 -c24 ./sw4.org hmr3.in.lassen
flux run -N1 -n4 -x -g1 -c24 ./sw4.org hmr3.in.lassen
flux run -N1 -n4 -x -g1 -c24 ./sw4.org hmr3.in.lassen

flux run -N1 -n4 -x -g1 -c21 ./sw4.org hmr3.in.lassen
flux run -N1 -n4 -x -g1 -c21 ./sw4.org hmr3.in.lassen
flux run -N1 -n4 -x -g1 -c21 ./sw4.org hmr3.in.lassen
flux run -N1 -n4 -x -g1 -c21 ./sw4.org hmr3.in.lassen

unset MPIBIND_RESTRICT
flux run -N1 -n4 -x -g1 -c24 ./sw4.org hmr3.in.lassen
flux run -N1 -n4 -x -g1 -c24 ./sw4.org hmr3.in.lassen
flux run -N1 -n4 -x -g1 -c24 ./sw4.org hmr3.in.lassen
flux run -N1 -n4 -x -g1 -c24 ./sw4.org hmr3.in.lassen

flux run -N1 -n4 -x -g1 -c21 ./sw4.org hmr3.in.lassen
flux run -N1 -n4 -x -g1 -c21 ./sw4.org hmr3.in.lassen
flux run -N1 -n4 -x -g1 -c21 ./sw4.org hmr3.in.lassen
flux run -N1 -n4 -x -g1 -c21 ./sw4.org hmr3.in.lassen


