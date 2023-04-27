#!/bin/sh
#Submit using flux batch <filename> or sbatch <filename> with flux_wrappers
#flux: -N1 -n64
#flux: --job-name=SW4
#flux: --output='SW4.{{id}}.out'
#flux: --error='SW4.{{id}}.err'
#flux: -t 20 
export MPICH_GPU_SUPPORT_ENABLED=1
flux run -n 8 -c 8 -g 1 -o gpu-affinity=per-task -o cpu-affinity=per-task  ./sw4 hmr3.in.lassen
flux run -n 8 -c 8 -g 1 -o gpu-affinity=per-task -o cpu-affinity=per-task  ./sw4 h.in.lassen
diff hayward-att-h270-ref-result/sta1.txt hayward-att-h270-ref-result/sta1.txt.last
diff hayward-att-h200-ref-result/sta1.txt hayward-att-h200-ref-result/sta1.txt.last

