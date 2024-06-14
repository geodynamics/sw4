#!/usr/bin/python3
import sys
import os
import flux
from flux.job import JobspecV1

def main():
    if len(sys.argv)!=3:
        print("Usage sw4.py <nodes> <input_file>")
        return 1
    handle = flux.Flux()
    nodes=int(sys.argv[1])
    print("Running sw4 on ",nodes, " using ",nodes*4," ranks with input file ",sys.argv[2])
    jobspec = JobspecV1.from_command(
        command=["./sw4",sys.argv[2]], num_tasks=nodes*4, num_nodes=nodes,
    )
    jobspec.cwd = os.getcwd()
    jobspec.exclusive=1
    jobspec.t=480
    jobspec.c=24
    #jobspec.setattr(thp,always)
    
    os.environ["MPICH_GPU_SUPPORT_ENABLED"]="1"
    os.environ["HSA_XNACK"]="1"
    jobspec.environment = dict(os.environ)
    
    jobspec.stdout="SW4.{{id}}.out"
    jobspec.stderr="SW4.{{id}}.err"
    print("Submitted ",flux.job.submit(handle, jobspec))


if __name__=="__main__":
    main()
