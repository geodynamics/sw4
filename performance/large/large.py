#!/usr/bin/python3
import sys
import os
import flux
from flux.job import JobspecV1

def main():

    handle = flux.Flux()
    jobspec = JobspecV1.from_command(
        command=["./sw4","hmr3XL.in"], num_tasks=496, num_nodes=124,
    )
    jobspec.cwd = os.getcwd()
    jobspec.exclusive=1
    jobspec.t=480
    jobspec.c=24
    jobspec.environment = dict(os.environ)
    jobspec.stdout="LargeTest.out"
    jobspec.stderr="LargeTest.err"
    print(flux.job.submit(handle, jobspec))


if __name__=="__main__":
    main()
