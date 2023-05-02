#!/usr/bin/python3
import sys
import os
import flux
from flux.job import JobspecV1

def main():

    handle = flux.Flux()
    jobspec = JobspecV1.from_command(
        command=["hostname"], num_tasks=4, num_nodes=1,cores_per_task=2
    )
    jobspec.cwd = os.getcwd()
    jobspec.environment = dict(os.environ)
    jobspec.stdout="Test.{{id}}.out"
    jobspec.stderr="Test.{{id}}.err"
    print(flux.job.submit(handle, jobspec))


if __name__=="__main__":
    main()
