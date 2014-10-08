#!/usr/bin/env python

# Arguments:
# check_results.py <filename> <errInf tolerance> <errL2 tolerance> <solInf tolerance>
# This assumes displacement variables and attenuation variables have the same tolerance

import sys

def run_checks(checks):
    fail = False
    for elem in checks:
        pass_check = (checks[elem][0] > checks[elem][1])
        print("%(check)s    tolerance: %(tol)f   calculated: %(val)f    pass: %(pass)s" %
                {"check": elem, "tol": checks[elem][0], "val": checks[elem][1], "pass": pass_check})
        if not pass_check: fail = True
    return fail

if len(sys.argv) != 5: exit(1)

fp = open(sys.argv[1])
fail = False
tols = [float(sys.argv[i+2]) for i in range(3)]
contents = fp.readlines()
for i in range(len(contents)):
    line = contents[i]
    if line.find("Atten") == 0 or line.find("Disp") == 0:
        num_line = contents[i+1]
        vals = [float(x) for x in num_line.split()]
        checks = {"errInf": [tols[0], vals[0]], "errL2": [tols[1], vals[1]], "solInf": [tols[2], 1.0-vals[2]]}
        fail = fail or run_checks(checks)

fp.close()


if fail: exit(1)
else: exit(0)

