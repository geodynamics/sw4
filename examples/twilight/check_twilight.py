#!/usr/bin/env python

# Arguments:
# check_twilight.py <filename> <errInf tolerance> <errL2 tolerance> <solInf tolerance>

import sys

if len(sys.argv) != 5: exit(1)

fp = open(sys.argv[1])
vals = [float(x) for x in fp.readlines()[-1].split()]
tols = [float(sys.argv[i+2]) for i in range(3)]
fp.close()

fail = False
checks = {"errInf": [tols[0], vals[0]], "errL2": [tols[1], vals[1]], "solInf": [tols[2], vals[2]]}
for elem in checks:
    pass_check = (checks[elem][0] > checks[elem][1])
    print("%(check)s    tolerance: %(tol)f   calculated: %(val)f    pass: %(pass)s" %
            {"check": elem, "tol": checks[elem][0], "val": checks[elem][1], "pass": pass_check})
    if not pass_check: fail = True

if fail: exit(1)
else: exit(0)

