#!/usr/bin/env python

# Arguments:
# check_results.py <err_tol|compare> <base_filename> <filename> <errInf tolerance> <errL2 tolerance> <solInf tolerance>
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

def main():
    if len(sys.argv) != 7: exit(1)

    if sys.argv[1] == "err_tol":
        # base_filename (argument 2) is ignored
        fp = open(sys.argv[3])
        fail = False
        tols = [float(sys.argv[i+4]) for i in range(3)]
        contents = fp.readlines()
        for i in range(len(contents)):
            line = contents[i]
            if line.find("Atten") == 0 or line.find("Disp") == 0:
                num_line = contents[i+1]
                vals = [float(x) for x in num_line.split()]
                checks = {"errInf": [tols[0], vals[0]], "errL2": [tols[1], vals[1]], "solInf": [tols[2], 1.0-vals[2]]}
                fail = fail or run_checks(checks)
        fp.close()
    elif sys.argv[1] == "compare":
        fail = False
        errInfTol = float(sys.argv[4])
        errL2Tol = float(sys.argv[5])
        base_file = open(sys.argv[2])
        test_file = open(sys.argv[3])
        base_data = base_file.readlines()
        test_data = test_file.readlines()
        base_file.close()
        test_file.close()
        if len(base_data) != len(test_data):
            print("ERROR: "+sys.argv[2]+" and "+sys.argv[3]+" are different lengths.")
            fail = True
        for i in range(min(len(test_data),len(base_data))):
            base_line = [float(x) for x in base_data[i].split()[1:3]]
            test_line = [float(x) for x in test_data[i].split()[1:3]]
            try:
                t0 = test_line[0]
                t1 = test_line[1]
                b0 = base_line[0]
                b1 = base_line[1]
                re0 = re1 = 0
                if t0 != 0 or b0 != 0: re0 = abs(b0-t0)/max(t0,b0)
                if t1 != 0 or b1 != 0: re1 = abs(b1-t1)/max(t1,b1)
                if re0 > errInfTol:
                    print("ERROR: Line %(line)d Tolerance: %(tol)f Actual: %(actual)f"%{"line":i, "tol": errInfTol, "actual": re0})
                    fail = True
                if re1 > errL2Tol:
                    print("ERROR: Line %(line)d Tolerance: %(tol)f Actual: %(actual)f"%{"line":i, "tol": errL2Tol, "actual": re1})
                    fail = True
            except:
                fail = True
    else:
        print("Unknown test type: "+sys.argv[1])
        fail = True

    if fail: exit(1)
    exit(0)

if __name__ == "__main__":
    main()

