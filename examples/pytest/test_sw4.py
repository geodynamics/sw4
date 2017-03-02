# Arguments:
# This assumes displacement variables and attenuation variables have the same tolerance

import os, sys

#------------------------------------------------
def run_checks(checks):
    fail = False
    for elem in checks:
        pass_check = (checks[elem][0] > checks[elem][1])
        print("%(check)s    tolerance: %(tol)f   calculated: %(val)f    pass: %(pass)s" %
                {"check": elem, "tol": checks[elem][0], "val": checks[elem][1], "pass": pass_check})
        if not pass_check: fail = True
    return fail

#------------------------------------------------
def compare_to_ref(base_file_name, test_file_name, errInfTol, errL2Tol):

    success = True

    base_file = open(base_file_name)
    test_file = open(test_file_name)
    base_data = base_file.readlines()
    test_data = test_file.readlines()

    # tmp
    #print('base_data:', base_data)
    #print('test_data:', test_data)

    base_file.close()
    test_file.close()
    if len(base_data) != len(test_data):
        print("ERROR: " + base_file_name + " and " + test_file_name + " are different lengths.")
        return False

    # for twilight tests, compare the 3 numbers on the last line of each file
    base_numbers = base_data[-1:][0]
    test_numbers = test_data[-1:][0]
    base_data = base_numbers.split() #base_data is a list of strings
    test_data = test_numbers.split()
    #print('base_data=', base_data)
    #print('test_data=', test_data)

    base_num = [float(x) for x in base_data] #base_num is a list of floats
    test_num = [float(x) for x in test_data]
    #print('base_num=', base_num);
    #print('test_num=', test_num);

    try:
        t0 = test_num[0]
        t1 = test_num[1]
        b0 = base_num[0]
        b1 = base_num[1]
        re0 = re1 = 0
        if t0 != 0 or b0 != 0: re0 = abs(b0-t0)/max(t0,b0)
        if t1 != 0 or b1 != 0: re1 = abs(b1-t1)/max(t1,b1)
        # tmp
        #print('t0=', t0, 'b0=', b0, 'rel err=', re0);
        #print('t1=', t1, 'b1=', b1, 'rel err=', re1);
        # end tmp
        if re0 > errInfTol:
            print("ERROR: Linf tolerance: %(tol)f Actual: %(actual)f"%{"tol": errInfTol, "actual": re0})
            success = False
        if re1 > errL2Tol:
            print("ERROR: L2 tolerance: %(tol)f Actual: %(actual)f"%{"tol": errL2Tol, "actual": re1})
            success = False
    except:
        success = False

    return success

#------------------------------------------------
def main_test():
    sep = '/'
    pytest_dir = os.getcwd()
    print('pytest_dir =', pytest_dir)
    pytest_dir_list = pytest_dir.split(sep)
    sw4_base_list = pytest_dir_list[:-2] # discard the last 2 sub-directories (examples/pytest)

    sw4_base_dir = sep.join(sw4_base_list)
    examples_dir = sw4_base_dir + '/examples'
    optimize_dir =  sw4_base_dir + '/optimize'

    print('sw4_base_dir =', sw4_base_dir)
    print('examples_dir =', examples_dir)
    print('optimize_dir =', optimize_dir)          
    
    sw4_exe = optimize_dir + '/sw4'
    #print('sw4-exe = ', sw4_exe)
    mpirun_base = 'mpirun'
    num_proc = 4
    mpirun_cmd = mpirun_base + ' -np ' + str(num_proc)
    print('mpirun_cmd = ', mpirun_cmd)
    sw4_mpi_run = mpirun_cmd + ' ' + sw4_exe
    #print('sw4_mpi_run = ', sw4_mpi_run)

    test_num=0
    
    # run all tests in the twilight directory
    test_dir = 'twilight'
    #for test_dir in ['twilight', 'lamb']: # this won't work because the base_cases are different in 'twilight' and 'lamb'
    
    #make a local test directory
    if not os.path.exists(test_dir):
        os.mkdir(test_dir)

    os.chdir(test_dir) # change to the new local directory

    result_file = 'TwilightErr.txt'

    #base_case = 'flat-twi'
    for base_case in ['flat-twi', 'gauss-twi']:
        # put all test cases in a list and add a for loop for each test_dir
        for i in [1,2,3]:
            test_num = test_num+1
        
            case_dir = base_case + '-' + str(i)
            test_case = case_dir + '.in'

            reference_dir = pytest_dir + '/reference' + sep + test_dir + sep + case_dir
            #print('reference_dir=', reference_dir)
    
            sw4_input_file = examples_dir + sep + test_dir + sep + test_case
            #print('sw4_input_file = ', sw4_input_file)

            sw4_stdout_file = case_dir + '.out'
            
            local_dir = pytest_dir + sep + test_dir
            #print('local_dir = ', local_dir)
    
            # pipe stdout and stderr to a temporary file
            run_cmd = sw4_mpi_run + ' ' + sw4_input_file + ' >& ' + sw4_stdout_file
            print('Test #', test_num, 'run_cmd:', run_cmd)

            # run sw4
            run_dir = os.getcwd()
            #print('Running sw4 from directory:', run_dir)
            status = os.system(run_cmd)
            print('Test #', test_num, 'execution status = ', status)
            print('Test #', test_num, 'output dirs: local case_dir =', case_dir, 'reference_dir =', reference_dir)

            # compare output
            success = compare_to_ref(reference_dir + sep + result_file, case_dir + sep + result_file, 1e-5, 1e-5)
            print('Test #', test_num, 'compare_to_ref returned:', success)

        # end for all base_cases

    # end for all cases in the test_dir
    
    os.chdir('..') # change back to the parent directory

#------------------------------------------------
if __name__ == "__main__":
   main_test()

