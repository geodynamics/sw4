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
def compare_to_ref(base_file_name, test_file_name, errTol):

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
        for jj in range(len(base_data)):
            t0 = test_num[jj]
            b0 = base_num[jj]
            re0 = 0
            if abs(b0) > 0: 
                re0 = abs(b0-t0)/abs(b0)
            else:
                re0 = abs(b0-t0)

            # tmp
            print('col=', jj, 'test=', t0, 'base=', b0, 'err=', re0);
            # end tmp
            if re0 > errTol:
                print('ERROR: compare_to_ref, base_data=', base_data, 'test_data=', test_data)
                success = False
            # end if
        # end for
    except:
        success = False

    return success

#------------------------------------------------
def guess_mpi_cmd():
    node_name = os.uname().nodename
    sys_name = os.uname().sysname
    print('node_name=', node_name)

    if 'quartz' in node_name:
        mpirun_cmd="srun -ppdebug -n 36"
    elif 'fourier' in node_name:
        mpirun_cmd="mpirun -np 4"
    # add more machine names here
    elif 'Linux' in sys_name:
        mpirun_cmd="mpirun -np 1"
    else:
        #default mpi command
        mpirun_cmd="mpirun -np 1"

    return mpirun_cmd

#------------------------------------------------
def main_test(testing_level=0):
    sep = '/'
    pytest_dir = os.getcwd()
    print('pytest_dir =', pytest_dir)
    pytest_dir_list = pytest_dir.split(sep)
    sw4_base_list = pytest_dir_list[:-1] # discard the last sub-directory (pytest)

    sw4_base_dir = sep.join(sw4_base_list)
    examples_dir = sw4_base_dir + '/examples'
    optimize_dir =  sw4_base_dir + '/optimize'

    print('sw4_base_dir =', sw4_base_dir)
    print('examples_dir =', examples_dir)
    print('optimize_dir =', optimize_dir)          
    
    sw4_exe = optimize_dir + '/sw4'
    #print('sw4-exe = ', sw4_exe)

    # guess the mpi run command from the uname info
    mpirun_cmd=guess_mpi_cmd()

    print('mpirun_cmd = ', mpirun_cmd)
    sw4_mpi_run = mpirun_cmd + ' ' + sw4_exe
    #print('sw4_mpi_run = ', sw4_mpi_run)

    num_test=0
    num_pass=0
    num_fail=0

    all_dirs = ['twilight', 'twilight', 'lamb']
    all_cases = ['flat-twi', 'gauss-twi', 'lamb']
    all_results =['TwilightErr.txt', 'TwilightErr.txt', 'LambErr.txt']

    # make num_meshes depend on the testing level
    print("Testing level=", testing_level)
    if testing_level == 0:
        num_meshes =[2, 2, 1]
    elif testing_level == 1:
        num_meshes =[3, 3, 2]
    
    # run all tests
    for qq in range(len(all_dirs)):
    
        test_dir = all_dirs[qq]
        base_case = all_cases[qq]
        result_file = all_results[qq]

        #make a local test directory
        if not os.path.exists(test_dir):
            os.mkdir(test_dir)

        os.chdir(test_dir) # change to the new local directory

        # put all test cases in a list and add a for loop for each test_dir
        for ii in range(num_meshes[qq]):
            num_test = num_test+1
        
            case_dir = base_case + '-' + str(ii+1)
            test_case = case_dir + '.in'
            print('Starting test #', num_test, 'in directory:', test_dir, 'with input file:', test_case)

            reference_dir = pytest_dir + '/reference' + sep + test_dir + sep + case_dir
            #print('reference_dir=', reference_dir)
    
            sw4_input_file = examples_dir + sep + test_dir + sep + test_case
            #print('sw4_input_file = ', sw4_input_file)

            sw4_stdout_file = case_dir + '.out'
            
            local_dir = pytest_dir + sep + test_dir
            #print('local_dir = ', local_dir)
    
            # pipe stdout and stderr to a temporary file
            run_cmd = sw4_mpi_run + ' ' + sw4_input_file + ' >& ' + sw4_stdout_file

            # run sw4
            run_dir = os.getcwd()
            #print('Running sw4 from directory:', run_dir)
            status = os.system(run_cmd)
            if status!=0:
                print('WARNING: Test #', num_test, 'returned non-zero exe status=', status)
                print('       run_cmd=', run_cmd)
            #print('Test #', num_test, 'output dirs: local case_dir =', case_dir, 'reference_dir =', reference_dir)

            # compare output
            success = compare_to_ref(reference_dir + sep + result_file, case_dir + sep + result_file, 1e-5)
            if success:        
                print('Test #', num_test, 'PASSED')
                num_pass += 1
            else:
                print('Test #', num_test, 'FAILED')
                num_fail += 1
            
        # end for qq in all_dirs[qq]

        os.chdir('..') # change back to the parent directory

    # end for all cases in the test_dir
    print('Out of', num_test, 'tests,', num_pass, 'PASSED and', num_fail, 'Failed')
    
#------------------------------------------------
if __name__ == "__main__":
   main_test()

