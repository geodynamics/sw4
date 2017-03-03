#!/usr/bin/env python3

# Arguments:
# This assumes displacement variables and attenuation variables have the same tolerance

import os, sys, argparse

#----(Currently not used)--------------------------------------------
def run_checks(checks):
    fail = False
    for elem in checks:
        pass_check = (checks[elem][0] > checks[elem][1])
        print("%(check)s    tolerance: %(tol)f   calculated: %(val)f    pass: %(pass)s" %
                {"check": elem, "tol": checks[elem][0], "val": checks[elem][1], "pass": pass_check})
        if not pass_check: fail = True
    return fail

#------------------------------------------------
def compare_last_line(base_file_name, test_file_name, errTol, absErrLimit, verbose):

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
        print("WARNING: " + base_file_name + " and " + test_file_name + " are different lengths.")


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
            if abs(b0) > absErrLimit: 
                re0 = abs(b0-t0)/abs(b0)
            else:
                re0 = abs(b0-t0)

            if verbose or re0 > errTol:
                print('INFO: compare_last_line: col=', jj, 'test=', t0, 'base=', b0, 'err=', re0);

            if re0 > errTol:
                print('ERROR: compare_last_line: err=', re0, '> tolerance=', errTol)
                success = False
            # end if
        # end for
    except:
        success = False

    return success

#------------------------------------------------
def guess_mpi_cmd(mpi_tasks, verbose):
    node_name = os.uname().nodename
    sys_name = os.uname().sysname
    if verbose: print('node_name=', node_name)

    if 'quartz' in node_name:
        if mpi_tasks<=0: mpi_tasks = 36
        mpirun_cmd="srun -ppdebug -n " + str(mpi_tasks)
    elif 'fourier' in node_name:
        if mpi_tasks<=0: mpi_tasks = 4
        mpirun_cmd="mpirun -np " + str(mpi_tasks)
    # add more machine names here
    elif 'Linux' in sys_name:
        if mpi_tasks<=0: mpi_tasks = 1
        mpirun_cmd="mpirun -np " + str(mpi_tasks)
    else:
        #default mpi command
        if mpi_tasks<=0: mpi_tasks = 1
        mpirun_cmd="mpirun -np " + str(mpi_tasks)

    if verbose: print('mpirun_cmd = ', mpirun_cmd)

    return mpirun_cmd

#------------------------------------------------
def main_test(testing_level=0, mpi_tasks=0, verbose=False):
    sep = '/'
    pytest_dir = os.getcwd()
    pytest_dir_list = pytest_dir.split(sep)
    sw4_base_list = pytest_dir_list[:-1] # discard the last sub-directory (pytest)

    sw4_base_dir = sep.join(sw4_base_list)
    optimize_dir =  sw4_base_dir + '/optimize'
    reference_dir = pytest_dir + '/reference' 

    if verbose: print('pytest_dir =', pytest_dir)
    if verbose: print('sw4_base_dir =', sw4_base_dir)
    if verbose: print('optimize_dir =', optimize_dir)          
    if verbose: print('reference_dir =', reference_dir)          
    
    sw4_exe = optimize_dir + '/sw4'
    #print('sw4-exe = ', sw4_exe)

    # guess the mpi run command from the uname info
    mpirun_cmd=guess_mpi_cmd(mpi_tasks, verbose)

    sw4_mpi_run = mpirun_cmd + ' ' + sw4_exe
    #print('sw4_mpi_run = ', sw4_mpi_run)

    num_test=0
    num_pass=0
    num_fail=0

    all_dirs = ['pointsource', 'twilight', 'twilight', 'lamb']
    all_cases = ['pointsource-sg', 'flat-twi', 'gauss-twi', 'lamb']
    all_results =['PointSourceErr.txt', 'TwilightErr.txt', 'TwilightErr.txt', 'LambErr.txt']
    num_meshes =[1, 2, 2, 1] # default number of meshes for level 0

    # add more tests for higher values of the testing level
    if verbose: print("Testing level=", testing_level)
    if testing_level == 1:
        num_meshes =[2, 3, 3, 2]
    elif testing_level == 2:
        num_meshes =[3, 3, 3, 3]
    
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
            if verbose: 
                print('Starting test #', num_test, 'in directory:', test_dir, 'with input file:', test_case)

            sw4_input_file = reference_dir + sep + test_dir + sep + test_case
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
                print('ERROR: Test', test_case, ': sw4 returned non-zero exit status=', status)
                print('       run_cmd=', run_cmd)
                print('Test #', num_test, "Input file:", test_case, 'FAILED')
                num_fail += 1
                continue #skip to next test

            ref_result = reference_dir + sep + test_dir + sep + case_dir + sep + result_file
            #print('Test #', num_test, 'output dirs: local case_dir =', case_dir, 'ref_result =', ref_result)

            # compare output
            success = compare_last_line(ref_result , case_dir + sep + result_file, 1e-5, 1e-10, verbose)
            if success:        
                print('Test #', num_test, "Input file:", test_case, 'PASSED')
                num_pass += 1
            else:
                print('Test #', num_test, "Input file:", test_case, 'FAILED')
                num_fail += 1
            
        # end for qq in all_dirs[qq]

        os.chdir('..') # change back to the parent directory

    # end for all cases in the test_dir
    print('Out of', num_test, 'tests,', num_pass, 'PASSED and', num_fail, 'Failed')
    
#------------------------------------------------
if __name__ == "__main__":
    # default arguments
    testing_level=0
    verbose=False
    mpi_tasks=0 # machine dependent default

    parser=argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
    parser.add_argument("-l", "--level", type=int, choices=[0, 1, 2], 
                        help="testing level")
    parser.add_argument("-m", "--mpitasks", type=int, help="number of mpi tasks")
    args = parser.parse_args()
    if args.verbose:
        #print("verbose mode enabled")
        verbose=True
    if args.level:
        #print("Testing level specified=", args.level)
        testing_level=args.level
    if args.mpitasks:
        #print("MPI-tasks specified=", args.mpitasks)
        if args.mpitasks > 0: mpi_tasks=args.mpitasks

    main_test(testing_level, mpi_tasks, verbose)

