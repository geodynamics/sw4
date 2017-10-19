#!/usr/bin/env python3

# Arguments:
# -h: help, -v: verbose mode -l testing level, -m mpi-tasks, -d sw4-exe-dir

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
def compare_one_line(base_file_name, test_file_name, errTol, absErrLimit, lineNum, verbose):

    success = True

    base_file = open(base_file_name)
    test_file = open(test_file_name)
    base_data = base_file.readlines() #reads all lines in the file
    test_data = test_file.readlines()

    # tmp
    #print('base_data file:', base_data)
    #print('test_data file:', test_data)

    base_file.close()
    test_file.close()
    if len(base_data) != len(test_data):
        print("WARNING: " + base_file_name + " and " + test_file_name + " are different lengths.")

    # for twilight tests, compare the 3 numbers on the last line of each file
    base_line = base_data[lineNum]
    test_line = test_data[lineNum]
    #print('base_line=', base_line)
    #print('test_line=', test_line)

    base_data = base_line.split() #base_data is a list of strings
    test_data = test_line.split()
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
                print('INFO: compare_one_line: col=', jj, 'test=', t0, 'base=', b0, 'err=', re0);

            if re0 > errTol:
                print('ERROR: compare_one_line: err=', re0, '> tolerance=', errTol)
                success = False
            # end if
        # end for
    except:
        success = False

    base_file.close()
    test_file.close()
    
    return success

#------------------------------------------------
def guess_mpi_cmd(mpi_tasks, verbose):
    if verbose: print('os.uname=', os.uname())
    node_name = os.uname()[1]
    if verbose: print('node_name=', node_name)
    sys_name = os.uname()[0]
    if verbose: print('sys_name=', sys_name)

    if 'quartz' in node_name:
        if mpi_tasks<=0: mpi_tasks = 36
        mpirun_cmd="srun -ppdebug -n " + str(mpi_tasks)
    elif 'cab' in node_name:
        if mpi_tasks<=0: mpi_tasks = 16
        mpirun_cmd="srun -ppdebug -n " + str(mpi_tasks)
    elif 'nid' in node_name:
        # all KNL nodes on cori have a node name starting with 'nid'
        if mpi_tasks<=0: mpi_tasks = 64
        mpirun_cmd="srun -c 4 --cpu_bind=cores -n " + str(mpi_tasks)
    elif 'fourier' in node_name:
        if mpi_tasks<=0: mpi_tasks = 4
        mpirun_cmd="mpirun -np " + str(mpi_tasks)
    elif 'ray' in node_name:
        if mpi_tasks<=0: mpi_tasks = 4
        mpirun_cmd="mpirun -np " + str(mpi_tasks) 
        if mpi_tasks > 1: mpirun_cmd += " mpibind"
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
def main_test(sw4lite_exe_dir="optimize_cuda_ray", testing_level=0, mpi_tasks=0, verbose=False):
    assert sys.version_info >= (3,2) # named tuples in Python version >=3.3
    sep = '/'
    pytest_dir = os.getcwd()
    pytest_dir_list = pytest_dir.split(sep)
    sw4lite_base_list = pytest_dir_list[:-1] # discard the last sub-directory (pytest)

    sw4lite_base_dir = sep.join(sw4lite_base_list)
    optimize_dir =  sw4lite_base_dir + sep + sw4lite_exe_dir
    reference_dir = pytest_dir + '/reference'

    # make sure the directories are there
    if not os.path.isdir(sw4lite_base_dir):
        print("ERROR: directory", sw4_base_dir, "does not exists")
        return False
    if not os.path.isdir(optimize_dir):
        print("ERROR: directory", optimize_dir, "does not exists (HINT: use -d 'sw4_exe_dir')")
        return False
    if not os.path.isdir(reference_dir):
        print("ERROR: directory", reference_dir, "does not exists")
        return False
    
    if verbose: print('pytest_dir =', pytest_dir)
    if verbose: print('sw4_base_dir =', sw4lite_base_dir)
    if verbose: print('optimize_dir =', optimize_dir)          
    if verbose: print('reference_dir =', reference_dir)          
    
    sw4lite_exe = optimize_dir + '/sw4lite'
    #print('sw4-exe = ', sw4_exe)

    # make sure sw4 is present in the optimize dir
    if not os.path.isfile(sw4lite_exe):
        print("ERROR: the file", sw4lite_exe, "does not exists (DID YOU FORGET TO BUILD SW4?)")
        return False

    # guess the mpi run command from the uname info
    mpirun_cmd=guess_mpi_cmd(mpi_tasks, verbose)

    sw4lite_mpi_run = mpirun_cmd + ' ' + sw4lite_exe
    #print('sw4_mpi_run = ', sw4_mpi_run)

    num_test=0
    num_pass=0
    num_fail=0

    num_meshes =[1, 1, 1] # default number of meshes for level 0

# Run in pytest/test_dir
    test_dirs   =['pointsource','loh1','loh2']
# Read input file from pytest/reference/test_dir/input_file
    input_files =['pointsource.in','LOH.1-h100.in','LOH.2-h100.in']
# Computed results files in pytest/test_dir/result_file
# Find reference results in pytest/reference/test_dir/result_file
    result_files=['pointsource-h0p04/PointSourceErr.txt','LOH1-h100/sta10.txt','LOH2-h100/sta10.txt']

    print("Running all tests for level", testing_level, "...")
    # run all tests
    for qq in range(len(test_dirs)):
        #make a local test directory
        if not os.path.exists(test_dirs[qq]):
            os.mkdir(test_dirs[qq])

        os.chdir(test_dirs[qq]) # change to the new local directory

        num_test = num_test+1

        if verbose: 
            print('Starting test in directory:', test_dirs[qq], 'with input file:', input_files[qq])

        sw4_input_file = reference_dir + sep + test_dirs[qq] + sep + input_files[qq]
        # print('sw4_input_file = ', sw4_input_file)

        sw4_stdout_file = test_dirs[qq] + str(mpi_tasks) + '.out'
            
        # pipe stdout and stderr to a temporary file
        run_cmd = sw4lite_mpi_run + ' ' + sw4_input_file + ' >& ' + sw4_stdout_file
        #print('run_cmd = ', run_cmd)

        # run sw4lite
        run_dir = os.getcwd()
        #print('Running sw4lite from directory:', run_dir)
        status = os.system(run_cmd)
        if status!=0:
            print('ERROR: Test in dir:', test_dirs[qq], 'Input file:', input_files[qq], ': sw4 returned non-zero exit status=', status, 'aborting test')
            print('run_cmd=', run_cmd)
            print("DID YOU USE THE CORRECT SW4LITE EXECUTABLE? (SPECIFY DIRECTORY WITH -d OPTION)")
            return False # bail out

        ref_result = reference_dir + sep + test_dirs[qq] + sep + result_files[qq]
        #print('ref result in ',ref_result)

            # compare output (always compare the last line)
        success = compare_one_line(ref_result , result_files[qq], 1e-5, 1e-10, -1, verbose)
        if success:        
            print('Test #', num_test, "Input file:", input_files[qq], 'PASSED')
            num_pass += 1
        else:
           print('Test #', num_test, "Input file:", input_files[qq], 'FAILED')
           num_fail += 1
        os.chdir('..') # change back to the parent directory

    # end for all cases in the test_dir
    print('Out of', num_test, 'tests,', num_fail, 'failed and ', num_pass, 'passed')
    # normal termination
    return True
    
#------------------------------------------------
if __name__ == "__main__":
    assert sys.version_info >= (3,2) # named tuples in Python version >=3.3
    # default arguments
    testing_level=0
    verbose=False
    mpi_tasks=0 # machine dependent default

    parser=argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
    parser.add_argument("-l", "--level", type=int, choices=[0, 1, 2], 
                        help="testing level")
    parser.add_argument("-m", "--mpitasks", type=int, help="number of mpi tasks")
    parser.add_argument("-d", "--sw4_exe_dir", help="name of directory for sw4 executable", default="optimize")
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
    if args.sw4_exe_dir:
        #print("sw4_exe_dir specified=", args.sw4_exe_dir)
        sw4_exe_dir=args.sw4_exe_dir

    if not main_test(sw4_exe_dir, testing_level, mpi_tasks, verbose):
        print("test_sw4 was unsuccessful")

