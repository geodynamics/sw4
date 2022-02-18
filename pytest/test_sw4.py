#!/usr/bin/env python3

# Arguments:
# -h: help, -v: verbose mode -l testing level, -m mpi-tasks, -d sw4-exe-dir -t omp-threads

import os, sys, argparse, subprocess

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
def compare_energy(test_file_name, errTol, verbose):

    success = True;

    f = open(test_file_name)

    Lnum=0;
    for line in f:
        Lnum = Lnum+1;
        if Lnum <= 3:
            thisEnergy = float(line);
            thirdEnergy = thisEnergy;
            if verbose and Lnum == 3:
                print('INFO: compare_energy: Ref. Energy =', thirdEnergy);

        else:
            prevEnergy = thisEnergy;
            thisEnergy = float(line);
# relative change from previous energy
            diff0 = (thisEnergy - prevEnergy)/thirdEnergy;

            if diff0 > errTol:
                print('ERROR: compare_energy: line =', Lnum, 'prev =', prevEnergy, 'this =', thisEnergy, 'rel diff =', diff0, '> tolerance=', errTol);
                success = False
                break            
            # end if
        # end if
    # end for
    if (Lnum < 4):
        print("ERROR: compare_energy: Less than 4 lines in the energy log!");
        success=False;
    elif verbose:
        print('INFO: compare_energy: line =', Lnum, 'prev =', prevEnergy, 'this =', thisEnergy, 'rel diff =', diff0);

    f.close()
    
    return success

#------------------------------------------------
def guess_mpi_cmd(mpi_tasks, omp_threads, cpu_allocation, verbose):
    if verbose: print('os.uname=', os.uname())
    node_name = os.uname()[1]
    if verbose: print('node_name=', node_name)
    sys_name = os.uname()[0]
    if verbose: print('sys_name=', sys_name)
    nersc_sys = os.getenv('NERSC_HOST')


    if 'quartz' in node_name:
        if omp_threads<=0: omp_threads=2;
        if mpi_tasks<=0: mpi_tasks = int(36/omp_threads)
        # the following setting is needed to combine h5py and subprocess.run on LC
        os.environ["PSM2_DEVICES"] = ""
        if cpu_allocation == "":
           mpirun_cmd="srun -ppdebug " + " -n " + str(mpi_tasks) + " -c " + str(omp_threads)
        else:
           mpirun_cmd="srun -ppdebug " + " -A " + cpu_allocation + " -n " + str(mpi_tasks) + " -c " + str(omp_threads)
    elif 'cab' in node_name:
        if omp_threads<=0: omp_threads=2;
        if mpi_tasks<=0: mpi_tasks = int(16/omp_threads)
        mpirun_cmd="srun -ppdebug -n " + str(mpi_tasks) + " -c " + str(omp_threads)
    elif 'nid' in node_name: # the cori knl nodes are called nid
        if omp_threads<=0: omp_threads=4;
        if mpi_tasks<=0: mpi_tasks = int(64/omp_threads) # for KNL nodes, use 64 hardware cores per node
        sw_threads = 4*omp_threads # Cori uses hyperthreading by default
        if mpi_tasks<=0: mpi_tasks = int(32/omp_threads) # for Haswell nodes
        sw_threads = omp_threads 
        mpirun_cmd="srun --cpu_bind=cores -n " + str(mpi_tasks) + " -c " + str(sw_threads)
        # For Perlmutter
        if nersc_sys == 'perlmutter':
            mpirun_cmd="srun -N 1 -n 4 -c 32 --gpu-bind=single:1 --cpu-bind=cores  --gpus-per-task=1 "
    elif 'fourier' in node_name:
        if omp_threads<=0: omp_threads=1;
        if mpi_tasks<=0: mpi_tasks = 4
        mpirun_cmd="mpirun -np " + str(mpi_tasks)
    elif 'batch' in node_name: # for summit
        if omp_threads<=0: omp_threads=7;
        if mpi_tasks<=0: mpi_tasks = 6
        mpirun_cmd="jsrun -a1 -c7 -r6 -l CPU-CPU -d packed -b packed:7 -n " + str(mpi_tasks)
        mpirun_cmd="jsrun -a1 -c7 -g1 -l CPU-CPU -d packed -b packed:7 -M -gpu -n " + str(mpi_tasks)
    elif 'rzansel' in node_name:
        os.environ["PSM2_DEVICES"] = ""
        if mpi_tasks<=0: mpi_tasks = 4
        mpirun_cmd="lrun -T4 "
    elif 'rznevada' in node_name:
        os.environ["PSM2_DEVICES"] = ""
        if mpi_tasks<=0: mpi_tasks = 4
        mpirun_cmd="srun -N2 -n2 "
    elif 'redwood' in node_name:
        os.environ["PSM2_DEVICES"] = ""
        if mpi_tasks<=0: mpi_tasks = 4
        mpirun_cmd="srun -p amdMI100 -n1 "
    elif 'login' in node_name:
        #os.environ["PSM2_DEVICES"] = ""
        #os.environ["ROMIO_FSTYPE_FORCE"]="ufs:"
        os.environ["MPICH_GPU_SUPPORT_ENABLED"]="1"
        if mpi_tasks<=0: mpi_tasks = 4
        mpirun_cmd="srun -N 1 -A GEO130_crusher -t 30 -n 8 -c8 --gpus-per-task=1 --gpu-bind=closest -p batch "
    elif 'lassen' in node_name:
        os.environ["PSM2_DEVICES"] = ""
        if mpi_tasks<=0: mpi_tasks = 4
        mpirun_cmd="lrun -T4 "
        #mpirun_cmd="jsrun -g4 -c40 -a4 -n" + str(mpi_tasks) # Simulate Summit runs with -g4
    # add more machine names here
    elif 'Linux' in sys_name:
        if omp_threads<=0: omp_threads=1;
        if mpi_tasks<=0: mpi_tasks = 1
        mpirun_cmd="srun -n 4 "
    elif 'redwood' in sys_name:
        mpirun_cmd="srun -N 1 -n 4 -c 32 -p amdMI100 ./bind "
    elif 'Darwin' in sys_name:
        if omp_threads<=0: omp_threads=1;
        if mpi_tasks<=0: mpi_tasks = 4
        mpirun_cmd="mpirun -np " + str(mpi_tasks)
    else:
        #default mpi command
        if omp_threads<=0: omp_threads=1;
        if mpi_tasks<=0: mpi_tasks = 1
        mpirun_cmd="mpirun -np " + str(mpi_tasks)

    return mpirun_cmd

#------------------------------------------------
def main_test(sw4_exe_dir="optimize_mp", pytest_dir ="none", testing_level=0, mpi_tasks=0, omp_threads=0, cpu_allocation="", verbose=False, usehdf5=False):

    assert sys.version_info >= (3,5) # named tuples in Python version >=3.3

    sep = '/'
    if pytest_dir == "none":
        pytest_dir = os.getcwd()

    pytest_dir_list = pytest_dir.split(sep)
    sw4_base_list = pytest_dir_list[:-1] # discard the last sub-directory (pytest)

    sw4_base_dir = sep.join(sw4_base_list)
    optimize_dir =  sw4_base_dir + sep + sw4_exe_dir
    reference_dir = pytest_dir + '/reference'

    # make sure the directories are there
    if not os.path.isdir(sw4_base_dir):
        print("ERROR: directory", sw4_base_dir, "does not exists")
        return False
    if not os.path.isdir(optimize_dir):
        print("ERROR: directory", optimize_dir, "does not exists (HINT: use -d 'sw4_exe_dir' or -p 'pytest_dir')")
        return False
    if not os.path.isdir(reference_dir):
        print("ERROR: directory", reference_dir, "does not exists")
        return False
    
    if verbose: print('pytest_dir =', pytest_dir)
    if verbose: print('sw4_base_dir =', sw4_base_dir)
    if verbose: print('optimize_dir =', optimize_dir)          
    if verbose: print('reference_dir =', reference_dir)          
    
    sw4_exe = optimize_dir + '/sw4'

    #print('sw4-exe = ', sw4_exe)

    # make sure sw4 is present in the optimize dir
    if not os.path.isfile(sw4_exe):
        print("ERROR: the file", sw4_exe, "does not exists (DID YOU FORGET TO BUILD SW4?)")
        return False

    # guess the mpi run command from the uname info
    mpirun_cmd=guess_mpi_cmd(mpi_tasks, omp_threads, cpu_allocation, verbose)

    sw4_mpi_run = mpirun_cmd + ' ' + sw4_exe
    if (omp_threads>0):
        os.putenv("OMP_NUM_THREADS", str(omp_threads))

    if verbose: print('sw4_mpi_run = ', sw4_mpi_run)

    num_test=0
    num_skip=0
    num_pass=0
    num_fail=0

    all_dirs = ['energy', 'energy', 'energy', 'energy',
                'meshrefine', 'meshrefine', 'meshrefine',
                'attenuation', 'attenuation', 'pointsource',
                'twilight', 'twilight', 'lamb',
                'curvimeshrefine', 'curvimeshrefine', 'curvimeshrefine','curvimeshrefine',
                'hdf5', 'hdf5', 'hdf5']

    all_cases = ['energy-nomr-2nd', 'energy-mr-4th', 'energy-mr-sg-order2', 'energy-mr-sg-order4', 'refine-el', 'refine-att', 'refine-att-2nd', 'tw-att', 'tw-topo-att', 'pointsource-sg', 'flat-twi', 'gauss-twi', 'lamb','gausshill-el','gauss-sg-mr','energy','gausshill-att','loh1-h100-mr-hdf5','loh1-h100-mr-hdf5-sfile', 'loh1-h100-mr-restart-hdf5']
    
    all_results =['energy.log', 'energy.log', 'energy.log', 'energy.log', 'TwilightErr.txt', 'TwilightErr.txt', 'TwilightErr.txt', 'TwilightErr.txt', 'TwilightErr.txt', 'PointSourceErr.txt', 'TwilightErr.txt', 'TwilightErr.txt', 'LambErr.txt','TwilightErr.txt','TwilightErr.txt','energy.log','TwilightErr.txt','hdf5.log','hdf5-sfile.log', 'hdf5-restart.log']
    num_meshes =[1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1] # default number of meshes for level 0

    # add more tests for higher values of the testing level
    if testing_level == 1:
        num_meshes =[1, 1, 1, 1, 2, 2, 2, 3, 2, 2, 3, 3, 2, 2, 2, 1, 2, 1, 1, 1]
    elif testing_level == 2:
        num_meshes =[1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 1, 3, 1, 1, 1]

    
    print("Running all tests for level", testing_level, "...")
    # run all tests
    for qq in range(len(all_dirs)):
   
        test_dir = all_dirs[qq]
        base_case = all_cases[qq]
        result_file = all_results[qq]

        # skip HDF5 test if specified
        if ('hdf5' not in base_case and usehdf5 != 4 and usehdf5 != 0) or (base_case == 'loh1-h100-mr-hdf5' and usehdf5 != 4 and usehdf5 != 1) or (base_case == 'loh1-h100-mr-hdf5-sfile' and usehdf5 != 4 and usehdf5 != 2) or (base_case == 'loh1-h100-mr-restart-hdf5' and usehdf5 != 4 and usehdf5 != 3):
            num_skip = num_skip+1
            num_test = num_test+1
            print('Test #', num_test, "Input file:", base_case+'-'+str(1)+'.in', 'SKIPPED')
        else:
 
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

                # sw4_stdout_file = case_dir + '.out'
                
                local_dir = pytest_dir + sep + test_dir
                #print('local_dir = ', local_dir)
                # pipe stdout and stderr to a temporary file
                run_cmd = mpirun_cmd.split() + [
                    sw4_exe,
                    sw4_input_file
                ]

                sw4_stdout_file = open(case_dir + '.out', 'wt')
                sw4_stderr_file = open(case_dir + '.err', 'wt')

                # pipe stdout and stderr to a temporary file
#                run_cmd = sw4_mpi_run + ' ' + sw4_input_file + ' >& ' + sw4_stdout_file

                # run sw4
                run_dir = os.getcwd()
                if verbose:
                    print('Running sw4 from directory:', run_dir)
                    print('run_cmd=', run_cmd)

# assign     OMP_NUM_THREADS
                status = subprocess.run(
                    run_cmd,
                    stdout=sw4_stdout_file,
                    stderr=sw4_stderr_file
                )

                sw4_stdout_file.close()
                sw4_stderr_file.close()

                # status = os.system(run_cmd)
                # if status!=0:
                #     print('ERROR: Test', test_case, ': sw4 returned non-zero exit status=', status, 'aborting test')
                #     print('run_cmd=', run_cmd)
                #     print("DID YOU USE THE CORRECT SW4 EXECUTABLE? (SPECIFY DIRECTORY WITH -d OPTION)")
                #     return False # bail out

                if status.returncode!=0:
                    print('ERROR: Test', test_case, ': sw4 returned non-zero exit status=', status.returncode, 'aborting test')
                    print('run_cmd=', run_cmd)
                    print("DID YOU USE THE CORRECT SW4 EXECUTABLE? (SPECIFY DIRECTORY WITH -d OPTION)")
                    return False # bail out


                ref_result = reference_dir + sep + test_dir + sep + case_dir + sep + result_file
                #print('Test #', num_test, 'output dirs: local case_dir =', case_dir, 'ref_result =', ref_result)


                if result_file == 'hdf5.log':
                    import verify_hdf5
                    success = True
                    sw4_stdout_file_tmp = open(case_dir + '.out', 'r')
                    for line in sw4_stdout_file_tmp:
                        if "Error" in line or "not compiled with HDF5" in line or "without sw4 compiled with HDF5" in line:
                            success = False
                            print('  SW4 is not compiled with HDF5 library!')
                            break;
                    sw4_stdout_file_tmp.close()
                    if success == True:
                        success = verify_hdf5.verify(pytest_dir, 1e-5)
                elif result_file == 'hdf5-sfile.log' or result_file == 'hdf5-restart.log':
                    import verify_hdf5
                    success = True
                    sw4_stdout_file_tmp = open(case_dir + '.out', 'r')
                    for line in sw4_stdout_file_tmp:
                        if "Error" in line or "not compiled with HDF5" in line or "without sw4 compiled with HDF5" in line:
                            success = False
                            print('  SW4 is not compiled with HDF5 library!')
                            break;
                    sw4_stdout_file_tmp.close()
                    if success == True:
                        success = verify_hdf5.verify_sac_image(pytest_dir, 1e-5)
                elif result_file == 'energy.log':
                    success = compare_energy(case_dir + sep + result_file, 1e-10, verbose)
                else:
                    # compare output with reference result (always compare the last line)
                    success = compare_one_line(ref_result , case_dir + sep + result_file, 1e-5, 1e-10, -1, verbose)
                    if success and 'attenuation' in test_dir:
                        # also compare the 3rd last line in the files
                        success = compare_one_line(ref_result , case_dir + sep + result_file, 1e-5, 1e-10, -3, verbose)

                if success:        
                    print('Test #', num_test, "Input file:", test_case, 'PASSED')
                    num_pass += 1
                else:
                    print('Test #', num_test, "Input file:", test_case, 'FAILED')
                    num_fail += 1
                
            os.chdir('..') # change back to the parent directory
            # end for ii in range(num_meshes[qq]):

        # end for qq in all_dirs[qq]

    # end for all cases in the test_dir
    print('Out of', num_test, 'tests,', num_fail, 'failed,', num_pass, 'passed, and', num_skip, 'skipped')
    # normal termination
    return True
    
#------------------------------------------------
if __name__ == "__main__":
    assert sys.version_info >= (3,5) # named tuples in Python version >=3.3
    # default arguments
    testing_level=0
    verbose=False
    usehdf5=4
    mpi_tasks=0 # machine dependent default
    omp_threads=0 #no threading by default
    cpu_allocation=""

    parser=argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
    parser.add_argument("-l", "--level", type=int, choices=[0, 1, 2], 
                        help="testing level")
    parser.add_argument("-m", "--mpitasks", type=int, help="number of mpi tasks")
    parser.add_argument("-t", "--ompthreads", type=int, help="number of omp threads per task")
    parser.add_argument("-d", "--sw4_exe_dir", help="name of directory for sw4 executable", default="optimize_mp")
    parser.add_argument("-p", "--pytest_dir", help="full path to the directory of pytest (/path/sw4/pytest)", default="none")
    parser.add_argument("-u", "--usehdf5", type=int, choices=[0, 1, 2, 3, 4], help="run HDF5 tests with, 0 no HDF5, 1 first HDF5 case, 2 second case,..., 4 all cases")
    parser.add_argument("-A","--cpu_allocation", help="name of cpu bank/allocation",default="")
    args = parser.parse_args()

    if args.usehdf5 == 0:
        usehdf5=0
    elif args.usehdf5 == 1:
        usehdf5=1
    elif args.usehdf5 == 2:
        usehdf5=2
    elif args.usehdf5 == 3:
        usehdf5=3
    elif args.usehdf5 == 4:
        usehdf5=4

    if args.verbose:
        #print("verbose mode enabled")
        verbose=True
    if args.level:
        #print("Testing level specified=", args.level)
        testing_level=args.level
    if args.mpitasks:
        #print("MPI-tasks specified=", args.mpitasks)
        if args.mpitasks > 0: mpi_tasks=args.mpitasks
    if args.ompthreads:
        #print("OMP-threads specified=", args.ompthreads)
        if args.ompthreads > 0: omp_threads=args.ompthreads
    if args.pytest_dir:
        #print("sw4_exe specified=", args.sw4_exe)
        pytest_dir=args.pytest_dir
    if args.sw4_exe_dir:
        #print("sw4_exe_dir specified=", args.sw4_exe_dir)
        sw4_exe_dir=args.sw4_exe_dir
    if args.cpu_allocation:
        #print("cpu_allocation specified=", args.cpu_allocation)
        cpu_allocation=args.cpu_allocation

    if not main_test(sw4_exe_dir, pytest_dir, testing_level, mpi_tasks, omp_threads, cpu_allocation, verbose, usehdf5):
        print("test_sw4 was unsuccessful")
        sys.exit(-1)
