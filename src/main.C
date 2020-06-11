//  SW4 LICENSE
// # ----------------------------------------------------------------------
// # SW4 - Seismic Waves, 4th order
// # ----------------------------------------------------------------------
// # Copyright (c) 2013, Lawrence Livermore National Security, LLC.
// # Produced at the Lawrence Livermore National Laboratory.
// #
// # Written by:
// # N. Anders Petersson (petersson1@llnl.gov)
// # Bjorn Sjogreen      (sjogreen2@llnl.gov)
// #
// # LLNL-CODE-643337
// #
// # All rights reserved.
// #
// # This file is part of SW4, Version: 1.0
// #
// # Please also read LICENCE.txt, which contains "Our Notice and GNU General
// Public License"
// #
// # This program is free software; you can redistribute it and/or modify
// # it under the terms of the GNU General Public License (as published by
// # the Free Software Foundation) version 2, dated June 1991.
// #
// # This program is distributed in the hope that it will be useful, but
// # WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
// # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
// # conditions of the GNU General Public License for more details.
// #
// # You should have received a copy of the GNU General Public License
// # along with this program; if not, write to the Free Software
// # Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA
//#include "mpi.h"

#include <mpi.h>
#include <omp.h>

#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "EW.h"
#include "Mspace.h"
#include "version.h"
#ifdef ENABLE_CUDA
#include "cuda_profiler_api.h"
#include "nvToolsExtCuda.h"
#endif

#if defined(SW4_SIGNAL_CHECKPOINT)
//
// Currently no way to get the singnal to all processes without killing the job
//
#include <csignal>
extern volatile std::sig_atomic_t signal_status;

void signal_handler(int signal) {
  signal_status = signal;
  std::cout << " RECEIVED SIGNAL " << signal << "\n";
}
#endif

using namespace std;

void usage(string thereason) {
  cout << endl
       << "sw4 - Summation by parts 4th order forward seismic wave propagator"
       << endl
       << endl
       << "Usage: sw4 [-v] file.in" << endl
       << "\t -v:      prints out the version info" << endl
       << "\t file.in: an input file" << endl
       << endl
       << "Reason for message: " << thereason << endl;
}

int main(int argc, char **argv) {
  // PROFILER_STOP;
  // cudaProfilerStop();
  int myRank = 0, nProcs = 0;
  string fileName;
  // bool checkmode = false;

  stringstream reason;
  presetGPUID();
  // Initialize MPI...
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

#if defined(SW4_SIGNAL_CHECKPOINT)
  std::signal(SIGUSR1, signal_handler);
#endif

#ifdef SW4_USE_UMPIRE
  umpire::ResourceManager &rma = umpire::ResourceManager::getInstance();
  auto allocator = rma.getAllocator("UM");
  // auto device_allocator = rma.getAllocator("DEVICE");

  const size_t pool_size = static_cast<size_t>(13) * 1024 * 1024 * 1024;

  auto pref_allocator = rma.makeAllocator<umpire::strategy::AllocationAdvisor>(
      "preferred_location_device", allocator, "PREFERRED_LOCATION",
      global_variables.device);

  auto pooled_allocator =
      rma.makeAllocator<umpire::strategy::DynamicPool, true>(
          string("UM_pool"), allocator, pool_size);

  const size_t pool_size_small = static_cast<size_t>(250) * 1024 * 1024;

  // This is a temporary workaround to the issue of Umpire always using device 0
  // for cudaMemAdvises using AllocationAdvisor.
  // if (global_variables.num_devices==1){

  // auto pooled_allocator_small =static_cast<size_t>(250)*1024*1024;
  auto pooled_allocator_small =
      rma.makeAllocator<umpire::strategy::DynamicPool, true>(
          string("UM_pool_temps"), pref_allocator, pool_size_small);

  const size_t object_pool_size = static_cast<size_t>(1024) * 1024 * 1024;

  // rma.makeAllocator<umpire::strategy::MonotonicAllocationStrategy,false>(string("UM_object_pool"),
  //					   object_pool_size,allocator);

  auto pooled_allocator_objects =
      rma.makeAllocator<umpire::strategy::DynamicPool, false>(
          string("UM_object_pool"), allocator, object_pool_size);

#ifdef SW4_MASS_PREFETCH
  std::cout << "Mass prefetch operational\n";
  global_variables.massprefetch.push_back(std::make_tuple(
      static_cast<char *>(pooled_allocator.allocate(1)), pool_size));
  global_variables.massprefetch.push_back(
      std::make_tuple(static_cast<char *>(pooled_allocator_small.allocate(1)),
                      pool_size_small));
  global_variables.massprefetch.push_back(
      std::make_tuple(static_cast<char *>(pooled_allocator_objects.allocate(1)),
                      object_pool_size));
#endif
  // rma.makeAllocator<umpire::strategy::MixedPool,false>(string("UM_object_pool"),
  //							   allocator,object_pool_size);

  // rma.makeAllocator<umpire::strategy::FixedPool,false>(string("UM_object_pool"),
  // 						     allocator,object_pool_size);

  // } else {
  //   auto pooled_allocator_small =
  //     rma.makeAllocator<umpire::strategy::DynamicPool,true>(string("UM_pool_temps"),
  //  							   pref_allocator,pool_size_small);
  // }

  // auto pooled_allocator2 =
  //   rma.makeAllocator<umpire::strategy::DynamicPool,false>(string("UM_pool_temps"),
  //                                                   allocator);
#endif

#ifdef ENABLE_TAU
  TAU_PROFILE_INIT(argc, argv);
#endif

  int status = 0;

  // mpi2 adds on four more args...  [-p4pg, dir, -p4wd, dir]
  int mpi2args = 4;

  if (argc != 2 && argc != 3 && argc != (2 + mpi2args) &&
      argc != (3 + mpi2args)) {
    reason << "Wrong number of args (1-2), not: " << argc - 1 << endl;

    if (myRank == 0) {
      for (int i = 0; i < argc; ++i)
        cout << "Argv[" << i << "] = " << argv[i] << endl;

      usage(reason.str());
    }
#ifdef USE_MAGMA
  if (magma_finalize()!=MAGMA_SUCCESS){
    std::cerr<<"ERROR MAGMA INIT FAILED\n";
    abort();
  }
#endif
    // Stop MPI
    MPI_Finalize();
    return 1;
  } else if (strcmp(argv[1], "-v") == 0) {
    if (myRank == 0) cout << ewversion::getVersionInfo() << endl;
    // Stop MPI
    MPI_Finalize();
    return status;
  } else if (argc == 1) {
    reason << "ERROR: ****No input file specified!" << endl;
    for (int i = 0; i < argc; ++i)
      reason << "Argv[" << i << "] = " << argv[i] << endl;

    if (myRank == 0) usage(reason.str());
    // Stop MPI
    MPI_Finalize();
    return 1;
  }

  else
    fileName = argv[1];

  if (myRank == 0) {
    cout << ewversion::getVersionInfo() << endl;
    cout << "Input file: " << fileName << endl;
  }

  MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

  cout.precision(8);
  // use sci format: 1.2345e-6
  cout << std::scientific;

  // Save the source description here
  vector<vector<Source *> > GlobalSources;
  // Save the time series here
  vector<vector<TimeSeries *> > GlobalTimeSeries;

// make a new simulation object by reading the input file 'fileName'
// nvtxRangePushA("outer");
#if defined(SW4_EXCEPTIONS)
  try {
#endif
    EW simulation(fileName, GlobalSources, GlobalTimeSeries);

    // nvtxRangePop();
    // PROFILER_STOP;
    if (!simulation.wasParsingSuccessful()) {
      if (myRank == 0) {
        cout << "Error: there were problems parsing the input file" << endl;
      }
      status = 1;
    } else {
      // get the simulation object ready for time-stepping
      // nvtxRangePushA("setup");
      simulation.setupRun(GlobalSources);
      // nvtxRangePop();

      if (!simulation.isInitialized()) {
        if (myRank == 0) {
          cout << "Error: simulation object not ready for time stepping"
               << endl;
        }
        status = 1;
      } else {
        if (myRank == 0) {
          int nth = 1;
#ifndef SW4_NOOMP
#pragma omp parallel
          {
            if (omp_get_thread_num() == 0) {
              nth = omp_get_num_threads();
            }
          }
#endif
          if (nth == 1) {
            if (nProcs > 1)
              cout << "Running sw4 on " << nProcs << " processors..." << endl;
            else
              cout << "Running sw4 on " << nProcs << " processor..." << endl;
          } else {
            if (nProcs > 1)
              // Assume same number of threads for each MPI-task.
              cout << "Running sw4 on " << nProcs << " processors, using "
                   << nth << " threads/processor..." << endl;
            else
              cout << "Running sw4 on " << nProcs << " processor, using " << nth
                   << " threads..." << endl;
          }
          if (simulation.m_croutines)
            cout << "   Using C routines." << endl;
          else
            cout << "   Using fortran routines." << endl;
          cout << "Writing output to directory: " << simulation.getPath()
               << endl;
        }
        // run the simulation
        simulation.solve(GlobalSources[0], GlobalTimeSeries[0], 0);

        // save all time series

        double myWriteTime = 0.0, allWriteTime;

        for (int ts = 0; ts < GlobalTimeSeries[0].size(); ts++) {
          GlobalTimeSeries[0][ts]->writeFile();
#ifdef USE_HDF5
          myWriteTime += GlobalTimeSeries[0][ts]->getWriteTime();
          if (ts == GlobalTimeSeries[0].size() - 1) {
            GlobalTimeSeries[0][ts]->closeHDF5File();

            MPI_Reduce(&myWriteTime, &allWriteTime, 1, MPI_DOUBLE, MPI_MAX, 0,
                       MPI_COMM_WORLD);
            if (myRank == 0)
              cout << "  ==> Max wallclock time to write time-series data is "
                   << allWriteTime << " seconds." << endl;
          }
#endif
        }

        if (myRank == 0) {
          cout << "============================================================"
               << endl
               << " program sw4 finished! " << endl
               << "============================================================"
               << endl;
        }

        status = 0;
      }
    }
#if defined(SW4_EXCEPTIONS)
  } catch (int e) {
    printf("Exception %d in EW init\n", e);
    MPI_Finalize();
    return status;
  }
#endif
  if (status == 1) {
    cout << "============================================================"
         << endl
         << "The execution on proc " << myRank << " was UNSUCCESSFUL." << endl
         << "============================================================"
         << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  print_hwm();

#ifdef USE_MAGMA
if (magma_finalize()!=MAGMA_SUCCESS){
  std::cerr<<"ERROR MAGMA FINALIZE FAILED\n";
  abort();
 } //else std::cout<<"MAGMA FINALIZE SUCCESSFULL\n"<<std::flush;
#endif
// Stop MPI
MPI_Finalize();
//std::cout<<"MPI_Finalize done\n"<<std::flush;
  return status;
}  // end of main
