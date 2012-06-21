#include "mpi.h"

#include "EW.h"

#include <cstring>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <mpi.h>
#include "version.h"

using namespace std;

void usage(string thereason)
{
  cout << endl
       << "sw4 - Summation by parts 4th order forward seismic wave propagator"  << endl << endl
       << "Usage: sw4 [-v] file.in" << endl
       << "\t -v:      prints out the version info" << endl
       << "\t file.in: an input file" << endl << endl
       << "Reason for message: " << thereason << endl;
}

int
main(int argc, char **argv)
{
  int myRank = 0, nProcs = 0;
  string fileName;
  bool checkmode = false;

  stringstream reason;

  // Initialize MPI...
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

#ifdef ENABLE_TAU
   TAU_PROFILE_INIT(argc, argv);
#endif

  int status = 0;

  // mpi2 adds on four more args...  [-p4pg, dir, -p4wd, dir]
  int mpi2args = 4;

  if (argc != 2 && argc != 3 && argc != (2+mpi2args) && argc != (3+mpi2args) )
    {
      reason << "Wrong number of args (1-2), not: " << argc-1 << endl; 
      
      if (myRank == 0)
      {
	for (int i = 0; i < argc; ++i)
	  cout << "Argv[" << i << "] = " << argv[i] << endl;

	usage(reason.str());
      }
      
// Stop MPI
      MPI_Finalize();
      return 1;
    }
  else if (strcmp(argv[1],"-v") == 0 )
  {
     if (myRank == 0)
        cout << ewversion::getVersionInfo() << endl;
// Stop MPI
     MPI_Finalize();
     return status;
  }
  else
     if (argc == 1) 
     {
        reason  << "ERROR: ****No input file specified!" << endl;
        for (int i = 0; i < argc; ++i)
           reason << "Argv[" << i << "] = " << argv[i] << endl;
        
        if (myRank == 0) usage(reason.str());
// Stop MPI
	MPI_Finalize();
        return 1;
     }

  else
    fileName = argv[1];

  if (myRank == 0) 
  {
    cout << ewversion::getVersionInfo() << endl;
    cout << "Input file: " << fileName << endl;
  }
  
  MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

  cout.precision(6);
  
// Save the source description here
  vector<Source*> GlobalSources; 
// Save the time series here
  vector<TimeSeries*> GlobalTimeSeries;

// make a new simulation object by reading the input file 'fileName'
  EW simulation(fileName, GlobalSources, GlobalTimeSeries);

  if (!simulation.wasParsingSuccessful())
  {
    if (myRank == 0)
    {
      cout << "Error: there were problems parsing the input file" << endl;
    }
    status=1;
  }
  else
  {
// get the simulation object ready for time-stepping
    simulation.setupRun( );
    simulation.preprocessSources( GlobalSources );

    if (!simulation.isInitialized())
    { 
      if (myRank == 0)
      {
	cout << "Error: simulation object not ready for time stepping" << endl;
      }
      status=1;
    }
    else
    {
      if (myRank == 0)
      {
	cout << "Running sw4 on " <<  nProcs << " processors..." << endl
	     << "Writing output to directory: " 
	     << simulation.getOutputPath() << endl;
      }
// run the simulation
      simulation.solve( GlobalSources, GlobalTimeSeries );

// save all time series
      
      for (int ts=0; ts<GlobalTimeSeries.size(); ts++)
      {
	GlobalTimeSeries[ts]->writeFile();
      }

      if( myRank == 0 )
      {
	cout << "============================================================" << endl
	     << " program sw4 finished! " << endl
	     << "============================================================" << endl;
      }

      status = 0;
    }
  }
  
  if( status == 1 )
    cout  << "============================================================" << endl
	  << "The execution on proc " << myRank << " was unsuccessful." << endl
	  << "============================================================" << endl;

// Stop MPI
  MPI_Finalize();
  return status;
} // end of main
