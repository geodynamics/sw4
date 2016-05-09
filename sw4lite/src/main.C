#include <mpi.h>
#include <string>
#include <iostream>
using namespace std;
#include "EW.h"

int main( int argc, char** argv )
{
   MPI_Init(&argc, &argv);
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   string filename;
   if( argc <= 1 )
   {
      if( myRank == 0 )
      {
	 cout  << "ERROR: ****No input file specified!" << endl;
	 for (int i = 0; i < argc; ++i)
	    cout << "Argv[" << i << "] = " << argv[i] << endl;
      }
      MPI_Finalize();
      return 1;
   }
   else
   {
      filename = argv[1];
      EW simulation(filename);
      MPI_Finalize();
      return 0;
   }
}
