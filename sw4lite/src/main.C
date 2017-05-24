#include <mpi.h>
#include <string>
#include <iostream>
using namespace std;
#include "EW.h"

int main( int argc, char** argv )
{
   //MPI_Init(&argc, &argv);
   int myRank;
   double  time_start, time_end;
#ifndef THREADED_MPI
   MPI_Init(&argc, &argv);
#else
   int provided;
   MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
   if (provided!=MPI_THREAD_MULTIPLE) std::cerr<<"ERROR:: MPI_THREAD_MULTIPLE is not supported\n";
   else std::cerr<<"Using Threaded MPI with MPI_THREAD_MULTIPLE\n";
#endif
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   MPI_Barrier(MPI_COMM_WORLD);
   time_start = MPI_Wtime();

   string filename;
   if( argc <= 1 )
   {
      if( myRank == 0 )
      {
	 cout  << "ERROR: ****No input file specified!" << endl;
	 for (int i = 0; i < argc; ++i)
	    cout << "Argv[" << i << "] = " << argv[i] << endl;
      }
      //MPI_Finalize();
      //return 1;
   }
   else
   {
      filename = argv[1];
      EW simulation(filename);
      //MPI_Finalize();
      //return 0;
   }

   MPI_Barrier(MPI_COMM_WORLD);
   time_end = MPI_Wtime();
   if(myRank == 0) cout <<  " Total running time: " << time_end - time_start << endl;

   MPI_Finalize();

   return 0;
}
