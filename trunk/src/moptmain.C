#include "EW.h"
#include "DataPatches.h"
#include <cstring>
#include "version.h"

void usage(string thereason)
{
  cout << endl
       << "sbp4mopt - Summation by parts 4th order inverse seismic wave solver"  << endl << endl
       << "Usage: sbp4mopt [-v] file.in" << endl
       << "\t -v:      prints out the version info" << endl
       << "\t file.in: an input file" << endl << endl
       << "Reason for message: " << thereason << endl;
}

void compute_f_and_df( EW& simulation, double xs[11], int nmpar, double* xm,
		       vector<Source*>& GlobalSources,
		       vector<TimeSeries*>& GlobalTimeSeries,
		       vector<TimeSeries*>& GlobalObservations, int varcase,
		       double& f, double dfs[11], double* dfm, int myrank )

//-----------------------------------------------------------------------
// Compute misfit and its gradient.
//
// Input: simulation - Simulation object
//        xs         - Vector of source parameters
//        nmpar      - Number of material parameters
//        xm         - Vector of material parameters
//        GlobalSources - The single source object
//        GlobalTimeSeries   - TimeSeries objects, number of objects and
//                    locations should agree with the GlobalObservations vector.
//        GlobalObservations - The observed data at receivers.
//        varcase  - 0 means 11 free parameters, 1 means assume fixed frequency
//                   2 means assume frequency and t0 fixed.
//
// Output: GlobalTimeSeries - The solution of the forward problem at the stations.
//         f                - The misfit.
//         dfs              - Gradient wrt to the source of misfit.
//         dfm              - Gradient wrt to the material of misfit.
//-----------------------------------------------------------------------
{
   vector<Source*> src(1);
   src[0] = GlobalSources[0]->copy(" ");
   src[0]->set_parameters( xs );

// Translate one-dimensional parameter vector xm to material data (rho,mu,lambda)
   int ng = simulation.mNumberOfGrids;
   vector<Sarray> rho(ng), mu(ng), lambda(ng);
   simulation.parameters_to_material( nmpar, xm, rho, mu, lambda );

// Run forward problem with guessed source, upred_saved,ucorr_saved are allocated
// inside solve_allpars. U and Um are final time solutions, to be used as 'initial' data
// when reconstructing U backwards.
   vector<DataPatches*> upred_saved(ng), ucorr_saved(ng);
   vector<Sarray> U(ng), Um(ng);
   simulation.solve_allpars( src, rho, mu, lambda, GlobalTimeSeries, U, Um, upred_saved, ucorr_saved );

// Compute misfit, 'diffs' will hold the source for the adjoint problem
   vector<TimeSeries*> diffs;
   for( int m=0 ; m < GlobalTimeSeries.size() ; m++ )
   {
      TimeSeries *elem = GlobalTimeSeries[m]->copy( &simulation, "diffsrc" );
      diffs.push_back(elem);
   }
   f = 0;
   double dshift, ddshift;
   for( int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
      f += GlobalTimeSeries[m]->misfit( *GlobalObservations[m], diffs[m], dshift, ddshift );

   double mftmp = f;
   MPI_Allreduce(&mftmp,&f,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   if( myrank == 0 )
   {
      cout.precision(16);  
      cout << " Misfit is = " << f << endl;
   }
// Get gradient by solving the adjoint problem:
   simulation.solve_backward_allpars( src, rho, mu, lambda,  diffs, U, Um, upred_saved, ucorr_saved, dfs, nmpar, dfm );
   if( varcase == 1 )
      dfs[10] = 0;
   else if( varcase == 2 )
      dfs[10] = dfs[9] = 0;
   
// diffs no longer needed, give back memory
   for( unsigned int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
      delete diffs[m];
   diffs.clear();

// Give back memory
   for( unsigned int g=0 ; g < ng ; g++ )
   {
      delete upred_saved[g];
      delete ucorr_saved[g];
   }
   delete src[0];
}

//-----------------------------------------------------------------------
int start_minv( int argc, char **argv, string& input_file, int& myRank,
		int& nProcs )
{
  stringstream reason;

  // Initialize MPI...
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

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
     return 2;
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
    input_file = argv[1];

  if (myRank == 0) 
  {
    cout << ewversion::getVersionInfo() << endl;
    cout << "Input file: " << input_file << endl;
  }
  MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
  return 0;
}

//-----------------------------------------------------------------------
int main(int argc, char **argv)
{
  string fileName;
  int myRank, nProcs;
  int status = start_minv( argc, argv, fileName, myRank, nProcs );
  
  if( status == 0 )
  {
// Save the source description here
     vector<Source*> GlobalSources; 
// Save the time series here
     vector<TimeSeries*> GlobalTimeSeries;
     vector<TimeSeries*> GlobalObservations;

// make a new simulation object by reading the input file 'fileName'
     EW simulation(fileName, GlobalSources, GlobalObservations, true );

     if (!simulation.wasParsingSuccessful())
     {
	if (myRank == 0)
	   cout << "Error: there were problems parsing the input file" << endl;
	status = 1;
     }
     else
     {
// get the simulation object ready for time-stepping
	simulation.setupRun( );
	simulation.preprocessSources( GlobalSources );
	if (!simulation.isInitialized())
	{ 
	   if (myRank == 0)
	      cout << "Error: simulation object not ready for time stepping" << endl;
	   status=1;
	}
	else if( GlobalSources.size() != 1 )
	{
	   if (myRank == 0)
	      cout << "Source optmization only implemented for a single source" << endl;
	}
	else
	{
// Successful initialization

//	   simulation.setQuiet(true);
	   simulation.setQuiet(false);
// Make observations aware of the utc reference time, if set.
// Filter observed data if required
	   for( int m = 0; m < GlobalObservations.size(); m++ )
	   {
	      //	      simulation.set_utcref( *GlobalObservations[m] );
	      if( simulation.m_prefilter_sources && simulation.m_filter_observations )
	      {
		 GlobalObservations[m]->filter_data( simulation.m_filterobs_ptr );
		 GlobalObservations[m]->writeFile( "_fi" );
	      }
	   }

//  First copy observations to GlobalTimeSeries, later, solve will insert 
//  the simulation time step and start time into GlobalTimeSeries.
	   for( int m = 0; m < GlobalObservations.size(); m++ )
	   {
	      string newname = "_out";
	      TimeSeries *elem = GlobalObservations[m]->copy( &simulation, newname, true );
	      GlobalTimeSeries.push_back(elem);
	   }

	   if (myRank == 0)
	   {
	      cout << "Running sw4mopt on " <<  nProcs << " processors..." << endl
		   << "Writing output to directory: " 
		   << simulation.getOutputPath() << endl;
	   }
	
// Default initial guess, the input source, stored in GlobalSources[0]
	   double xs[11];
	   GlobalSources[0]->get_parameters(xs);

// Perturb material if gradient testing
           simulation.perturb_mtrl();
// Default initial guess material = the material given in the input file
           int nmpar;
	   double* xm=NULL;
	   simulation.get_nr_of_material_parameters( nmpar );
	   //           cout << "nmpar " << nmpar << endl;
           if( nmpar > 0 )           
	   {
	      xm = new double[nmpar];
	      simulation.get_material_parameter( nmpar, xm );
	   }

	   if( myRank == 0 )
	   {
	      cout << "Initial source guess : \n";
	      cout << "   x0 = " << xs[0] << " y0 = " << xs[1] << " z0 = " << xs[2] <<endl;
	      cout << "  mxx = " << xs[3] << " mxy= " << xs[4] << " mxz= " << xs[5] << endl;
	      cout << "  myy = " << xs[6] << " myz= " << xs[7] << " mzz= " << xs[8] << endl;
	      cout << "   t0 = " << xs[9] << " freq = " << xs[10] << endl;
	   }

// figure out how many parameters we need
	   int maxit, maxrestart, varcase=0, stepselection=0;
	   bool dolinesearch, fletcher_reeves=true;
	   double tolerance;
	   //	   simulation.get_cgparameters( maxit, maxrestart, tolerance, fletcher_reeves, stepselection,
	   //				     dolinesearch, varcase );
           varcase = 0;
	   int nvar=11;
	   if( varcase == 1 )
	      nvar = 10;
	   else if( varcase == 2 )
	      nvar = 9;

           double f, dfs[11];
	   double* dfm = NULL;
           if( nmpar > 0 )
	      dfm = new double[nmpar];

	   compute_f_and_df( simulation, xs, nmpar, xm, GlobalSources, GlobalTimeSeries,
			     GlobalObservations, varcase, f, dfs, dfm, myRank );

	   if( myRank == 0 )
	   {
	      cout << "============================================================" << endl
		   << " sw4opt ( Material/Source estimation solver) finished! " << endl
		   << "============================================================" << endl;
	   }
	}
     }
  }
  else if( status == 1 )
     cout  << "============================================================" << endl
	   << "The execution on proc " << myRank << " was unsuccessful." << endl
	   << "============================================================" << endl;
  if( status == 2 )
     status = 0;
// Stop MPI
  MPI_Finalize();
  return 0;
  // Note: Always return 0, to avoid having one error message per process from LC slurmd.
  //  return status;
} 

//-----------------------------------------------------------------------

   
