#include "EW.h"
#include "DataPatches.h"
#include <cstring>
#include "version.h"

#include <fcntl.h>
#include <unistd.h>

void lbfgs( EW& simulation, int ns, double xs[11], double sf[11], int nmpar,
	    double* xm, double* sfm,
	    vector<Source*>& GlobalSources,
	    vector<TimeSeries*>& GlobalTimeSeries,
	    vector<TimeSeries*> & GlobalObservations,
	    int m, int myRank );


void usage(string thereason)
{
  cout << endl
       << "sbp4mopt - Summation by parts 4th order inverse seismic wave solver"  << endl << endl
       << "Usage: sbp4mopt [-v] file.in" << endl
       << "\t -v:      prints out the version info" << endl
       << "\t file.in: an input file" << endl << endl
       << "Reason for message: " << thereason << endl;
}



//-----------------------------------------------------------------------
void compute_f( EW& simulation, int ns, double xs[11], int nm, double* xm,
		vector<Source*>& GlobalSources,
		vector<TimeSeries*>& GlobalTimeSeries,
		vector<TimeSeries*>& GlobalObservations,
		double& mf )
//-----------------------------------------------------------------------
// Compute misfit.
//
// Input: simulation - Simulation object
//        x          - Vector of source parameters
//        GlobalSources - The single source object
//        GlobalTimeSeries   - TimeSeries objects, number of objects and
//                    locations should agree with the GlobalObservations vector.
//        GlobalObservations - The observed data at receivers.
//
// Output: GlobalTimeSeries - The solution of the forward problem at the stations.
//         mf               - The misfit.
//-----------------------------------------------------------------------
{
   vector<Source*> src(1);
   src[0] = GlobalSources[0]->copy(" ");
   src[0]->set_parameters( xs );

// Translate one-dimensional parameter vector xm to material data (rho,mu,lambda)
   int ng = simulation.mNumberOfGrids;
   vector<Sarray> rho(ng), mu(ng), lambda(ng);
   simulation.parameters_to_material( nm, xm, rho, mu, lambda );

// Run forward problem with guessed source, upred_saved,ucorr_saved are allocated
// inside solve_allpars. U and Um are final time solutions, to be used as 'initial' data
// when reconstructing U backwards.
   vector<DataPatches*> upred_saved(ng), ucorr_saved(ng);
   vector<Sarray> U(ng), Um(ng);
   simulation.solve_allpars( src, rho, mu, lambda, GlobalTimeSeries, U, Um, upred_saved, ucorr_saved );

// Compute misfit,
   mf = 0;
   double dshift, ddshift, dd1shift;
   for( int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
      mf += GlobalTimeSeries[m]->misfit( *GlobalObservations[m], NULL, dshift, ddshift, dd1shift );
   double mftmp = mf;
   MPI_Allreduce(&mftmp,&mf,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

// Give back memory
   for( unsigned int g=0 ; g < ng ; g++ )
   {
      delete upred_saved[g];
      delete ucorr_saved[g];
   }
   delete src[0];
}

//-----------------------------------------------------------------------
void compute_f_and_df( EW& simulation, int ns, double xs[11], int nmpar, double* xm,
		       vector<Source*>& GlobalSources,
		       vector<TimeSeries*>& GlobalTimeSeries,
		       vector<TimeSeries*>& GlobalObservations, 
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
   double dshift, ddshift, dd1shift;
   for( int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
      f += GlobalTimeSeries[m]->misfit( *GlobalObservations[m], diffs[m], dshift, ddshift, dd1shift );

   double mftmp = f;
   MPI_Allreduce(&mftmp,&f,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   if( myrank == 0 )
   {
      cout.precision(16);  
      cout << " Misfit is = " << f << endl;
   }
// Get gradient by solving the adjoint problem:
   simulation.solve_backward_allpars( src, rho, mu, lambda,  diffs, U, Um, upred_saved, ucorr_saved, dfs, nmpar, dfm );
   for( int i=ns ; i < 11 ; i++ )
      dfs[i] = 0;
   
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
void restrict( int active[6], int wind[6], double* xm, double* xmi )
{
   int ni = (wind[1]-wind[0]+1);
   int nj = (wind[3]-wind[2]+1);
   int nia= (active[1]-active[0]+1);
   int nja= (active[3]-active[2]+1);

   for( int k=wind[4] ; k<=wind[5] ; k++ )
      for( int j=wind[2] ; j<=wind[3] ; j++ )
	 for( int i=wind[0] ; i<=wind[1] ; i++ )
	 {
            size_t indi = (i-wind[0]) + ni*(j-wind[2]) + ni*nj*(k-wind[4]);
	    size_t ind  = (i-active[0]) + nia*(j-active[2]) + nia*nja*(k-active[4]);
	    xmi[3*indi]   = xm[3*ind];
	    xmi[3*indi+1] = xm[3*ind+1];
	    xmi[3*indi+2] = xm[3*ind+2];
	 }
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
	simulation.setupRun( GlobalSources );
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
	   bool dolinesearch, fletcher_reeves=true, testing;
	   double tolerance;
	   simulation.get_cgparameters( maxit, maxrestart, tolerance, fletcher_reeves, stepselection,
					dolinesearch, varcase, testing );
           varcase = 4;
	   int nvar=11;
           if( varcase == 0 )
	      nvar = 11;
	   else if( varcase == 1 )
	      nvar = 10;
	   else if( varcase == 2 )
	      nvar = 9;
           else if( varcase == 4 )
	      nvar = 0;
	   else
	      if( myRank == 0 )
		 cout << "Error: varcase = " << varcase << " not defined for sw4mopt" << endl;
	   
           double f, dfs[11], sf[11];
	   double* dfm = NULL;
           double* sfm = NULL;

	   // Will not use source inversion, for the first tests, scale to 1.
           for( int i=0 ; i<11 ;i++)
	      sf[i]=1;
	   
           if( nmpar > 0 )
	   {
	      dfm = new double[nmpar];
              sfm = new double[nmpar];
	   }

	   // Try nondimensional units for first tests.
           for( int i=0 ; i<nmpar ;i++)
	      sfm[i]=1;

	   // New try with setting some scale factors
	   //	   simulation.get_scale_factors( nmpar, sfm );

           if( simulation.m_opttest == 2 )
	   {
	      compute_f_and_df( simulation, nvar, xs, nmpar, xm, GlobalSources, GlobalTimeSeries,
	   		        GlobalObservations, f, dfs, dfm, myRank );

	      double* dfmp= new double[nmpar]; 
              double h =1e-6;
              int grid = 0;
              int active[6], activeg[6];
	      activeg[0] = simulation.m_iStartActGlobal[grid];
	      activeg[1] = simulation.m_iEndActGlobal[grid];
	      activeg[2] = simulation.m_jStartActGlobal[grid];
	      activeg[3] = simulation.m_jEndActGlobal[grid];
	      activeg[4] = simulation.m_kStartActGlobal[grid];
	      activeg[5] = simulation.m_kEndActGlobal[grid];
	      active[0] = simulation.m_iStartAct[grid];
	      active[1] = simulation.m_iEndAct[grid];
	      active[2] = simulation.m_jStartAct[grid];
	      active[3] = simulation.m_jEndAct[grid];
	      active[4] = simulation.m_kStartAct[grid];
	      active[5] = simulation.m_kEndAct[grid];
              int wind[6];
	      // wind is intersection of 'active' and 'interior of proc.'
              for( int i=0 ; i < 6  ; i++ )
		 wind[i] = active[i];
	      if( active[1] >= active[0] )
	      {
		 if( wind[0] < simulation.m_iStartInt[grid] )
		    wind[0] = simulation.m_iStartInt[grid];
		 if( wind[1] > simulation.m_iEndInt[grid] )
		    wind[1] = simulation.m_iEndInt[grid];
		 if( wind[0] > wind[1] )
		    wind[1] = wind[0]-1;
	      }
	      if( active[3] >= active[2] )
	      {
		 if( wind[2] < simulation.m_jStartInt[grid] )
		    wind[2] = simulation.m_jStartInt[grid];
		 if( wind[3] > simulation.m_jEndInt[grid] )
		    wind[3] = simulation.m_jEndInt[grid];
		 if( wind[2] > wind[3] )
		    wind[3] = wind[2]-1;
	      }
	      if( active[5] >= active[4] )
	      {
		 if( wind[4] < simulation.m_kStartInt[grid] )
		    wind[4] = simulation.m_kStartInt[grid];
		 if( wind[5] > simulation.m_kEndInt[grid] )
		    wind[5] = simulation.m_kEndInt[grid];
		 if( wind[4] > wind[5] )
		    wind[5] = wind[4]-1;
	      }
	 
              int start[3]   ={wind[0]-activeg[0],wind[2]-activeg[2],wind[4]-activeg[4]};
	      int locsize[3] ={wind[1]-wind[0]+1,wind[3]-wind[2]+1,wind[5]-wind[4]+1};
	      int globsize[3]={activeg[1]-activeg[0]+1,activeg[3]-activeg[2]+1,activeg[5]-activeg[4]+1};
              
	      //              cout << myRank << " " << activeg[2] << " " << activeg[3] << " " <<  active[2] << " " << active[3] <<  " " << wind[2] << " " << wind[3] << endl;
     //              cout << myRank << " " << simulation .m_jStart[grid] << " " << simulation.m_jEnd[grid] << endl;
	      //              int start[3]   ={active[0]-activeg[0],active[2]-activeg[2],active[4]-activeg[4]};
	      //	      int locsize[3] ={active[1]-active[0]+1,active[3]-active[2]+1,active[5]-active[4]+1};
	      //	      int globsize[3]={activeg[1]-activeg[0]+1,activeg[3]-activeg[2]+1,activeg[5]-activeg[4]+1};
              int nptsbuf = 1000000;
              Parallel_IO* pio = new Parallel_IO(1,1,globsize,locsize,start,nptsbuf,0);
              int fid;
              if( myRank == 0 )
	      {
		 // create file, write header
		 fid = open( "hessian.bin", O_CREAT | O_TRUNC | O_WRONLY, 0660 ); 
		 if (fid == -1 )
		 {
		    VERIFY2(0, "ERROR: error opening file hessian.bin for writing header");
                    exit(-1);
		 }
                 int prec = 8;
		 size_t nr=write(fid,&prec,sizeof(int));
		 int npatch=1;
		 nr=write(fid,&npatch,sizeof(int));
                 int nc = 3;
		 nr=write(fid,&nc,sizeof(int));
		 double gz=simulation.mGridSize[0];
		 nr=write(fid,&gz,sizeof(double));
		 int dims[6]={activeg[0],globsize[0],activeg[2],globsize[1],activeg[4],globsize[2]};
		 nr=write(fid,dims,6*sizeof(int));
	      }
              else
		 fid = open( "hessian.bin", O_WRONLY );
              size_t offset = sizeof(double)+9*sizeof(int);
              size_t colsize = globsize[0]*((size_t)globsize[1])*globsize[2]*3;

	      // Verify IO
              int npts = (wind[1]-wind[0]+1)*(wind[3]-wind[2]+1)*(wind[5]-wind[4]+1);
              double* xmi;
              if( npts > 0 )
                 xmi = new double[3*npts];
	      //	          simulation.get_material_parameter( nmpar, xm );
	      //		  if( npts > 0 )
	      //		     restrict( active, wind, xm, xmi );
	      //			  pio->write_array( &fid, 3, xmi, offset, "double");
	      //			  close(fid);
	      //			  exit(0);
              for( int kper = activeg[4] ; kper <= activeg[5] ; kper++ )
		 for( int jper = activeg[2] ; jper <= activeg[3] ; jper++ )
		    for( int iper = activeg[0] ; iper <= activeg[1] ; iper++ )
                       for( int var = 0 ; var < 3 ; var++ )
		       {
		       // perturb material 
			  simulation.perturb_mtrl(iper,jper,kper,h,grid,var);
			  simulation.get_material_parameter( nmpar, xm );
			  compute_f_and_df( simulation, nvar, xs, nmpar, xm, GlobalSources, GlobalTimeSeries,
					    GlobalObservations, f, dfs, dfmp, myRank );
			  for( int p= 0 ; p < nmpar ; p++ )
			     dfmp[p] = (dfmp[p]-dfm[p])/h;
			  // Save Hessian column
			  restrict( active, wind, dfmp, xmi );
			  pio->write_array( &fid, 3, xmi, offset, "double");
			  offset += colsize*sizeof(double);
		       // restore material 
			  simulation.perturb_mtrl(iper,jper,kper,-h,grid,var);
		       }
	      close(fid);
	   }
           else
	   {
	      int method, bfgs_m;
	      simulation.get_optmethod( method, bfgs_m );
	   
	      lbfgs( simulation, nvar, xs, sf, nmpar, xm, sfm, GlobalSources, GlobalTimeSeries,
		     GlobalObservations, bfgs_m, myRank );
	   }

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

   
