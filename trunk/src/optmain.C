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
       << "sbp4opt - Summation by parts 4th order inverse seismic wave solver"  << endl << endl
       << "Usage: sbp4opt [-v] file.in" << endl
       << "\t -v:      prints out the version info" << endl
       << "\t file.in: an input file" << endl << endl
       << "Reason for message: " << thereason << endl;
}

void guess_source( EW & simulation, vector<Source*>& sources, vector<TimeSeries*>& timeseries,
		   vector<TimeSeries*>& observations, double* xv, int myRank );
extern "C" { void linsolvelu_( int*, double*, double*, double* ); }

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
  vector<TimeSeries*> GlobalObservations;

// make a new simulation object by reading the input file 'fileName'
  EW simulation(fileName, GlobalSources, GlobalObservations, true );

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

//  First copy observations to GlobalTimeSeries, and 
//  then setupRun will insert the simulation time step
//  and start time into GlobalTimeSeries.
     for( int m = 0; m < GlobalObservations.size(); m++ )
     {
	TimeSeries *elem = GlobalObservations[m]->copy( &simulation, "tscopy" );
	GlobalTimeSeries.push_back(elem);
     }

// get the simulation object ready for time-stepping
     simulation.setupRun( GlobalSources );

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
	cout << "Running sbp4opt on " <<  nProcs << " processors..." << endl
	     << "Writing output to directory: " 
	     << simulation.getOutputPath() << endl;
      }

// Source initial guess 
      if( GlobalSources.size() != 1 )
      {
	 cout << "Source optmization only implemented for a single source" << endl;
      }
      else
      {
         double xv[11];
	 guess_source( simulation, GlobalSources, GlobalTimeSeries, GlobalObservations, xv, myRank );
         if( myRank == 0 )
	 {
	    cout << "Initial source guess : \n";
	    cout << "   x0 = " << xv[0] << " y0 = " << xv[1] << " z0 = " << xv[2] <<endl;
	    cout << "  mxx = " << xv[3] << " mxy= " << xv[4] << " mxz= " << xv[5] << endl;
	    cout << "  myy = " << xv[6] << " myz= " << xv[7] << " mzz= " << xv[8] << endl;
	    cout << "   t0 = " << xv[9] << " freq = " << xv[10] << endl;
	 }
      }

// Run forward problem with guessed source
      simulation.setupRun( GlobalSources );
      simulation.solve( GlobalSources, GlobalTimeSeries );

// Compute misfit
      double mf = 0;
      for( int m=0 ; m < GlobalTimeSeries.size() ; m++ )
	 mf += GlobalTimeSeries[m]->misfit( *GlobalObservations[m] );
      double mftmp=mf;
      MPI_Allreduce(&mftmp,&mf,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      mf *= 0.5;
      if( myRank == 0 )
	 cout << "Misfit = " << mf << endl;

// For gradient, compute backwards problem:

// save all time series

      for (int ts=0; ts<GlobalTimeSeries.size(); ts++)
      {
	GlobalTimeSeries[ts]->writeFile();
      }

      if( myRank == 0 )
      {
	cout << "============================================================" << endl
	     << " sbp4f ( Summation by parts 4th order forward wave propagator) finished! " << endl
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

//-----------------------------------------------------------------------
void guess_source( EW &  simulation, vector<Source*>& sources, vector<TimeSeries*>& timeseries,
		   vector<TimeSeries*>& observations, double* xv, int myRank )
{
   // 1. Guess position and time offset
   int n = observations.size();
   double* dist = new double[n];
   double* xr   = new double[n];
   double* yr   = new double[n];
   double* zr   = new double[n];
   double cp, cs;
   simulation.average_speeds(cp,cs);
   cout << "average speeds " << cp << " " << cs << endl;

   for( int s= 0 ; s < observations.size() ; s++ )
      if( observations[s]->myPoint() )
      {
	 dist[s] = cp*observations[s]->arrival_time( 1e-6 );
	 xr[s] = observations[s]->getX();
	 yr[s] = observations[s]->getY();
	 zr[s] = observations[s]->getZ();
      }

   // Least squares with initial guess from 'source' object
   double d0 = cp*sources[0]->getOffset();
   double x0 = sources[0]->getX0();
   double y0 = sources[0]->getY0();
   double z0 = sources[0]->getZ0();
   double a[16], b[4], dx[4];
   int it=0, maxit=20;
   double err0=1, tol=1e-12, err=1;
   while( err > tol && it < maxit )
   {
      for( int i=0 ; i<16 ;i++ )
	 a[i] = 0;
      b[0]=b[1]=b[2]=b[3]=0;
      for( int s = 0 ; s < n ; s++ )
      {
         if( observations[s]->myPoint() )
	 {
	    a[0] += 4*(x0-xr[s])*(x0-xr[s]);
	    a[1] += 4*(x0-xr[s])*(y0-yr[s]);
	    a[2] += 4*(x0-xr[s])*(z0-zr[s]);
	    a[3] -= 4*(x0-xr[s])*(d0-dist[s]);

	    a[5] += 4*(y0-yr[s])*(y0-yr[s]);
	    a[6] += 4*(y0-yr[s])*(z0-zr[s]);
	    a[7] -= 4*(y0-yr[s])*(d0-dist[s]);

	    a[10] += 4*(z0-zr[s])*(z0-zr[s]);
	    a[11] -= 4*(z0-zr[s])*(d0-dist[s]);

	    a[15] += 4*(d0-dist[s])*(d0-dist[s]);

	    double nrm=(x0-xr[s])*(x0-xr[s])+(y0-yr[s])*(y0-yr[s])+
	       (z0-zr[s])*(z0-zr[s])-(d0-dist[s])*(d0-dist[s]);
	    b[0]  += 2*(x0-xr[s])*nrm;
	    b[1]  += 2*(y0-yr[s])*nrm;
	    b[2]  += 2*(z0-zr[s])*nrm;
	    b[3]  -= 2*(d0-dist[s])*nrm;
	 }
      }
 // Assemble matrix from all processors
      double ain[16], bin[4];
      for( int i = 0 ; i < 16 ; i++ )
	 ain[i] = a[i];
      for( int i = 0 ; i < 4 ; i++ )
	 bin[i] = b[i];

      MPI_Allreduce( ain, a, 16, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
      MPI_Allreduce( bin, b,  4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
      
      // symmetric matrix
      a[4]  = a[1];
      a[8]  = a[2];
      a[12] = a[3];
      a[9]  = a[6];
      a[13] = a[7];
      a[14] = a[11];
      int four=4;
      linsolvelu_( &four, a, b, dx );
      err = 0;
      for( int i=0 ; i < 4 ; i++ )
	 err = err > abs(dx[i]) ? err : abs(dx[i]);
      if( it == 0 )
	 err0 = err;
      err /= err0;
      //      cout <<  "Gauss-Newton iteration " << it << " " << err << endl;
      x0 -= dx[0];
      y0 -= dx[1];
      z0 -= dx[2];
      d0 -= dx[3];
      it++;
   }
   delete[] dist;
   delete[] xr;
   delete[] yr;
   delete[] zr;
   if( err > tol )
   {
      if( myRank == 0 )
	 cout << "Error: No convergence in Gauss-Newton iteration for initial guess\n";
   }
   double t0 = d0/cp;
   if( myRank == 0 )
   {
      cout << "Guessed position of source " << x0 << " " << y0 << " " << z0 << endl;
      cout << "Guessed time offset " << t0 << "..will not be used." << endl;
   }
   // 2. Guess t0 and freq. Take from given source for now.
   t0  = sources[0]->getOffset();
   double freq= sources[0]->getFrequency();
   //   x0 = sources[0]->getX0();
   //   y0 = sources[0]->getY0();
   //   z0 = sources[0]->getZ0();

   // 3. Tensor components guess, assume moment tensor source
   timeDep tfunc = sources[0]->getTfunc();
   double amp    = sources[0]->getAmplitude();
   double mxx, mxy, mxz, myy, myz, mzz;
   if( sources[0]->isMomentSource() )
   {
      vector<TimeSeries*> tsxx(n), tsxy(n), tsxz(n), tsyy(n), tsyz(n), tszz(n);
      for( int s=0 ; s < n ; s++ )
      {
         tsxx[s] = timeseries[s]->copy(&simulation, "tsxxcopy" );
         tsxy[s] = timeseries[s]->copy(&simulation, "tsxycopy" );
         tsxz[s] = timeseries[s]->copy(&simulation, "tsxzcopy" );
         tsyy[s] = timeseries[s]->copy(&simulation, "tsyycopy" );
         tsyz[s] = timeseries[s]->copy(&simulation, "tsyzcopy" );
         tszz[s] = timeseries[s]->copy(&simulation, "tszzcopy" );
      }

      Source* onesrc = new Source(&simulation, amp, freq, t0, x0, y0, z0, 1, 0, 0, 0, 0, 0, tfunc, "xx");
      vector<Source*> src(1);
      src[0] = onesrc;
      simulation.setupRun( src );
      simulation.solve( src, tsxx );
      delete onesrc;

      onesrc = new Source(&simulation, amp, freq, t0, x0, y0, z0, 0, 1, 0, 0, 0, 0, tfunc, "xy");
      src[0] = onesrc;
      //      simulation.setupRun( src );
      simulation.solve( src, tsxy );
      delete onesrc;

      onesrc = new Source(&simulation, amp, freq, t0, x0, y0, z0, 0, 0, 1, 0, 0, 0, tfunc, "xz");
      src[0] = onesrc;
      //      simulation.setupRun( src );
      simulation.solve( src, tsxz );
      delete onesrc;

      onesrc = new Source(&simulation, amp, freq, t0, x0, y0, z0, 0, 0, 0, 1, 0, 0, tfunc, "yy");
      src[0] = onesrc;
      //      simulation.setupRun( src );
      simulation.solve( src, tsyy );
      delete onesrc;

      onesrc = new Source(&simulation, amp, freq, t0, x0, y0, z0, 0, 0, 0, 0, 1, 0, tfunc, "yz");
      src[0] = onesrc;
      //      simulation.setupRun( src );
      simulation.solve( src, tsyz );
      delete onesrc;

      onesrc = new Source(&simulation, amp, freq, t0, x0, y0, z0, 0, 0, 0, 0, 0, 1, tfunc, "zz");
      src[0] = onesrc;
      //      simulation.setupRun( src );
      simulation.solve( src, tszz );
      delete onesrc;

      // Solve linear least squares problem for the moment components
      double a[36], b[6], x[6];
      for( int i=0 ; i < 36 ; i++ )
	 a[i] = 0;
      for( int i=0 ; i < 6 ; i++ )
	 b[i] = 0;

      for( int s=0 ; s < n ; s++ )
      {
	 TimeSeries* tsobs = tsxx[s]->copy(&simulation, "tsobscopy");
	 tsobs->interpolate( *observations[s] );

	 // assemble the matrix
	 if( tsxx[s]->myPoint() )
	 {
	    int nsteps = tsobs->getNsteps();
	    double** tsxxp = tsxx[s]->getRecordingArray();
	    double** tsxyp = tsxy[s]->getRecordingArray();
	    double** tsxzp = tsxz[s]->getRecordingArray();
	    double** tsyyp = tsyy[s]->getRecordingArray();
	    double** tsyzp = tsyz[s]->getRecordingArray();
	    double** tszzp = tszz[s]->getRecordingArray();
	    double** tsobsp = tsobs->getRecordingArray();
	    
#define amat(i,j) a[i-1+6*(j-1)]
	    for( int i=0 ; i < nsteps ; i++ )
	    {
	       amat(1,1) += (tsxxp[0][i]*tsxxp[0][i]+tsxxp[1][i]*tsxxp[1][i]+tsxxp[2][i]*tsxxp[2][i]);
	       amat(2,2) += (tsxyp[0][i]*tsxyp[0][i]+tsxyp[1][i]*tsxyp[1][i]+tsxyp[2][i]*tsxyp[2][i]);
	       amat(3,3) += (tsxzp[0][i]*tsxzp[0][i]+tsxzp[1][i]*tsxzp[1][i]+tsxzp[2][i]*tsxzp[2][i]);
	       amat(4,4) += (tsyyp[0][i]*tsyyp[0][i]+tsyyp[1][i]*tsyyp[1][i]+tsyyp[2][i]*tsyyp[2][i]);
	       amat(5,5) += (tsyzp[0][i]*tsyzp[0][i]+tsyzp[1][i]*tsyzp[1][i]+tsyzp[2][i]*tsyzp[2][i]);
	       amat(6,6) += (tszzp[0][i]*tszzp[0][i]+tszzp[1][i]*tszzp[1][i]+tszzp[2][i]*tszzp[2][i]);
	       amat(1,2) += (tsxxp[0][i]*tsxyp[0][i]+tsxxp[1][i]*tsxyp[1][i]+tsxxp[2][i]*tsxyp[2][i]);
	       amat(1,3) += (tsxxp[0][i]*tsxzp[0][i]+tsxxp[1][i]*tsxzp[1][i]+tsxxp[2][i]*tsxzp[2][i]);
	       amat(1,4) += (tsxxp[0][i]*tsyyp[0][i]+tsxxp[1][i]*tsyyp[1][i]+tsxxp[2][i]*tsyyp[2][i]);
	       amat(1,5) += (tsxxp[0][i]*tsyzp[0][i]+tsxxp[1][i]*tsyzp[1][i]+tsxxp[2][i]*tsyzp[2][i]);
	       amat(1,6) += (tsxxp[0][i]*tszzp[0][i]+tsxxp[1][i]*tszzp[1][i]+tsxxp[2][i]*tszzp[2][i]);
	       amat(2,3) += (tsxyp[0][i]*tsxzp[0][i]+tsxyp[1][i]*tsxzp[1][i]+tsxyp[2][i]*tsxzp[2][i]);
	       amat(2,4) += (tsxyp[0][i]*tsyyp[0][i]+tsxyp[1][i]*tsyyp[1][i]+tsxyp[2][i]*tsyyp[2][i]);
	       amat(2,5) += (tsxyp[0][i]*tsyzp[0][i]+tsxyp[1][i]*tsyzp[1][i]+tsxyp[2][i]*tsyzp[2][i]);
	       amat(2,6) += (tsxyp[0][i]*tszzp[0][i]+tsxyp[1][i]*tszzp[1][i]+tsxyp[2][i]*tszzp[2][i]);
	       amat(3,4) += (tsxzp[0][i]*tsyyp[0][i]+tsxzp[1][i]*tsyyp[1][i]+tsxzp[2][i]*tsyyp[2][i]);
	       amat(3,5) += (tsxzp[0][i]*tsyzp[0][i]+tsxzp[1][i]*tsyzp[1][i]+tsxzp[2][i]*tsyzp[2][i]);
	       amat(3,6) += (tsxzp[0][i]*tszzp[0][i]+tsxzp[1][i]*tszzp[1][i]+tsxzp[2][i]*tszzp[2][i]);
	       amat(4,5) += (tsyyp[0][i]*tsyzp[0][i]+tsyyp[1][i]*tsyzp[1][i]+tsyyp[2][i]*tsyzp[2][i]);
	       amat(4,6) += (tsyyp[0][i]*tszzp[0][i]+tsyyp[1][i]*tszzp[1][i]+tsyyp[2][i]*tszzp[2][i]);
	       amat(5,6) += (tsyzp[0][i]*tszzp[0][i]+tsyzp[1][i]*tszzp[1][i]+tsyzp[2][i]*tszzp[2][i]);
	       b[0] += tsxxp[0][i]*tsobsp[0][i]+tsxxp[1][i]*tsobsp[1][i]+tsxxp[2][i]*tsobsp[2][i];
	       b[1] += tsxyp[0][i]*tsobsp[0][i]+tsxyp[1][i]*tsobsp[1][i]+tsxyp[2][i]*tsobsp[2][i];
	       b[2] += tsxzp[0][i]*tsobsp[0][i]+tsxzp[1][i]*tsobsp[1][i]+tsxzp[2][i]*tsobsp[2][i];
	       b[3] += tsyyp[0][i]*tsobsp[0][i]+tsyyp[1][i]*tsobsp[1][i]+tsyyp[2][i]*tsobsp[2][i];
	       b[4] += tsyzp[0][i]*tsobsp[0][i]+tsyzp[1][i]*tsobsp[1][i]+tsyzp[2][i]*tsobsp[2][i];
	       b[5] += tszzp[0][i]*tsobsp[0][i]+tszzp[1][i]*tsobsp[1][i]+tszzp[2][i]*tsobsp[2][i];
	    }
	 }	    
         delete tsobs;
      }
      double ain[36], bin[6];
      for( int i = 0 ; i < 36 ; i++ )
	 ain[i] = a[i];
      for( int i = 0 ; i < 6 ; i++ )
	 bin[i] = b[i];
      MPI_Allreduce( ain, a, 36, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
      MPI_Allreduce( bin, b,  6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
      // symmetric matrix
      amat(2,1) = amat(1,2);
      amat(3,1) = amat(1,3);
      amat(4,1) = amat(1,4);
      amat(5,1) = amat(1,5);
      amat(6,1) = amat(1,6);
      amat(3,2) = amat(2,3);
      amat(4,2) = amat(2,4);
      amat(5,2) = amat(2,5);
      amat(6,2) = amat(2,6);
      amat(4,3) = amat(3,4);
      amat(5,3) = amat(3,5);
      amat(6,3) = amat(3,6);
      amat(5,4) = amat(4,5);
      amat(6,4) = amat(4,6);
      amat(6,5) = amat(5,6);
      int six=6;
      //      if( myRank == 0 )
      //      {
      //	 cout << "Tensor system matrix : " << endl;
      //	 for( int i=1 ; i <= 6 ; i++ )
      //	    cout << " " << amat(i,1)<<" " << amat(i,2) <<" " << amat(i,3) <<" " << amat(i,4) <<" " << amat(i,5) <<" " << amat(i,6) << endl;
      //	 cout << " ... and right hand side : " << endl;
      //	 for( int i=1 ; i <= 6 ; i++ )
      //	    cout << " " << b[i-1] << endl;
      //	 cout << " ... amplitude = " << amp << endl;
      //      }
      linsolvelu_( &six, a, b, x );
      mxx = x[0]/amp;
      mxy = x[1]/amp;
      mxz = x[2]/amp;
      myy = x[3]/amp;
      myz = x[4]/amp;
      mzz = x[5]/amp;
#undef amat      
      if( myRank == 0 )
      {
	 cout << "Moment tensor guess: mxx = " << mxx << " mxy= " << mxy << " mxz= " << mxz <<endl;
	 cout << "                     myy = " << myy << " myz= " << myz << " mzz= " << mzz <<endl;
      }
      for( int s=0 ; s < timeseries.size() ; s++ )
      {
         delete tsxx[s];
	 delete tsxy[s];
	 delete tsxz[s];
	 delete tsyy[s];
	 delete tsyz[s];
	 delete tszz[s];
      }
   }
   else
   {
      cout << "Initial guess for point sources not implemented \n";
   }

   // Construct guessed source object and return in sources vector, instead of initial source
   Source *sguess = new Source( &simulation, amp, freq, t0, x0, y0, z0,
			       mxx, mxy, mxz, myy, myz, mzz, tfunc, "srcguess");
   
   delete sources[0];
   sources[0] = sguess;
   // Store in array for future use in CG iteration
   xv[0] = x0;
   xv[1] = y0;
   xv[2] = z0;
   xv[3] = mxx;
   xv[4] = mxy;
   xv[5] = mxz;
   xv[6] = myy;
   xv[7] = myz;
   xv[8] = mzz;
   xv[9] = t0;
   xv[10] = freq;
}
