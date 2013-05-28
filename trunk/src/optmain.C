#include "mpi.h"

#include "EW.h"
#include "ConvParOutput.h"

#include <fcntl.h>

#include <cstring>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <mpi.h>
#include "version.h"

#include "F77_FUNC.h"

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

void guess_source_position_layer( EW & simulation, vector<Source*>& sources, 
			    vector<TimeSeries*>& observations,
			    double& x0, double& y0, double& z0, double& t0, int myRank );

void guess_source_position( EW & simulation, vector<Source*>& sources, 
			    vector<TimeSeries*>& timeseries, vector<TimeSeries*>& observations,
			    double& x0, double& y0, double& z0, double& t0, int myRank );

void guess_source_t0freq( EW & simulation, vector<Source*>& sources,
			  double tstart, double& t0, double& freq, int myRank );

void guess_source_moments( EW & simulation, vector<Source*>& sources, vector<TimeSeries*>& timeseries,
			   vector<TimeSeries*>& observations, double* xv, int myRank );

void guess_shift( int n, double xs[n], vector<TimeSeries*>& timeseries,
		  vector<TimeSeries*>& observations, double dt );

extern "C" { void F77_FUNC(linsolvelu,LINSOLVEU)( int*, double*, double*, double* );
             void F77_FUNC(dgesv,DGESV)( int*, int*, double*, int*, int*, double*, int*, int*);
             void F77_FUNC(dgels,DGELS)( char*, int*, int*, int*, double*, int*, double*, int*,
					 double*, int*, int*); }


//-----------------------------------------------------------------------
void testsourced2( EW & simulation, vector<Source*>& GlobalSources )
// Test d^2(src)/dp_m dp_k implementation.
{
   double gradient[11], hess[121];
   simulation.testsourcediff( GlobalSources, gradient, hess );
   for( int i=0 ; i < 11 ; i++ )
      for( int j=0 ; j<i ; j++ )
	 hess[i+11*j] = hess[j+11*i];

   cout << "======================================================================="<< endl;
   vector<Source*> persrc(1);
   for( int testcomp = 0 ; testcomp < 11 ; testcomp++ )
   {
      double h = 0.00001;
      persrc[0] = GlobalSources[0]->copy( " " );
      persrc[0]->perturb( h, testcomp );
      double gradp[11], hessp[121];
      simulation.testsourcediff( persrc, gradp, hessp );
      cout << " hessian, colum, "  << testcomp << " numerically computed " << endl;
      for( int m=0 ; m < 11 ; m++ )
	 cout << " " << -hess[m+testcomp*11] << "   " << -(gradp[m]-gradient[m])/h << endl;
      delete persrc[0];
   }
}
//-----------------------------------------------------------------------
void test_gradient( EW& simulation, vector<Source*>& GlobalSources,
		    vector<TimeSeries*>& GlobalTimeSeries,
		    vector<TimeSeries*> & GlobalObservations, int myRank,
		    double sf[11] )
{
// Run forward problem with guessed source
   simulation.solve( GlobalSources, GlobalTimeSeries );
   //      for( int m=0 ; m < GlobalTimeSeries.size(); m++ )
   //         GlobalTimeSeries[m]->writeFile( "_1" );

// Compute misfit, 'diffs' will hold the source for the adjoint problem
   vector<TimeSeries*> diffs;

// Save all time series
//   for (int ts=0; ts<GlobalTimeSeries.size(); ts++)
//      GlobalTimeSeries[ts]->writeFile();

// Create and initialize 'diffs'
   char str[10];
   for( int m=0 ; m < GlobalTimeSeries.size() ; m++ )
   {
      snprintf(str,10,"%i",m);
      string name = "diffsrc";
      name.append(str);
      //      TimeSeries *elem = GlobalTimeSeries[m]->copy( &simulation, "diffsrc" );
      TimeSeries *elem = GlobalTimeSeries[m]->copy( &simulation, name );
      diffs.push_back(elem);
   }

// Assign misfit time series to 'diffs'
   double mf = 0;
   double* dshiftloc  = new double[GlobalObservations.size()];
   double* ddshiftloc = new double[GlobalObservations.size()];
   for( int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
   {
      double dd1shift;
      mf += GlobalTimeSeries[m]->misfit( *GlobalObservations[m], diffs[m], dshiftloc[m], ddshiftloc[m], dd1shift );
   }
   double mftmp = mf;
   MPI_Allreduce(&mftmp,&mf,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   if( myRank == 0 )
      cout << "Misfit = " << mf << endl;
   double *dshift = new double[GlobalObservations.size()];
   double *ddshift = new double[GlobalObservations.size()];
   MPI_Allreduce(dshiftloc,dshift,GlobalObservations.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   MPI_Allreduce(ddshiftloc,ddshift,GlobalObservations.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   delete[] dshiftloc;
   delete[] ddshiftloc;
   if( myRank == 0 )
   {
      cout << "Observation gradient = " << endl;
      for( int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
	 cout << " " << dshift[m] << endl;
   }

   // Test Gradient computation
   // Get gradient by computing the backwards problem:
   double gradient[11], hess[121];
   simulation.solve_backward( GlobalSources, diffs, gradient, hess );
   if( myRank == 0 )
   {
      cout << "Gradient, by adjoint equation = " << endl;
      for( int i = 0 ; i < 11 ; i++ )
	 cout << "   " << gradient[i] << endl;

      FILE *fd = fopen("gradient-adj.txt","w");
      for( int i= 0 ; i < 11 ; i++ )
      {
	 fprintf( fd, " %12.5g ", gradient[i] );
	 fprintf( fd, "\n" );
      }
      fclose(fd);
      fd = fopen("gradient-adj.bin","w");
      fwrite(gradient,sizeof(double),11,fd);
      fclose(fd);
   }

   // Validate shift derivatives by comparing numerical gradient
   for( int testcomp = 0 ; testcomp < GlobalTimeSeries.size() ; testcomp++ )
   {
      double h = 2.98e-8*sqrt(mf);
      GlobalObservations[testcomp]->add_shift( h );
      double mfp = 0, dshift, ddshift, dd1shift;
      for( int m=0 ; m < GlobalTimeSeries.size() ; m++ )
	 mfp += GlobalTimeSeries[m]->misfit( *GlobalObservations[m], NULL, dshift, ddshift, dd1shift );
      mftmp=mfp;
      MPI_Allreduce(&mftmp,&mfp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      if( myRank == 0 )
      {
	 cout << "Misfit of perturbed problem = " << mfp << endl;
	 cout << "Difference = " << mfp-mf << endl;
	 cout << "Observation " << testcomp+1 << " Numerical derivative = " << (mfp-mf)/h << endl;
      }
      GlobalObservations[testcomp]->add_shift( -h );
   }

   // Find parameter sizes, used for approximate gradient by difference quotients.
   double xv[11], gnum[11];
   GlobalSources[0]->get_parameters(xv);
   double lscale = sqrt((xv[0]*xv[0]+xv[1]*xv[1]+xv[2]*xv[2])/3.0);
   double mscale = sqrt((xv[3]*xv[3]+xv[4]*xv[4]+xv[5]*xv[5]+xv[6]*xv[6]+xv[7]*xv[7]+xv[8]*xv[8])/6.0);
   double tscale = abs(xv[9])+1;
   double fscale = abs(xv[10])+1;
   double sizes[11]={lscale,lscale,lscale,mscale,mscale,mscale,mscale,mscale,mscale,tscale,fscale};
   //   for( int c=0; c < 11 ; c++ )
   //      cout << "size" << c << "= " << sizes[c] << endl;

   for( int testcomp = 0 ; testcomp < 11 ; testcomp++ )
   {
// Validation, compute derivative by forward solve
      GlobalSources[0]->set_derivative(testcomp);
      simulation.solve( GlobalSources, GlobalTimeSeries );
      double dxidp=0;
      for( int m=0 ; m < GlobalTimeSeries.size(); m++ )
	 dxidp += GlobalTimeSeries[m]->product( *diffs[m] );
      double dxidptmp = dxidp;
      MPI_Allreduce(&dxidptmp,&dxidp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);      
      if( myRank == 0 )
	 cout << "Component " << testcomp+1 << " Forward solve derivative = " << dxidp << endl;
      GlobalSources[0]->set_noderivative();         

// Validation, compute numerical gradient:
//      double h = 1e-8*sizes[testcomp];
      double h = 2.98e-8*sqrt(mf)*sf[testcomp];
      //      cout << "h = " << h << endl;
//      double h = 0.00;
      vector<Source*> persrc(1);
      persrc[0] = GlobalSources[0]->copy( " " );
      //      cout << "unper " << *GlobalSources[0] << endl;
      persrc[0]->perturb( h, testcomp );
	 //		 GlobalSources[0]->perturb(h,testcomp);
      //		 simulation.preprocessSources( persrc );
      //      cout << " per " << *persrc[0] << endl;
      simulation.solve( persrc, GlobalTimeSeries );
      //      for( int m=0 ; m < GlobalTimeSeries.size(); m++ )
      //         GlobalTimeSeries[m]->writeFile( "_2" );

      double mfp = 0, dshift, ddshift, dd1shift;
      for( int m=0 ; m < GlobalTimeSeries.size() ; m++ )
	 mfp += GlobalTimeSeries[m]->misfit( *GlobalObservations[m], NULL, dshift, ddshift, dd1shift );
      mftmp=mfp;
      MPI_Allreduce(&mftmp,&mfp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      if( myRank == 0 )
      {
	 cout << "Misfit of perturbed problem = " << mfp << endl;
	 cout << "Difference = " << mfp-mf << endl;
	 cout << "Component " << testcomp+1 << " Numerical derivative = " << (mfp-mf)/h << endl;
         gnum[testcomp] = (mfp-mf)/h;
      }
      delete persrc[0];
   }

   if( myRank == 0 )
   {
      FILE *fd = fopen("gradient-num.txt","w");
      for( int i= 0 ; i < 11 ; i++ )
      {
	 fprintf( fd, " %12.5g ", gnum[i] );
	 fprintf( fd, "\n" );
      }
      fclose(fd);
      fd = fopen("gradient-num.bin","w");
      fwrite(gnum,sizeof(double),11,fd);
      fclose(fd);
   }
}

//-----------------------------------------------------------------------
void test_hessian(  EW& simulation, vector<Source*>& GlobalSources,
		    vector<TimeSeries*>& GlobalTimeSeries,
		    vector<TimeSeries*> & GlobalObservations, int myRank,
                    double sf[11] )
{
   // Test Hessian computation
   // Run forward problem with guessed source
   simulation.solve( GlobalSources, GlobalTimeSeries );
   // Compute misfit, 'diffs' will hold the source for the adjoint problem
   vector<TimeSeries*> diffs;
   for( int m=0 ; m < GlobalTimeSeries.size() ; m++ )
   {
      TimeSeries *elem = GlobalTimeSeries[m]->copy( &simulation, "diffsrc" );
      diffs.push_back(elem);
   }
   double mf = 0, dshift, ddshift, dd1shift;
   for( int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
      mf += GlobalTimeSeries[m]->misfit( *GlobalObservations[m], diffs[m], dshift, ddshift, dd1shift );
   double mftmp = mf;
   MPI_Allreduce(&mftmp,&mf,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   if( myRank == 0 )
      cout << "Misfit = " << mf << endl;

   //   for( int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
   //      diffs[m]->writeFile("_1");

// Get Hessian by solving the backwards problem:
   bool skipthis = false;
   double gradient[11], hess[121], hnum[121], gnum[11];
   simulation.solve_backward( GlobalSources, diffs, gradient, hess );

// Assemble the first part of hessian matrix
   vector<vector<TimeSeries*> > dudp(11);
   if( !skipthis )
   {
      for( int comp = 0 ; comp < 11 ; comp++ )
      {
         if( myRank == 0 )
	    cout << "Solving forward problem no. " << comp+1 << endl;
	 for( int m=0 ; m < GlobalTimeSeries.size() ; m++ )
	 {
	    TimeSeries *elem = GlobalTimeSeries[m]->copy( &simulation, "dudpsrc" );
	    dudp[comp].push_back(elem);
	 }
	 GlobalSources[0]->set_derivative(comp);
	 simulation.solve( GlobalSources, dudp[comp] );
      }
      double hess1[121];
      for( int m= 0 ; m < 121 ; m++ )
	 hess1[m] = 0;
      for( int m=0 ; m < GlobalTimeSeries.size() ; m++ )
      {
	 for( int comp = 0 ; comp < 11 ; comp++ )
	    for(int comp1 = comp; comp1 < 11 ; comp1++ )
	       hess1[comp+11*comp1] += dudp[comp][m]->product_wgh(*dudp[comp1][m]);
      }
      double hess1tmp[121];
      for( int m=0 ; m < 121 ; m++ )
	 hess1tmp[m] = hess1[m];
      MPI_Allreduce(hess1tmp,hess1,121,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

      for( int comp = 0 ; comp < 11 ; comp++ )
	 for(int comp1 = 0 ; comp1 < comp; comp1++ )
	    hess1[comp+11*comp1] = hess1[comp1+11*comp];


	      // Add to second part of Hessian
      for( int comp = 0 ; comp < 11 ; comp++ )
	 for(int comp1 = 0 ; comp1 < 11; comp1++ )
	    hess[comp+11*comp1] += hess1[comp+11*comp1];

      if( myRank == 0 )
      {
	 cout << "Hessian, by adjoint equation = " << endl;
	 for( int i = 0 ; i < 11 ; i++ )
	 {
	    for( int j= 0 ; j<11 ; j++ )
	       cout << "   " << hess[i+11*j] ;
	    cout << endl;
	 }

	 FILE *fd = fopen("hessian-adj.txt","w");
	 for( int i= 0 ; i < 11 ; i++ )
	 {
	    for( int j= 0 ; j < 11 ; j++ )
	       fprintf( fd, " %12.5g ", hess[i+11*j] );
	    fprintf( fd, "\n" );
	 }
	 fclose(fd);
	 fd = fopen("hessian-adj.bin","w");
	 fwrite(hess,sizeof(double),121,fd);
	 fclose(fd);

	 fd = fopen("gradient-adj.txt","w");
	 for( int i= 0 ; i < 11 ; i++ )
	 {
	    fprintf( fd, " %12.5g ", gradient[i] );
	    fprintf( fd, "\n" );
	 }
	 fclose(fd);
	 fd = fopen("gradient-adj.bin","w");
	 fwrite(gradient,sizeof(double),11,fd);
	 fclose(fd);
      }
   }
   // Find parameter sizes, used for approximate hessian by difference quotients.
   //   double lscale = sqrt((gradient[0]*gradient[0]+gradient[1]*gradient[1]+gradient[2]*gradient[2])/3.0);
   //   double mscale = sqrt((gradient[3]*gradient[3]+gradient[4]*gradient[4]+gradient[5]*gradient[5]+
   //			 gradient[6]*gradient[6]+gradient[7]*gradient[7]+gradient[8]*gradient[8])/6.0);
   //   double tscale = abs(gradient[9])+1;
   //   double fscale = abs(gradient[10])+1;

   double xv[11];
   GlobalSources[0]->get_parameters(xv);
   double lscale = sqrt((xv[0]*xv[0]+xv[1]*xv[1]+xv[2]*xv[2])/3.0);
   double mscale = sqrt((xv[3]*xv[3]+xv[4]*xv[4]+xv[5]*xv[5]+xv[6]*xv[6]+xv[7]*xv[7]+xv[8]*xv[8])/6.0);
   double tscale = abs(xv[9])+1;
   double fscale = abs(xv[10])+1;
   double sizes[11]={lscale,lscale,lscale,mscale,mscale,mscale,mscale,mscale,mscale,tscale,fscale};

   // Validate by numerical differentiation
   vector<Source*> persrc(1);
   GlobalSources[0]->set_noderivative();
   for( int testcomp = 0 ; testcomp < 11 ; testcomp++ )
   {
      //      double h = 1e-6*sizes[testcomp];
      double h = 2.98e-8*sqrt(mf)*sf[testcomp];
      //            double h=0;
      persrc[0] = GlobalSources[0]->copy( " " );
      persrc[0]->perturb( h, testcomp );

// Run forward problem with perturbed source
      simulation.solve( persrc, GlobalTimeSeries );

// Compute misfit, 'diffs' will hold the source for the adjoint problem
      diffs.clear();
      for( int m=0 ; m < GlobalTimeSeries.size() ; m++ )
      {
	 TimeSeries *elem = GlobalTimeSeries[m]->copy( &simulation, "diffsrc" );
	 diffs.push_back(elem);
      }
      double mfp = 0, dshift, ddshift, dd1shift;
      for( int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
	 mfp += GlobalTimeSeries[m]->misfit( *GlobalObservations[m], diffs[m], dshift, ddshift, dd1shift );
      double mftmp = mfp;
      MPI_Allreduce( &mftmp, &mfp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

      //      cout << "testcomp = " << testcomp << " misfit p " << mfp << " h= " << h << endl;
      //      for( int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
      //	 diffs[m]->writeFile("_2");

      double gradp[11];
      simulation.solve_backward( persrc, diffs, gradp, hess );
      if( myRank == 0 )
      {
	 cout << "Hessian, col. no " << testcomp+1 << " , by numerical derivative = " << endl;
	 for( int i = 0 ; i < 11 ; i++ )
	 {
	    cout << "   " << (gradp[i]-gradient[i])/h << endl;
            hnum[i+11*testcomp] = (gradp[i]-gradient[i])/h;
	 }
         gnum[testcomp] = (mfp-mf)/h;
      }
      delete persrc[0];
   }
   if( myRank == 0 )
   {
      FILE *fd = fopen("hessian-num.txt","w");
      for( int i= 0 ; i < 11 ; i++ )
      {
	 for( int j= 0 ; j < 11 ; j++ )
	    fprintf( fd, " %12.5g ", hnum[i+11*j] );
	 fprintf( fd, "\n" );
      }
      fclose(fd);
      fd = fopen("hessian-num.bin","w");
      fwrite(hnum,sizeof(double),121,fd);
      fclose(fd);

      fd = fopen("gradient-num.txt","w");
      for( int i= 0 ; i < 11 ; i++ )
      {
	 fprintf( fd, " %12.5g ", gnum[i] );
	 fprintf( fd, "\n" );
      }
      fclose(fd);
      fd = fopen("gradient-num.bin","w");
      fwrite(gnum,sizeof(double),11,fd);
      fclose(fd);
   }
}

//-----------------------------------------------------------------------
void compute_scalefactors(  EW& simulation, vector<Source*>& GlobalSources,
			    vector<TimeSeries*>& GlobalTimeSeries,
			    vector<TimeSeries*> & GlobalObservations, int myRank,
			    int varcase, int nvar, double* sf )
{
   // Compute Hessian for input source:
   bool write_hessian = true;
   if( !simulation.compute_sf() )
      simulation.get_scalefactors(sf);
   if( simulation.compute_sf() || varcase == 3 )
   {
   // Run forward problem 
      if(myRank == 0 )
	 cout << endl << "*** Computing scale factors..." << endl << 
	    "Solving forward problem for evaluating misfit" << endl;
      simulation.solve( GlobalSources, GlobalTimeSeries );
   // Compute misfit, 'diffs' will hold the source for the adjoint problem
      vector<TimeSeries*> diffs;
      for( int m=0 ; m < GlobalTimeSeries.size() ; m++ )
      {
	 TimeSeries *elem = GlobalTimeSeries[m]->copy( &simulation, "diffsrc" );
	 diffs.push_back(elem);
      }

      double mf = 0;
      double* dshiftloc  = new double[GlobalTimeSeries.size()];
      double* ddshiftloc = new double[GlobalTimeSeries.size()];
      double* dd1shiftloc = new double[GlobalTimeSeries.size()];
      for( int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
	 mf += GlobalTimeSeries[m]->misfit( *GlobalObservations[m], diffs[m], dshiftloc[m], ddshiftloc[m], dd1shiftloc[m] );

      if( varcase == 3 )
      {
	 double* dshift  = new double[GlobalTimeSeries.size()];
	 double* ddshift = new double[GlobalTimeSeries.size()];
	 double* dd1shift = new double[GlobalTimeSeries.size()];
	 MPI_Allreduce( dshiftloc, dshift, GlobalTimeSeries.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	 MPI_Allreduce( ddshiftloc, ddshift, GlobalTimeSeries.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	 MPI_Allreduce( dd1shiftloc, dd1shift, GlobalTimeSeries.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	 //         if( myRank == 0 )
	 //	    for( int m=0 ; m < GlobalTimeSeries.size() ; m++ )
	 //	    {
	 //	       cout << "scale factors " << m << " : dshift= " << dshift[m] << " ddshift= " << ddshift[m] << " dd1shift= " << dd1shift[m] << endl;
	 //	    }
	 for( int m=0 ; m < GlobalTimeSeries.size() ; m++ )
	    if( ddshift[m] > 0 )
	       sf[11+m] = 1/sqrt(ddshift[m]);
	    else
	       sf[11+m] = 1/sqrt(dd1shift[m]);
	 //         sf[11] = 0.0127969;
	 //	 sf[12] =  0.00759026;
	 //         sf[13] = 0.00510912;
	 //         sf[14] =0.0103381;
	 //	 sf[15] = 0.0122798;
	 delete[] dshift;
	 delete[] ddshift;
	 delete[] dd1shift;
      }
      delete[] dshiftloc;
      delete[] ddshiftloc;
      delete[] dd1shiftloc;

      if( simulation.compute_sf() )
      {
// Get second part of Hessian by solving the backwards problem:
	 double gradient[11], hess[121];
	 if( myRank == 0 )
	    cout << "Solving one backwards problem to calculate H2..." << endl;
	 simulation.solve_backward( GlobalSources, diffs, gradient, hess );

// Assemble the first part of hessian matrix by solving 11 forward problems
	 vector<vector<TimeSeries*> > dudp(11);
	 if( myRank == 0 )
	    cout << "**** Starting the calculation of H1..." << endl;
	 for( int comp = 0 ; comp < nvar ; comp++ )
	 {
	    if( myRank == 0 )
	       cout << "Solving forward problem no. " << comp+1 << endl;
	    for( int m=0 ; m < GlobalTimeSeries.size() ; m++ )
	    {
	       TimeSeries *elem = GlobalTimeSeries[m]->copy( &simulation, "dudpsrc" );
	       dudp[comp].push_back(elem);
	    }
	    GlobalSources[0]->set_derivative(comp);
	    simulation.solve( GlobalSources, dudp[comp] );
	 }
	 GlobalSources[0]->set_noderivative();
	 double hess1[121];
	 for( int m= 0 ; m < 121 ; m++ )
	    hess1[m] = 0;
	 for( int m=0 ; m < GlobalTimeSeries.size() ; m++ )
	 {
	    for( int comp = 0 ; comp < nvar ; comp++ )
	       for(int comp1 = comp; comp1 < nvar ; comp1++ )
		  hess1[comp+11*comp1] += dudp[comp][m]->product_wgh(*dudp[comp1][m]);
	 }
	 double hess1tmp[121];
	 for( int m=0 ; m < 121 ; m++ )
	    hess1tmp[m] = hess1[m];
	 MPI_Allreduce(hess1tmp,hess1,121,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

	 // Fill in symmetric part
	 for( int comp = 0 ; comp < nvar ; comp++ )
	    for(int comp1 = 0 ; comp1 < comp; comp1++ )
	       hess1[comp+11*comp1] = hess1[comp1+11*comp];

   // Add to second part of Hessian
	 for( int comp = 0 ; comp < nvar ; comp++ )
	    for(int comp1 = 0 ; comp1 < nvar; comp1++ )
	       hess[comp+11*comp1] += hess1[comp+11*comp1];


   // Check if diagonal elements are positive
	 bool diagpositive = true;
	 for( int comp = 0 ; comp < nvar ; comp++ )
	    if( hess[comp+11*comp] <= 0 )
	       diagpositive = false;

	 if( !diagpositive )
	 {
	    if( myRank == 0 )
	    {
	       cout << "*** The full Hessian has negative diagonal elements, " << endl;
	       cout << "*** scaling factors will use only part H1 of the Hessian" << endl;
	    }
	    for( int comp= 0 ; comp < nvar ; comp++ )
	       if( hess1[comp+11*comp] != 0 )
		  sf[comp]=1/sqrt(hess1[comp+11*comp]);
	       else
		  sf[comp]=1;
	 }
	 else
	 {
	    for( int comp= 0 ; comp < nvar ; comp++ )
	       sf[comp]=1/sqrt(hess[comp+11*comp]);
	 }

// define sf for the remaining elements of sf[11]
	 for( int comp = nvar; comp < 11; comp++ )
	    sf[comp]=1;

	 if( myRank == 0 && write_hessian )
	 {
	    const string hfiletxt = simulation.getOutputPath() + "hessian0.txt";
	    FILE *fd = fopen(hfiletxt.c_str(),"w");
      //      FILE *fd = fopen("hessian0.txt","w");
	    for( int i= 0 ; i < 11 ; i++ )
	    {
	       for( int j= 0 ; j < 11 ; j++ )
		  fprintf( fd, " %12.5g ", hess[i+11*j] );
	       fprintf( fd, "\n" );
	    }
	    fclose(fd);
	    const string hfilebin = simulation.getOutputPath() + "hessian0.bin";
	    fd = fopen(hfilebin.c_str(),"w");
	 //      fd = fopen("hessian0.bin","w");
	    fwrite(hess,sizeof(double),121,fd);
	    fclose(fd);
	 }

// dudp no longer needed, give back memory
	 for( unsigned int comp = 0 ; comp < nvar ; comp++ )
	 {
	    for( unsigned int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
	       delete dudp[comp][m];
	    dudp[comp].clear();
	 }
	 dudp.clear();
      }
// diffs no longer needed, give back memory
      for( unsigned int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
	 delete diffs[m];
      diffs.clear();
   }
}

//-----------------------------------------------------------------------
void get_testmatrix( double a[121], double b[11] )
{
   int matnr = 2;
   if( matnr == 1 )
   // Matrix, condition number 42.7
   {
a[0]= 3.7995;
a[1]= 1.1685;
a[2]= 1.2067;
a[3]= 0.70703;
a[4]= 1.3233;
a[5]= 0.82453;
a[6]= 1.0939;
a[7]= 1.6855;
a[8]= 1.0282;
a[9]= 1.3346;
a[10]= 1.2185;
a[11]= 1.1685;
a[12]= 3.8464;
a[13]= 0.28999;
a[14]= 0.60157;
a[15]= 1.2261;
a[16]= 0.55725;
a[17]= 1.4118;
a[18]= 1.136;
a[19]= 0.49832;
a[20]= 1.2101;
a[21]= 0.83279;
a[22]= 1.2067;
a[23]= 0.28999;
a[24]= 3.6891;
a[25]= 0.66593;
a[26]= 0.27825;
a[27]= 0.79807;
a[28]= 1.2069;
a[29]= 1.1285;
a[30]= 1.3184;
a[31]= 1.5865;
a[32]= 0.94359;
a[33]= 0.70703;
a[34]= 0.60157;
a[35]= 0.66593;
a[36]= 3.9984;
a[37]= 0.51509;
a[38]= 0.72423;
a[39]= 1.1635;
a[40]= 0.62092;
a[41]= 0.53648;
a[42]= 1.6473;
a[43]= 1.5052;
a[44]= 1.3233;
a[45]= 1.2261;
a[46]= 0.27825;
a[47]= 0.51509;
a[48]= 2.3784;
a[49]= 1.2237;
a[50]= 0.65149;
a[51]= 0.28021;
a[52]= 1.2922;
a[53]= 1.0603;
a[54]= 1.3863;
a[55]= 0.82453;
a[56]= 0.55725;
a[57]= 0.79807;
a[58]= 0.72423;
a[59]= 1.2237;
a[60]= 2.9778;
a[61]= 1.455;
a[62]= 1.0325;
a[63]= 1.1104;
a[64]= 0.83746;
a[65]= 0.4189;
a[66]= 1.0939;
a[67]= 1.4118;
a[68]= 1.2069;
a[69]= 1.1635;
a[70]= 0.65149;
a[71]= 1.455;
a[72]= 3.371;
a[73]= 1.0505;
a[74]= 0.68197;
a[75]= 0.72974;
a[76]= 1.4663;
a[77]= 1.6855;
a[78]= 1.136;
a[79]= 1.1285;
a[80]= 0.62092;
a[81]= 0.28021;
a[82]= 1.0325;
a[83]= 1.0505;
a[84]= 3.7767;
a[85]= 1.2186;
a[86]= 0.81181;
a[87]= 1.1149;
a[88]= 1.0282;
a[89]= 0.49832;
a[90]= 1.3184;
a[91]= 0.53648;
a[92]= 1.2922;
a[93]= 1.1104;
a[94]= 0.68197;
a[95]= 1.2186;
a[96]= 2.1314;
a[97]= 1.3262;
a[98]= 1.0394;
a[99]= 1.3346;
a[100]= 1.2101;
a[101]= 1.5865;
a[102]= 1.6473;
a[103]= 1.0603;
a[104]= 0.83746;
a[105]= 0.72974;
a[106]= 0.81181;
a[107]= 1.3262;
a[108]= 3.584;
a[109]= 0.61031;
a[110]= 1.2185;
a[111]= 0.83279;
a[112]= 0.94359;
a[113]= 1.5052;
a[114]= 1.3863;
a[115]= 0.4189;
a[116]= 1.4663;
a[117]= 1.1149;
a[118]= 1.0394;
a[119]= 0.61031;
a[120]= 3.8237;

// RHS
b[0]= 0.96878;
b[1]= 0.029845;
b[2]= 0.16637;
b[3]= 0.65361;
b[4]= 0.081468;
b[5]= 0.85381;
b[6]= 0.21996;
b[7]= 0.54435;
b[8]= 0.64101;
b[9]= 0.44907;
b[10]= 0.56633;

// solution
//x[0]= 0.47476;
//x[1]= 0.18175;
//x[2]= -0.3066;
//x[3]= 0.037958;
//x[4]= -1.066;
//x[5]= 0.6123;
//x[6]= -0.31587;
//x[7]= -0.36386;
//x[8]= 0.64496;
//x[9]= 0.02004;
//x[10]= 0.38605;
   }
   else if( matnr == 2 )
   {
      // LOH1 Hessian, condition number 4.73e38
a[0]= 0.0073804;
a[1]= 7.4506e-18;
a[2]= 8.2107e-17;
a[3]= -1.4517e-32;
a[4]= -1.0401e-31;
a[5]= -1.169e-34;
a[6]= 1.3703e-32;
a[7]= 2.1565e-18;
a[8]= 1.4943e-33;
a[9]= 6.2096e-13;
a[10]= -5.3488e-14;
a[11]= 7.4506e-18;
a[12]= 0.0073804;
a[13]= 7.6634e-17;
a[14]= 1.5086e-32;
a[15]= -1.0321e-31;
a[16]= 2.1565e-18;
a[17]= -1.4345e-32;
a[18]= 4.012e-33;
a[19]= -2.7255e-34;
a[20]= 6.3043e-13;
a[21]= -5.4247e-14;
a[22]= 8.2107e-17;
a[23]= 7.6634e-17;
a[24]= 0.0023531;
a[25]= 4.0842e-35;
a[26]= -1.6347e-18;
a[27]= -1.6785e-32;
a[28]= 2.1021e-35;
a[29]= -1.4197e-32;
a[30]= -2.7949e-34;
a[31]= 9.2974;
a[32]= -0.26103;
a[33]= -1.4517e-32;
a[34]= 1.5086e-32;
a[35]= 4.0842e-35;
a[36]= 1.3677e-33;
a[37]= -6.5483e-50;
a[38]= -4.2764e-49;
a[39]= -3.0916e-34;
a[40]= 1.3297e-49;
a[41]= -1.1475e-33;
a[42]= 6.2794e-31;
a[43]= 3.8317e-32;
a[44]= -1.0401e-31;
a[45]= -1.0321e-31;
a[46]= -1.6347e-18;
a[47]= -6.5483e-50;
a[48]= 3.6923e-33;
a[49]= 2.0848e-49;
a[50]= -3.2073e-50;
a[51]= -2.1382e-49;
a[52]= 6.4146e-50;
a[53]= -1.5005e-18;
a[54]= 7.1005e-16;
a[55]= -1.169e-34;
a[56]= 2.1565e-18;
a[57]= -1.6785e-32;
a[58]= -4.2764e-49;
a[59]= 2.0848e-49;
a[60]= 6.0437e-33;
a[61]= -1.0691e-49;
a[62]= 5.1317e-49;
a[63]= 2.5659e-49;
a[64]= -9.6205e-29;
a[65]= -3.3339e-30;
a[66]= 1.3703e-32;
a[67]= -1.4345e-32;
a[68]= 2.1021e-35;
a[69]= -3.0916e-34;
a[70]= -3.2073e-50;
a[71]= -1.0691e-49;
a[72]= 1.3677e-33;
a[73]= -9.5718e-50;
a[74]= -1.1475e-33;
a[75]= -9.0213e-31;
a[76]= -8.1078e-32;
a[77]= 2.1565e-18;
a[78]= 4.012e-33;
a[79]= -1.4197e-32;
a[80]= 1.3297e-49;
a[81]= -2.1382e-49;
a[82]= 5.1317e-49;
a[83]= -9.5718e-50;
a[84]= 6.0437e-33;
a[85]= -2.4589e-49;
a[86]= -1.0475e-28;
a[87]= -3.2624e-30;
a[88]= 1.4943e-33;
a[89]= -2.7255e-34;
a[90]= -2.7949e-34;
a[91]= -1.1475e-33;
a[92]= 6.4146e-50;
a[93]= 2.5659e-49;
a[94]= -1.1475e-33;
a[95]= -2.4589e-49;
a[96]= 2.9241e-33;
a[97]= 8.3791e-32;
a[98]= 5.0386e-32;
a[99]= 6.2096e-13;
a[100]= 6.3043e-13;
a[101]= 9.2974;
a[102]= 6.2794e-31;
a[103]= -1.5005e-18;
a[104]= -9.6205e-29;
a[105]= -9.0213e-31;
a[106]= -1.0475e-28;
a[107]= 8.3791e-32;
a[108]= 153360;
a[109]= -0.031911;
a[110]= -5.3488e-14;
a[111]= -5.4247e-14;
a[112]= -0.26103;
a[113]= 3.8317e-32;
a[114]= 7.1005e-16;
a[115]= -3.3339e-30;
a[116]= -8.1078e-32;
a[117]= -3.2624e-30;
a[118]= 5.0386e-32;
a[119]= -0.031911;
a[120]= 256.91;

// RHS
      b[0]= 0;
      b[1]= 0;
      b[2]= 0;
      b[3]= 0;
      b[4]= 0;
      b[5]= 0;
      b[6]= 0;
      b[7]= 0;
      b[8]= 0;
      b[9]= 0;
      b[10]= 0;
   }
}

//-----------------------------------------------------------------------
void compute_f( EW& simulation, int n, double x[n], vector<Source*>& GlobalSources,
		vector<TimeSeries*>& GlobalTimeSeries,
		vector<TimeSeries*>& GlobalObservations, double& mf, bool testing=false )
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
   //   bool testing = true;
   if( testing )
   {
      double a[121],b[11];
      get_testmatrix( a, b );
      double xax=0, bx=0;
      
      for( int j=0 ; j<11 ; j++ )
      {
	 //         cout << "adiag(" << j+1 << ") = " << a[j+11*j] << endl;
	 for( int i=0 ; i < 11 ; i++ )
	    xax += x[i]*a[i+11*j]*x[j];
         bx += x[j]*b[j];
      }
      mf = 0.5*xax-bx;
      return;
   }
   vector<Source*> src(1);
   src[0] = GlobalSources[0]->copy( " " );
   src[0]->set_parameters( x );

// Run forward problem with guessed source
   simulation.solve( src, GlobalTimeSeries );

// Compute misfit,
   mf = 0;
   double dshift, ddshift, dd1shift;
   for( int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
   {
      if( 11+m<n )
	 GlobalObservations[m]->add_shift(  x[11+m] );
      mf += GlobalTimeSeries[m]->misfit( *GlobalObservations[m], NULL, dshift, ddshift, dd1shift );
      if( 11+m<n )
	 GlobalObservations[m]->add_shift( -x[11+m] );
   }
   double mftmp = mf;
   MPI_Allreduce(&mftmp,&mf,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   delete src[0];
}

//-----------------------------------------------------------------------
void compute_f_and_df( EW& simulation, int n, double x[n],
		       vector<Source*>& GlobalSources,
		       vector<TimeSeries*>& GlobalTimeSeries,
		       vector<TimeSeries*>& GlobalObservations, int varcase,
		       double& f, double df[n], double ddf[121], bool testing=false )
//-----------------------------------------------------------------------
// Compute misfit and its gradient.
//
// Input: simulation - Simulation object
//        x          - Vector of source parameters
//        GlobalSources - The single source object
//        GlobalTimeSeries   - TimeSeries objects, number of objects and
//                    locations should agree with the GlobalObservations vector.
//        GlobalObservations - The observed data at receivers.
//        varcase  - 0 means 11 free parameters, 1 means assume fixed frequency
//                   2 means assume frequency and t0 fixed.
//
// Output: GlobalTimeSeries - The solution of the forward problem at the stations.
//         f                - The misfit.
//         df               - Gradient of misfit.
//         ddf              - The part of the Hessian that is associated 
//                            with the backward solve.
//-----------------------------------------------------------------------
{
   //   bool testing = true;
   if( testing )
   {
      double a[121], b[11];
      get_testmatrix( a, b );
      double xax=0, bx=0;
      for( int i=0 ; i < 11 ; i++ )
      {
	 df[i] = 0;
	 for( int j=0 ; j<11 ; j++ )
	 {
	    xax += x[i]*a[i+11*j]*x[j];
	    df[i] += a[i+11*j]*x[j];
	 }
	 df[i] -= b[i];
         bx += x[i]*b[i];
      }
      f = 0.5*xax-bx;
      //      for( int j=0 ; j<11 ; j++ )
      //	 for( int i=0 ; i < 11 ; i++ )
      //	    cout << "a(" << i+1 << "," << j+1 << ")="<< a[i+11*j] << ";"<<endl;
      return;
   }
   vector<Source*> src(1);
   src[0] = GlobalSources[0]->copy( " " );
   src[0]->set_parameters( x );

// Run forward problem with guessed source
   simulation.solve( src, GlobalTimeSeries );
// Compute misfit, 'diffs' will hold the source for the adjoint problem
   vector<TimeSeries*> diffs;
   for( int m=0 ; m < GlobalTimeSeries.size() ; m++ )
   {
      TimeSeries *elem = GlobalTimeSeries[m]->copy( &simulation, "diffsrc" );
      diffs.push_back(elem);
   }
   if( varcase == 3 )
   {
      f = 0;
      double* dshiftloc  = new double[GlobalObservations.size()];
      for( int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
	 dshiftloc[m] = 0;

      double ddshift, dd1shift;
      for( int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
      {
	 if( m+11 < n )
	    GlobalObservations[m]->add_shift( x[11+m] );
	 f += GlobalTimeSeries[m]->misfit( *GlobalObservations[m], diffs[m], dshiftloc[m], ddshift, dd1shift );
	 if( m+11 < n )
	    GlobalObservations[m]->add_shift( -x[11+m] );
      }
      MPI_Allreduce(dshiftloc,&df[11],GlobalObservations.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      delete[] dshiftloc;
   }
   else
   {
      f = 0;
      double dshift, ddshift, dd1shift;
      for( int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
	 f += GlobalTimeSeries[m]->misfit( *GlobalObservations[m], diffs[m], dshift, ddshift, dd1shift );
   }
   double mftmp = f;
   MPI_Allreduce(&mftmp,&f,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

   // Get gradient by solving the adjoint problem:
   simulation.solve_backward( src, diffs, df, ddf );

   if( varcase == 1 )
      df[10] = 0;
   else if( varcase == 2 || varcase == 3 )
      df[10] = df[9] = 0;

// diffs no longer needed, give back memory
   for( unsigned int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
      delete diffs[m];
   diffs.clear();
   delete src[0];
}

//-----------------------------------------------------------------------
void compute_dtd2fd( EW& simulation, int n, double x[n], vector<Source*>& GlobalSources,
		     vector<TimeSeries*>& GlobalTimeSeries,
		     vector<TimeSeries*>& GlobalObservations, double d[n],
		     int varcase, double& dtHd, bool ddf_isdefined,
		     double ddf[121], int myRank, bool testing=false )
//-----------------------------------------------------------------------
// Compute d'*H*d, where H is the Hessian.
//
// Input: simulation - Simulation object
//        x          - Vector of source parameters
//        GlobalSources - The single source object
//        GlobalTimeSeries   - TimeSeries objects, number of objects and
//                    locations should agree with the GlobalObservations vector.
//        GlobalObservations - The observed data at receivers.
//        d          - Vector to multiply the Hessian from left and right.
//        varcase  - 0 means 11 free parameters, 
//                   1 means assume fixed frequency
//                   2 means assume frequency and t0 fixed.
//                   3 means frequency and t0 fixed, + extra shift parameters
//        ddf_isdefined - 'true' means that the matrix ddf on input contains the part 
//                        of the Hessian that is associated with the backward solve.
//                        'false' means that ddf has not been computed. This routine
//                        will compute it, increasing the computational cost.
//
// Output: GlobalTimeSeries - The solution of the forward problem with a gradient source.
//         dtHd             - The product d'*H*d.
//
//-----------------------------------------------------------------------

{
   //   bool testing = true;
   if( testing )
   {
      double a[121],b[11];
      get_testmatrix( a, b );
      dtHd=0;
      for( int j=0 ; j < 11 ; j++ )
	 for( int i=0 ; i < 11 ; i++ )
	    dtHd += d[i]*a[i+11*j]*d[j];
      return;
   }

   vector<Source*> src(1);
   src[0] = GlobalSources[0]->copy( " ");
   src[0]->set_parameters( x );
   if( !ddf_isdefined )
   {
      simulation.solve( src, GlobalTimeSeries );

// Compute sources for the adjoint problem
      vector<TimeSeries*> diffs;
      for( int m=0 ; m < GlobalTimeSeries.size() ; m++ )
      {
	 TimeSeries *elem = GlobalTimeSeries[m]->copy( &simulation, "diffsrc" );
	 diffs.push_back(elem);
      }
      double dshift, ddshift, dd1shift;
      for( int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
	 GlobalTimeSeries[m]->misfit( *GlobalObservations[m], diffs[m], dshift, ddshift, dd1shift );

// Backward solve for second part of Hessian
      double df[11];
      src[0]->set_noderivative();
      simulation.solve_backward( src, diffs, df, ddf );

// diffs no longer needed, give back memory
      for( unsigned int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
	 delete diffs[m];
      diffs.clear();
   }

   double dtHdshift = 0;
   if( varcase == 3 )
   {
      // Here GlobalTimeSeries contains the current solution
      double* ddshiftloc  = new double[GlobalObservations.size()];
      for( int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
	 ddshiftloc[m] = 0;
      double dshift, dd1shift;
      for( int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
      {
	 if( m+11 < n )
	    GlobalObservations[m]->add_shift( x[11+m] );
	 GlobalTimeSeries[m]->misfit( *GlobalObservations[m], NULL, dshift, ddshiftloc[m], dd1shift );
	 if( m+11 < n )
	    GlobalObservations[m]->add_shift( -x[11+m] );
      }
      double *ddshift = new double[GlobalObservations.size()];
      MPI_Allreduce(ddshiftloc,ddshift,GlobalObservations.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      delete[] ddshiftloc;
      for( int m=0 ; m < GlobalTimeSeries.size() ; m++ )
	 dtHdshift += d[11+m]*ddshift[m]*d[11+m];
      delete[] ddshift;
   }

   // Forward solve for first part of Hessian
   src[0]->set_dirderivative( d );
   simulation.solve( src, GlobalTimeSeries );
   dtHd = 0;
   for( int m=0 ; m < GlobalTimeSeries.size() ; m++ )
      dtHd += GlobalTimeSeries[m]->product_wgh( *GlobalTimeSeries[m] );
   double dtHdtmp = dtHd;
   MPI_Allreduce(&dtHdtmp,&dtHd,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

   double d2h = 0;
   if( varcase == 1 )
      for( int i=0 ; i < 10 ; i++ )
	 for( int j=0 ; j < 10 ; j++ )
	    d2h += ddf[i+11*j]*d[i]*d[j];
   else if( varcase == 2 || varcase == 3 )
      for( int i=0 ; i < 9 ; i++ )
	 for( int j=0 ; j < 9 ; j++ )
	    d2h += ddf[i+11*j]*d[i]*d[j];
   else 
      for( int i=0 ; i < 11 ; i++ )
	 for( int j=0 ; j < 11 ; j++ )
	    d2h += ddf[i+11*j]*d[i]*d[j];

   //   if( myRank == 0 )
   //   {
   //      cout << "Hessian first part = " << dtHd << endl;
   //      cout << "Hessian second part= " << d2h << endl;
   //   }
   dtHd += d2h;
   dtHd += dtHdshift;

   delete src[0];
}

//-----------------------------------------------------------------------
void linesearch( EW& simulation, vector<Source*>& GlobalSources,
		 vector<TimeSeries*>& GlobalTimeSeries, vector<TimeSeries*>& GlobalObservations,
		 int nvar, int n, double x[n], double f, double df[n], double p[n],
		 double cgstep, double maxstep, double steptol, double xnew[n],
		 double& fnew, double sf[n], int myRank, int& retcode, int& nstep_reductions, bool testing )

//-----------------------------------------------------------------------
// Line seach by backtracking for CG.
// Algorithm A6.3.1 in book by Dennis and Schnabel, adapted for CG.
// 
// CG update without scaling is x^{k+1} = x^k + alpha_k p^k
//
// Scaled update: x^{k+1} = x^k + lambda*lambda_{max}*p^k/|| D p^k ||
//  where D is the diagonal matrix of scale factors.
// lambda is a normalized step length, where lambda = 1 corresponds 
// to the maximum allowed step. 
//
//  The maximum step is restricted by lambda_{max}. In unscaled units 
//  the max step becomes is the minimum of alpha_k and lambda_{max}/|| D p^k ||. 
//
//  The scaled update gives:
//   || Dx^{k+1} - Dx^k || = lambda*lambda_{max} 
//  Thus lambda_{max} is the maximum change of x in the scaled units (lambda=1).
//
//  The minimum step length is set to approximately lambda = steptol/lambda_{max}
//  so that 
//   || Dx^{k+1} - Dx^k || = steptol
//  in maximum norm, i.e., steptol is the smallest allowable change in scaled units.
//
// Input: simulation - Simulation object
//        GlobalSources - The single source object
//        GlobalTimeSeries   - TimeSeries objects, number of objects and
//                    locations should agree with the GlobalObservations vector.
//        GlobalObservations - The observed data at receivers.
//        x          - Vector of source parameters, current iterate
//        f          - Misfit at x.
//        df         - Gradient at x.
//        p          - Conjugate gradient vector = search direction
//        cgstep     - Original estimate of step length, alpha_k above.
//        maxstep    - Maximum allowed step lambda_{max} above. 
//        steptol    - Minimum allowed step, for the scaled variables.
//        sf         - Scale factors.
//        myRank     - ID of this processor.
//
// Output: xnew    - New iterate.
//         fnew    - Misfit at new iterate.
//         retcode - Return code, indicating 0 - success, 1 - failure
//
//-----------------------------------------------------------------------
{
   bool maxtaken = false;
   double alpha = 1e-6, cglen=0, ang=0;
   vector<double> lambdasave, fcnsave;

   for( int i=0 ; i < n ; i++ )
      xnew[i] = x[i];

   retcode = 2;
   for( int i=0 ; i < n ; i++ )
   {
      cglen += p[i]*p[i]/(sf[i]*sf[i]);
      ang   += p[i]*df[i];
   }
   if( ang >= 0 )
   {
      if( myRank == 0 )
      {
	 cout << "LINESEARCH: Warning, direction is not a descent direction" << endl;
         cout << "   switching to gradient direction search" << endl;
      }
      cglen = 0;
      for( int i=0 ; i < n ; i++ )
      {
	 p[i] = -df[i]/(sf[i]*sf[i]);
	 cglen += p[i]*p[i]/(sf[i]*sf[i]);
      }
   }
   cglen = sqrt(cglen);
   // p is scaled such that lambda=1 corresponds to the 
   // cg step x^{k+1} = x^k + cgstep*p^k unless maxstep is set lower.
   // i.e., maxstep is the maximum allowable relative change in x, ||x^{k+1}-x^k||/typx < maxstep
      if( cglen*cgstep > maxstep )
      {
         for( int i=0 ; i < n ; i++ )
            p[i] = p[i]*(maxstep/cglen);
	 cglen = maxstep;
      }
      else
      {
         for( int i=0; i < n ; i++ )
            p[i] = cgstep*p[i];
         cglen = cgstep*cglen;
      }
      // minlambda is scaled such that x^{k+1}-x^k = lambda*p^k with lambda=minlambda implies
      //            || x^{k+1}-x^k ||_\infty < minlambda (termination criterium)
      double initslope = 0;
      for( int i=0; i < n ; i++ )
         initslope += df[i]*p[i];

      double rellength = 0;

      for( int i=0; i < nvar ; i++ )
      {
         double rlocal;
         if( fabs(x[i]) > sf[i] )
	    rlocal = fabs(p[i])/fabs(x[i]);
	 else
	    rlocal = fabs(p[i])/sf[i];
	 if( rlocal > rellength )
	    rellength = rlocal;
      }
      for( int i=11; i < n ; i++ )
      {
         double rlocal;
         if( fabs(x[i]) > sf[i] )
	    rlocal = fabs(p[i])/fabs(x[i]);
	 else
	    rlocal = fabs(p[i])/sf[i];
	 if( rlocal > rellength )
	    rellength = rlocal;
      }
      double minlambda = steptol/rellength;
      double lambda = 1, fnewprev, lambdaprev;
      nstep_reductions = 0;
      while( retcode == 2 )
      {
	 for( int i=0; i < n ; i++ )
	    xnew[i] = x[i] + lambda*p[i];

	 // prevent z from becoming negative
	 //         if( !testing && (xnew[2] < 0 && p[2] != 0) )
	 //	 {
	 //            lambda = -x[2]/p[2];
	 //	    for( int i=0; i < n ; i++ )
	 //	       xnew[i] = x[i] + lambda*p[i];
	 //	 }	 
	 compute_f( simulation, n, xnew, GlobalSources, GlobalTimeSeries, GlobalObservations, fnew, testing );
         lambdasave.push_back(lambda);
         fcnsave.push_back(fnew);
         if( fnew < f + alpha*lambda*initslope )
	 {
	    // Found satisfactory step
	    retcode = 0;
	    if( lambda == 1 && (cglen > 0.99*maxstep) )
	       maxtaken = true;
	    return;
	 }
	 else if( lambda<minlambda )
	 {
	    // Could not find satisfactory step
	    ang = 0;
	    for( int i=0 ; i < n ; i++ )
	       ang += p[i]*df[i];

            if( myRank == 0 )
	    {
	       cout << "LINESEARCH: no satsifactory step found \n";
	       cout << "cg-alpha = " << cgstep << endl;
	       cout << "scprod = " << ang << endl;
	       cout << "search direction = " << endl;
	       for( int i=0 ; i < n ; i++ )
		  cout << " " << p[i] << endl;
	       cout << "gradient direction = " << endl;
	       for( int i=0 ; i < n ; i++ )
		  cout << " " << df[i] << endl;
	       cout << "maxstep   = " << maxstep << endl;
	       cout << "minlambda = " << minlambda << endl;
	       cout << "initslope = " << initslope << endl;
	       cout << " f = " << f << " fnew = " << fnew <<  " lambda = " << lambda ;
	       cout << " fnew-(f+alpha*lambda*initslope) = " << fnew-(f+alpha*lambda*initslope) << endl;
               cout << "lambda and f(lambda) tried = " << endl;
	       for( int i=0 ; i < lambdasave.size() ; i++ )
                  cout << lambdasave[i] << " " << fcnsave[i] << endl;
	       cout << "x = " << endl;
	       for( int i=0 ; i < n ; i++ )
		  cout << " " << x[i] << endl;
	    }
	    retcode = 1;
	    // Take a full step
	    for( int i=0 ; i < n ; i++ )
	       xnew[i] = x[i] + p[i];

	    // Compute return value for fnew
	    //            if( !testing && (xnew[2] < 0 && p[2] != 0) )
	    //	    {
	    //               lambda = -x[2]/p[2];
	    //	       for( int i=0 ; i < n ; i++ )
	    //		  xnew[i] = x[i] + lambda*p[i];
	    //	    }
	    compute_f( simulation, n, xnew, GlobalSources, GlobalTimeSeries, GlobalObservations, fnew, testing );
	    return;
	 }
	 else
	 {
	    // Reduce step size
            double ltemp;
	    if( lambda == 1 )
	       ltemp = -initslope/(2*(fnew-f-initslope));
	    else
	    {
               double r1 =(fnew-f-lambda*initslope)/(lambda-lambdaprev);
               double r2 =(fnewprev-f-lambdaprev*initslope)/(lambda-lambdaprev);
               double a = 1/(lambda*lambda)*r1 - 1/(lambdaprev*lambdaprev)*r2;
               double b = -lambdaprev/(lambda*lambda)*r1 +
		                       lambda/(lambdaprev*lambdaprev)*r2;
               double disc = b*b-3*a*initslope;
               if( a == 0 )
                  ltemp = -initslope/(2*b);
               else
                  ltemp = (-b+sqrt(disc))/(3*a);
               if( ltemp> 0.5*lambda )
                  ltemp = 0.5*lambda;
	    }
	    lambdaprev = lambda;
	    fnewprev = fnew;
	    if( ltemp < 0.1*lambda )
	       lambda = 0.1*lambda;
	    else
	       lambda = ltemp;
	    nstep_reductions++;
	 }
      }
}

//-----------------------------------------------------------------------
void cg( EW& simulation, int n, double x[n], double sf[n],
	 vector<Source*>& GlobalSources,
	 vector<TimeSeries*>& GlobalTimeSeries,
	 vector<TimeSeries*> & GlobalObservations,
	 int myRank )
//-----------------------------------------------------------------------
// Conjugated gradient minimization of misfit function.
// 
// Input: simulation - simulation object
//        x          - Parameter vector, used as initial guess.
//        sf         - Scaling factors.
//        GlobalSource - Object representing the unique source.
//        GlobalTimeSeries   - TimeSeries objects, number of objects and
//                    locations should agree with the GlobalObservations vector.
//        GlobalObservations - The observed data at receivers.
//        myRank      - ID of this processor.
//
// Output: x - The minimum is returned.
//
// Note: Configuration of the c-g algorithm is done by a call 
//       to get_cgparameters of the simulation object.
//
//-----------------------------------------------------------------------
{
   int nvar, j, k;
   bool done = false;
   bool dolinesearch = true;
   int maxit, maxrestart, varcase=0, stepselection=0;
   bool fletcher_reeves=true;

   // Do not allow source to rise closer than this to the surface:
   double tolerance;
   double f, d2f[121], rnorm;
   double* df, *d, *da, *xa, *dfp;
   bool testing=false, stepfileio=true;
   int verbose=3, ind=0, nreductions=0;
   int* iconvdata;
   double* convdata, *paradata;


   simulation.get_cgparameters( maxit, maxrestart, tolerance, fletcher_reeves, stepselection,
				dolinesearch, varcase, testing );
   if( maxrestart == 0 )
      return;

   df  = new double[n];
   d   = new double[n];
   da  = new double[n];
   xa  = new double[n];
   dfp = new double[n];

   if( varcase == 0 )
   {
      nvar  = 11;
      maxit = 11;
   }
   if( varcase == 1 )
   {
      nvar = 10;
      maxit= 10;
   }
   else if( varcase == 2 )
   {
      nvar = 9;
      maxit= 9;
   }
   else if( varcase == 3 )
   {
      nvar  = 9;
      maxit = 9 + GlobalObservations.size();
   }

   if( testing )
   {
      nvar  = 11;
      maxit = 11;
   }

   ConvParOutput outinfo( simulation, n, nvar, myRank, verbose, true );

   if( myRank == 0 && verbose > 0 )
      cout << "Begin CG iteration by evaluating initial misfit and gradient..." << endl;

  // Compute gradient and misfit
   compute_f_and_df( simulation, n, x, GlobalSources, GlobalTimeSeries, GlobalObservations,
		     varcase, f, df, d2f, testing );

 // Compute norm of gradient
   rnorm = 0;
   for( int i=0 ; i < n ; i++ )
      rnorm = rnorm > fabs(df[i])*sf[i] ? rnorm : fabs(df[i])*sf[i];

   outinfo.print_dfmsg(f,df,sf);
   outinfo.save_step(f,df,sf,x,rnorm,0.0,0,0,0);

   j = 1;
   while( j <= maxrestart && !done )
   {
      for( int i=0 ; i < n ; i++ )
         d[i] = -df[i]*sf[i]*sf[i];

      // Compute norm of gradient
      //      rnorm = 0;
      //      for( int i=0 ; i < n ; i++ )
      //	 rnorm = rnorm > fabs(df[i])*sf[i] ? rnorm : fabs(df[i])*sf[i];

      int k = 1;
      while( k <= maxit && rnorm > tolerance )
      {  
	 double alpha, fp;
	 if( stepselection == 0 )
	 {
	    double den = 0, dnrm=0, xnrm=0;
	    for( int i= 0 ; i < nvar ; i++ )
	    {
	       den  += d[i]*df[i];
	       dnrm += d[i]*d[i];
	       xnrm += x[i]*x[i];
	    }
	    for( int i= 11 ; i < n ; i++ )
	    {
	       den  += d[i]*df[i];
	       dnrm += d[i]*d[i];
	       xnrm += x[i]*x[i];
	    }
            double stlim = 0.1*sqrt(xnrm/dnrm);
            double h=-2*f/den;
            if( h > stlim )
	       h = stlim;
	    for( int i=0; i<n ; i++ )
	       xa[i] = x[i] + h*d[i];
            if( myRank == 0 && verbose > 2 )
	       cout << "Step length computation a) " << endl;
            compute_f( simulation, n, xa, GlobalSources, GlobalTimeSeries, GlobalObservations, fp, testing );
	    alpha = h*f/(fp+f);
	 }
	 else
	 {
            double dtHd;
            if( myRank == 0 && verbose > 2 )
	       cout << "Step length computation b) " << endl;

	    compute_dtd2fd( simulation, n, x, GlobalSources, GlobalTimeSeries, GlobalObservations, d, varcase, dtHd,
			    true, d2f, myRank, testing );
            outinfo.print_scalar( dtHd, "dtHd= " ,3 );
            outinfo.print_vector( d, "d= " ,3 );

            alpha = 0;
            for( int i=0; i<n ; i++ )
	       alpha += df[i]*d[i];
	    alpha = -alpha/dtHd;
	 }
         if( dolinesearch )
         {
	    for( int i=0 ; i < n ; i++ )
	       da[i] = d[i];
            int retcode;
            if( myRank == 0 && verbose > 2 )
	       cout << "Line search.. " << endl;
	    linesearch( simulation, GlobalSources, GlobalTimeSeries, GlobalObservations, nvar, n,
			x, f, df, da, fabs(alpha), 10.0, tolerance*0.01, xa, fp, sf, myRank, retcode,
			nreductions, testing );
            if( myRank == 0 && verbose > 2 )
	       cout << " .. return code "  << retcode << " misfit changed from " << f << " to " << fp << endl;
	 }
         else
	 {
	    for( int i=0 ; i < n ; i++ )
	       xa[i] = x[i] + alpha*d[i];
	 }         

	 // xa is now the new iterate
	 double dxnorm = 0;
	 for( int i=0 ; i < nvar ; i++ )
	 {
            double locnorm = fabs(x[i]-xa[i]);
	    if( fabs(xa[i])> sf[i] )
	       locnorm /= fabs(xa[i]);
	    else
	       locnorm /= sf[i];
            if( locnorm > dxnorm )
	       dxnorm = locnorm;
	    x[i] = xa[i];
	 }
	 for( int i=11 ; i < n ; i++ )
	 {
            double locnorm = fabs(x[i]-xa[i]);
	    if( fabs(xa[i])> sf[i] )
	       locnorm /= fabs(xa[i]);
	    else
	       locnorm /= sf[i];
            if( locnorm > dxnorm )
	       dxnorm = locnorm;
	    x[i] = xa[i];
	 }
	 compute_f_and_df( simulation, n, x, GlobalSources, GlobalTimeSeries, GlobalObservations, varcase, f, dfp,
			   d2f, testing );

// Compute norm of new gradient
         rnorm = 0;
	 for( int i=0 ; i < n ; i++ )
	    if( fabs(dfp[i])*sf[i] > rnorm )
	       rnorm = fabs(dfp[i])*sf[i];

// Save the time series after each sub-iteration.
         for( int ts=0 ; ts < GlobalTimeSeries.size() ; ts++ )
	    GlobalTimeSeries[ts]->writeFile();

// Print some information on standard output
         outinfo.print_xmsg(x,rnorm,dxnorm,j,k);
         outinfo.print_dfmsg(f,dfp,sf);

// Save convergence and parameters to file
         outinfo.save_step(f,dfp,sf,x,rnorm,dxnorm,nreductions,j,k);

         double beta;
         if( k < maxit )
	 {
            double num=0, mix=0, den=0;
	    for( int i=0 ; i < n ; i++ )
	    {
               num +=dfp[i]*dfp[i]*sf[i]*sf[i];
	       den +=  df[i]*df[i]*sf[i]*sf[i];
	       mix += df[i]*dfp[i]*sf[i]*sf[i];
	    }
	    if( fletcher_reeves )
	       beta = num/den;
	    else
	    {
	       beta = (num-mix)/den;
               if( beta <= 0 )
	       {
                  k = maxit;
		  beta = 0;
                  if( myRank == 0 && verbose > 2 )
		     cout << "P-R restart because beta <= 0 " << endl;
	       }
	    }
	 }
	 else
	    beta = 0;
         for( int i=0 ; i < n ; i++ )
	 {
	    d[i] = -sf[i]*sf[i]*dfp[i] + beta*d[i];
            df[i] = dfp[i];
	 }
         k++;
      }
      j++;
      done = rnorm < tolerance;
   }
   outinfo.finish();

   delete[] df;
   delete[] d;
   delete[] da;
   delete[] xa;
   delete[] dfp;
}

//-----------------------------------------------------------------------
void lbfgs( EW& simulation, int n, double x[n], double sf[n], vector<Source*>& GlobalSources,
	    vector<TimeSeries*>& GlobalTimeSeries,
	    vector<TimeSeries*> & GlobalObservations,
	    int m, int myRank )
//-----------------------------------------------------------------------
// l-BFGS minimization of misfit function.
// 
// Input: simulation - simulation object
//        x          - Parameter vector, used as initial guess.
//        sf         - Scaling factors.
//        GlobalSource - Object representing the unique source.
//        GlobalTimeSeries   - TimeSeries objects, number of objects and
//                    locations should agree with the GlobalObservations vector.
//        GlobalObservations - The observed data at receivers.
//        m           - Number of previous vectors to save.
//        myRank      - ID of this processor.
//
// Output: x - The minimum is returned.
//
// Note: Configuration of the algorithm is done by a call 
//       to get_cgparameters of the simulation object.
//
//-----------------------------------------------------------------------
{
   int nvar, j, k;
   bool done = false;
   bool dolinesearch = true;
   bool fletcher_reeves=true;
   int maxit, maxrestart, varcase=0, stepselection=0, verbose=2;

   // Do not allow source to rise closer than this to the surface:
   double tolerance;
   double f,  d2f[121], rnorm;
   double* df, *d, *da, *xa, *dfp, *dx;
   bool testing=false, hscale=true;
   int nreductions = 0;
   
   simulation.get_cgparameters( maxit, maxrestart, tolerance, fletcher_reeves, stepselection,
				dolinesearch, varcase, testing );
   if( maxrestart == 0 )
      return;

   df  = new double[n];
   d   = new double[n];
   da  = new double[n];
   xa  = new double[n];
   dfp = new double[n];
   dx  = new double[n];

   if( varcase == 0 )
      nvar = 11;
   if( varcase == 1 )
      nvar = 10;
   else if( varcase == 2 )
      nvar = 9;
   else if( varcase == 3 )
      nvar = 9;

   // Initial guess for testing
   if( testing )
      nvar  = 11;

   ConvParOutput outinfo( simulation, n, nvar, myRank, verbose, false );

   if( myRank == 0 )
      cout << "Begin L-BFGS iteration by evaluating initial misfit and gradient..." << endl;

   compute_f_and_df( simulation, n, x, GlobalSources, GlobalTimeSeries, GlobalObservations,
		     varcase, f, df, d2f, testing );
   rnorm = 0;
   for( int i = 0 ; i < n ; i++ )
      rnorm = rnorm > fabs(df[i])*sf[i] ? rnorm : fabs(df[i])*sf[i];

   outinfo.print_dfmsg(f,df,sf);
   outinfo.save_step(f,df,sf,x,rnorm,0.0,0,0);

   int it = 1;
   // s and y stores the m vectors
   double* s = new double[m*n]; 
   double* y = new double[m*n]; 
   double* rho = new double[m];
   double* al = new double[m];
   double* dftemp = new double[n];

   // kf points to the first vector
   int kf = 0;

   while( it <= maxrestart && !done )
   {
      // perform the two-loop recursion (Alg 7.4) to compute search direction d=-H*df
      int nv = it-1 < m ? it-1 : m ;
      if( nv == 0 )
	 for( int i=0 ; i < n ; i++ )
	    d[i] = -sf[i]*sf[i]*df[i];
      else
      {
	 for( int j=0 ; j < n ; j++ )
	    dftemp[j] = df[j];
	 
         int vi = kf+nv-1;
	 if( vi > m-1 )
	    vi = vi-m;
	 for( int i=nv ; i >= 1 ; i-- )
	 {
            double scprod1 = 0, scprod2=0;
	    for( int j=0 ; j < n ; j++ )
	    {
	       scprod1 += s[j+n*vi]*y[j+n*vi];
	       scprod2 += s[j+n*vi]*dftemp[j];
	    }
	    rho[vi] = 1/scprod1;
	    al[vi] = scprod2*rho[vi];
	    for( int j=0 ; j < n ; j++ )
	       dftemp[j] -= al[vi]*y[j+n*vi];
	    vi--;
	    if( vi == -1 )
	       vi = nv-1;
	 }

	 // H0*df
         if( hscale || nv == 0 )
	 {
	    // Use scale factors for diagonal H0
	    for( int  i=0 ; i < n ; i++ )
	       d[i] = sf[i]*sf[i]*dftemp[i];
	 }
	 else
	 {
	    // Use diagonal H0 given by formula (7.20) 
	    double gamn = 0, gamd=0;
            vi = kf+nv-1;
	    if( vi > m-1 )
	       vi = vi-m;
	    for( int i=0 ; i < n ; i++ )
	    {
	       gamn += s[i+n*vi]*y[i+n*vi];
	       gamd += y[i+n*vi]*y[i+n*vi];
	    }
	    gamn = gamn/gamd;
	    for( int  i=0 ; i < n ; i++ )
	       d[i] = gamn*dftemp[i];
	 }
	 
         vi = kf;
	 for( int i=1 ; i <= nv ; i++ )
	 {
            double scprod1 = 0;
	    for( int j=0 ; j < n ; j++ )
	       scprod1 += y[j+n*vi]*d[j];
	    for( int j=0 ; j < n ; j++ )
	       d[j] += (al[vi]-rho[vi]*scprod1)*s[j+n*vi];
            vi++;
	    if( vi > nv-1 )
	       vi = 0;
	 }
	 for( int i=0 ; i < n ; i++ )
	    d[i] = -d[i];
      }
   // Done with two-loop recursion (Alg 7.4) 

   // Next do line search, or update with constant step length
      double alpha = 1.0;
      if( dolinesearch )
      {
	 for( int i=0 ; i < n ; i++ )
	    da[i] = d[i];
	 int retcode;
         double fp;
	 if( myRank == 0 && verbose > 2 )
	    cout << "Line search.. " << endl;
	 linesearch( simulation, GlobalSources, GlobalTimeSeries, GlobalObservations, nvar, n,
		     x, f, df, da, fabs(alpha), 10.0, tolerance*0.01, xa, fp, sf, myRank, retcode, nreductions, testing );
	 if( myRank == 0 && verbose > 2 )
	    cout << " .. return code "  << retcode << " misfit changed from " << f << " to " << fp << endl;
      }
      else
	 for( int i=0 ; i < n ; i++ )
	    xa[i] = x[i] + alpha*d[i];
      
   // xa is now the new iterate
      double dxnorm = 0;
      for( int i=0 ; i < nvar ; i++ )
      {
	 double locnorm = fabs(x[i]-xa[i]);
	 if( fabs(xa[i])> sf[i] )
	    locnorm /= fabs(xa[i]);
	 else
	    locnorm /= sf[i];
	 if( locnorm > dxnorm )
	    dxnorm = locnorm;
	 dx[i] = xa[i] - x[i];
	 x[i]  = xa[i];
      }
      for( int i=nvar ; i < 11 ; i++ )
	 dx[i] = 0;
      for( int i=11 ; i < n ; i++ )
      {
	 double locnorm = fabs(x[i]-xa[i]);
	 if( fabs(xa[i])> sf[i] )
	    locnorm /= fabs(xa[i]);
	 else
	    locnorm /= sf[i];
	 if( locnorm > dxnorm )
	    dxnorm = locnorm;
	 dx[i] = xa[i] - x[i];
	 x[i]  = xa[i];
      }

      compute_f_and_df( simulation, n, x, GlobalSources, GlobalTimeSeries, GlobalObservations, varcase,
			f, dfp, d2f, testing );
      rnorm = 0;
      for( int i=0 ; i < n ; i++ )
	 if( fabs(dfp[i])*sf[i] > rnorm )
	    rnorm = fabs(dfp[i])*sf[i];

// Save the time series after each sub-iteration.
      for( int ts=0 ; ts < GlobalTimeSeries.size() ; ts++ )
	 GlobalTimeSeries[ts]->writeFile();

// Print some information on standard output
      outinfo.print_xmsg( x, rnorm, dxnorm, it );
      outinfo.print_dfmsg( f, dfp, sf );

// Save convergence and parameters to file
      outinfo.save_step( f, dfp, sf, x, rnorm, dxnorm, nreductions, it );

   // Update s and y vectors
      if( it > m )
      {
	 for( int i=0 ; i < n ; i++ )
	 {
	    s[i+n*kf] = dx[i];
	    y[i+n*kf] = dfp[i]-df[i];
	 }
	 kf++;
	 if( kf > m-1 )
	    kf = 0;
      }
      else
      {
	 for( int i=0 ; i < n ; i++ )
	 {
	    s[i+n*(it-1)] = dx[i];
	    y[i+n*(it-1)] = dfp[i]-df[i];
	 }
      }

   // Advance one iteration
      it++;
      done = rnorm < tolerance;
      for( int i=0 ; i < n ; i++ )
	 df[i] = dfp[i];
   }
   outinfo.finish();
   delete[] s;
   delete[] y;
   delete[] dftemp;
   delete[] al;
   delete[] rho;
   delete[] df;
   delete[] d;
   delete[] da;
   delete[] xa;
   delete[] dfp;
   delete[] dx;
}

//-----------------------------------------------------------------------
void bfgs( EW& simulation, int n, double x[n], double sf[n],
	   vector<Source*>& GlobalSources,
	   vector<TimeSeries*>& GlobalTimeSeries,
	   vector<TimeSeries*> & GlobalObservations,
	   int myRank )

//-----------------------------------------------------------------------
// BFGS minimization of misfit function.
// 
// Input: simulation - simulation object
//        x          - Parameter vector, used as initial guess.
//        sf         - Scaling factors.
//        GlobalSource - Object representing the unique source.
//        GlobalTimeSeries   - TimeSeries objects, number of objects and
//                    locations should agree with the GlobalObservations vector.
//        GlobalObservations - The observed data at receivers.
//        myRank      - ID of this processor.
//
// Output: x - The minimum is returned.
//
// Note: Configuration of the algorithm is done by a call 
//       to get_cgparameters of the simulation object.
//
//-----------------------------------------------------------------------
{
   int nvar;
   bool done = false;
   bool dolinesearch = true;
   int maxit, maxrestart, varcase=0, stepselection=0, verbose=2;
   bool fletcher_reeves=true;

   // Do not allow source to rise closer than this to the surface:
   double tolerance;
   double f, d2f[121], rnorm;
   double* df, *d, *hinv, *da, *xa, *dfp, *y, *z;
   bool testing=false;
   int nreductions=0;

   simulation.get_cgparameters( maxit, maxrestart, tolerance, fletcher_reeves, stepselection,
				dolinesearch, varcase, testing );
   if( maxrestart == 0 )
      return;

   if( varcase == 0 )
      nvar = 11;
   else if( varcase == 1 )
      nvar = 10;
   else if( varcase == 2 )
      nvar = 9;
   else if( varcase == 3 )
      nvar = 9;

   if( testing )
      nvar = 11;

   df  = new double[n];
   d   = new double[n];
   da  = new double[n];
   xa  = new double[n];
   dfp = new double[n];
   y   = new double[n];
   z   = new double[n];
   hinv= new double[n*n];

   ConvParOutput outinfo( simulation, n, nvar, myRank, verbose, false );

   if( myRank == 0 )
      cout << "Begin BFGS iteration by evaluating initial misfit and gradient..." << endl;
   compute_f_and_df( simulation, n, x, GlobalSources, GlobalTimeSeries, GlobalObservations,
		     varcase, f, df, d2f, testing );

   rnorm = 0;
   for( int i = 0 ; i < n ; i++ )
      rnorm = rnorm > fabs(df[i])*sf[i] ? rnorm : fabs(df[i])*sf[i];

   outinfo.print_dfmsg(f,df,sf);
   outinfo.save_step(f,df,sf,x,rnorm,0.0,0,0);

   // Initial guess inverse Hessian, diagonal
   for( int i=0 ; i < n*n; i++ )
      hinv[i] = 0;
   for( int i=0 ; i < n ; i++ )
      hinv[i+n*i] = sf[i]*sf[i];
   
   int it = 1;
   while( it <= maxrestart && !done )
   {
      // Compute step direction
      for( int i=0 ; i < n ; i++ )
      {
         d[i] = 0;
         for( int j=0 ; j<n ; j++ )
	    d[i] -= hinv[i+n*j]*df[j];
      }

      // Compute step length and update x
      double alpha = 1.0;
      if( dolinesearch )
      {
	 for( int i=0 ; i < n ; i++ )
	    da[i] = d[i];
	 int retcode;
	 if( myRank == 0 && verbose > 2 )
	    cout << "Line search.. " << endl;
         double fp;
	 linesearch( simulation, GlobalSources, GlobalTimeSeries, GlobalObservations, nvar, n,
		     x, f, df, da, fabs(alpha), 10.0, tolerance*0.01, xa, fp, sf, myRank, retcode, nreductions, testing );
	 if( myRank == 0 && verbose > 2 )
	    cout << " .. return code "  << retcode << " misfit changed from " << f << " to " << fp << endl;
      }
      else
      {
	 for( int i=0 ; i < n ; i++ )
	    xa[i] = x[i] + alpha*d[i];
      }         

      double dxnorm = 0;
      for( int i=0 ; i < nvar ; i++ )
      {
	 double locnorm = fabs(x[i]-xa[i]);
	 if( fabs(xa[i])> sf[i] )
	    locnorm /= fabs(xa[i]);
	 else
	    locnorm /= sf[i];
	 if( locnorm > dxnorm )
	    dxnorm = locnorm;
         d[i] = xa[i]-x[i];
	 x[i] = xa[i];
      }
      for( int i=11 ; i < n ; i++ )
      {
	 double locnorm = fabs(x[i]-xa[i]);
	 if( fabs(xa[i])> sf[i] )
	    locnorm /= fabs(xa[i]);
	 else
	    locnorm /= sf[i];
	 if( locnorm > dxnorm )
	    dxnorm = locnorm;
         d[i] = xa[i]-x[i];
	 x[i] = xa[i];
      }
      compute_f_and_df( simulation, n, x, GlobalSources, GlobalTimeSeries, GlobalObservations, varcase,
			f, dfp, d2f, testing );
      rnorm = 0;
      for( int i=0 ; i < n ; i++ )
      {
	 if( fabs(dfp[i])*sf[i] > rnorm )
	    rnorm = fabs(dfp[i])*sf[i];
         y[i]  = dfp[i] - df[i];
         df[i] = dfp[i];
      }
      // Save the time series after each sub-iteration.
      for( int ts=0 ; ts < GlobalTimeSeries.size() ; ts++ )
	 GlobalTimeSeries[ts]->writeFile();
 
// Print some information on standard output
      outinfo.print_xmsg( x, rnorm, dxnorm, it );
      outinfo.print_dfmsg( f, dfp, sf );

// Save convergence and parameters to file
      outinfo.save_step( f, dfp, sf, x, rnorm, dxnorm, nreductions, it );

  // Update inverse Hessian
      // z = G*y, sc1=y'*s, sc2=y'*z
      double sc1=0, sc2=0;
      for( int i=0 ; i < n ; i++ )
      {
         z[i] = 0;
	 for( int j=0 ; j < n ; j++ )
	    z[i] += hinv[i+n*j]*y[j];
         sc1 += y[i]*d[i];
	 sc2 += y[i]*z[i];
      }
      // Perform the update:      
      double a1=1/sc1;
      double a2=sc2*a1*a1+a1;
      for( int j=0 ; j < n ; j++ )
	 for( int i=0 ; i < n ; i++ )
	    hinv[i+n*j] += -a1*(d[i]*z[j]+z[i]*d[j])+a2*d[i]*d[j];

      // Move on to the next iteration
      it++;
      done = rnorm < tolerance;
   }
   outinfo.finish();

   delete[] df;
   delete[] d;
   delete[] da;
   delete[] xa;
   delete[] dfp;
   delete[] y;
   delete[] z;
   delete[] hinv;
}

//-----------------------------------------------------------------------
void steepest_descent( EW& simulation, int n, double x[n], double sf[n],
		       vector<Source*>& GlobalSources,
		       vector<TimeSeries*>& GlobalTimeSeries,
		       vector<TimeSeries*> & GlobalObservations,
		       int myRank )
   
//-----------------------------------------------------------------------
// Steepest descent minimization of misfit function.
// 
// Input: simulation - simulation object
//        x          - Parameter vector, used as initial guess.
//        sf         - Scaling factors.
//        GlobalSource - Object representing the unique source.
//        GlobalTimeSeries   - TimeSeries objects, number of objects and
//                    locations should agree with the GlobalObservations vector.
//        GlobalObservations - The observed data at receivers.
//        myRank      - ID of this processor.
//
// Output: x - The minimum is returned.
//
// Note: Configuration of the c-g algorithm is done by a call 
//       to get_cgparameters of the simulation object.
//
//-----------------------------------------------------------------------
{
   int nvar;
   bool done = false;
   bool dolinesearch = true;
   int maxit, maxrestart, varcase=0, stepselection=0,verbose=2;
   bool fletcher_reeves=true;

   double tolerance;
   double f, d2f[121], rnorm;
   double* df, *d, *da, *xa;
   bool testing=false;
   int nreductions=0;

   simulation.get_cgparameters( maxit, maxrestart, tolerance, fletcher_reeves, stepselection,
				dolinesearch, varcase, testing );
   if( maxrestart == 0 )
      return;
   
   df  = new double[n];
   d   = new double[n];
   da  = new double[n];
   xa  = new double[n];

   if( varcase == 0 )
      nvar = 11;
   else if( varcase == 1 )
      nvar = 10;
   else if( varcase == 2 )
      nvar = 9;
   else if( varcase == 3 )
      nvar = 9;

// Initial guess for testing
   if( testing )
      nvar = 11;

   ConvParOutput outinfo( simulation, n, nvar, myRank, verbose, false );

   if( myRank == 0 && verbose > 0 )
      cout << "Begin steepest descent iteration by evaluating initial misfit and gradient..." << endl;
   compute_f_and_df( simulation, n, x, GlobalSources, GlobalTimeSeries, GlobalObservations,
		     varcase, f, df, d2f, testing );

// Compute norm of gradient
   rnorm = 0;
   for( int i=0 ; i < n ; i++ )
      rnorm = rnorm > fabs(df[i])*sf[i] ? rnorm : fabs(df[i])*sf[i];

   outinfo.print_dfmsg(f,df,sf);
   outinfo.save_step(f,df,sf,x,rnorm,0.0,0,0);

   int it = 0;
   while( it <= maxrestart && !done )
   {
      for( int i=0 ; i < n ; i++ )
         d[i] = -df[i]*sf[i]*sf[i];
      double alpha, fp;

// Compute initial step length
      if( stepselection == 0 )
      {
         double den = 0, dnrm=0, xnrm=0;
	 for( int i= 0 ; i < n ; i++ )
	 {
	    den  += d[i]*df[i];
	    dnrm += d[i]*d[i];
	    xnrm += x[i]*x[i];
	 }
	 double stlim = 0.1*sqrt(xnrm/dnrm);
	 double h=-2*f/den;
	 if( h > stlim )
	    h = stlim;
	 for( int i=0; i<n ; i++ )
	    xa[i] = x[i] + h*d[i];
	 if( myRank == 0 && verbose > 2 )
	       cout << "Step length computation a) " << endl;
	 compute_f( simulation, n, xa, GlobalSources, GlobalTimeSeries, GlobalObservations, fp, testing );
	 alpha = h*f/(fp+f);
      }
      else
      {
	 double dtHd;
	 if( myRank == 0 && verbose > 2 )
	    cout << "Step length computation b) " << endl;
	 compute_dtd2fd( simulation, n, x, GlobalSources, GlobalTimeSeries, GlobalObservations, d, varcase, dtHd,
			    true, d2f, myRank, testing );
	 outinfo.print_scalar( dtHd, "dtHd= " ,3 );
	 outinfo.print_vector( d, "d= " ,3 );
	 alpha = 0;
	 for( int i=0; i<n ; i++ )
	    alpha += df[i]*d[i];
	 alpha = -alpha/dtHd;
      }

// Determine step length and update x.  xa is the new iterate.
      if( dolinesearch )
      {
	 for( int i=0 ; i < n ; i++ )
	    da[i] = d[i];
	 int retcode;
	 if( myRank == 0 && verbose > 2 )
	    cout << "Line search.. " << endl;
	 linesearch( simulation, GlobalSources, GlobalTimeSeries, GlobalObservations, nvar, n,
		     x, f, df, da, fabs(alpha), 100.0, tolerance*0.01, xa, fp, sf, myRank, retcode, nreductions, testing );
	 if( myRank == 0 && verbose > 2 )
	    cout << " .. return code "  << retcode << " misfit changed from " << f << " to " << fp << endl;
      }
      else
      {
	 for( int i=0 ; i < n ; i++ )
	    xa[i] = x[i] + alpha*d[i];
      }         

 // Compute the norm of x_k-x_{k+1}, and update x := xa
      double dxnorm = 0;
      for( int i=0 ; i < nvar ; i++ )
      {
	 double locnorm = fabs(x[i]-xa[i]);
	 if( fabs(xa[i])> sf[i] )
	    locnorm /= fabs(xa[i]);
	 else
	    locnorm /= sf[i];
	 if( locnorm > dxnorm )
	    dxnorm = locnorm;
	 x[i]=xa[i];
      }
      for( int i=11 ; i < n ; i++ )
      {
	 double locnorm = fabs(x[i]-xa[i]);
	 if( fabs(xa[i])> sf[i] )
	    locnorm /= fabs(xa[i]);
	 else
	    locnorm /= sf[i];
	 if( locnorm > dxnorm )
	    dxnorm = locnorm;
	 x[i]=xa[i];
      }
      compute_f_and_df( simulation, n, x, GlobalSources, GlobalTimeSeries, GlobalObservations, varcase, f, df, d2f,
			testing );

// Compute the norm of the gradient
      rnorm = 0;
      for( int i=0 ; i < n ; i++ )
	 if( fabs(df[i])*sf[i] > rnorm )
	    rnorm = fabs(df[i])*sf[i];

// Save the time series after each iteration.
      for( int ts=0 ; ts < GlobalTimeSeries.size() ; ts++ )
	 GlobalTimeSeries[ts]->writeFile();

// Print some information on standard output
      outinfo.print_xmsg(x,rnorm,dxnorm,it);
      outinfo.print_dfmsg(f,df,sf);

// Save convergence and parameters to file
      outinfo.save_step(f,df,sf,x,rnorm,dxnorm,nreductions,it);

// Move on to the next iteration
      it++;
      done = rnorm < tolerance;
   }
   outinfo.finish();
   delete[] df;
   delete[] d;
   delete[] da;
   delete[] xa;
}

//-----------------------------------------------------------------------
int main(int argc, char **argv)
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
  vector<TimeSeries*> GlobalTimeSeries; // GlobalTimeSeries holds the synthetic motions
  vector<TimeSeries*> GlobalObservations; // GlobalObservations holds the observed motions

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

// get the simulation object ready for time-stepping
     simulation.setupRun( GlobalSources );
     //     simulation.preprocessSources( GlobalSources );  // AP: preprocessSources called for GlobalSources
     if (!simulation.isInitialized())
     { 
	if (myRank == 0)
	{
	   cout << "Error: simulation object not ready for time stepping" << endl;
	}
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

// run the forward solver in stealth mode
       simulation.setQuiet(true);

     // Filter observed data if required
	for( int m = 0; m < GlobalObservations.size(); m++ )
	{
	   GlobalObservations[m]->set_utc_to_simulation_utc();
	   if( simulation.m_prefilter_sources && simulation.m_filter_observations )
	   {
	      GlobalObservations[m]->filter_data( simulation.m_filterobs_ptr );
     //	      GlobalObservations[m]->writeFile( "_fd" ); // AP changed the extension to _fd = filtered data
	   }
	}

//  First copy observations to GlobalTimeSeries, later, solve will insert 
//  the simulation time step and start time into GlobalTimeSeries.
	for( int m = 0; m < GlobalObservations.size(); m++ )
	{
	   string newname = "_out";
	   TimeSeries *elem = GlobalObservations[m]->copy( &simulation, newname, true );
	   // Reset station utc reference point to simulation utc reference point,
	   // reasonable since GlobalTimeSeries will be used with the forward solver.
           elem->set_utc_to_simulation_utc();
	   GlobalTimeSeries.push_back(elem);
	}

	if (myRank == 0)
	{
	   cout << "Running sw4opt on " <<  nProcs << " processors..." << endl
		<< "Writing output to directory: " 
		<< simulation.getOutputPath() << endl;
	}
	
// Variables needed:
	int maxit, maxrestart, varcase=0, stepselection=0;
	bool dolinesearch, fletcher_reeves=true, testing=false;
	double tolerance;
	simulation.get_cgparameters( maxit, maxrestart, tolerance, fletcher_reeves, stepselection,
				     dolinesearch, varcase, testing );
	// varcase: 0-solve for 11 source parameters
	//          1-solve for 10 source parameters (no frequency),
	//          2-solve for 9 source parameters, (no t0, no frequency)
	//          3-solve for 9 source parameters + one shift parameter for each observation
	// The scalefactor array (sf) and the parameter array (xv) are allocated to always
	// have 11 source parameters. 
	double* xv, *sf;
        int n;
	if( varcase != 3 )
	{
	   xv = new double[11];
           sf = new double[11];
           n  = 11;
	}
	else
	{
	   sf = new double[11+GlobalObservations.size()];
	   xv = new double[11+GlobalObservations.size()];
           n  = 11+GlobalObservations.size();
	}
     
// Initial guess
        bool output_initial_seismograms = false;

//   Default guess, the input source, stored in GlobalSources[0]
        bool guesspos, guesst0fr, guessmom, guessshifts;
        simulation.compute_guess( guesspos, guesst0fr, guessmom, guessshifts, output_initial_seismograms );
	GlobalSources[0]->get_parameters(xv);

	// Compute one forward solve, and guess shift. Use input source parameters for the synthetics.
        if( output_initial_seismograms || (varcase==3 && guessshifts) )
	{
           simulation.print_utc();
           vector<TimeSeries*> localTimeSeries;
	   for( int m = 0; m < GlobalObservations.size(); m++ )
	   {
	      string newname = "_ini";
	      TimeSeries *elem = GlobalObservations[m]->copy( &simulation, newname, true );
	      elem->set_utc_to_simulation_utc();
	      localTimeSeries.push_back(elem);
	      if (simulation.getVerbosity()>=2)
		 GlobalObservations[m]->print_timeinfo();
	   }
	   if (myRank == 0)
	   {
	     cout << "Solving a forward problem to compute initial seismograms..." << endl;
	   }
	   simulation.solve( GlobalSources, localTimeSeries );
           if( output_initial_seismograms )
	   {
	      for (int ts=0; ts<localTimeSeries.size(); ts++)
	      {
		 localTimeSeries[ts]->writeFile();
		 if (simulation.getVerbosity()>=2)
		    localTimeSeries[ts]->print_timeinfo();
	      }
	   }
	   if( varcase == 3 && guessshifts )
	   {
              double dt = simulation.getTimeStep();
              guess_shift( n, xv, localTimeSeries, GlobalObservations, dt );
		 if( myRank == 0 )
		 {
		    cout  << "Computed initial time shifts... "<<endl;
		    for (int ts=0; ts<localTimeSeries.size(); ts++)
		      cout << "station '" << localTimeSeries[ts]->getStationName() << "' t-shift = " << xv[11+ts] << endl;
		 }
	   }
	   for( int ts=0 ; ts<localTimeSeries.size();ts++)
	      delete localTimeSeries[ts];
	} // end if( output_initial_seismograms || (varcase==3 && guessshifts) )

//  Set initial shifts to zero if not estimated.
        if( varcase == 3 && !guessshifts )
	   for( int i=11 ; i < n ; i++ )
	      xv[i] = 0;

// Output shifted observation data
        for( int m=0 ; m < GlobalObservations.size() ; m++ )
	{
           if( varcase == 3 )
	      GlobalObservations[m]->add_shift( xv[11+m]);
           GlobalObservations[m]->writeFile( "_fd" ); 
	}
 	   
	if( guesspos )
	{
	   // Guess source position, and perhaps t0
	   double tstart;
	   //	   guess_source_position( simulation, GlobalSources, GlobalTimeSeries, GlobalObservations,
	   //				  xv[0], xv[1], xv[2], tstart, myRank );
	   guess_source_position_layer( simulation, GlobalSources, GlobalObservations,
				  xv[0], xv[1], xv[2], tstart, myRank );
	   if( guesst0fr )
	      // Guess t0
	      guess_source_t0freq( simulation, GlobalSources,
				   tstart, xv[9], xv[10], myRank );
	   // Update source with guess(es)
	   GlobalSources[0]->set_parameters( xv );
	}
	if( guessmom )
	// Estimate the moment tensor components, will update both GlobalSources[0] and xv
	   guess_source_moments( simulation, GlobalSources, GlobalTimeSeries, GlobalObservations,
				 xv, myRank );

	if( myRank == 0 )
	{
	   cout << "Initial source guess : \n";
	   cout << "   x0 = " << xv[0] << " y0 = " << xv[1] << " z0 = " << xv[2] <<endl;
	   cout << "  mxx = " << xv[3] << " mxy= " << xv[4] << " mxz= " << xv[5] << endl;
	   cout << "  myy = " << xv[6] << " myz= " << xv[7] << " mzz= " << xv[8] << endl;
	   cout << "   t0 = " << xv[9] << " freq = " << xv[10] << endl;
	}

	// nvar is number of source parameters
	int nvar;
        if( varcase == 0 )
	   nvar = 11;
	else if( varcase == 1 )
	   nvar = 10;
	else if( varcase == 2 )
	   nvar = 9;
	else if( varcase == 3 )
	   nvar = 9;

	compute_scalefactors( simulation, GlobalSources, GlobalTimeSeries, GlobalObservations, myRank,
			      varcase, nvar, sf );
	if( myRank == 0 )
	{
           cout << "scalefactors x0=" << sf[0] << " y0=" << sf[1] << " z0=" << sf[2] << " Mxx=" << sf[3]
		<< " Mxy=" << sf[4] << " Mxz=" << sf[5] << " Myy=" << sf[6]  << " Myz=" << sf[7]
		<< " Mzz=" << sf[8] << " t0=" << sf[9] << " freq=" << sf[10] << endl;
           if( varcase == 3 )
	      for( int m= 0 ; m < GlobalTimeSeries.size() ; m++ )
		 cout << "sf-obs" << m << "= " << sf[11+m]<<endl;
	}

	// Shift back, the optimization code assumes that GlobalObservations are unshifted.
        if( varcase == 3 )
	   for( int m=0 ; m < GlobalObservations.size() ; m++ )
	      GlobalObservations[m]->add_shift(-xv[11+m]);

        if( simulation.m_opttest == 1 )
	   testsourced2( simulation, GlobalSources );
        else if( simulation.m_opttest == 2 )
	   test_gradient( simulation, GlobalSources, GlobalTimeSeries, GlobalObservations, myRank, sf );
	else if( simulation.m_opttest == 3 )
	   test_hessian( simulation, GlobalSources, GlobalTimeSeries, GlobalObservations, myRank, sf );
        else if( simulation.m_opttest == 4 )
           for( int i=1 ; i < 16 ;  i++ )
	   {
	      // Function to minimize, plotted along one direction.
	      //	      double mzz= (-1+2.0*i/30.0)*1e18;
	      //              double la = 0.0066*(-1+2.0*i/15);
              double la = -i/1500000.0;
              double f;
	      //              double p[11]= {0.0090176, 0.0090176, 0.0143213, 1.07286e+12, 2.44711e+12, -6.93935e+10,
	      //                             1.07286e+12, -6.93935e+10, 4.12262e+09, 4.60064e-06, 1.0537e-06};
	      double p[11] = {5.2552e-03, 1.0995e-02, -1.8408e-03, 3.5146e-03, 1.5270e-03, -2.0704e-04, 7.3706e-04, -4.1574e-05, 3.7204e-07, -2.9325e-04, 3.9289e-02};
              double yv[11];
	      for( int j=0; j < 11 ; j++ )
		 yv[j] = xv[j] + la*p[j];
              compute_f( simulation, 11, yv, GlobalSources, GlobalTimeSeries, GlobalObservations, f );
              if( myRank == 0 )
		 //		 cout << mzz << " " << f << endl;
		 cout << la << " " << f << endl;
	   }
        else if( simulation.m_opttest == 5 )
	{
	   // Function surface
           int ni = 50, nj=50;
	   // Standard
           xv[0] = 15000;
	   xv[1] = 15000;
           xv[2] = 2000;
	   // Offset
	   //	   xv[0] = 13000;
	   //	   xv[1] = 13000;
	   //	   xv[2] =  2000;
	   xv[3] = 0;
           xv[4] = 1e18;
	   xv[5] = xv[6] = xv[7] = xv[8] = 0;
	   xv[9] = 1.45;
	   xv[10] = 6.0;
           double xmin = 10000, ymin=10000, xmax=20000, ymax=20000;
           double zmin = 200, zmax=8200;
	   double* fsurf = new double[ni*nj];
           for( int j= 0 ; j < nj ; j++ )
	      for( int i= 0 ; i < ni ; i++ )
	      {
                 xv[0] = xmin + i*(xmax-xmin)/(ni-1);
                 xv[1] = ymin + j*(ymax-ymin)/(nj-1);
		 //                 xv[2] = zmin + j*(zmax-zmin)/(nj-1);
		 compute_f( simulation, 11, xv, GlobalSources, GlobalTimeSeries, GlobalObservations, fsurf[i+ni*j] );
	      }
           if( myRank == 0 )
	   {
              int fd=open("funcsurf.bin", O_CREAT | O_TRUNC | O_WRONLY , 0660 );
	      int nr = write(fd,&ni,sizeof(int));
	      write(fd,&nj,sizeof(int));
              write(fd,&xmin,sizeof(double));
              write(fd,&ymin,sizeof(double));
              write(fd,&xmax,sizeof(double));
              write(fd,&ymax,sizeof(double));
	      write(fd,fsurf,sizeof(double)*ni*nj);
	      close(fd);
	   }
	}
	else
	{
	   // run optimization
           int method;
	   int bfgs_m;
	   simulation.get_optmethod( method, bfgs_m );
	   if( method == 1 )
	      cg( simulation, n, xv, sf, GlobalSources, GlobalTimeSeries,
		  GlobalObservations, myRank );
           else if( method == 2 )
	      lbfgs( simulation, n, xv, sf, GlobalSources, GlobalTimeSeries,
		     GlobalObservations, bfgs_m, myRank );
	   else if( method == 3 )
	      bfgs( simulation, n, xv, sf, GlobalSources, GlobalTimeSeries,
		    GlobalObservations, myRank );
	   else if( method == 4 )
	      steepest_descent( simulation, n, xv, sf, GlobalSources, GlobalTimeSeries,
				GlobalObservations, myRank );
// Save all time series
//	   for (int ts=0; ts<GlobalTimeSeries.size(); ts++)
//	   {
//	      GlobalTimeSeries[ts]->writeFile();
//	   }
// Done, remove memory
	   for (int ts=0; ts<GlobalTimeSeries.size(); ts++)
	      delete GlobalTimeSeries[ts];
	   for (int ts=0; ts<GlobalObservations.size(); ts++)
	      delete GlobalObservations[ts];
	}

	if( myRank == 0 )
	{
	   cout << "============================================================" << endl
		<< " sw4opt ( Source estimation solver) finished! " << endl
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
void exact_guess_source_position( double xr[4], double yr[4], double dist[4],
				  double& x0, double& y0, double& z0, double& d0, int& fail )
//
// Assumes that zr=0 at all stations
//
{
   double mat[9], x[3], rhs[3];
   int j;
   for( int i=0 ; i<3 ; i++ )
   {
      j = 0;
      mat[i+3*j] = xr[0]-xr[i+1];
      j = 1;
      mat[i+3*j] = yr[0]-yr[i+1];
      j = 2;
      mat[i+3*j] = -(dist[0]-dist[i+1]);
      rhs[i] = 0.5*(xr[0]*xr[0]+yr[0]*yr[0]-dist[0]*dist[0]-(x[i+1]*x[i+1]+yr[i+1]*yr[i+1]-dist[i+1]*dist[i+1]));
   }
   int three=3, one=1, ipiv[3], info;
   F77_FUNC(dgesv,DGESV)( &three, &one, mat, &three, ipiv, rhs, &three, &info );
   if( info != 0 )
   {
      fail = info;
      return;
   }
   x0 = rhs[0];
   y0 = rhs[1];
   d0 = rhs[2];
   z0 = - (x0-xr[0])*(x0-xr[0]) - (y0-yr[0])*(y0-yr[0]) + (d0-dist[0])*(d0-dist[0]);
   if( z0 < 0 )
   {
      fail = 4;
      //      z0   = 0;
   }
   else
   {
      fail = 0;
      z0   = sqrt(z0);
   }
}

//-----------------------------------------------------------------------
void one_ray( double xr, double yr, double zr, double xsrc, double ysrc,
	      double zsrc, vector<double>& cpv, vector<double>& zv,
	      int& n, double& t, int ifail[2], int myRank )
{
   ifail[0] = 0;
   ifail[1] = 0;
   if( zsrc < zv[0] )
   {
      // Homogeneous material
      t = sqrt( (xr-xsrc)*(xr-xsrc)+(yr-ysrc)*(yr-ysrc)+(zr-zsrc)*(zr-zsrc))/cpv[0];
      n = 1;
   }
   else
   {
      // find layer of source
      n = zv.size();
      while( zsrc < zv[n-1] && (n>0) )
	 n--;
      double* xpts = new double[n+1];
      double* ypts = new double[n+1];
      // initial guess
      //      cout << "n= " << n << endl;
      for( int i=1 ; i <= n ; i++ )
      {
	 xpts[i] = xr + i*(xsrc-xr)/(n+1);
	 ypts[i] = yr + i*(ysrc-yr)/(n+1);
      }
      //            if( myRank == 0 )
      //            {
      //	       cout << " zsrc= " << zsrc << endl;
      //      	 cout << "guess in ray n= " << n << endl; 
      //               for( int i=0 ; i < n ; i++ )
      //		  cout << " cp = " << cpv[i] << " z=" << zv[i] << endl;
      //               for( int i=1 ; i <= n ; i++ )
      //            cout << " x,y= " << xpts[i] << " " << ypts[i] << endl;
      //            }
      // Fixed point iteration
      int maxit = 200, it=0;
      double tol=1e-10, er=1;
      double *b = new double[n];
      double *c = new double[n];
      double *rhs = new double[2*n];
      double *a = new double[(n)*(n)];
      double *x = new double[n];
      double *y = new double[n];
      double *d = new double[n+1];
      int* ipiv = new int[n];

      while( (er > tol) && (it<maxit) )
      {
	 for( int i=0 ; i < n ; i++ )
	    b[i]=c[i]=x[i]=y[i]=0;
	 for( int i=0 ; i < 2*n ; i++ )
	    rhs[i]=0;
	 for( int i=0 ; i < n*n ; i++ )
	    a[i]=0;
	 d[0] = cpv[0]*sqrt( (xpts[1]-xr)*(xpts[1]-xr) + (ypts[1]-yr)*(ypts[1]-yr)+(zv[0]-zr)*(zv[0]-zr));
         for( int i=2 ; i <= n ; i++ )
	    d[i-1] = cpv[i-1]*sqrt( (xpts[i]-xpts[i-1])*(xpts[i]-xpts[i-1])+
			     (ypts[i]-ypts[i-1])*(ypts[i]-ypts[i-1])+ (zv[i-1]-zv[i-2])*(zv[i-1]-zv[i-2]));
	 d[n]= cpv[n]*sqrt( (xsrc-xpts[n])*(xsrc-xpts[n])+(ysrc-ypts[n])*(ysrc-ypts[n])+
			 (zsrc-zv[n-1])*(zsrc-zv[n-1]));
	 if( n > 1 )
	 {
	    //            b[0]   = xr/d[0];
	    //	    b[n-1] = xsrc/d[n];
	    //            c[0]   = yr/d[0];
	    //	    c[n-1] = ysrc/d[n];
            rhs[0  +n*0] =  xr/d[0];
	    rhs[n-1+n*0] = xsrc/d[n];
	    rhs[0  +n*1] = yr/d[0];
	    rhs[n-1+n*1] = ysrc/d[n];
	 }
         else
	 {
	    //            b[0] = xr/d[0] + xsrc/d[n];
	    //            c[0] = yr/d[0] + ysrc/d[n];
	    rhs[0] = xr/d[0]+xsrc/d[n];
	    rhs[1] = yr/d[0]+ysrc/d[n];
	 }
	 for( int i=1 ; i <= n-1 ; i++ )
	 {
            a[i-1+n*(i-1)] = 1/d[i-1]+1/d[i];
            a[i-1+n*i]     = -1/d[i];
	    a[i+n*(i-1)]   = -1/d[i];
	 }
	 a[n-1+n*(n-1)] = 1/d[n-1]+1/d[n];
         int two=2, info;
         F77_FUNC(dgesv,DGESV)( &n, &two, a, &n, ipiv, rhs, &n, &info );
	 //	 F77_FUNC(linsolvelu,LINSOLVEU)( &n, a, b, x );
	 //	 F77_FUNC(linsolvelu,LINSOLVEU)( &n, a, c, y );
         if( info != 0 )
	 {
	    cout <<"Error solving system in one_ray. info= " << info << endl;
            ifail[0] = 1;
	    ifail[1] = info;
	    return;
	 }
         for( int i=0 ; i < n ; i++ )
	 {
            x[i] = rhs[i+n*0];
	    y[i] = rhs[i+n*1];
	 }
	 er = 0;
	 for( int i=1 ; i <= n ; i++ )
	 {
            if( er < fabs(xpts[i]-x[i-1]) )
	       er = fabs(xpts[i]-x[i-1]);
	    if( er < fabs(ypts[i]-y[i-1]) )
	       er = fabs(ypts[i]-y[i-1]);
            xpts[i] = x[i-1];
	    ypts[i] = y[i-1];
	 }
	 //         if( myRank == 0 )
	 //         cout << "ray err = " << er << endl;
         it++;
      }
      if( er > tol )
      {
	 ifail[0] = 2;
         if( myRank == 0 )
	    cout << "WARNING: no convergence in one_ray err= "<< er << " tol= " << tol << endl;
      }
      //      if( myRank == 0 )
      //	 for( int i=1; i <= n ; i++ )
      //	    cout << " x1 = " << xpts[i]  << " y1 = "  << ypts[i] << endl;
      t = 0;
      for( int i=0 ; i <= n ;i++ )
	 t += d[i]/(cpv[i]*cpv[i]);
      delete[] a, b, c, d, x, y, xpts, ypts, ipiv, rhs;
   }
}

//-----------------------------------------------------------------------
void guess_source_position_layer( EW& simulation, vector<Source*>& sources, 
				  vector<TimeSeries*>& observations,
				  double& x0, double& y0, double& z0, double& t0, int myRank )
{

   int no = observations.size();
   double* trl = new double[no];
   double* xrl = new double[no];
   double* yrl = new double[no];
   double* zrl = new double[no];

   vector<double> cpv;
   vector<double> zv;
   simulation.layered_speeds(cpv,zv);
   //   // LOH1
   //   cpv.push_back(4000);
   //   cpv.push_back(6000);
   //   zv.push_back(1000);
   for( int s= 0 ; s < observations.size() ; s++ )
      if( observations[s]->myPoint() )
      {
	 trl[s] = observations[s]->arrival_time( 1e-6 );
	 xrl[s] = observations[s]->getX();
	 yrl[s] = observations[s]->getY();
	 zrl[s] = observations[s]->getZ();
      }
      else
	 trl[s] = xrl[s] = yrl[s] = zrl[s] = 0;


// Collect station info across the machine:
   double* tr = new double[no];
   double* xr = new double[no];
   double* yr = new double[no];
   double* zr = new double[no];

   MPI_Allreduce( trl, tr, no, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   MPI_Allreduce( xrl, xr, no, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   MPI_Allreduce( yrl, yr, no, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   MPI_Allreduce( zrl, zr, no, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

   delete[] trl, xrl, yrl, zrl;

   double x, y, t, h;
   int n;

   // Initial guess for t0, manually determined arrivals (for geyser's eq.):
   t0 = 14.1;
   tr[0] = 34;
   tr[1] = 33;
   tr[2] = 24;
   tr[3] = 32.3;
   tr[4] = 21;
   tr[5] = 27;
   tr[6] = 35.5;

   // Initial guess for position
   x0 = sources[0]->getX0();
   y0 = sources[0]->getY0();
   z0 = sources[0]->getZ0();

// Gauss-Newton to solve for the source position
   h  = 1e-6*sqrt(x0*x0+y0*y0+z0*z0);
   int maxit = 500, it=0;
   double tol=1e-6, er=1;
   
   double* jac = new double[4*no];
   double* b   = new double[no];

   char job='N';
   int four=4, one=1, lwork=16, info, fail[2];   
   double diagelement=1e-8;
   double* work = new double[lwork];

   double* jactjac=new double[16];
   double* jactb = new double[4];
   int* ipiv=new int[4];

   if( myRank == 0 )
   {
      cout << "stations " << endl;
      for( int s=0 ; s < no ; s++ )
	 cout << " x = " << xr[s] << " y= " << yr[s] << " z= " << zr[s] << " t= " << tr[s] << endl;
   }
   while( (er > tol ) && (it < maxit) )
   {
      for( int r=0 ; r < no ; r++ )
      {
         one_ray( xr[r], yr[r], zr[r], x0, y0, z0, cpv, zv, n, t, fail, myRank );
	 b[r] = t + t0 - tr[r];
         double tl = t;
	 one_ray( xr[r], yr[r], zr[r], x0+h, y0, z0, cpv, zv, n, t, fail, myRank );
	 jac[r+no*0] = (t-tl)/h;
	 one_ray( xr[r], yr[r], zr[r], x0, y0+h, z0, cpv, zv, n, t, fail, myRank );
	 jac[r+no*1] = (t-tl)/h;
	 one_ray( xr[r], yr[r], zr[r], x0, y0, z0+h, cpv, zv, n, t, fail, myRank );
	 jac[r+no*2] = (t-tl)/h;
         jac[r+no*3] = 1;
      }
// Solve in least squares sense
//      if(myRank==0 )
//      {
//         cout << "matrix and rhs " << endl;
//	 for( int r=0 ; r < no ; r++ )
//	    cout << jac[r+no*0] << " " << jac[r+no*1] << " " << jac[r+no*2] << " " << jac[r+no*3] << " rhs=" << b[r] << endl;// 
//      }
      // Solve QR directly
      //      F77_FUNC(dgels,DGELS)( &job, &no, &four, &one, jac, &no, b, &no, work, &lwork, &info );
      //      if( info != 0 )
      //      {
      //	 cout <<"Error solving system in guess_source_position_layer. info= " << info << endl;
      //	 return;
      //      }

//  Form the normal equations, and add diagonal stabilizing element
      // 1. Form Jac'*Jac
      for( int j=0 ; j < 4 ; j++ )
	 for( int i=0 ; i < 4 ; i++ )
	 {
            jactjac[i+4*j] = 0;
	    for( int k=0 ; k < no ; k++ )
	       jactjac[i+4*j] += jac[k+no*j]*jac[k+no*i];
	 }
      // 2. Add diagonal stabilizer, and form Jac'*b
      for( int i= 0 ; i < 4 ; i++ )
      {
         jactjac[i+4*i] += diagelement;
         jactb[i] = 0;
	 for( int k=0 ; k < no ; k++ )
	    jactb[i] += jac[k+no*i]*b[k];
      }
      int four=4,one=1;

      //      if(myRank==0 )
      //      {
      //         cout << "normal eq matrix and rhs " << endl;
      //	 for( int r=0 ; r < 4 ; r++ )
      //	    cout << jactjac[r+4*0] << " " << jactjac[r+4*1] << " " << jactjac[r+4*2] << " " << jactjac[r+4*3] << " rhs=" << jactb[r] << endl; 
      //      }

      F77_FUNC(dgesv,DGESV)( &four, &one, jactjac, &four, ipiv, jactb, &four, &info );      
      if( info != 0 )
      {
	 cout <<"Error solving system in guess_source_position_layer. info= " << info << endl;
	 return;
      }
      for( int k=0 ; k < 4 ;k++ )
	 b[k] = jactb[k];

      er = 0;
      for( int i= 0 ; i<4 ; i++ )
	 if( er < fabs(b[i]) )
	    er = fabs(b[i]);

      if( myRank == 0 )
      {
	 cout << "err = " << er << " x0=" << x0 << " y0=" << y0 << " z0=" << z0 << " t0=" << t0 << endl;
         cout << "      dx0 = " << b[0] << " dy0 = " << b[1] << " dz0 = " << b[2] << " dt0 = " << b[3] << endl;
      }
      x0 = x0 - b[0];
      y0 = y0 - b[1];
      z0 = z0 - b[2];
      t0 = t0 - b[3];
      it++;
   }
   delete[] jac, b, work, xr, yr, zr, tr;
   delete[] jactjac, jactb, ipiv;
   if( myRank == 0 )
   {
      cout << "Guessed position of source " << x0 << " " << y0 << " " << z0 << endl;
      cout << "Guessed time offset " << t0 << endl;
   }
}

//-----------------------------------------------------------------------
void guess_source_position( EW &  simulation, vector<Source*>& sources, 
			    vector<TimeSeries*>& timeseries, vector<TimeSeries*>& observations,
			    double& x0, double& y0, double& z0, double& t0, int myRank )
{
   int n = observations.size();
   double* dist = new double[n];
   double* xr   = new double[n];
   double* yr   = new double[n];
   double* zr   = new double[n];
   double cp, cs;
   simulation.average_speeds(cp,cs);

   //   if( myRank == 0 )
   //   {
   //
   //      for( int i=0 ; i < zv.size() ; i++ )
   //         cout << " c" << i+1 << " = " << cpv[i] << " z" << i+1 << " = " << zv[i] << endl;
   //      cout << " c" << zv.size()+1 << " = " << cpv[cpv.size()-1] << endl;

   //      double xr=2.0, yr=0.0, zr=0.0;
   //      double xsrc=-0.5, ysrc=1.4, zsrc=5.0;
   //      double tt;
   //      int nn, ifail[2];
   //      vector<double> zs(3), cps(4);
   //      zs[0] = 1;
   //      zs[1] = 2;
   //      zs[2] = 4;
   //      cps[0] = 0.2;
   //      cps[1] = 2;
   //      cps[2] = 8;
   //      cps[3] = 10;
   //      one_ray( xr, yr, zr, xsrc, ysrc, zsrc, cps, zs, nn, tt, ifail );
   //      cout << " tmin = "  << tt << " n= " << nn << endl;
   //   }
   //   cp = 6000;
   //   cs = 3464;
   if( myRank == 0 )
      cout << "average speeds " << cp << " " << cs << endl;

   for( int s= 0 ; s < observations.size() ; s++ )
      if( observations[s]->myPoint() )
      {
	 dist[s] = cp*observations[s]->arrival_time( 1e-3 );
	 xr[s] = observations[s]->getX();
	 yr[s] = observations[s]->getY();
	 zr[s] = observations[s]->getZ();
//        cout << "s = " << s << " xr, yr, zr, dist = " << xr[s] << " " << yr[s] << " " << zr[s] << " " << dist[s] << endl;
      }
      else
	 dist[s] = xr[s] = yr[s] = zr[s] = 0;
// Compute initial guess for the non-linear problem by solving the problem with four stations exactly.
//
// Start with the four first stations. Collect them across the machine:
   double xrg[4],yrg[4],drg[4];
   for( int s=0 ; s < 4 ; s++ )
   {
      MPI_Allreduce( &xr[s],   &xrg[s], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
      MPI_Allreduce( &yr[s],   &yrg[s], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
      MPI_Allreduce( &dist[s], &drg[s], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   }
   int fail;
   double d0;
   exact_guess_source_position( xrg, yrg, drg, x0, y0, z0, d0, fail );

   if( fail != 0 )
   {
      if( myRank == 0 )
      {
         cout << "fail in initial guess. Fail =  " << fail << endl;
	 for( int s=0 ; s < 4 ; s++ )
	    cout << " xr= " << xrg[s] << " yr= " << yrg[s] << " dr= " << drg[s] << endl;
         if( fail == 4 )
	    cout << " z0^2 = " << z0 << endl;
      }
      x0 = sources[0]->getX0();
      y0 = sources[0]->getY0();
      z0 = sources[0]->getZ0();
      d0 = 0;
   }      

   
   //   // Least squares with initial guess from 'source' object
   //   x0 = sources[0]->getX0();
   //   y0 = sources[0]->getY0();
   //   z0 = sources[0]->getZ0();
   //   // Initial guess for d0:
   //   double d0 = 0;

   //   double sdsq=0,ressq=0, dsum=0;
   //   for( int s= 0 ; s < observations.size() ; s++ )
   //      if( observations[s]->myPoint() )
   //      {
   //         dsum  += dist[s];
   //         sdsq  += dist[s]*dist[s];
   //	 ressq += (x0-xr[s])*(x0-xr[s])+(y0-yr[s])*(y0-yr[s])+(z0-zr[s])*(z0-zr[s]);
   //      }
   //   double dsumtmp=dsum, sdsqtmp=sdsq, ressqtmp=ressq;
   //   MPI_Allreduce( &dsumtmp, &dsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   //   MPI_Allreduce( &sdsqtmp, &sdsq, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   //   MPI_Allreduce( &ressqtmp, &ressq,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   //   int nobs = observations.size();
   //   d0 = dsum/nobs + sqrt((dsum/nobs)*(dsum/nobs)-sdsq/nobs+ressq/nobs);
   if( myRank == 0 )
      cout << "initial guess x0, y0, z0, d0 " << x0 << " " << y0 << " " << z0 << " " << d0 << endl;

   double a[16], b[4], dx[4];
   int it=0, maxit=20;
   double err0=1, tol=1e-8, err=1;
   while( err > tol && it < maxit )
   {
      double res = 0;
      for( int i=0 ; i<16 ;i++ )
	 a[i] = 0;
      b[0]=b[1]=b[2]=b[3]=0;
      for( int s = 0 ; s < n ; s++ )
      {
         if( observations[s]->myPoint() )
	 {
	    double nrm=(x0-xr[s])*(x0-xr[s])+(y0-yr[s])*(y0-yr[s])+
	       (z0-zr[s])*(z0-zr[s])-(d0-dist[s])*(d0-dist[s]);
	    //            double nrm = 0;

	    a[0] += 4*(x0-xr[s])*(x0-xr[s]) + 2*nrm;
	    a[1] += 4*(x0-xr[s])*(y0-yr[s]);
	    a[2] += 4*(x0-xr[s])*(z0-zr[s]);
	    a[3] -= 4*(x0-xr[s])*(d0-dist[s]);

	    a[5] += 4*(y0-yr[s])*(y0-yr[s]) + 2*nrm;
	    a[6] += 4*(y0-yr[s])*(z0-zr[s]);
	    a[7] -= 4*(y0-yr[s])*(d0-dist[s]);

	    a[10] += 4*(z0-zr[s])*(z0-zr[s]) + 2*nrm;
	    a[11] -= 4*(z0-zr[s])*(d0-dist[s]);

	    a[15] += 4*(d0-dist[s])*(d0-dist[s]) -2*nrm;

	    //	    nrm=(x0-xr[s])*(x0-xr[s])+(y0-yr[s])*(y0-yr[s])+
	    //	       (z0-zr[s])*(z0-zr[s])-(d0-dist[s])*(d0-dist[s]);
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
      F77_FUNC(linsolvelu,LINSOLVEU)( &four, a, b, dx );
      err = 0;
      for( int i=0 ; i < 4 ; i++ )
	 err = err > abs(dx[i]) ? err : abs(dx[i]);
      if( it == 0 )
	 err0 = err;
      err /= err0;
      if( myRank == 0 )
      {
	 cout <<  "Newton iteration " << it << " " << err << endl;
         cout << "x0 y0 z0 d0 " << x0 << " " << y0 << " " << z0 << " " << d0 << endl;
      }
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
	 cout << "Error: No convergence in Newton iteration for initial guess\n";
   }
   // Both z0 and -z0 can be roots, since all stations have z=0.
   if( z0 < 0 )
      z0 = -z0;
   t0 = d0/cp;
   if( myRank == 0 )
   {
      cout << "Guessed position of source " << x0 << " " << y0 << " " << z0 << endl;
      //      cout << "Guessed time offset " << t0 << "..will not be used." << endl;
   }
}

//-----------------------------------------------------------------------
void guess_source_t0freq( EW &  simulation, vector<Source*>& sources,
			  double tstart, double& t0, double& freq, int myRank )
{
   //-----------------------------------------------------------------------
   // Guess t0 and freq in the source. Frequency is taken from the input source,
   // need to come up with a method for estimating it. t0 is the sum of the 
   // event time and the rise time, should be modified to be different for
   // different source functions, e.g., for a Gaussian:
   // t0 = t0 + sqrt(-2*log(sqrt(2*pi)*1e-6/freq))/freq;
   //
   // 
   // Input: simulation  - Simulation object
   //        sources     - Initial guess source.
   //        tstart      - Estimated start time of event
   //        myRank      - This processor's ID in the MPI communicator.
   // Output: t0   - Estimated t0.
   //         freq - Estimated frequency
   //-----------------------------------------------------------------------

   freq= sources[0]->getFrequency();

   // Add Approximate rise time to start of event.
   t0 = 5/freq + tstart;
   if( myRank == 0 )
   {
      cout << "Guessed t0 and freq " << t0 << " and " << freq << endl;
   }
}

//-----------------------------------------------------------------------
void guess_source_moments( EW &  simulation, vector<Source*>& sources, vector<TimeSeries*>& timeseries,
			   vector<TimeSeries*>& observations, double* xv, int myRank )
{
   //-----------------------------------------------------------------------
   // Compute moment components by least square minimization of the full waveform misfit.
   // The solution of this quadratic minimization problem, is found by direct solution
   // of the normal equations.
   // 
   // Input: simulation  - Simulation object
   //        sources     - Initial guess source, everything except the moment components will be used.
   //        timeseries  - Used as template for creating simulation output time series.
   //        observations- Observed data.
   //        myRank      - This processor's ID in the MPI communicator.
   // Output: sources    - The computed moment components will be returned in the input source.
   //         xv         - The computed moment components will also be inserted into the parameter vector, xv.
   //-----------------------------------------------------------------------

   bool debug=false;
   timeDep tfunc = sources[0]->getTfunc();
   double x0 = sources[0]->getX0();
   double y0 = sources[0]->getY0();
   double z0 = sources[0]->getZ0();
   double t0 = sources[0]->getOffset();
   double freq= sources[0]->getFrequency();
   double mxx, mxy, mxz, myy, myz, mzz;

   if( sources[0]->isMomentSource() )
   {
      int n = observations.size();
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

      if( myRank == 0 )
	cout << endl << "*** Estimating the moment tensor..." << endl << "Solving forward problem for Mxx=1" << endl;
      Source* onesrc = sources[0]->copy("xx");
      onesrc->setMoments(1,0,0,0,0,0);

      vector<Source*> src(1);
      src[0] = onesrc;

      //      simulation.preprocessSources( src );  // AP: preprocessSources has already been called for sources[0]
      simulation.solve( src, tsxx );
      delete onesrc;

      if( myRank == 0 )
	cout << "Solving forward problem for Mxy=1" << endl;
      onesrc = sources[0]->copy("xy");
      onesrc->setMoments(0,1,0,0,0,0);
      src[0] = onesrc;

      //      simulation.preprocessSources( src );// AP: preprocessSources has already been called for sources[0]
      simulation.solve( src, tsxy );
      delete onesrc;

      if( myRank == 0 )
	cout << "Solving forward problem for Mxz=1" << endl;
      onesrc = sources[0]->copy("xz");
      onesrc->setMoments(0,0,1,0,0,0);
      src[0] = onesrc;

      //      simulation.preprocessSources( src );// AP: preprocessSources has already been called for sources[0]
      simulation.solve( src, tsxz );
      delete onesrc;

      if( myRank == 0 )
	cout << "Solving forward problem for Myy=1" << endl;
      onesrc = sources[0]->copy("yy");
      onesrc->setMoments(0,0,0,1,0,0);
      src[0] = onesrc;

      //      simulation.preprocessSources( src );// AP: preprocessSources has already been called for sources[0]
      simulation.solve( src, tsyy );
      delete onesrc;

      if( myRank == 0 )
	cout << "Solving forward problem for Myz=1" << endl;
      onesrc = sources[0]->copy("yz");
      onesrc->setMoments(0,0,0,0,1,0);
      src[0] = onesrc;

      //      simulation.preprocessSources( src );// AP: preprocessSources has already been called for sources[0]
      simulation.solve( src, tsyz );
      delete onesrc;

      if( myRank == 0 )
	cout << "Solving forward problem for Mzz=1" << endl;
      onesrc = sources[0]->copy("zz");
      onesrc->setMoments(0,0,0,0,0,1);
      src[0] = onesrc;

      //      simulation.preprocessSources( src );// AP: preprocessSources has already been called for sources[0]
      simulation.solve( src, tszz );
      delete onesrc;

      double a[36], b[6], x[6];
      for( int i=0 ; i < 36 ; i++ )
	 a[i] = 0;
      for( int i=0 ; i < 6 ; i++ )
	 b[i] = 0;

      for( int s=0 ; s < n ; s++ )
      {
	 TimeSeries* tsobs = tsxx[s]->copy(&simulation, "tsobscopy");
	 tsobs->interpolate( *observations[s] );

 // Assemble the matrix
	 if( tsxx[s]->myPoint() )
	 {
#define amat(i,j) a[i-1+6*(j-1)]
            amat(1,1) += tsxx[s]->product_wgh( *tsxx[s] );
            amat(2,2) += tsxy[s]->product_wgh( *tsxy[s] );
            amat(3,3) += tsxz[s]->product_wgh( *tsxz[s] );
            amat(4,4) += tsyy[s]->product_wgh( *tsyy[s] );
            amat(5,5) += tsyz[s]->product_wgh( *tsyz[s] );
            amat(6,6) += tszz[s]->product_wgh( *tszz[s] );
            amat(1,2) += tsxx[s]->product_wgh( *tsxy[s] );
	    amat(1,3) += tsxx[s]->product_wgh( *tsxz[s] );
	    amat(1,4) += tsxx[s]->product_wgh( *tsyy[s] );
	    amat(1,5) += tsxx[s]->product_wgh( *tsyz[s] );
	    amat(1,6) += tsxx[s]->product_wgh( *tszz[s] );
	    amat(2,3) += tsxy[s]->product_wgh( *tsxz[s] );
	    amat(2,4) += tsxy[s]->product_wgh( *tsyy[s] );
	    amat(2,5) += tsxy[s]->product_wgh( *tsyz[s] );
	    amat(2,6) += tsxy[s]->product_wgh( *tszz[s] );
	    amat(3,4) += tsxz[s]->product_wgh( *tsyy[s] );
	    amat(3,5) += tsxz[s]->product_wgh( *tsyz[s] );
	    amat(3,6) += tsxz[s]->product_wgh( *tszz[s] );
	    amat(4,5) += tsyy[s]->product_wgh( *tsyz[s] );
	    amat(4,6) += tsyy[s]->product_wgh( *tszz[s] );
	    amat(5,6) += tsyz[s]->product_wgh( *tszz[s] );
            b[0] += tsxx[s]->product_wgh( *tsobs );
            b[1] += tsxy[s]->product_wgh( *tsobs );
            b[2] += tsxz[s]->product_wgh( *tsobs );
            b[3] += tsyy[s]->product_wgh( *tsobs );
            b[4] += tsyz[s]->product_wgh( *tsobs );
            b[5] += tszz[s]->product_wgh( *tsobs );
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

      if( myRank == 0 && debug )
      {
	 cout << "Tensor system matrix : " << endl;
	 for( int i=1 ; i <= 6 ; i++ )
	    cout << " " << amat(i,1)<<" " << amat(i,2) <<" " << amat(i,3) <<" " << amat(i,4) <<
	                     " " << amat(i,5) <<" " << amat(i,6) << endl;
	 cout << " ... and right hand side : " << endl;
	 for( int i=1 ; i <= 6 ; i++ )
	    cout << " " << b[i-1] << endl;
      }
      int six=6;
      F77_FUNC(linsolvelu,LINSOLVEU)( &six, a, b, x );
      mxx = x[0];
      mxy = x[1];
      mxz = x[2];
      myy = x[3];
      myz = x[4];
      mzz = x[5];
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
   Source* sguess = sources[0]->copy("srcguess");
   sguess->setMoments(mxx,mxy,mxz,myy,myz,mzz);
   xv[3] = mxx;
   xv[4] = mxy;
   xv[5] = mxz;
   xv[6] = myy;
   xv[7] = myz;
   xv[8] = mzz;

   delete sources[0];
   sources[0] = sguess;
}

//-----------------------------------------------------------------------
void guess_shift( int n, double xs[n], vector<TimeSeries*>& timeseries,
		  vector<TimeSeries*>& observations, double dt )
{
   bool l2misfit = false;
   for( int m = 0 ; m < observations.size() ; m++ )
   {
      int nst = 0;
      nst = timeseries[m]->getNsteps();
      int nsttmp = nst;
      MPI_Allreduce( &nsttmp, &nst, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
      int smin, smax;
      if( nst % 2 == 0 )
	 smin = -nst/2;
      else
	 smin = -(nst+1)/2;
      smax = smin+nst-1;
      double mfmin=1e38;      
      double shift = 0;
      for( int s=smin ; s <= smax ; s++ )
      {
	 observations[m]->add_shift(s*dt);
	 double dshift, ddshift, dd1shift;
	 double mf=0;
         if( l2misfit )
	    mf = timeseries[m]->misfit( *observations[m], NULL, dshift, ddshift, dd1shift );
	 else
	    mf = timeseries[m]->misfit2( *observations[m] );
         double mftmp = mf;
	 MPI_Allreduce( &mftmp, &mf, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	 if( mf < mfmin )
	 {
	    mfmin = mf;
	    shift = s*dt;
	 }
	 observations[m]->add_shift(-s*dt);
      }
      xs[11+m] = shift;
   }
}
