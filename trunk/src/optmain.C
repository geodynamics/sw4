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
      persrc[0] = GlobalSources[0]->copy( &simulation );
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
		    vector<TimeSeries*> & GlobalObservations, int myRank )
{
// Run forward problem with guessed source
   simulation.solve( GlobalSources, GlobalTimeSeries );
// Compute misfit, 'diffs' will hold the source for the adjoint problem
   vector<TimeSeries*> diffs;
   for( int m=0 ; m < GlobalTimeSeries.size() ; m++ )
   {
      TimeSeries *elem = GlobalTimeSeries[m]->copy( &simulation, "diffsrc" );
      diffs.push_back(elem);
   }

   double mf = 0;
   for( int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
      mf += GlobalTimeSeries[m]->misfit( *GlobalObservations[m], diffs[m] );
   double mftmp = mf;
   MPI_Allreduce(&mftmp,&mf,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   if( myRank == 0 )
      cout << "Misfit = " << mf << endl;

   // Test Gradient computation
   // Get gradient by computing the backwards problem:
   double gradient[11], hess[121];
   simulation.solve_backward( GlobalSources, diffs, gradient, hess );
   if( myRank == 0 )
   {
      cout << "Gradient, by adjoint equation = " << endl;
      for( int i = 0 ; i < 11 ; i++ )
	 cout << "   " << gradient[i] << endl;
   }
   for( int testcomp = 0 ; testcomp < 11 ; testcomp++ )
   {
// Validation, compute derivative by forward solve
      GlobalSources[0]->set_derivative(testcomp);
      simulation.solve( GlobalSources, GlobalTimeSeries );
      double dxidp=0;
      for( int m=0 ; m < GlobalTimeSeries.size(); m++ )
	 dxidp += GlobalTimeSeries[m]->product( *diffs[m] );
      cout << "Component " << testcomp+1 << " Forward solve derivative = " << dxidp << endl;
      GlobalSources[0]->set_noderivative();         

// Validation, compute numerical gradient:
      double h = 0.00001;
      vector<Source*> persrc(1);
      persrc[0] = GlobalSources[0]->copy( &simulation );
      persrc[0]->perturb( h, testcomp );
	 //		 GlobalSources[0]->perturb(h,testcomp);
	 //		 simulation.preprocessSources( GlobalSources );
      simulation.solve( persrc, GlobalTimeSeries );
      double mfp = 0;
      for( int m=0 ; m < GlobalTimeSeries.size() ; m++ )
	 mfp += GlobalTimeSeries[m]->misfit( *GlobalObservations[m], NULL );
      mftmp=mfp;
      MPI_Allreduce(&mftmp,&mfp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      if( myRank == 0 )
      {
	 cout << "Misfit of perturbed problem = " << mfp << endl;
	 cout << "Difference = " << mfp-mf << endl;
	 cout << "Component " << testcomp+1 << " Numerical derivative = " << (mfp-mf)/h << endl;
      }
      delete persrc[0];
   }
}

//-----------------------------------------------------------------------
void test_hessian(  EW& simulation, vector<Source*>& GlobalSources,
		    vector<TimeSeries*>& GlobalTimeSeries,
		    vector<TimeSeries*> & GlobalObservations, int myRank )
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

   double mf = 0;
   for( int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
      mf += GlobalTimeSeries[m]->misfit( *GlobalObservations[m], diffs[m] );
   double mftmp = mf;
   MPI_Allreduce(&mftmp,&mf,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   if( myRank == 0 )
      cout << "Misfit = " << mf << endl;


// Get Hessian by solving the backwards problem:
   bool skipthis = false;
   double gradient[11], hess[121];
   simulation.solve_backward( GlobalSources, diffs, gradient, hess );

// Assemble the first part of hessian matrix
   vector<vector<TimeSeries*> > dudp(11);
   if( !skipthis )
   {
      for( int comp = 0 ; comp < 11 ; comp++ )
      {
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
      }
   }
   // Validate by numerical differentiation
   double h = 0.0001;
   vector<Source*> persrc(1);
   GlobalSources[0]->set_noderivative();
   for( int testcomp = 0 ; testcomp < 11 ; testcomp++ )
   {
      persrc[0] = GlobalSources[0]->copy( &simulation );
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
      for( int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
	 GlobalTimeSeries[m]->misfit( *GlobalObservations[m], diffs[m] );

		 //		 GlobalSources[0]->perturb(h,testcomp);
		 //		 simulation.preprocessSources( GlobalSources );
      double gradp[11];
      simulation.solve_backward( persrc, diffs, gradp, hess );
      if( myRank == 0 )
      {
	 cout << "Hessian, col. no " << testcomp+1 << " , by numerical derivative = " << endl;
	 for( int i = 0 ; i < 11 ; i++ )
	    cout << "   " << (gradp[i]-gradient[i])/h << endl;
	 
      }
      delete persrc[0];
   }
}

//-----------------------------------------------------------------------
void compute_scalefactors(  EW& simulation, vector<Source*>& GlobalSources,
			    vector<TimeSeries*>& GlobalTimeSeries,
			    vector<TimeSeries*> & GlobalObservations, int myRank,
			    double sf[11] )
{
   // Compute Hessian for input source:

   // Run forward problem 
   simulation.solve( GlobalSources, GlobalTimeSeries );
   // Compute misfit, 'diffs' will hold the source for the adjoint problem
   vector<TimeSeries*> diffs;
   for( int m=0 ; m < GlobalTimeSeries.size() ; m++ )
   {
      TimeSeries *elem = GlobalTimeSeries[m]->copy( &simulation, "diffsrc" );
      diffs.push_back(elem);
   }

   double mf = 0;
   for( int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
      mf += GlobalTimeSeries[m]->misfit( *GlobalObservations[m], diffs[m] );
   //   double mftmp = mf;
   //   MPI_Allreduce(&mftmp,&mf,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   //   if( myRank == 0 )
   //      cout << "Misfit = " << mf << endl;


// Get second part of Hessian by solving the backwards problem:
   double gradient[11], hess[121];
   simulation.solve_backward( GlobalSources, diffs, gradient, hess );

// Assemble the first part of hessian matrix by solving 11 forward problems
   vector<vector<TimeSeries*> > dudp(11);
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

   // Fill in symmetric part
   for( int comp = 0 ; comp < 11 ; comp++ )
      for(int comp1 = 0 ; comp1 < comp; comp1++ )
	 hess1[comp+11*comp1] = hess1[comp1+11*comp];

   // Add to second part of Hessian
   for( int comp = 0 ; comp < 11 ; comp++ )
      for(int comp1 = 0 ; comp1 < 11; comp1++ )
	 hess[comp+11*comp1] += hess1[comp+11*comp1];


   // Check if diagonal elements are positive
   bool diagpositive = true;
   for( int comp = 0 ; comp < 11 ; comp++ )
      if( hess[comp+11*comp] <= 0 )
	 diagpositive = false;

   if( !diagpositive )
   {
      if( myRank == 0 )
      {
	 cout << "Hessian have negative diagonal elements, " << endl;
         cout << " scaling factors will use only the first part of the Hessian" << endl;
      }
      for( int comp= 0 ; comp < 11 ; comp++ )
	 sf[comp]=1/sqrt(hess1[comp+11*comp]);
   }
   else
   {
      for( int comp= 0 ; comp < 11 ; comp++ )
	 sf[comp]=1/sqrt(hess[comp+11*comp]);
   }

// diffs no longer needed, give back memory
   for( unsigned int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
      delete diffs[m];
   diffs.clear();
   
// dudp no longer needed, give back memory
   for( unsigned int comp = 0 ; comp < 11 ; comp++ )
   {
      for( unsigned int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
	 delete dudp[comp][m];
      dudp[comp].clear();
   }
   dudp.clear();
}
void get_testmatrix( double a[121], double b[11] )
{
   // Matrix
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

//-----------------------------------------------------------------------
void compute_f( EW& simulation, double x[11], vector<Source*>& GlobalSources,
		vector<TimeSeries*>& GlobalTimeSeries,
		vector<TimeSeries*>& GlobalObservations, double& mf )
// Compute misfit
{
   bool testing = false;
   if( testing )
   {
      double a[121],b[11];
      get_testmatrix( a, b );
      double xax=0, bx=0;
      
      for( int j=0 ; j<11 ; j++ )
      {
	 for( int i=0 ; i < 11 ; i++ )
	    xax += x[i]*a[i+11*j]*x[j];
         bx += x[j]*b[j];
      }
      mf = 0.5*xax-bx;
      return;
   }
   vector<Source*> src(1);
   src[0] = GlobalSources[0]->copy( &simulation );
   src[0]->set_parameters( x );

// Run forward problem with guessed source
   simulation.solve( src, GlobalTimeSeries );
// Compute misfit,
   mf = 0;
   for( int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
      mf += GlobalTimeSeries[m]->misfit( *GlobalObservations[m], NULL );
   double mftmp = mf;
   MPI_Allreduce(&mftmp,&mf,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   delete src[0];
}

//-----------------------------------------------------------------------
void compute_fanddf( EW& simulation, double x[11], vector<Source*>& GlobalSources,
		     vector<TimeSeries*>& GlobalTimeSeries,
		     vector<TimeSeries*>& GlobalObservations, double& f, double df[11] )
   // Compute misfit and its gradient
{
   bool testing = false;
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
      return;
   }
   vector<Source*> src(1);
   src[0] = GlobalSources[0]->copy( &simulation );
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
   f = 0;
   for( int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
      f += GlobalTimeSeries[m]->misfit( *GlobalObservations[m], diffs[m] );
   double mftmp = f;
   MPI_Allreduce(&mftmp,&f,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

   // Get gradient by solving the adjoint problem:
   double hess[121];
   simulation.solve_backward( src, diffs, df, hess );

// diffs no longer needed, give back memory
   for( unsigned int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
      delete diffs[m];
   diffs.clear();

   delete src[0];
}

//-----------------------------------------------------------------------
void compute_dtd2fd( EW& simulation, double x[11], vector<Source*>& GlobalSources,
		     vector<TimeSeries*>& GlobalTimeSeries,
		     vector<TimeSeries*>& GlobalObservations, double d[11], double& dtHd )
   // Compute d'*H*d
{
   bool testing = false;
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
   src[0] = GlobalSources[0]->copy( &simulation );
   src[0]->set_parameters( x );

   simulation.solve( src, GlobalTimeSeries );

// Compute sources for the adjoint problem
   vector<TimeSeries*> diffs;
   for( int m=0 ; m < GlobalTimeSeries.size() ; m++ )
   {
      TimeSeries *elem = GlobalTimeSeries[m]->copy( &simulation, "diffsrc" );
      diffs.push_back(elem);
   }
   for( int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
      GlobalTimeSeries[m]->misfit( *GlobalObservations[m], diffs[m] );

   // Forward solve for first part of Hessian
   src[0]->set_dirderivative( d );
   simulation.solve( src, GlobalTimeSeries );
   dtHd = 0;
   for( int m=0 ; m < GlobalTimeSeries.size() ; m++ )
      dtHd += GlobalTimeSeries[m]->product_wgh( *GlobalTimeSeries[m] );
   double dtHdtmp = dtHd;
   MPI_Allreduce(&dtHdtmp,&dtHd,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

   // Backward solve for second part of Hessian
   double df[11], hess[121];
   simulation.solve_backward( src, diffs, df, hess );
   double d2h = 0;
   for( int i=0 ; i < 11 ; i++ )
      for( int j=0 ; j < 11 ; j++ )
	 d2h += hess[i+11*j]*d[i]*d[j];

   dtHd += d2h;
   delete src[0];
}

//-----------------------------------------------------------------------
void linesearch( EW& simulation, vector<Source*>& GlobalSources,
		 vector<TimeSeries*>& GlobalTimeSeries, vector<TimeSeries*>& GlobalObservations,
		 double x[11], double f, double df[11], double p[11],
		 double cgstep, double maxstep, double steptol, double xnew[11],
		 double fnew, double sf[11], int myRank, int& retcode )
{
   bool maxtaken = false;
   double alpha = 1e-6, cglen=0, ang=0;

   retcode = 2;
   for( int i=0 ; i < 11 ; i++ )
   {
      cglen += p[i]*p[i]/(sf[i]*sf[i]);
      ang += p[i]*df[i];
   }
   if( ang >= 0 )
   {
      if( myRank == 0 )
      {
	 cout << "LINESEARCH: Warning, direction is not a descent direction" << endl;
         cout << "   switching to gradient direction search" << endl;
      }
      cglen = 0;
      for( int i=0 ; i < 11 ; i++ )
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
         for( int i=0 ; i < 11 ; i++ )
            p[i] = p[i]*(maxstep/cglen);
	 cglen = maxstep;
      }
      else
      {
         for( int i=0; i < 11 ; i++ )
            p[i] = cgstep*p[i];
         cglen = cgstep*cglen;
      }
      // minlambda is scaled such that x^{k+1}-x^k = lambda*p^k with lambda=minlambda implies
      //            || x^{k+1}-x^k ||_\infty < minlambda (termination criterium)
      double initslope = 0;
      double rellength = 0;
      for( int i=0; i < 11 ; i++ )
      {
         initslope = initslope + df[i]*p[i];
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
      while( retcode == 2 )
      {
	 for( int i=0; i < 11 ; i++ )
	    xnew[i] = x[i] + lambda*p[i];
	 compute_f( simulation, xnew, GlobalSources, GlobalTimeSeries, GlobalObservations, fnew );
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
	    for( int i=0 ; i < 11 ; i++ )
	       ang += p[i]*df[i];
            cout << "LINESEARCH: no satsifactory step found \n";
	    cout << "cg-alpha = " << cgstep << endl;
	    cout << "scprod = " << ang << endl;
	    cout << "search direction = " << endl;
	    for( int i=0 ; i < 11 ; i++ )
	       cout << " " << p[i] << endl;
	    cout << "gradient direction = " << endl;
	    for( int i=0 ; i < 11 ; i++ )
	       cout << " " << df[i] << endl;
	    cout << "maxstep   = " << maxstep << endl;
	    cout << "minlambda = " << minlambda << endl;
	    cout << "initslope = " << initslope << endl;
	    cout << " f = " << f << " fnew = " << fnew <<  " lambda = " << lambda ;
	    cout << " fnew-(f+alpha*lambda*initslope) = " << fnew-(f+alpha*lambda*initslope) << endl;
	    for( int i=0 ; i < 11 ; i++ )
	       cout << " " << x[i] << endl;
	    retcode = 1;
	    // Take a full step
	    for( int i=0 ; i < 11 ; i++ )
	       xnew[i] = x[i] + p[i];
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
	 }
      }
}

//-----------------------------------------------------------------------
void cg( EW& simulation, double x[11], double sf[11], vector<Source*>& GlobalSources,
	 vector<TimeSeries*>& GlobalTimeSeries,
	 vector<TimeSeries*> & GlobalObservations,
	 int myRank, bool fletcher_reeves, int stepselection )
{
   const int n = 11;
   int j, k;
   bool done = false;
   bool dolinesearch = true;
   int maxit, maxrestart;
   double tolerance;
   double f, df[11], d[11], d2f[121], da[11], xa[11], dfp[11], rnorm, errnorm;
   //   cout << "Enter CG " << endl;
   simulation.get_cgparameters( maxit, maxrestart, tolerance );
   if( maxit > n )
      maxit = n;
   //   cout << "maxit = " << maxit << " maxrestart= " << maxrestart << " tolerance= " << tolerance << endl;
   compute_fanddf( simulation, x, GlobalSources, GlobalTimeSeries, GlobalObservations, f, df );
   cout << "Begin iteration misfit= "  << f << endl;
   cout << " gradient = " ;
   for( int i=0 ; i < n ; i++ )
   {
      cout << df[i] << " ";
      if( i==5 )
	 cout << endl << "      " ;
   }
   cout << endl;

   j = 1;
   while( j <= maxrestart && !done )
   {
      for( int i=0 ; i < n ; i++ )
         d[i] = -df[i]*sf[i]*sf[i];
      rnorm = 0;
      for( int i=0 ; i < n ; i++ )
	 rnorm = rnorm > fabs(df[i])*sf[i] ? rnorm : fabs(df[i])*sf[i];
      errnorm = rnorm;
      int k = 1;
      while( k <= maxit && errnorm > tolerance )
      {  
         double den = 0, dnrm=0, xnrm=0;
	 for( int i= 0 ; i < n ; i++ )
	 {
	    den  += d[i]*df[i];
	    dnrm += d[i]*d[i];
	    xnrm += x[i]*x[i];
	 }
	 double alpha, fp;
	 if( stepselection == 0 )
	 {
            double h=-2*f/den;
	    for( int i=0; i<n ; i++ )
	       xa[i] = x[i] + h*d[i];
            compute_f( simulation, xa, GlobalSources, GlobalTimeSeries, GlobalObservations, fp );
	    alpha = h*f/(fp+f);
	 }
	 else
	 {
            double dtHd;
	    compute_dtd2fd( simulation, x, GlobalSources, GlobalTimeSeries, GlobalObservations, d, dtHd );
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
	    linesearch( simulation, GlobalSources, GlobalTimeSeries, GlobalObservations,
			x, f, df, da, fabs(alpha), 1.0, tolerance*0.001, xa, fp, sf, myRank, retcode );
	 }
         else
	    for( int i=0 ; i < n ; i++ )
	       xa[i] = x[i] + alpha*d[i];

	 // xa is now the new iterate
	 double dxnorm = 0;
	 for( int i=0 ; i < n ; i++ )
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

	 compute_fanddf( simulation, x, GlobalSources, GlobalTimeSeries, GlobalObservations, f, dfp );
         rnorm = 0;
	 for( int i=0 ; i < n ; i++ )
	    if( fabs(dfp[i])*sf[i] > rnorm )
	       rnorm = fabs(dfp[i])*sf[i];

         cout << "-----------------------------------------------------------------------" << endl;
         cout << " it=" << j << "," << k << " dfnorm= " << rnorm << " dxnorm= " << dxnorm << endl;
         cout << " new x = " ;
	 for( int i=0 ; i < n ; i++ )
	 {
	    cout << x[i] << " ";
	    if( i==5 )
	       cout << endl << "      " ;
	 }
	 cout << endl;
         cout << " Misfit= "  << f << endl;
         cout << " scaled gradient = " ;
	 for( int i=0 ; i < n ; i++ )
	 {
	    cout << dfp[i]*sf[i] << " ";
	    if( i==5 )
	       cout << endl << "      " ;
	 }
	 cout << endl;

	 errnorm = rnorm;
         double beta;
         if( k < maxit )
	 {
            double num=0, mix=0;
            den = 0;
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
                  if( myRank == 0 )
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
        char str[10];
        snprintf(str,10,"%i",m);
	string newname = "tscopy";
        newname.append(str);
	TimeSeries *elem = GlobalObservations[m]->copy( &simulation, newname );
	GlobalTimeSeries.push_back(elem);
     }

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
     else if( GlobalSources.size() != 1 )
     {
	if (myRank == 0)
	   cout << "Source optmization only implemented for a single source" << endl;
     }
     else
     {
// Successful initialization

	if (myRank == 0)
	{
	   cout << "Running sbp4opt on " <<  nProcs << " processors..." << endl
		<< "Writing output to directory: " 
		<< simulation.getOutputPath() << endl;
	}
	
// Variables needed:
//        bool guess = false; // auto or from input file
	int testcase  = 4; // job to do
	//	int scale_factor = 1; // Method for computing scale factors, could be from input file
	//        int maxit, maxrestarts;
	//	double tolerance;
	double xv[11], sf[11];

// Initial guess
	if( simulation.compute_guess() )
	   guess_source( simulation, GlobalSources, GlobalTimeSeries, GlobalObservations, xv, myRank );
        else
	   GlobalSources[0]->get_parameters(xv);
	
	if( myRank == 0 )
	{
	   cout << "Initial source guess : \n";
	   cout << "   x0 = " << xv[0] << " y0 = " << xv[1] << " z0 = " << xv[2] <<endl;
	   cout << "  mxx = " << xv[3] << " mxy= " << xv[4] << " mxz= " << xv[5] << endl;
	   cout << "  myy = " << xv[6] << " myz= " << xv[7] << " mzz= " << xv[8] << endl;
	   cout << "   t0 = " << xv[9] << " freq = " << xv[10] << endl;
	}

	if( simulation.compute_sf() )
           compute_scalefactors( simulation, GlobalSources, GlobalTimeSeries, GlobalObservations, myRank, sf );
	else
           simulation.get_scalefactors(sf);
	      

	if( myRank == 0 )
	{
	   cout << " Scaling factors : " << endl;
	   for( int comp = 0 ; comp < 11 ; comp++ )
	   {
	      cout << " sf" << comp+1 << " = " << sf[comp] ;
	      cout << endl;
	   }
	}

        if( testcase == 1 )
	   testsourced2( simulation, GlobalSources );
        else if( testcase == 2 )
	   test_gradient( simulation, GlobalSources, GlobalTimeSeries, GlobalObservations, myRank );
	else if( testcase == 3 )
	   test_hessian( simulation, GlobalSources, GlobalTimeSeries, GlobalObservations, myRank );
	else
	{
	   // run optimization
           bool fletcher_reeves = true;
           int stepselection = 0 ; // 0 - C.Tape choice, 1 - uTHu	   
	   cg( simulation, xv, sf, GlobalSources, GlobalTimeSeries,
	       GlobalObservations, myRank, fletcher_reeves, stepselection );
	}
// Save all time series
//	for (int ts=0; ts<GlobalTimeSeries.size(); ts++)
//	{
//	   GlobalTimeSeries[ts]->writeFile();
//	}

	if( myRank == 0 )
	{
	   cout << "============================================================" << endl
		<< " sbp4f ( Summation by parts 4th order inverse wave solver) finished! " << endl
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
      simulation.setupRun( );
      simulation.preprocessSources( src );
      simulation.solve( src, tsxx );
      delete onesrc;

      onesrc = new Source(&simulation, amp, freq, t0, x0, y0, z0, 0, 1, 0, 0, 0, 0, tfunc, "xy");
      src[0] = onesrc;

      simulation.preprocessSources( src );
      simulation.solve( src, tsxy );
      delete onesrc;

      onesrc = new Source(&simulation, amp, freq, t0, x0, y0, z0, 0, 0, 1, 0, 0, 0, tfunc, "xz");
      src[0] = onesrc;

      simulation.preprocessSources( src );
      simulation.solve( src, tsxz );
      delete onesrc;

      onesrc = new Source(&simulation, amp, freq, t0, x0, y0, z0, 0, 0, 0, 1, 0, 0, tfunc, "yy");
      src[0] = onesrc;

      simulation.preprocessSources( src );
      simulation.solve( src, tsyy );
      delete onesrc;

      onesrc = new Source(&simulation, amp, freq, t0, x0, y0, z0, 0, 0, 0, 0, 1, 0, tfunc, "yz");
      src[0] = onesrc;

      simulation.preprocessSources( src );
      simulation.solve( src, tsyz );
      delete onesrc;

      onesrc = new Source(&simulation, amp, freq, t0, x0, y0, z0, 0, 0, 0, 0, 0, 1, tfunc, "zz");
      src[0] = onesrc;

      simulation.preprocessSources( src );
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
