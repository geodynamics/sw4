#include "sw4.h"

void get_data( double x, double y, double z, double& u, double& v, double& w, 
               double& mu, double& lambda, double& rho );

void fg(double x,double y,double z, double eqs[3] );

void generate_ghgrid( int ib, int ie, int jb, int je, int kb, int ke,
		      double h, double* x, double* y, double* z, double topo_zmax );

void rhs4sgcurv( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
		 float_sw4* a_u, float_sw4* a_mu, float_sw4* a_lambda, float_sw4* a_met,
		 float_sw4* a_jac, float_sw4* a_lu, int* onesided, float_sw4* a_acof,
		 float_sw4* a_bope, float_sw4* a_ghcof, float_sw4* a_strx, float_sw4* a_stry );

void rhs4sgcurv_rev( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
		 float_sw4* a_u, float_sw4* a_mu, float_sw4* a_lambda, float_sw4* a_met,
		 float_sw4* a_jac, float_sw4* a_lu, int* onesided, float_sw4* a_acof,
		 float_sw4* a_bope, float_sw4* a_ghcof, float_sw4* a_strx, float_sw4* a_stry );

int metric( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
	    float_sw4* x, float_sw4* y, float_sw4* z, float_sw4* met, float_sw4* jac );
int metric_rev( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
	    float_sw4* x, float_sw4* y, float_sw4* z, float_sw4* met, float_sw4* jac );

void addsgd4fort( int ifirst, int ilast, int jfirst, int jlast,
		  int kfirst, int klast,
		  float_sw4* a_up, float_sw4* a_u, float_sw4* a_um, float_sw4* a_rho,
		  float_sw4* a_dcx,  float_sw4* a_dcy,  float_sw4* a_dcz,
		  float_sw4* a_strx, float_sw4* a_stry, float_sw4* a_strz,
		  float_sw4* a_cox,  float_sw4* a_coy,  float_sw4* a_coz,
		  float_sw4 beta );
void addsgd4fort_unrl( int ifirst, int ilast, int jfirst, int jlast,
		  int kfirst, int klast,
		  float_sw4* a_up, float_sw4* a_u, float_sw4* a_um, float_sw4* a_rho,
		  float_sw4* a_dcx,  float_sw4* a_dcy,  float_sw4* a_dcz,
		  float_sw4* a_strx, float_sw4* a_stry, float_sw4* a_strz,
		  float_sw4* a_cox,  float_sw4* a_coy,  float_sw4* a_coz,
		  float_sw4 beta );
void addsgd4fort_rev( int ifirst, int ilast, int jfirst, int jlast,
		  int kfirst, int klast,
		  float_sw4* a_up, float_sw4* a_u, float_sw4* a_um, float_sw4* a_rho,
		  float_sw4* a_dcx,  float_sw4* a_dcy,  float_sw4* a_dcz,
		  float_sw4* a_strx, float_sw4* a_stry, float_sw4* a_strz,
		  float_sw4* a_cox,  float_sw4* a_coy,  float_sw4* a_coz,
		  float_sw4 beta );

extern "C" {
   void rhs4th3fortsgstr_( int*, int*, int*, int*, int*, int*, int*, int*,
			   double*, double*, double*, double*, double*, double*, double*,
			   double*, double*, double*, double*, char* );
   void rhs4sg( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast, int nk,
	     int* onesided, double* acof, double* bope, double* ghcof, float_sw4* a_lu,
	     float_sw4* a_u, float_sw4* a_mu, float_sw4* a_lambda, float_sw4 h,
		float_sw4* a_strx, float_sw4* a_stry, float_sw4* a_strz );
		//	     float_sw4* lu1, float_sw4* lu2, float_sw4* lu3 );
   void rhs4sg_rev( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast, int nk,
	     int* onesided, double* acof, double* bope, double* ghcof, float_sw4* a_lu,
	     float_sw4* a_u, float_sw4* a_mu, float_sw4* a_lambda, float_sw4 h,
	     float_sw4* a_strx, float_sw4* a_stry, float_sw4* a_strz  );
   double gettimec_();
   //   double gettimec2();
   void varcoeffs4_( double*, double* );
   void wavepropbop_4_( double*, double*, double*, double*, double*, double*, double* );
   void bopext4th_( double*, double* );
   void curvilinear4sg_( int*, int*, int*, int*, int*, int*, double*, double*, double*,
			 double*, double*, double*, int*, double*, double*, double*,
			 double*, double*, char* );
   void metric_( int*, int*, int*, int*, int*, int*, double*, double*, double*, double*,
		 double*, int* );
   void addsgd4_( double*, double*, double*, double*, double*, double*, double*, double*, double*,
		  double*, double*, double*, double*, double*, double*,
		  int*, int*, int*, int*, int*, int*, double* );
}

#include <cstring>
#include <cstdlib>
#include <cmath>

#include <fstream>
#include <iostream>

using namespace std;

int main( int argc, char** argv )
{
   int nb = 2;
   //   int ib = 5, jb=3, kb=-3;
   //   int ib=0, jb=0, kb=0;
   int ib=-1, jb=-1, kb=-1;
   int ni=25, nj=25, nk=25;
   int a=1;
   int onesided[6]={0,0,0,0,0,0};
   bool append=true;
   int cartesian = 1;
   bool supergrid = false;
   float_sw4 dt =0.01;

   while( a < argc )
   {
      if( strcmp(argv[a],"-n")==0 )
      {
	 ni = atoi(argv[a+1]);
	 nj = ni;
	 nk = ni;
	 a += 2;
      }
      else if( strcmp(argv[a],"-ni")==0 )
      {
	 ni = atoi(argv[a+1]);
	 a += 2;
      }
      else if( strcmp(argv[a],"-nj")==0 )
      {
	 nj = atoi(argv[a+1]);
	 a += 2;
      }
      else if( strcmp(argv[a],"-nk")==0 )
      {
	 nk = atoi(argv[a+1]);
	 a += 2;
      }
      else if( strcmp(argv[a],"-osu")== 0 )
      {
	 onesided[4] = 1;
	 a++;
      }
      else if( strcmp(argv[a],"-osl")== 0 )
      {
	 onesided[5] = 1;
	 a++;
      }
      else if( strcmp(argv[a],"-newfile")== 0 )
      {
	 append = false;
	 a++;
      }
      else if( strcmp(argv[a],"-curvilinear")== 0 )
      {
	 cartesian = 0;
	 a++;
      }
      else if( strcmp(argv[a],"-supergrid")== 0 )
      {
	 supergrid = true;
	 cartesian = 1;
	 a++;
      }
      else
      {
	 cout << "Unknown option " << argv[a] << endl;
	 a++;
      }
   }

   double h = 1.0/(ni-1);
   int ie = ni+ib-1, je=nj+jb-1, ke=nk+kb-1;
  
   double* u      = new double[ni*nj*nk*3];
   double* urev   = new double[ni*nj*nk*3];
   double* eqs    = new double[ni*nj*nk*3];
   double* lu     = new double[ni*nj*nk*3];
   double* lurev  = new double[ni*nj*nk*3];
   double* lu2    = new double[ni*nj*nk*3];
   double* mu     = new double[ni*nj*nk];
   double* lambda = new double[ni*nj*nk];
   double* rho    = new double[ni*nj*nk];
   double* strx   = new double[ni];
   double* stry   = new double[nj];
   double* strz   = new double[nk];
   double* met    = new double[ni*nj*nk*4];
   double* metrev = new double[ni*nj*nk*4];
   double* jac    = new double[ni*nj*nk];

   double* dcx   = new double[ni];
   double* dcy   = new double[nj];
   double* dcz   = new double[nk];
   double* cox   = new double[ni];
   double* coy   = new double[nj];
   double* coz   = new double[nk];
   //   double* lucmp1 = new double[ni*nj*nk];
   //   double* lucmp2 = new double[ni*nj*nk];
   //   double* lucmp3 = new double[ni*nj*nk];

 // Populate the arrays
   size_t npts = ni*nj*nk;
   for( int k=0 ; k< nk ; k++ )
      for( int j=0 ; j< nj ; j++ )
	 for( int i=0 ; i< ni ; i++ )
	 {
	    int ind = i+ni*j+ni*nj*k;
	    get_data( i*h, j*h, k*h, u[3*ind], u[3*ind+1], u[3*ind+2],
		      mu[ind], lambda[ind], rho[ind]);
	    // reversed index on u
	    get_data( i*h, j*h, k*h, urev[ind], urev[ind+npts], urev[ind+2*npts],
		      mu[ind], lambda[ind], rho[ind]);
	    fg( i*h, j*h, k*h, &eqs[3*ind] );
	    //	    lu[3*ind]  =1e38;
	    //	    lu[3*ind+1]=1e38;
	    //	    lu[3*ind+2]=1e38;
	    //	    lu2[3*ind]  =1e30;
	    //	    lu2[3*ind+1]=1e30;
	    //	    lu2[3*ind+2]=1e30;
	    //	    lurev[ind]  =1e25;
	    //	    lurev[ind+npts]=1e25;
	    //	    lurev[ind+2*npts]=1e25;
	 }
// Supergrid info
   double beta = 0.0001;
   for( int i=0 ; i< ni ; i++ )
   {
      strx[i] = 1;
      dcx[i]= 0.02;
      cox[i]= 0.01;
   }
   for( int j=0 ; j< nj ; j++ )
   {
      stry[j] = 1;
      dcy[j]= 0.02;
      coy[j]= 0.01;
   }
   for( int k=0 ; k< nk ; k++ )
   {
      strz[k] = 1;
      dcz[k]= 0.02;
      coz[k]= 0.01;
   }

 // SBP boundary operator definition
   double acof[384], bope[48], ghcof[6], bop[24];
   double m_iop[5], m_iop2[5], m_bop2[24], gh2, m_hnorm[4], m_sbop[5];
   varcoeffs4_( acof, ghcof );
   wavepropbop_4_(m_iop, m_iop2, bop, m_bop2, &gh2, m_hnorm, m_sbop );
   bopext4th_( bop, bope );

   if( !cartesian )
   {
 // setup grid, metric mapping and jacobian
      double* x = new double[ni*nj*nk];
      double* y = new double[ni*nj*nk];
      double* z = new double[ni*nj*nk];
      generate_ghgrid( ib, ie, jb, je, kb, ke, h, x, y, z, 1.0 );
      int ierr;
      metric_( &ib, &ie, &jb, &je, &kb, &ke, x, y, z, met, jac, &ierr );
      if( ierr != 0 )
      {
	 cout << "ERROR, metric returned ierr= " << ierr << endl;
	 exit(-1);
      }
      ierr = metric_rev( ib, ie, jb, je, kb, ke, x, y, z, metrev, jac );
      if( ierr != 0 )
      {
	 cout << "ERROR, metric_rev returned ierr= " << ierr << endl;
	 exit(-1);
      }
   }

   int nkupbndry = ke-2;


// Run and time C and fortran routines
   double tc[10], tf[10], tcr[10];
   int nopsint, nopsbnd;
   if( supergrid )
   {
      for( size_t i=0 ; i < 3*ni*nj*nk ; i++ )
      {
	 lu[i] = u[i]+0.2;
	 lu2[i] = u[i]-0.3;
      }
      nopsint = 375;
      nopsbnd = 0;
      for( int s=0 ; s < 10 ; s++ )
      {
	 tc[s] = gettimec_();
	 addsgd4fort( ib, ie, jb, je, kb, ke, u, lu, lu2, rho, dcx, dcy, dcz,
		      strx, stry, strz, cox, coy, coz, beta ); 
	 tc[s] = gettimec_()-tc[s];

	 tcr[s] = gettimec_();
	 addsgd4fort_rev( ib, ie, jb, je, kb, ke, urev, lu, lu2, rho, dcx, dcy, dcz,
		      strx, stry, strz, cox, coy, coz, beta ); 
	 tcr[s] = gettimec_()-tcr[s];
	 tf[s] = gettimec_();
	 //	 double t1 = gettimec2();
	 addsgd4_( &dt, &h, u, lu, lu2, rho, dcx, dcy, dcz,
		   strx, stry, strz, cox, coy, coz,
		   &ib, &ie, &jb, &je, &kb, &ke, &beta );  
	 //	 double t2=gettimec2();
	 //	 cout << t1 << " " << t2 << " diff = " << t2-t1 << endl;
	 tf[s] = gettimec_() - tf[s];

      }
   }
   else
   {
      if( cartesian )
      {
	 nopsint = 666;
	 nopsbnd = 1244;
	 for( int s=0 ; s < 10 ; s++ )
	 {
	    tc[s] = gettimec_();
	    rhs4sg( ib, ie, jb, je, kb, ke, nkupbndry, onesided, acof, bope, ghcof,
		    lu, u, mu, lambda, h, strx, stry, strz );
      //	      lu, u, mu, lambda, h, strx, stry, strz, lucmp1, lucmp2, lucmp3 );
	    tc[s] = gettimec_()-tc[s];
      //      for( int i=0 ; i < npts ; i++ )
      //      {
      //	 lu[3*i]=lucmp1[i];
      //	 lu[3*i+1]=lucmp2[i];
      //	 lu[3*i+2]=lucmp3[i];
      //      }
	    tcr[s] = gettimec_();
	    rhs4sg_rev( ib, ie, jb, je, kb, ke, nkupbndry, onesided, acof, bope, ghcof,
			lurev, urev, mu, lambda, h, strx, stry, strz );
	    tcr[s] = gettimec_()-tcr[s];
	    char op='=';
	    tf[s] = gettimec_();
	    rhs4th3fortsgstr_( &ib, &ie, &jb, &je, &kb, &ke, 
			       &nkupbndry, onesided, acof, bope, ghcof,
			       lu2, u, mu, lambda, &h, strx, stry, strz, &op );
	    tf[s] = gettimec_() - tf[s];
	 }
      }
      else
      {
	 nopsint = 2126;
	 nopsbnd = 6049;
	 for( int s=0 ; s < 10 ; s++ )
	 {
	    tc[s] = gettimec_();
	    rhs4sgcurv( ib, ie, jb, je, kb, ke, u, mu, lambda, met, jac, lu,
			onesided, acof, bope, ghcof, strx, stry );
	    tc[s] = gettimec_()-tc[s];
	    tcr[s] = gettimec_();
	    rhs4sgcurv_rev( ib, ie, jb, je, kb, ke, urev, mu, lambda, metrev, jac, lurev,
			    onesided, acof, bope, ghcof, strx, stry );
	    tcr[s] = gettimec_()-tcr[s];
	    char op='=';
	    tf[s] = gettimec_();
	    curvilinear4sg_( &ib, &ie, &jb, &je, &kb, &ke, 
			     u, mu, lambda, met, jac, lu2, onesided,
			     acof, bope, ghcof, strx, stry, &op );
	    tf[s] = gettimec_() - tf[s];
	 }
      }
 // Compute approximation error and make a consistency check that both routines give the same result:
      double er[3]={0,0,0};
      for( int k=nb ; k< nk-nb ; k++ )
	 for( int j=nb ; j< nj-nb ; j++ )
	    for( int i=nb ; i< ni-nb ; i++ )
	    {
	       int ind = i+ni*j+ni*nj*k;
	       for( int m= 0 ; m<3 ; m++ )
	       {
		  double err = eqs[3*ind+m]*rho[ind]-lu[3*ind+m];
		  if( fabs(err) > er[m] )
		     er[m] = fabs(err);
		  if( fabs(lu[3*ind+m]-lu2[3*ind+m])> 1e-7 )
		  {
		     cout << " component " << m << " at i= " << i << " j= " << j << " k= " << k << " : " << endl;
		     cout << "         lu(fortran) = " << lu2[3*ind+m] << " lu(C) = "  << lu[3*ind+m] << endl;
		  }
		  if( fabs(lurev[ind+m*npts]-lu2[3*ind+m])> 1e-7 )
		  {
		     cout << " component " << m << " at i= " << i << " j= " << j << " k= " << k << " : " << endl;
		     cout << "         lu(fortran) = " << lu2[3*ind+m] << " lu(Crev) = "  << lurev[ind+m*npts] << endl;
		  }
	       }
	    }
      cout << "Errors " << er[0] << " " << er[1] << " " << er[2] << endl;
   }
// Output timing results:
   cout <<  "Time C code and fortran code :" << endl;
   double tcavg=0, tfavg=0, travg=0;
   for( int s=0 ; s < 10 ; s++ )
   {
      cout << "C " << tc[s] << " sec. , C(rev)  " <<  tcr[s] << " sec., Fortran   " << tf[s] << " sec. " << endl;
      if( s > 1 )
      {
	 tcavg += tc[s];
	 tfavg += tf[s];
	 travg += tcr[s];
      }
   }
   tcavg /= 8;
   tfavg /= 8;
   travg /= 8;
   int totpts = (ni-4)*(nj-4)*(nk-4);
   double cmflop, fmflop, crflop;
   cmflop = nopsint*(totpts/tcavg)/1e9;
   crflop = nopsint*(totpts/travg)/1e9;
   fmflop = nopsint*(totpts/tfavg)/1e9;
   cout << " C code, average (8 last)       " << tcavg << " sec = " << cmflop << " Gflops" << endl;
   cout << " C code (rev), average (8 last) " << travg << " sec = " << crflop << " Gflops" << endl;
   cout << " Fortran code, average (8 last) " << tfavg << " sec = " << fmflop << " Gflops" << endl;

// Save timing numbers to file for plotting purpose
   ofstream* utfil;
   if( append )
      utfil = new ofstream("result.dat",ofstream::app);
   else
      utfil = new ofstream("result.dat",ofstream::trunc);
   *utfil << ni << " " << nj << " " << nk << " " << tcavg << " " << tfavg << " " << travg << " " << cmflop << " "
	  << fmflop << " " << crflop << endl;

   utfil->close();
}

void get_data( double x, double y, double z, double& u, double& v, double& w, 
               double& mu, double& lambda, double& rho )
{
lambda = cos(x)*pow(sin(3*y),2)*cos(z);
mu     = sin(3*x)*sin(y)*sin(z);
rho    = x*x*x+1+y*y+z*z;
u      = cos(x*x)*sin(y*x)*z*z;
v      = sin(x)*cos(y*y)*sin(z);
w      = cos(x*y)*sin(z*y);
}


void fg(double x,double y,double z, double eqs[3] )
{
  double t1;
  double t10;
  double t109;
  double t11;
  double t112;
  double t113;
  double t117;
  double t119;
  double t120;
  double t126;
  double t127;
  double t13;
  double t135;
  double t14;
  double t148;
  double t149;
  double t150;
  double t151;
  double t152;
  double t159;
  double t16;
  double t17;
  double t173;
  double t175;
  double t18;
  double t187;
  double t19;
  double t190;
  double t2;
  double t20;
  double t21;
  double t22;
  double t24;
  double t25;
  double t26;
  double t27;
  double t29;
  double t3;
  double t32;
  double t33;
  double t36;
  double t37;
  double t38;
  double t39;
  double t4;
  double t41;
  double t5;
  double t50;
  double t51;
  double t56;
  double t59;
  double t61;
  double t64;
  double t69;
  double t70;
  double t71;
  double t72;
  double t75;
  double t77;
  double t78;
  double t79;
  double t8;
  double t85;
  double t88;
  double t9;
  double t95;
  double t99;
  {
    t1 = 3.0*x;
    t2 = cos(t1);
    t3 = sin(y);
    t4 = t2*t3;
    t5 = sin(z);
    t8 = sin(x);
    t9 = 3.0*y;
    t10 = sin(t9);
    t11 = t10*t10;
    t13 = cos(z);
    t14 = t8*t11*t13;
    t16 = x*x;
    t17 = sin(t16);
    t18 = t17*x;
    t19 = y*x;
    t20 = sin(t19);
    t21 = z*z;
    t22 = t20*t21;
    t24 = 2.0*t18*t22;
    t25 = cos(t16);
    t26 = cos(t19);
    t27 = t25*t26;
    t29 = t27*y*t21;
    t32 = sin(t1);
    t33 = t32*t3;
    t36 = cos(x);
    t37 = t36*t11;
    t38 = t37*t13;
    t39 = 2.0*t33*t5+t38;
    t41 = t25*t16*t22;
    t50 = t25*t20;
    t51 = y*y;
    t56 = t8*t8;
    t59 = sin(t51);
    t61 = t59*y*t5;
    t64 = t36*t36;
    t69 = z*y;
    t70 = cos(t69);
    t71 = t26*t70;
    t72 = t71*y;
    t75 = t20*t51*t70;
    t77 = cos(y);
    t78 = t32*t77;
    t79 = cos(t51);
    t85 = t5*(t36*t79*t5+t27*x*t21);
    t88 = y*t5;
    t95 = sin(t69);
    t99 = -t20*y*t95+2.0*t50*z;
    t109 = 1/(t16*x+1.0+t51+t21);
    eqs[0] = ((6.0*t4*t5-t14)*(-t24+t29)+t39*(-4.0*t41-2.0*t17*t20*t21-4.0*t18*
t26*y*t21-t50*t51*t21)+2.0*t56*t11*t13*t61-2.0*t64*t11*t13*t61-t14*t72-t38*t75+
t78*t85+t33*t5*(-2.0*t36*t59*t88-t41)+t33*t13*t99+t33*t5*(-t75+2.0*t50))*t109;
    t112 = t8*t79;
    t113 = t112*t5;
    t117 = 2.0*t17*t16*t26*t21;
    t119 = t50*t19*t21;
    t120 = t27*t21;
    t126 = t36*t10;
    t127 = cos(t9);
    t135 = t39*t8;
    t148 = t20*x;
    t149 = t70*y;
    t150 = t148*t149;
    t151 = t26*t95;
    t152 = t151*t69;
    t159 = -t148*t95+t71*z+t112*t13;
    eqs[1] = (3.0*t4*t85+t33*t5*(-t113-t117-t119+t120)-2.0*(2.0*t78*t5+6.0*t126
*t13*t127)*t8*t61-4.0*t135*t79*t51*t5-2.0*t135*t59*t5+6.0*t126*t13*(-t24+t29+
t72)*t127+t37*t13*(-t117-t119+t120-t150-t152+t71)+t33*t13*t159+t33*t5*(-t150-
t152+t71-t113))*t109;
    t173 = 4.0*t18*t20*z;
    t175 = 2.0*t27*t69;
    t187 = t8*t59;
    t190 = 2.0*t187*y*t13;
    eqs[2] = (3.0*t4*t5*t99+t33*t5*(-t26*t51*t95-t173+t175)+t78*t5*t159+t33*t5*
(-t26*t16*t95-2.0*t148*t70*z-t151*t21-t190)+(2.0*t33*t13-t37*t5)*t26*t149-t39*
t26*t95*t51-t37*t5*(-t24+t29-2.0*t187*t88)+t37*t13*(-t173+t175-t190))*t109;
    return;
  }
}
