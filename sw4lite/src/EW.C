#include "sw4.h"

#include "EW.h"

#include <sstream>
#include <fstream>

#ifdef SW4_OPENMP
#include <omp.h>
#endif

#include <mpi.h>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <algorithm>

#include "Source.h"
#include "GridPointSource.h"
#include "CheckPoint.h"
#include "MaterialBlock.h"
#include "TimeSeries.h"

#include "F77_FUNC.h"
#include "EWCuda.h"
#include "mynvtx.h"
#include "policies.h"
#include <chrono>
#include <thread>

void PrintPointerAttributes(void *ptr);
#ifdef CUDA_CODE
bool IsManaged(void *ptr);
#endif
#ifndef SW4_CROUTINES
extern "C" {
   void F77_FUNC(rhs4th3fortsgstr,RHS4TH3FORTSGSTR)( int*, int*, int*, int*, int*, int*, int*, int*, 
						     float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
						     float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, char* );

   void F77_FUNC(corrfort,CORRFORT)(int*, int*, int*, int*, int*, int*, 
				 float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4* );

   void F77_FUNC(dpdmtfort,DPDMTFORT)(int*, int* , int*, int*, int*, int*, 
				  float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4* );    

   void F77_FUNC(predfort,PREDFORT)(int*, int*, int*, int*, int*, int*, 
				 float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*);    

   void F77_FUNC(solerr3, SOLERR3)(int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
				   float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
				   float_sw4*, int*, int*, int*, int*, int*, int* );

   void F77_FUNC(bcfortsg, BCFORTSG)( int*, int*, int*, int*, int*, int*, 
				      int *, int*, int*, int*, float_sw4*, float_sw4*, boundaryConditionType*,
				      float_sw4 *, float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, 
				      float_sw4*, float_sw4*, float_sw4*, float_sw4*, 
				      float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4* );

   void F77_FUNC(wavepropbop_4, WAVEPROPBOP_4)(float_sw4 *, float_sw4 *, float_sw4 *, float_sw4 *, float_sw4 *, 
                                               float_sw4 *, float_sw4 *);

   void F77_FUNC(varcoeffs4,VARCOEFFS4)(float_sw4 *, float_sw4 *);

   void F77_FUNC(bopext4th,BOPEXT4TH)(float_sw4 *, float_sw4 *);

   void F77_FUNC(addsgd4,ADDSGD4) (float_sw4* dt, float_sw4 *h, float_sw4 *a_Up, float_sw4*a_U, float_sw4*a_Um, float_sw4* Rho,
				float_sw4 *sg_dc_x, float_sw4* sg_dc_y, float_sw4* sg_dc_z, float_sw4* sg_str_x, float_sw4* sg_str_y, float_sw4* sg_str_z,
				float_sw4* sg_corner_x, float_sw4* sg_corner_y, float_sw4* sg_corner_z,
				int *ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast, float_sw4* damping_coefficient );
    void F77_FUNC(addsgd4c,ADDSGD4C) ( float_sw4* dt, float_sw4* a_Up, float_sw4* a_U, float_sw4* a_Um, 
				       float_sw4* Rho, float_sw4* sg_dc_x, float_sw4* sg_dc_y,
				       float_sw4* sg_str_x, float_sw4* sg_str_y, float_sw4* mJ,
				       float_sw4* sg_corner_x, float_sw4* sg_corner_y,
				       int* ifirst, int* ilast, int* jfirst, int* jlast, int* kfirst, 
				       int* klast, float_sw4* damping_coefficient );

   void F77_FUNC(addsgd6,ADDSGD6) (float_sw4* dt, float_sw4 *h, float_sw4 *a_Up, float_sw4*a_U, float_sw4*a_Um, float_sw4* Rho,
				float_sw4 *sg_dc_x, float_sw4* sg_dc_y, float_sw4* sg_dc_z, float_sw4* sg_str_x, float_sw4* sg_str_y, float_sw4* sg_str_z,
				float_sw4* sg_corner_x, float_sw4* sg_corner_y, float_sw4* sg_corner_z,
				int *ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast, float_sw4* damping_coefficient );

    void F77_FUNC(addsgd6c,ADDSGD6C) ( float_sw4* dt, float_sw4* a_Up, float_sw4* a_U, float_sw4* a_Um, 
				       float_sw4* Rho, float_sw4* sg_dc_x, float_sw4* sg_dc_y,
				       float_sw4* sg_str_x, float_sw4* sg_str_y, float_sw4* mJ,
				       float_sw4* sg_corner_x, float_sw4* sg_corner_y,
				       int* ifirst, int* ilast, int* jfirst, int* jlast, int* kfirst, 
				       int* klast, float_sw4* damping_coefficient );

   void F77_FUNC(gridinfo,GRIDINFO) (int *ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast, 
				     float_sw4* met, float_sw4* jac, float_sw4* minj, float_sw4* maxj );
   void F77_FUNC(metric,METRIC)( int *ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast, 
				 float_sw4* x, float_sw4* y, float_sw4* z, float_sw4* met, float_sw4* jac, int* ierr );
   void F77_FUNC(metricexgh,METRICEXGH)( int *ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast, 
					 int* nx, int* ny, int* nz, float_sw4* x, float_sw4* y, float_sw4* z, float_sw4* met,
					 float_sw4* jac, int* order, float_sw4* sb, float_sw4* zmax, float_sw4* amp, 
					 float_sw4* xc, float_sw4* yc, float_sw4* xl, float_sw4* yl );
   void F77_FUNC(curvilinear4sg,CURVILINEAR4SG)( int *ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast, 
						 float_sw4* u, float_sw4* mu, float_sw4* lambda, float_sw4* met, float_sw4* jac,
						 float_sw4* lu, int* onesided, float_sw4* acof, float_sw4* bope, float_sw4* ghcof,
						 float_sw4* strx, float_sw4* stry, char* op );
   void F77_FUNC(freesurfcurvisg,FREESURFCURVISG)( int *ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast, 
						   int *nz, int* side, float_sw4* u, float_sw4* mu, float_sw4* lambda, float_sw4* met,
						   float_sw4* s, float_sw4* forcing, float_sw4* strx, float_sw4* stry );
   void F77_FUNC(dspev,DSPEV)(char & JOBZ, char & UPLO, int & N, double *AP, double *W, double *Z, int & LDZ, double *WORK, int & INFO);
}
#else
extern "C" 
{
   void F77_FUNC(dspev,DSPEV)(char & JOBZ, char & UPLO, int & N, double *AP, double *W, double *Z, int & LDZ, double *WORK, int & INFO);
}
void rhs4sg_rev( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
	     int nk, int* onesided, const float_sw4* a_acof, const float_sw4* a_bope, const float_sw4* a_ghcof,
	     float_sw4* a_lu, const float_sw4* a_u, const float_sw4* a_mu, const float_sw4* a_lambda, 
	     float_sw4 h, const float_sw4* a_strx, const float_sw4* a_stry, const float_sw4* a_strz  );

void rhs4sg( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
	     int nk, int* onesided, const float_sw4* a_acof, const float_sw4* a_bope, const float_sw4* a_ghcof,
	     float_sw4* a_lu, const float_sw4* a_u, const float_sw4* a_mu, const float_sw4* a_lambda, 
	     float_sw4 h, const float_sw4* a_strx, const float_sw4* a_stry, const float_sw4* a_strz  );

void rhs4sgcurv_rev( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
		     const float_sw4* a_u, const float_sw4* a_mu, const float_sw4* a_lambda, const float_sw4* a_met,
		     const float_sw4* a_jac, float_sw4* a_lu, const int* onesided, const float_sw4* acof,
		     const float_sw4* bope, const float_sw4* ghcof, const float_sw4* a_strx, const float_sw4* a_stry );

void rhs4sgcurv( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
		 float_sw4* a_u, float_sw4* a_mu, float_sw4* a_lambda, float_sw4* a_met,
		 float_sw4* a_jac, float_sw4* a_lu, int* onesided, float_sw4* acof,
		 float_sw4* bope, float_sw4* ghcof, float_sw4* a_strx, float_sw4* a_stry );
#endif

EW::EW( const string& filename ) :
   mCFL(1.3),
   mTstart(0.0),
   mTmax(0.0),
   mTimeIsSet(false),
   mNumberOfTimeSteps(-1),
   mPrintInterval(100),
   m_ghost_points(2),
   m_ext_ghost_points(2),
   m_ppadding(2),
   mVerbose(0),
   mQuiet(false),
   m_supergrid_damping_coefficient(0.02),
   m_sg_damping_order(4),
   m_sg_gp_thickness(20),
   m_use_supergrid(false),
   m_checkfornan(false),
   m_topography_exists(false),
   m_grid_interpolation_order(4),
   m_zetaBreak(0.95),
   m_point_source_test(false),
   mPath("./"),
   m_moment_test(false),
   m_pfs(false),
   m_nwriters(8),
   m_output_detailed_timing(false),
   m_save_trace(false),
   m_ndevice(0),
   m_corder(false),
   m_use_dg(false),
   mGeoAz(0.0),
   mLonOrigin(-118.0), 
   mLatOrigin(37.0), 
   mMetersPerDegree(111319.5),
   mMetersPerLongitude(87721.0),
   mConstMetersPerLongitude(false)
{
   m_gpu_blocksize[0] = 16;
   m_gpu_blocksize[1] = 16;
   m_gpu_blocksize[2] = 1;

   if( sizeof(float_sw4) == 4 )
      m_mpifloat = MPI_FLOAT;
   else if( sizeof(float_sw4) == 8 )
      m_mpifloat = MPI_DOUBLE;
   else
      CHECK_INPUT(false,"Error, could not identify float_sw4");

   MPI_Comm_rank( MPI_COMM_WORLD, &m_myrank );
   MPI_Comm_size( MPI_COMM_WORLD, &m_nprocs );

#ifdef CUDA_CODE
   m_iop = newmanaged(5+5+24+5+384+24+48+6+4);
   memset(m_iop,0,sizeof(float_sw4)*(5+5+24+5+384+24+48+6+4));
   m_iop2 = m_iop+5;
   m_bop2 = m_iop2+5;
   m_sbop = m_bop2+24;
   m_acof = m_sbop+5;
   m_bop = m_acof+384;
   m_bope = m_bop+24;
   m_ghcof = m_bope+48;
   m_hnorm = m_ghcof+6;
   
#endif


   m_restart_check_point = CheckPoint::nil;
   parseInputFile( filename );
   setupRun( );
   if(m_use_dg){
       timeStepLoopdGalerkin();
   }
   else
   {
       timesteploop( mU, mUm );
   }
    

}
EW::~EW(){
#ifdef CUDA_CODE
  delmanaged(m_iop);
  m_iop2=m_bop2=m_sbop=m_acof=m_bop=m_bope=m_ghcof=m_hnorm=NULL;
#endif
}
  
//-----------------------------------------------------------------------
int EW::computeEndGridPoint( float_sw4 maxval, float_sw4 h )
{
   const float_sw4 reltol = 1e-5;
   const float_sw4 abstol = 1e-12;
   float_sw4 fnpts = round(maxval/h+1);
   int npts;
   if( fabs((fnpts-1)*h-maxval) < reltol*fabs(maxval)+abstol )
      npts = static_cast<int>(fnpts);
   else
      npts = static_cast<int>(fnpts)+1;
   return npts;
}

//-----------------------------------------------------------------------
bool EW::startswith(const char begin[], char *line)
{
  int lenb = strlen(begin);

  // We ignore any preceeding whitespace
  while (strncmp(line, " ", 1) == 0 || strncmp(line, "\t", 1) == 0)
    line++;
    
  if (strncmp(begin, line, lenb) == 0)
     return true;
  else 
     return false;
}

//-----------------------------------------------------------------------
void EW::badOption(string name, char* option) const
{
   if (m_myrank == 0)
      cout << "\tWarning: ignoring " << name << " line option '" << option << "'" << endl;
}

//-----------------------------------------------------------------------
void EW::processGrid( char* buffer )
{
   float_sw4 x = 0.0, y=0.0, z=0.0, h=0.0;
   int nx=0, ny=0, nz=0;
   stringstream gridSetupErrStream;
   gridSetupErrStream << endl
		     << "----------------------------------------" << endl
		     << " Only five ways to setup grid: " << endl
		     << "  1. provide h and nx, ny, nz " << endl
		     << "  2. provide h and x, y, z " << endl
		     << "  3. provide x,y,z and nx " << endl
		     << "  4. provide x,y,z and ny " << endl
		     << "  5. provide x,y,z and nz " << endl
		     << "----------------------------------------" << endl
		     << endl;
   string gridSetupErr = gridSetupErrStream.str();

   char* token = strtok(buffer, " \t");

   token = strtok(NULL, " \t");
   string err = "ERROR in ProcessGrid: ";
   if( m_myrank == 0 )
      cout << endl << "* Processing the grid command..." << endl;
   while (token != NULL)
   {
     // while there are tokens in the string still
      if (startswith("#", token) || startswith(" ", buffer))
        // Ignore commented lines and lines with just a space.
	 break;
      if (startswith("ny=", token))
      {
	 token += 3;
	 CHECK_INPUT(atoi(token) > 0, 
		err << "ny is not a positive integer: " << token);
	 ny = atoi(token);
      }
      else if (startswith("nx=", token))
      {
	 token += 3;
	 CHECK_INPUT(atoi(token) > 0, 
		err << "nx is not a positive integer: " << token);
	 nx = atoi(token);
      }
      else if (startswith("nz=", token))
      {
	 token += 3;
	 CHECK_INPUT(atoi(token) >= 0, 
		err << "nz is not a positive integer: " << token);
	 nz = atoi(token);
      }
      else if (startswith("x=", token))
      {
	 token += 2;
	 CHECK_INPUT(atof(token) > 0.0, err << "x is not a positive float: " << token);
	 x = atof(token);
      }
      else if (startswith("y=", token))
      {
	 token += 2;
	 CHECK_INPUT(atof(token) >= 0.0, err << "y is negative: " << token);
	 y = atof(token);
      }
      else if (startswith("z=", token))
      {
	 token += 2; 
	 CHECK_INPUT(atof(token) > 0.0, err << "z is not a positive float: " << token);
	 z = atof(token);
      }
      else if (startswith("h=", token))
      {
	 token += 2;
	 CHECK_INPUT(atof(token) > 0.0, 
 	       err << "h is not a positive float: " << token);
	 h = atof(token);
      }
      else
      {
	 badOption("grid", token);
      }
      token = strtok(NULL, " \t");
   }
  
  //--------------------------------------------------------------------
  // There are only three ways to specify a grid.
  //--------------------------------------------------------------------
   if (h != 0.0)
   {
      if (nx > 0 || nz > 0 || ny > 0)
      {
      //----------------------------------------------------------------
      // 1.  nx, [ny], nz and h
      //----------------------------------------------------------------
	 CHECK_INPUT(nx && nz, gridSetupErr);
	 CHECK_INPUT(x == 0.0 && y == 0.0 && z == 0.0, gridSetupErr);
      }
      else
      {
      //--------------------------------------------------------------
      // 2.  x, [y], z and h
      //--------------------------------------------------------------
	 CHECK_INPUT(x > 0.0 && z > 0.0, gridSetupErr);
	 CHECK_INPUT(nx == 0 && ny == 0 && nz == 0, gridSetupErr);
      }
   }
   else
   {
    //--------------------------------------------------------------------
    // 3.  x, [y], z and nx|ny|nz
    //--------------------------------------------------------------------
      CHECK_INPUT(x > 0.0 && z > 0.0, gridSetupErr);
      CHECK_INPUT((nx > 0) + (ny > 0) + (nz > 0) == 1, gridSetupErr);
   }

   int nxprime, nyprime, nzprime;
   float_sw4 xprime, yprime, zprime;
   if (nx > 0 && h == 0.0)
   {
    // we set the number grid points in the x direction
    // so we'll compute the grid spacing from that.
      h = x / (nx-1);
      if (m_myrank == 0)
	   cout << "* Setting h to " << h << " from  x/(nx-1) (x=" << x << ", nx=" << nx << ")" << endl;
      
      nxprime = nx;
      nzprime = computeEndGridPoint(z, h);
      nyprime = computeEndGridPoint(y, h);
   }
   else if (ny > 0 && h == 0.0)
   {
    // set hte number of grid points from y direction and ny
      h = y/(ny-1);
      if (m_myrank == 0)
	   cout << "* Setting h to " << h << " from  y/(ny-1) (y=" << y << ", ny=" << ny << ")" << endl;
      nyprime = ny;
      nxprime = computeEndGridPoint(x, h);
      nzprime = computeEndGridPoint(z, h);
   }
   else if (nz > 0 && h == 0.0)
   {
    // set the number of grid points from z direction and nz
      h = z/(nz-1);
      if (m_myrank == 0)
	 cout << "* Setting h to " << h << " from  z/(nz-1) (z=" << z << ", nz=" << nz << ")" << endl;
      nzprime = nz;
      nxprime = computeEndGridPoint(x, h);
      nyprime = computeEndGridPoint(y, h);
   }
   else
   {
	//----------------------------------------------------
	// h was set by the user, so compute the appropriate
	// nx, ny, and nz or x, y, z.
	//----------------------------------------------------
      if (nx == 0 && x != 0.0)
	 nxprime = computeEndGridPoint(x, h);
      else if (nx != 0)
	 nxprime = nx;
      else
	 CHECK_INPUT(0, gridSetupErr);

      if (nz == 0 && z != 0.0)
	 nzprime = computeEndGridPoint(z, h);
      else if (nz != 0)
	 nzprime = nz;
      else
	 CHECK_INPUT(0, gridSetupErr);

      if (ny == 0 && y != 0.0)
	 nyprime = computeEndGridPoint(y, h);
      else if (ny != 0)
	 nyprime = ny;
      else
	 CHECK_INPUT(0, gridSetupErr);
   }
   
   if (m_myrank == 0 && mVerbose >=3)
      printf("**** Setting up the grid for a non-periodic problem\n");
    
   if (nxprime != nx && m_myrank == 0)
      cout << "* Setting nx to " << nxprime << " to be consistent with h=" << h << endl;
   if (nyprime != ny && m_myrank == 0)
      cout << "* Setting ny to " << nyprime << " to be consistent with h=" << h << endl;
   if (nzprime != nz && m_myrank == 0)
      cout << "* Setting nz to " << nzprime << " to be consistent with h=" << h << endl;

    // -------------------------------------------------------------
    // Now we adjust the geometry bounds based on the actual 
    // number of grid points used in each dimension.
    // -------------------------------------------------------------
   xprime = (nxprime-1)*h;
   zprime = (nzprime-1)*h;
   yprime = (nyprime-1)*h;
  
   float_sw4 eps = 1.e-9*sqrt(xprime*xprime+yprime*yprime+zprime*zprime);
  
   if (fabs(xprime-x) > eps && m_myrank == 0)
      cout << "* Changing x from " << x << " to " << xprime << " to be consistent with h=" << h << endl;
   if (fabs(zprime-z) > eps && m_myrank == 0)
      cout << "* Changing z from " << z << " to " << zprime << " to be consistent with h=" << h << endl;
   if (fabs(yprime-y) > eps && m_myrank == 0)
      cout << "* Changing y from " << y << " to " << yprime << " to be consistent with h=" << h << endl;

   m_nx_base = nxprime;
   m_ny_base = nyprime;
   m_nz_base = nzprime;
   m_h_base = h;
   m_global_xmax = xprime;
   m_global_ymax = yprime;
   m_global_zmax = zprime;
   m_global_zmin = 0; 
}
   
//-----------------------------------------------------------------------
void EW::processTime(char* buffer)
{
   float_sw4 t=0.0;
   int steps = -1;
   char* token = strtok(buffer, " \t");
   token = strtok(NULL, " \t");
   string err = "ERROR in processTime: ";
   while (token != NULL)
   {
      // while there are still tokens in the string
      if (startswith("#", token) || startswith(" ", buffer))
          // Ignore commented lines and lines with just a space.
          break;
      if (startswith("t=", token))
      {
	 token += 2; // skip t=
	 CHECK_INPUT(atof(token) >= 0.0, err << "t is not a positive float: " << token);
	 t = atof(token);
      }
      else if (startswith("steps=", token))
      {
	 token += 6; // skip steps=
	 CHECK_INPUT(atoi(token) >= 0, err << "steps is not a non-negative integer: " << token);
	 steps = atoi(token);
      }
      else
      {
	 badOption("time", token);
      }
      token = strtok(NULL, " \t");
   }
   CHECK_INPUT(!( (t > 0.0) && (steps >= 0) ),
          "Time Error: Cannot set both t and steps for time");
   if (t > 0.0)
   {
      mTmax   = t;
      mTstart = 0;
      mTimeIsSet = true;
   }
   else if (steps >= 0)
   {
      mTstart = 0;
      mNumberOfTimeSteps = steps;
      mTimeIsSet = false;
   }
     // Set UTC as current date
   time_t tsec;
   time( &tsec );
   struct tm *utctime = gmtime( &tsec );
   m_utc0[0] = utctime->tm_year+1900;
   m_utc0[1] = utctime->tm_mon+1;
   m_utc0[2] = utctime->tm_mday;
   m_utc0[3] = utctime->tm_hour;
   m_utc0[4] = utctime->tm_min;
   m_utc0[5] = utctime->tm_sec;
   m_utc0[6] = 0; //milliseconds not given by 'time', not needed here.
}

//-----------------------------------------------------------------------
void EW::processTopography(char * buffer )
{
   //
   // Note, m_topoFileName, m_topoExtFileName, m_maxIter, m_EFileResolution, m_QueryTyp could
   // have been declared local variables in EW::parseInputFile, and transfered as
   // procedure parameters to smoothTopography and getEfileInfo
   //
    char* token = strtok(buffer, " \t");
    CHECK_INPUT(strcmp("topography", token) == 0, 
 	    "ERROR: not a topography line...: " << token);
    string topoFile="surf.tp", style, fileName;
    bool needFileName=false, gotFileName=false;

    m_zetaBreak=0.95;
    m_grid_interpolation_order = 4;
    m_use_analytical_metric = false;

    token = strtok(NULL, " \t");

    while (token != NULL)
    {
      // while there are still tokens in the string 
       if (startswith("#", token) || startswith(" ", buffer))
        // Ignore commented lines and lines with just a space.
	  break;
       if (startswith("zmax=", token))
       {
	  token += 5; // skip logfile=
	  m_topo_zmax = atof(token);
       }
// //                       1234567890
       else if (startswith("order=", token))
       {
	  token += 6; // skip logfile=
	  m_grid_interpolation_order = atoi(token);
	  if (m_grid_interpolation_order < 2 || m_grid_interpolation_order > 7)
	  {
	     if (m_myrank == 0)
		cout << "order needs to be 2,3,4,5,6,or 7 not: " << m_grid_interpolation_order << endl;
	     MPI_Abort(MPI_COMM_WORLD, 1);
	  }
       
       }
//                          123456789
       else if( startswith("zetabreak=", token) ) // developer option: not documented in user's guide
       {
	  token += 10;
	  m_zetaBreak = atof(token);
	  CHECK_INPUT( m_zetaBreak > 0 && m_zetaBreak <= 1, "Error: zetabreak must be in [0,1], not " << m_zetaBreak);
       }
       //       else if (startswith("smooth=", token))
       //       {
       //	  token += 7; // skip smooth=
       //	  m_maxIter = atoi(token);
       //	  if (m_maxIter < 0 || m_maxIter > 1000)
       //	  {
       //	     if (m_myRank == 0)
       //		cout << "Number of smoothing iterations needs to be >=0 and <=1000, not: "<< m_maxIter << endl;
       //	     MPI_Abort(MPI_COMM_WORLD, 1);
       //	  }
       //       }
       else if( startswith("input=", token ) )
       {
	  token += 6;
	  style = token;
// new keyword: geographic, but keeping grid for backwards compatibility with WPP
//	  if (strcmp("grid", token) == 0 || strcmp("geographic", token) == 0)
//	  {
//	     m_topoInputStyle=GridFile;
//	     m_topography_exists=true;
//	     needFileName=true;
//	  }
//	  else if (strcmp("cartesian", token) == 0)
//	  {
//	     m_topoInputStyle=CartesianGrid;
//	     m_topography_exists=true;
//	     needFileName=true;
//	  }
//	  else if (strcmp("efile", token) == 0)
//	  {
//	     m_topoInputStyle=Efile;
//	     m_topography_exists=true;
//	     needFileName=true; // we require the file name to be given on the topography command line
//	  }
//	  else if (strcmp("rfile", token) == 0)
//	  {
//	     m_topoInputStyle=Rfile;
//	     m_topography_exists=true;
//	     needFileName=true; // we require the file name to be given on the topography command line
//	  }
//	  else if (strcmp("image", token) == 0)
//	  {
//	     m_topoInputStyle=TopoImage;
//	     m_topography_exists=true;
//	     needFileName=true; // we require the file name to be given on the topography command line
//	  }
	  if( strcmp("gaussian", token) == 0)
	     //	  else if (strcmp("gaussian", token) == 0)
	  {
	     m_topoInputStyle=GaussianHill;
	     m_topography_exists=true;
	  }
	  else
	  {
	     badOption("topography> input", token);
	  }
       }
       else if( startswith("file=", token ) )
       {
	  token += 5;
	  m_topoFileName = token;
	  gotFileName=true;
	//        if (m_myRank==0)
	// 	 cout << "read topo file name=" << m_topoFileName <<endl;
       }
       //#ifdef ENABLE_ETREE
       //       else if( startswith("etree=", token ) )
       //       {
       //	  token += 6;
       //	  m_topoFileName = token;
       //	  m_topoInputStyle=Efile;
       //	  m_topography_exists=true;
       //	  gotFileName=true;
       //       }
       //      else if (startswith("xetree=", token))
       //      {
       //	 token += 7; // skip xetree=
       //	 m_topoExtFileName = token;
       //      }
       //#endif
//                        12345678901
//       else if( startswith("resolution=", token ) )
//       {
//	  token += 11;
//	  m_EFileResolution = atof(token);
//	  CHECK_INPUT(m_EFileResolution>0.,"Resolution must be positive, not " << m_EFileResolution);
//       }
//                        123456789012
       else if( startswith("gaussianAmp=", token ) )
       {
	  token += 12;
	  m_GaussianAmp = atof(token);
       }
//                        123456789012
       else if( startswith("gaussianXc=", token ) )
       {
	  token += 11;
	  m_GaussianXc = atof(token);
       }
//                        123456789012
       else if( startswith("gaussianYc=", token ) )
       {
	  token += 11;
	  m_GaussianYc = atof(token);
       }
// //                        123456789012
       else if( startswith("gaussianLx=", token ) )
       {
	  token += 11;
	  m_GaussianLx = atof(token);
       }
//                        123456789012
       else if( startswith("gaussianLy=", token ) )
       {
	  token += 11;
	  m_GaussianLy = atof(token);
       }
       else if( startswith("analyticalMetric=", token ) )
       {
	  token += 17;
	  m_use_analytical_metric = strcmp(token,"1")==0 ||
	     strcmp(token,"true")==0 || strcmp(token,"yes")==0;
       }
       else
       {
	  badOption("topography", token);
       }
       token = strtok(NULL, " \t");
    }
    if (needFileName)
       CHECK_INPUT(gotFileName, 
		   "ERROR: no topography file name specified...: " << token);

    CHECK_INPUT(m_topoInputStyle == GaussianHill,
		"Topography style " << m_topoInputStyle << " not yet implemented " << endl);

    if( m_topoInputStyle != GaussianHill && m_use_analytical_metric )
    {
       m_use_analytical_metric = false;
       if( m_myrank == 0 )
	  cout << "Analytical metric only defined for Gaussian Hill topography" <<
	     " topography analyticalMetric option will be ignored " << endl;
    }
}

//-----------------------------------------------------------------------
void EW::processFileIO(char* buffer)
{
   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("fileio", token) == 0, "ERROR: not a fileio line...: " << token);
   token = strtok(NULL, " \t");
   string err = "FileIO Error: ";

   while (token != NULL)
    {
       if (startswith("#", token) || startswith(" ", buffer))
          break;
       if(startswith("path=", token)) {
          token += 5; // skip path=
	  mPath = token;
	  mPath += '/';
	  //          path = token;
       }
       else if (startswith("verbose=", token))
       {
          token += 8; // skip verbose=
          CHECK_INPUT(atoi(token) >= 0, err << "verbose must be non-negative, not: " << token);
          mVerbose = atoi(token);
       }
       else if (startswith("printcycle=", token))
       {
          token += 11; // skip printcycle=
          CHECK_INPUT(atoi(token) > -1,
	         err << "printcycle must be zero or greater, not: " << token);
          mPrintInterval = atoi(token);
       }
       else if (startswith("pfs=", token))
       {
          token += 4; // skip pfs=
          m_pfs = (atoi(token) == 1);
       }
//                          1234567890
       else if (startswith("nwriters=", token))
       {
          token += 9; // skip nwriters=
          CHECK_INPUT(atoi(token) > 0,
	         err << "nwriters must be positive, not: " << token);
          m_nwriters = atoi(token);
       }
       else
       {
          badOption("fileio", token);
       }
       token = strtok(NULL, " \t");
    }
}

//-----------------------------------------------------------------------
void EW::processCheckPoint(char* buffer)
{
   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("checkpoint", token) == 0, "ERROR: not a checkpoint line...: " << token);
   token = strtok(NULL, " \t");
   string err = "CheckPoint Error: ";
   int cycle=-1, cycleInterval=0;
   float_sw4 time=0.0, timeInterval=0.0;
   bool timingSet=false;
   string filePrefix = "restart";
   size_t bufsize=10000000;

   while (token != NULL)
    {
       if (startswith("#", token) || startswith(" ", buffer))
          break;
      if (startswith("time=", token) )
      {
	 token += 5; // skip time=
	 CHECK_INPUT( atof(token) >= 0., err << "time must be a non-negative number, not: " << token);
	 time = atof(token);
	 timingSet = true;
      }
      else if (startswith("timeInterval=", token) )
      {
	 token += 13; // skip timeInterval=
	 CHECK_INPUT( atof(token) >= 0., err<< "timeInterval must be a non-negative number, not: " << token);
	 timeInterval = atof(token);
	 timingSet = true;
      }
      else if (startswith("cycle=", token) )
      {
	 token += 6; // skip cycle=
	 CHECK_INPUT( atoi(token) >= 0., err << "cycle must be a non-negative integer, not: " << token);
	 cycle = atoi(token);
	 timingSet = true;
      }
      else if (startswith("cycleInterval=", token) )
      {
	 token += 14; // skip cycleInterval=
	 CHECK_INPUT( atoi(token) >= 0., err << "cycleInterval must be a non-negative integer, not: " << token);
	 cycleInterval = atoi(token);
	 timingSet = true;
      }
      else if (startswith("file=", token))
      {
	 token += 5; // skip file=
	 filePrefix = token;
      }
      else if (startswith("bufsize=", token))
      {
	 token += 8; // skip bufsize=
	 bufsize = atoi(token);
      }
      else
      {
	 badOption("checkpoint", token);
      }
      token = strtok(NULL, " \t");
   }
   CHECK_INPUT( timingSet, "Processing checkpoint command: " << 
		"at least one timing mechanism must be set: cycle, time, cycleInterval or timeInterval"  << endl );
   CheckPoint* chkpt = new CheckPoint( this, time, timeInterval, cycle, cycleInterval, filePrefix, bufsize ); 
   m_check_points.push_back(chkpt);
}

//-----------------------------------------------------------------------
void EW::processRestart(char* buffer)
{
   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("restart", token) == 0, "ERROR: not a restart line...: " << token);
   token = strtok(NULL, " \t");
   //   string err = "Restart Error: ";
   //   int cycle=-1, cycleInterval=0;
   //   double time=0.0, timeInterval=0.0;
   //   bool timingSet=false;
   string fileName;
   bool filenamegiven = false;
   size_t bufsize=10000000;

   while (token != NULL)
    {
       if (startswith("#", token) || startswith(" ", buffer))
          break;
      if (startswith("file=", token) )
      {
	 token += 5; // skip file=
	 fileName = token;
	 filenamegiven = true;
      }
      else if (startswith("bufsize=", token))
      {
	 token += 8; // skip bufsize=
	 bufsize = atoi(token);
      }
      else
      {
	 badOption("restart", token);
      }
      token = strtok(NULL, " \t");
   }
   CHECK_INPUT( filenamegiven, "Processing restart command: " << 
		"restart file name must be given"  << endl );
   CHECK_INPUT( m_restart_check_point == CheckPoint::nil, "Processing restart command: "<<
		" There can only be one restart file");
   m_restart_check_point = new CheckPoint( this, fileName, bufsize );

}

//-----------------------------------------------------------------------
void EW::processTestPointSource(char* buffer)
{
   char* token = strtok(buffer, " \t");
   token = strtok(NULL, " \t");
   float_sw4 cs = 1.0, rho=1.0, cp=sqrt(3.0);
   bool free_surface=false;
   while (token != NULL)
   {
      if (startswith("#", token) || startswith(" ", buffer))
	 break;
      if (startswith("cp=", token))
      {
	 token += 3; 
	 cp = atof(token);
      }
      else if (startswith("cs=", token))
      {
	 token += 3; 
	 cs = atof(token);
      }
      else if (startswith("rho=", token))
      {
	 token += 4; 
	 rho = atof(token);
      }
      else if (startswith("diractest=", token))
      {
	 token += 10; 
	 if( strcmp(token,"1")==0 || strcmp(token,"true")==0 )
	    m_moment_test = true;
      }
      else if (startswith("halfspace=", token))
      {
	 token += 10; 
	 free_surface = ( strcmp(token,"1")==0 || strcmp(token,"true")==0 );
      }
      else
      {
	 badOption("testpointsource", token);
      }
      token = strtok(NULL, " \t");
   }
   m_point_source_test = true;
   float_sw4 mu = rho*cs*cs;
   float_sw4 la = rho*cp*cp-2*mu;
   for( int g=0 ; g < mNumberOfGrids ; g++ )
   {
      mRho[g].set_value(rho);
      mMu[g].set_value(mu);
      mLambda[g].set_value(la);
   }
   for( int side=0 ; side < 6 ; side++ )
      mbcGlobalType[side]=bSuperGrid;
   if( free_surface )
      mbcGlobalType[4]=bStressFree;
}

//----------------------------------------------------------------------------
void EW::processSource( char* buffer )
{
   Source* sourcePtr;
   float_sw4 m0 = 1.0;
   float_sw4 t0=0.0, f0=1.0, freq=1.0;
  // Should be center of the grid
   float_sw4 x = 0.0, y = 0.0, z = 0.0;
   //  int i = 0, j = 0, k = 0;
   float_sw4 mxx=0.0, mxy=0.0, mxz=0.0, myy=0.0, myz=0.0, mzz=0.0;
   //  float_sw4 strike=0.0, dip=0.0, rake=0.0;
   float_sw4 fx=0.0, fy=0.0, fz=0.0;
   int isMomentType = -1;
  
   //  float_sw4 lat = 0.0, lon = 0.0, depth = 0.0;
   float_sw4 depth= 0.0;
   bool topodepth = false, depthSet=false, zSet=false;
  
   bool cartCoordSet = false;
   //  bool geoCoordSet = false;
   //  bool strikeDipRake = false;
   //  bool dfileset=false;
   //  bool sacbaseset = false;

   //   int ncyc = 0;
   //   bool ncyc_set = false;

   float_sw4* par=NULL;
   int* ipar=NULL;
   int npar=0, nipar=0;
   int ncyc = 5;

   timeDep tDep = iRickerInt;
   char formstring[100];
   //  char dfile[1000];
   strcpy(formstring, "Ricker");

   char* token = strtok(buffer, " \t");
   token = strtok(NULL, " \t");
   string err = "ERROR in ProcessSource: ";

   //  string cartAndGeoErr = "source command: Cannot set both a geographical (lat,lon) and cartesian coordinate (x,y)";
   string pointAndMomentErr = "source command: Cannot set both a point source and moment tensor formulation";
   while (token != NULL)
   {
      // while there are tokens in the string still
      if (startswith("#", token) || startswith(" ", buffer))
          // Ignore commented lines and lines with just a space.
          break;
      if (startswith("m0=", token) )
      {
	 token += 3; // skip m0=
	 CHECK_INPUT(atof(token) >= 0.0, 
                  err << "source command: scalar moment term must be positive, not: " << token);
	 m0 = atof(token);
      }
      else if (startswith("x=", token))
      {
         token += 2; // skip x=
         x = atof(token);
         cartCoordSet = true; 
      }
      else if (startswith("y=", token))
      {
         token += 2; // skip y=
         y = atof(token);
         cartCoordSet = true;
      }
      else if (startswith("z=", token))
      {
         token += 2; // skip z=
// with topography, the z-coordinate can have both signs!
         z = atof(token);
	 topodepth=false; // this is absolute depth
         zSet = true;
      }
      else if (startswith("depth=", token)) // this is the same as topodepth: different from WPP
      {
         token += 6; // skip depth=
         depth = atof(token);
	 topodepth = true;
         CHECK_INPUT(depth >= 0.0,
		     err << "source command: Depth below topography must be greater than or equal to zero");
	 depthSet=true;
      }
//                         1234567890
      else if (startswith("Mxx=", token) || startswith("mxx=", token))
      {
         CHECK_INPUT(isMomentType != 0, err << pointAndMomentErr);
         token += 4; // skip Mxx=
         mxx = atof(token);
         isMomentType = 1;
      }
      else if (startswith("Mxy=", token) || startswith("mxy=", token))
      {
         CHECK_INPUT(isMomentType != 0, err << pointAndMomentErr);
         token += 4; // skip Mxy=
         mxy = atof(token);
	  isMomentType = 1;
      }
      else if (startswith("Mxz=", token) || startswith("mxz=", token))
      {
         CHECK_INPUT(isMomentType != 0, err << pointAndMomentErr);
         token += 4; // skip Mxz=
         mxz = atof(token);
         isMomentType = 1;
      }
      else if (startswith("Myy=", token) || startswith("myy=", token))
      {
         CHECK_INPUT(isMomentType != 0, err << pointAndMomentErr);
         token += 4; // skip Myy=
         myy = atof(token);
         isMomentType = 1;
      }
      else if (startswith("Myz=", token) || startswith("myz=", token))
      {
         CHECK_INPUT(isMomentType != 0, err << pointAndMomentErr);
         token += 4; // skip Myz=
         myz = atof(token);
         isMomentType = 1;
      }
      else if (startswith("Mzz=", token) || startswith("mzz=", token))
      {
         CHECK_INPUT(isMomentType != 0, err << pointAndMomentErr);
         token += 4; // skip Mzz=
         mzz = atof(token);
         isMomentType = 1;
      }
      else if (startswith("Fz=", token) || startswith("fz=", token))
      {
         CHECK_INPUT(isMomentType != 1, err << pointAndMomentErr);
         token += 3; // skip Fz=
         fz = atof(token);
         isMomentType = 0;
      }
      else if (startswith("Fx=", token) || startswith("fx=", token))
      {
         CHECK_INPUT(isMomentType != 1, err << pointAndMomentErr);
         token += 3; // skip Fx=
         fx = atof(token);
         isMomentType = 0;
      }
      else if (startswith("Fy=", token) || startswith("fy=", token))
      {
         CHECK_INPUT(isMomentType != 1, err << pointAndMomentErr);
         token += 3; // skip Fy=
         fy = atof(token);
         isMomentType = 0;
      }
      else if (startswith("t0=", token))
      {
         token += 3; // skip t0=
         t0 = atof(token);
      }
      else if (startswith("freq=", token))
      {
         token += 5; // skip freq=
         freq = atof(token);
         CHECK_INPUT(freq > 0,
                 err << "source command: Frequency must be > 0");
      }
      else if (startswith("f0=", token))
      {
         CHECK_INPUT(isMomentType != 1,
                 err << "source command: Cannot set force amplitude for moment tensor terms");
	 token += strlen("f0=");
         f0 = atof(token);
      }
      else if (startswith("type=",token))
      {
         token += 5;
         strncpy(formstring, token,100);
         if (!strcmp("Ricker",formstring))
            tDep = iRicker;
         else if (!strcmp("Gaussian",formstring))
            tDep = iGaussian;
         else if (!strcmp("Ramp",formstring))
            tDep = iRamp;
         else if (!strcmp("Triangle",formstring))
            tDep = iTriangle;
         else if (!strcmp("Sawtooth",formstring))
            tDep = iSawtooth;
         else if (!strcmp("SmoothWave",formstring))
            tDep = iSmoothWave;
         else if (!strcmp("Erf",formstring) || !strcmp("GaussianInt",formstring) )
            tDep = iErf;
         else if (!strcmp("VerySmoothBump",formstring))
            tDep = iVerySmoothBump;
         else if (!strcmp("RickerInt",formstring) )
            tDep = iRickerInt;
         else if (!strcmp("Brune",formstring) )
	    tDep = iBrune;
         else if (!strcmp("BruneSmoothed",formstring) )
	    tDep = iBruneSmoothed;
         else if (!strcmp("DBrune",formstring) )
	    tDep = iDBrune;
         else if (!strcmp("GaussianWindow",formstring) )
	    tDep = iGaussianWindow;
         else if (!strcmp("Liu",formstring) )
	    tDep = iLiu;
         else if (!strcmp("Dirac",formstring) )
	    tDep = iDirac;
         else if (!strcmp("C6SmoothBump",formstring) )
	    tDep = iC6SmoothBump;
	 else
            if (m_myrank == 0)
	      cout << "unknown time function: " << formstring << endl << " using default RickerInt function." << endl;
      }
      else
      {
         badOption("source", token);
      }
      token = strtok(NULL, " \t");
   }

   CHECK_INPUT(depthSet || zSet,
	       err << "source command: depth, topodepth or z-coordinate must be specified");
   if (depthSet)
   {
      z = depth;
   }
   if (cartCoordSet)
   {
      float_sw4 xmin = 0.;
      float_sw4 ymin = 0.;
      float_sw4 zmin;

// only check the z>zmin when we have topography. For a flat free surface, we will remove sources too 
// close or above the surface in the call to mGlobalUniqueSources[i]->correct_Z_level()
      if (m_topography_exists) // topography command must be read before the source command
	 zmin = m_global_zmin;
      else
	 zmin = 0;

      if ( (m_topography_exists && (x < xmin || x > m_global_xmax || y < ymin || y > m_global_ymax )) ||
	   (!m_topography_exists && (x < xmin || x > m_global_xmax || y < ymin || y > m_global_ymax || 
				    z < zmin || z > m_global_zmax)) )
      {
	 stringstream sourceposerr;
	 sourceposerr << endl
		   << "***************************************************" << endl
		   << " FATAL ERROR:  Source positioned outside grid!  " << endl
		   << endl
		   << " Source Type: " << formstring << endl
		   << "              @ x=" << x 
		   << " y=" << y << " z=" << z << endl 
		   << endl;
	    
	 if ( x < xmin )
	    sourceposerr << " x is " << xmin - x << 
	  " meters away from min x (" << xmin << ")" << endl;
	 else if ( x > m_global_xmax)
	    sourceposerr << " x is " << x - m_global_xmax << 
	       " meters away from max x (" << m_global_xmax << ")" << endl;
	 if ( y < ymin )
	    sourceposerr << " y is " << ymin - y << 
	  " meters away from min y (" << ymin << ")" << endl;
	 else if ( y > m_global_ymax)
	    sourceposerr << " y is " << y - m_global_ymax << 
	  " meters away from max y (" << m_global_ymax << ")" << endl;
	 if ( z < zmin )
	    sourceposerr << " z is " << zmin - z << 
	  " meters away from min z (" << zmin << ")" << endl;
	 else if ( z > m_global_zmax)
	    sourceposerr << " z is " << z - m_global_zmax << 
	  " meters away from max z (" << m_global_zmax << ")" << endl;
	 sourceposerr << "***************************************************" << endl;
	 if (m_myrank == 0)
	    cout << sourceposerr.str();
	 MPI_Abort(MPI_COMM_WORLD, 1);
      }
   }
   if (isMomentType)
   {
    // Remove amplitude variable
      mxx *= m0;
      mxy *= m0;
      mxz *= m0;
      myy *= m0;
      myz *= m0;
      mzz *= m0;
    // these have global location since they will be used by all processors
      sourcePtr = new Source(this, freq, t0, x, y, z, mxx, mxy, mxz, myy, myz, mzz,
			     tDep, formstring, topodepth, ncyc, par, npar, ipar, nipar, false ); // false is correctStrengthForMu

      if (sourcePtr->ignore())
      {
	 delete sourcePtr;
      }
      else
      {
	 m_globalUniqueSources.push_back(sourcePtr);
      }
   }
   else // point forcing
   {
    // Remove amplitude variable
      fx *= f0;
      fy *= f0;
      fz *= f0;
      // global version (gets real coordinates)
      sourcePtr = new Source(this, freq, t0, x, y, z, fx, fy, fz, tDep, formstring, topodepth, ncyc,
			   par, npar, ipar, nipar, false ); // false is correctStrengthForMu

      //...and add it to the list of forcing terms
      if (sourcePtr->ignore())
      {
	 delete sourcePtr;
      }
      else
      {
	 m_globalUniqueSources.push_back(sourcePtr);
      }
   }
}

//-----------------------------------------------------------------------
void EW::processSuperGrid(char *buffer)
{
   char* token = strtok(buffer, " \t");
   token = strtok(NULL, " \t");
   int sg_thickness; // sg_transition;
   float_sw4 sg_coeff;
   bool thicknessSet=false, dampingCoeffSet=false; // , transitionSet=false
   while (token != NULL)
   {
      if (startswith("#", token) || startswith(" ", buffer))
        // Ignore commented lines and lines with just a space.
	 break;

//                  1234567890
      if (startswith("gp=", token)) // in number of grid sizes (different from WPP)
      {
	 token += 3;
	 sg_thickness = atoi(token);
	 CHECK_INPUT(sg_thickness>0, "The number of grid points in the supergrid damping layer must be positive, not: "<< sg_thickness);
	 thicknessSet = true;
      }
      else if (startswith("dc=", token))
      {
	 token += 3;
	 sg_coeff = atof(token);
	 CHECK_INPUT(sg_coeff>=0., "The supergrid damping coefficient must be non-negative, not: "<<sg_coeff);
	 dampingCoeffSet=true;
      }
      else
      {
	 badOption("supergrid", token);
      }
      token = strtok(NULL, " \t");
   } // end while token
  
   if (thicknessSet)
      m_sg_gp_thickness = sg_thickness;

   if (dampingCoeffSet)
      m_supergrid_damping_coefficient = sg_coeff;
   else if( m_sg_damping_order == 4 )
      m_supergrid_damping_coefficient = 0.02;
   else if( m_sg_damping_order == 6 )
      m_supergrid_damping_coefficient = 0.005;
}

//-----------------------------------------------------------------------
void EW::processDeveloper(char* buffer)
{
   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("developer", token) == 0, "ERROR: not a developer line...: " << token);
   token = strtok(NULL, " \t");
   while (token != NULL)
   {
     // while there are tokens in the string still
     if (startswith("#", token) || startswith(" ", buffer))
       // Ignore commented lines and lines with just a space.
        break;
     if( startswith("cfl=",token) )
     {
	token += 4;
	float_sw4 cfl = atof(token);
	CHECK_INPUT( cfl > 0, "Error negative CFL number");
	//	set_cflnumber( cfl );
	mCFL = cfl;
     }
     else if( startswith("checkfornan=",token) )
     {
	token += 12;
	m_checkfornan = strcmp(token,"1")==0 || strcmp(token,"on")==0 || strcmp(token,"yes")==0;
     }
     else if( startswith("reporttiming=",token) )
     {
	token += 13;
	m_output_detailed_timing = strcmp(token,"1")==0 || strcmp(token,"on")==0 || strcmp(token,"yes")==0;
     }
     else if( startswith("trace=",token) )
     {
	token += 6;
	m_save_trace = strcmp(token,"yes")==0
	   || strcmp(token,"1")==0 || strcmp(token,"on")==0;
     }
     else if( startswith("thblocki=",token) )
     {
	token += 9;
	m_gpu_blocksize[0] = atoi(token);
     }
     else if( startswith("thblockj=",token) )
     {
	token += 9;
	m_gpu_blocksize[1] = atoi(token);
     }
     else if( startswith("thblockk=",token) )
     {
	token += 9;
	m_gpu_blocksize[2] = atoi(token);
     }
     else if( startswith("corder=",token) )
     {
	token += 7;
	m_corder = strcmp(token,"yes")==0
	   || strcmp(token,"1")==0 || strcmp(token,"on")==0;
	Sarray::m_corder = m_corder;
#ifndef SW4_CROUTINES
	CHECK_INPUT(m_corder==0,"ERROR: developer option corder, must be zero when fortran routines are used");
#endif	
     }
     else
     {
       badOption("developer", token);
     }
     token = strtok(NULL, " \t");
   }
}

//------------------------------------------------------------------------
void EW::processMaterialBlock( char* buffer )
{
  float_sw4 vpgrad=0.0, vsgrad=0.0, rhograd=0.0;
  bool x1set=false, x2set=false, y1set=false, y2set=false, 
    z1set=false, z2set=false;

  float_sw4 x1=0.0, x2=0.0, y1=0.0, y2=0.0, z1=0.0, z2=0.0;


  string name = "Block";

  char* token = strtok(buffer, " \t");
  CHECK_INPUT(strcmp("block", token) == 0,
	      "ERROR: material block can be set by a block line, not: " << token);

  string err = token;
  err += " Error: ";

  token = strtok(NULL, " \t");

  float_sw4 vp=-1, vs=-1, rho=-1, qp=-1, qs=-1, freq=1;
  bool absDepth=false;

  while (token != NULL)
    {
      // while there are tokens in the string still
      if (startswith("#", token) || startswith(" ", buffer))
          // Ignore commented lines and lines with just a space.
	break;
// the xygrad keywords must occur before the corresponding xy keywords
      if (startswith("rhograd=", token))
      {
         token += 8; // skip rhograd=
         rhograd = atof(token);
      }
      else if (startswith("vpgrad=", token))
      {
         token += 7; // skip vpgrad=
         vpgrad = atof(token);
      }
      else if (startswith("vsgrad=", token))
      {
         token += 7; // skip vsgrad=
         vsgrad = atof(token);
      }
      else if (startswith("vp=", token) )
      {
         token += 3; // skip vp=
         vp = atof(token);
      }
      else if (startswith("vs=", token) )
      {
         token += 3; // skip vs=
         vs = atof(token);
      }
      else if (startswith("rho=", token))
      {
         token += 4; // skip rho=
         rho = atof(token);
      }
      else if (startswith("r=", token)) // superseded by rho=, but keep for backward compatibility
      {
         token += 2; // skip r=
         rho = atof(token);
      }
      else if (startswith("Qs=", token) || startswith("qs=",token) )
      {
         token += 3; // skip qs=
         qs = atof(token);
      }
      else if (startswith("Qp=", token) || startswith("qp=",token) )
      {
         token += 3; // skip qp=
         qp = atof(token);
      }
//                         1234567890
      else if (startswith("absdepth=", token) )
      {
	token += 9; // skip absdepth=
	absDepth = (bool) atoi(token);
      }
      else if (startswith("x1=", token))
      {
         token += 3; // skip x1=
         x1 = atof(token);
         x1set = true;
      }
      else if (startswith("x2=", token))
      {
         token += 3; // skip x2=
         x2 = atof(token);
         x2set = true;
      }
      else if (startswith("y1=", token))
      {
         token += 3; // skip y1=
         y1 = atof(token);
         y1set = true;
      }
      else if (startswith("y2=", token))
      {
         token += 3; // skip y2=
         y2 = atof(token);
         y2set = true;
      }
      else if (startswith("z1=", token))
      {
         token += 3; // skip z1=
         z1 = atof(token);
         z1set = true;
      }
      else if (startswith("z2=", token))
      {
         token += 3; // skip z2=
         z2 = atof(token);
         z2set = true;
      }
      else
      {
         badOption("block", token);
      }
      token = strtok(NULL, " \t");
    }
  // End parsing...
  

  // Set up a block on the EW object.

  if (x1set)
  {
     // CHECK_INPUT(x1 >= 0.,
     // 	     err << "x1 is less than the minimum x, " 
     // 	     << x1 << " < " << 0.);
     CHECK_INPUT(x1 <= m_global_xmax,
	     err << "x1 is greater than the maximum x, " 
	     << x1 << " > " << m_global_xmax);
  }
  else
    x1 = -m_global_xmax; //x1 = 0.;

  if (x2set)
  {
     CHECK_INPUT(x2 >= 0.,
             err << "x2 is less than the minimum x, " 
             << x2 << " < " << 0.);
     // CHECK_INPUT(x2 <= m_global_xmax,
     //         err << "x2 is greater than the maximum x, " 
     //         << x2 << " > " << m_global_xmax);
  }
  else
    x2 = 2.*m_global_xmax;//x2 = m_global_xmax;

  CHECK_INPUT( x2 >= x1, " (x1..x2), upper bound is smaller than lower bound");
  
  //--------------------------------------------------------
  // Set j bounds, goes with Y in WPP
  //--------------------------------------------------------
  if (y1set)
  {
     // CHECK_INPUT(y1 >= 0.,
     // 	     err << "y1 is less than the minimum y, " << y1 << " < " << 0.);

     CHECK_INPUT(y1 <= m_global_ymax,
		  err << "y1 is greater than the maximum y, " << y1 << " > " << m_global_ymax);
  }
  else
    y1 = -m_global_ymax;//y1 = 0.;
      
  if (y2set)
  {
     CHECK_INPUT(y2 >= 0.,
	     err << "y2 is less than the minimum y, " << y2 << " < " << 0.);
  }
  else
    y2 = 2.*m_global_ymax;//y2 = m_global_ymax;

  CHECK_INPUT( y2 >= y1, " (y1..y2), upper bound is smaller than lower bound");

  if (z1set)
  {
    // CHECK_INPUT(topographyExists() || z1 >= 0.,
    //         err << "z1 is less than the minimum z, " << z1 << " < " << 0.);
    CHECK_INPUT(z1 <= m_global_zmax, 
            err << "z1 is greater than the maximum z, " << z1 << " > " << m_global_zmax);
  }
  else
    z1 = m_global_zmin - (m_global_zmax-m_global_zmin);

  if (z2set)
  {
    CHECK_INPUT(topographyExists() || z2 >= 0.,
            err << "z2 is less than the minimum z, " << z2 << " < " << 0.);
    // CHECK_INPUT(z2 <= m_global_zmax,
    // 		err << "z2 is greater than the maximum z, " << z2 << " > " << m_global_zmax);
  }
  else
    z2 = m_global_zmax + (m_global_zmax-m_global_zmin);

  CHECK_INPUT( z2 >= z1, " (z1..z2), upper bound is smaller than lower bound");

  if( getVerbosity() >=2 &&  m_myrank == 0 )
     cout << name << " has bounds " << x1 << " " << x2 << " " << y1 << " "
	  << y2 << " " << z1 << " " << z2 << endl;

  CHECK_INPUT( vs > 0 && vp > 0 && rho > 0 , "Error in block " << name << " vp vs rho are   "
	       << vp << " " << vs << " " << rho );

  MaterialBlock* bl = new MaterialBlock( this ,rho, vs, vp, x1, x2, y1, y2, z1, z2, qs, qp, freq );
  bl->set_gradients( rhograd, vsgrad, vpgrad );
  bl->set_absoluteDepth( absDepth );
  m_mtrlblocks.push_back(bl);
}

//-----------------------------------------------------------------------
void EW::processReceiver(char* buffer )
{
  float_sw4 x=0.0, y=0.0, z=0.0;
  float_sw4 lat = 0.0, lon = 0.0, depth = 0.0;
  bool cartCoordSet = false, geoCoordSet = false;
  string fileName = "station";
  string staName = "station";
  bool staNameGiven=false;
  
  int writeEvery = 1000;

  bool topodepth = false;

  bool usgsformat = 0, sacformat=1; // default is to write sac files
  TimeSeries::receiverMode mode=TimeSeries::Displacement;

  char* token = strtok(buffer, " \t");
  bool nsew=false; 
  //int vel=0;

// tmp
//  cerr << "******************** INSIDE process receiver *********************" << endl;

  CHECK_INPUT(strcmp("rec", token) == 0 || strcmp("sac", token) == 0, "ERROR: not a rec line...: " << token);
  token = strtok(NULL, " \t");

  string err = "RECEIVER Error: ";

//* testing
  // if (proc_zero())
  //   cout << "start parsing of receiver command, token:" << token << "(end token)" << endl;

  while (token != NULL)
  {
     // while there are tokens in the string still
     //     cout << m_myRank << " token " << token <<"x"<<endl;

     if (startswith("#", token) || startswith(" ", buffer))
        // Ignore commented lines and lines with just a space.
        break;
     if (startswith("x=", token))
     {
        CHECK_INPUT(!geoCoordSet,
                err << "receiver command: Cannot set both a geographical (lat, lon) and a cartesian (x,y) coordinate");
        token += 2; // skip x=
        cartCoordSet = true;
        x = atof(token);
        CHECK_INPUT(x >= 0.0,
		    "receiver command: x must be greater than or equal to 0, not " << x);
        CHECK_INPUT(x <= m_global_xmax,
		    "receiver command: x must be less than or equal to xmax, not " << x);
     }
     else if (startswith("y=", token))
     {
        CHECK_INPUT(!geoCoordSet,
                err << "receiver command: Cannot set both a geographical (lat, lon) and a cartesian (x,y) coordinate");
        token += 2; // skip y=
        cartCoordSet = true;
        y = atof(token);
        CHECK_INPUT(y >= 0.0,
                "receiver command: y must be greater than or equal to 0, not " << y);
        CHECK_INPUT(y <= m_global_ymax,
		    "receiver command: y must be less than or equal to ymax, not " << y);
     }
     else if (startswith("lat=", token))
     {
        CHECK_INPUT(!cartCoordSet,
                err << "receiver command: Cannot set both a geographical (lat, lon) and a cartesian (x,y) coordinate");
        token += 4; // skip lat=
        lat = atof(token);
        CHECK_INPUT(lat >= -90.0,
                "receiver command: lat must be greater than or equal to -90 degrees, not " 
                << lat);
        CHECK_INPUT(lat <= 90.0,
                "receiver command: lat must be less than or equal to 90 degrees, not "
                << lat);
        geoCoordSet = true;
     }
     else if (startswith("lon=", token))
     {
        CHECK_INPUT(!cartCoordSet,
                err << "receiver command: Cannot set both a geographical (lat, lon) and a cartesian (x,y) coordinate");
        token += 4; // skip lon=
        lon = atof(token);
        CHECK_INPUT(lon >= -180.0,
                "receiver command: lon must be greater or equal to -180 degrees, not " 
                << lon);
        CHECK_INPUT(lon <= 180.0,
                "receiver command: lon must be less than or equal to 180 degrees, not "
                << lon);
        geoCoordSet = true;
     }
     else if (startswith("z=", token))
     {
       token += 2; // skip z=
       depth = z = atof(token);
       topodepth = false; // absolute depth (below mean sea level)
       CHECK_INPUT(z <= m_global_zmax,
		   "receiver command: z must be less than or equal to zmax, not " << z);
     }
     else if (startswith("depth=", token))
     {
        token += 6; // skip depth=
       z = depth = atof(token);
       topodepth = true; // by depth we here mean depth below topography
       CHECK_INPUT(depth >= 0.0,
	       err << "receiver command: depth must be greater than or equal to zero");
       CHECK_INPUT(depth <= m_global_zmax,
		   "receiver command: depth must be less than or equal to zmax, not " << depth);
     }
//                        1234567890
     else if (startswith("topodepth=", token))
     {
        token += 10; // skip topodepth=
       z = depth = atof(token);
       topodepth = true; // by depth we here mean depth below topography
       CHECK_INPUT(depth >= 0.0,
	       err << "receiver command: depth must be greater than or equal to zero");
       CHECK_INPUT(depth <= m_global_zmax,
		   "receiver command: depth must be less than or equal to zmax, not " << depth);
     }
     else if(startswith("file=", token))
     {
        token += 5; // skip file=
        fileName = token;
     }
     else if (startswith("sta=", token))
     {
        token += strlen("sta=");
        staName = token;
	staNameGiven=true;
     }
     else if( startswith("nsew=", token) )
     {
        token += strlen("nsew=");
        nsew = atoi(token) == 1;
     }
     // else if( startswith("velocity=", token) )
     // {
     //    token += strlen("velocity=");
     //    vel = atoi(token);
     // }
     else if (startswith("writeEvery=", token))
     {
       token += strlen("writeEvery=");
       writeEvery = atoi(token);
       CHECK_INPUT(writeEvery >= 0,
	       err << "sac command: writeEvery must be set to a non-negative integer, not: " << token);
     }
     else if( startswith("usgsformat=", token) )
     {
        token += strlen("usgsformat=");
        usgsformat = atoi(token);
     }
     else if( startswith("sacformat=", token) )
     {
        token += strlen("sacformat=");
        sacformat = atoi(token);
     }
     else if( startswith("variables=", token) )
     {
//* testing
       // if (proc_zero())
       // 	 printf("Inside rec command, before parsing 'variables=', token:'%s'(end token)\n", token);
       
       token += strlen("variables=");

//* testing
       // if (proc_zero())
       // 	 printf("Inside rec command, after parsing 'variables=', token:'%s'(end token)\n", token);

       if( strcmp("displacement",token)==0 )
       {
	 mode = TimeSeries::Displacement;
       }
       else if( strcmp("velocity",token)==0 )
       {
	 mode = TimeSeries::Velocity;
       }
       else if( strcmp("div",token)==0 )
       {
	 mode = TimeSeries::Div;
       }
       else if( strcmp("curl",token)==0 )
       {
	 mode = TimeSeries::Curl;
       }
       else if( strcmp("strains",token)==0 )
       {
	 mode = TimeSeries::Strains;
       }
       else if( strcmp("displacementgradient",token)==0 )
       {
	 mode = TimeSeries::DisplacementGradient;
       }
       else
       {
	 if (m_myrank == 0 )
	   cout << "receiver command: variables=" << token << " not understood" << endl
		<< "using default mode (displacement)" << endl << endl;
	 mode = TimeSeries::Displacement;
       }
       
     }
     else
     {
        badOption("receiver", token);
     }
     token = strtok(NULL, " \t");
//* testing
     // if (proc_zero())
     //   cout << "rec command: Bottom of while loop, token:" << token << "(end token)" << endl;
     
  }  
  //  cout << "end receiver " << m_myRank << endl;

  if (geoCoordSet)
  {
    computeCartesianCoord(x, y, lon, lat);
// check if (x,y) is within the computational domain
  }

  if (!staNameGiven)
    staName = fileName;

  bool inCurvilinear=false;
// we are in or above the curvilinear grid 
  if ( topographyExists() && z < m_zmin[mNumberOfCartesianGrids-1])
  {
    inCurvilinear = true;
  }
      
// check if (x,y,z) is not in the global bounding box
  if ( !( (inCurvilinear || z >= 0) && x>=0 && x<=m_global_xmax && y>=0 && y<=m_global_ymax))
  {
// The location of this station was outside the domain, so don't include it in the global list
    if (m_myrank == 0 && getVerbosity() > 0)
    {
      stringstream receivererr;
  
      receivererr << endl 
		  << "***************************************************" << endl
		  << " WARNING:  RECEIVER positioned outside grid!" << endl;
      receivererr << " No RECEIVER file will be generated for file = " << fileName << endl;
      if (geoCoordSet)
      {
	receivererr << " @ lon=" << lon << " lat=" << lat << " depth=" << depth << endl << endl;
      }
      else
      {
	receivererr << " @ x=" << x << " y=" << y << " z=" << z << endl << endl;
      }
      
      receivererr << "***************************************************" << endl;
      cerr << receivererr.str();
      cerr.flush();
    }
  }
  else
  {
    TimeSeries *ts_ptr = new TimeSeries(this, fileName, staName, mode, sacformat, usgsformat, x, y, depth, 
					topodepth, writeEvery, !nsew);
// include the receiver in the global list
    m_GlobalTimeSeries.push_back(ts_ptr);
  }
}

//-----------------------------------------------------------------------
void EW::defineDimensionsGXY( )
{
//
// Defines the number of grids and dimensions in the x- and y-directions,
// It also defines the parallel decomposition, which is only made in the x-y directions.
//
// The z-direction requires topography to be known before computing dimensions.
// x- and y-dimensions must be defined before the topography is read. 
// Hence, we have to 1. Define x and y dimensions, 
//                   2. Read the topography
//                   3. Define z dimensions.
   if (mVerbose && m_myrank == 0 )
      printf("defineDimensionsGXY: #ghost points=%i, #parallel padding points=%i\n", m_ghost_points, m_ppadding);

 // Grids are enumerated from bottom to the top, i.e, g=0 is at the bottom, and g=mNumberOfGrids-1 is at the top.
 // Note, this is oposite to the z-coordinate which is largest at the bottom and smallest at the top.
   if( m_nz_base > 1 && !m_topography_exists )
   {
      // Flat
      mNumberOfCartesianGrids = mNumberOfGrids = 1;
      m_is_curvilinear.push_back(false);
   }
   else if( m_nz_base > 1 && m_topography_exists )
   {
      // Curvilinear
      mNumberOfGrids = 2;
      mNumberOfCartesianGrids = 1;
      m_is_curvilinear.push_back(false);
      m_is_curvilinear.push_back(true);
   }
   else
      if( m_myrank == 0 )
	 cout << "ERROR in defineDimensionsXY, domain could not be defined" << endl;

// Compute parallel decomposition
   int nx_finest_w_ghost = m_nx_base+2*m_ghost_points;
   int ny_finest_w_ghost = m_ny_base+2*m_ghost_points;
   proc_decompose_2d( nx_finest_w_ghost, ny_finest_w_ghost, m_nprocs, m_nprocs_2d );
   int is_periodic[2]={0,0};

   MPI_Cart_create( MPI_COMM_WORLD, 2, m_nprocs_2d, is_periodic, true, &m_cartesian_communicator );
   //   int my_proc_coords[2];
   MPI_Cart_get( m_cartesian_communicator, 2, m_nprocs_2d, is_periodic, m_myrank_2d );
   MPI_Cart_shift( m_cartesian_communicator, 0, 1, m_neighbor, m_neighbor+1 );
   MPI_Cart_shift( m_cartesian_communicator, 1, 1, m_neighbor+2, m_neighbor+3 );

   if( m_myrank == 0 && mVerbose >= 3)
   {
     cout << " Grid distributed on " << m_nprocs << " processors " << endl;
     cout << " Finest grid size    " << nx_finest_w_ghost << " x " << ny_finest_w_ghost << endl;
     cout << " Processor array     " << m_nprocs_2d[0] << " x " << m_nprocs_2d[1] << endl;
   }
   int ifirst, ilast, jfirst, jlast;
   decomp1d( nx_finest_w_ghost, m_myrank_2d[0], m_nprocs_2d[0], ifirst, ilast );
   decomp1d( ny_finest_w_ghost, m_myrank_2d[1], m_nprocs_2d[1], jfirst, jlast );

   ifirst -= m_ghost_points;
   ilast  -= m_ghost_points;
   jfirst -= m_ghost_points;
   jlast  -= m_ghost_points;

   // Define dimension arrays
   mGridSize.resize(mNumberOfGrids);
   m_global_nx.resize(mNumberOfGrids);
   m_global_ny.resize(mNumberOfGrids);
   
   m_iStart.resize(mNumberOfGrids);
   m_iEnd.resize(mNumberOfGrids);
   m_jStart.resize(mNumberOfGrids);
   m_jEnd.resize(mNumberOfGrids);

   m_iStartInt.resize(mNumberOfGrids);
   m_iEndInt.resize(mNumberOfGrids);
   m_jStartInt.resize(mNumberOfGrids);
   m_jEndInt.resize(mNumberOfGrids);

 // Compute decomposition of x-y dimensions.
   for( int g = 0 ; g < mNumberOfGrids; g++ )
   {
      mGridSize[g]   = m_h_base;
      m_global_nx[g] = m_nx_base;
      m_global_ny[g] = m_ny_base;

 // save the local index bounds
      m_iStart[g] = ifirst;
      m_iEnd[g]   = ilast;
      m_jStart[g] = jfirst;
      m_jEnd[g]   = jlast;

 // local index bounds for interior points (= no ghost or parallel padding points)
      if (ifirst == 1-m_ghost_points)
	 m_iStartInt[g] = 1;
      else
	 m_iStartInt[g] = ifirst+m_ppadding;

      if (ilast == m_global_nx[g] + m_ghost_points)
	 m_iEndInt[g]   = m_global_nx[g];
      else
	 m_iEndInt[g]   = ilast - m_ppadding;

      if (jfirst == 1-m_ghost_points)
	 m_jStartInt[g] = 1;
      else
	 m_jStartInt[g] = jfirst+m_ppadding;

      if (jlast == m_global_ny[g] + m_ghost_points)
	 m_jEndInt[g]   = m_global_ny[g];
      else
	 m_jEndInt[g]   = jlast - m_ppadding;
   }

// Set up arrays of arrays.

// Materials
   mMu.resize(mNumberOfGrids);
   mLambda.resize(mNumberOfGrids);
   mRho.resize(mNumberOfGrids);

   // Super-grid data
   m_sg_dc_x.resize(mNumberOfGrids);
   m_sg_dc_y.resize(mNumberOfGrids);
   m_sg_dc_z.resize(mNumberOfGrids);
   m_sg_str_x.resize(mNumberOfGrids);
   m_sg_str_y.resize(mNumberOfGrids);
   m_sg_str_z.resize(mNumberOfGrids);
   m_sg_corner_x.resize(mNumberOfGrids);
   m_sg_corner_y.resize(mNumberOfGrids);
   m_sg_corner_z.resize(mNumberOfGrids);

   // Boundary information   
   m_onesided.resize(mNumberOfGrids);
   m_bcType.resize(mNumberOfGrids);

   // Default values
   for( int g= 0 ;g < mNumberOfGrids ; g++ )
   {
      m_onesided[g] = new int[6];
      m_bcType[g] = new boundaryConditionType[6];
      for( int side =0 ; side < 6 ; side++ )
      {
	 m_onesided[g][side] = 0;
	 m_bcType[g][side] = bProcessor;
      }
   }
}

//-----------------------------------------------------------------------
void EW::defineDimensionsZ()
{
   // Assumes that topography is known, and computes the z-direction
   // dimensions of arrays.

   // Compute average elevation 
   float_sw4 topo_avg=0;
   if( m_topography_exists )
   {
      float_sw4 tzmin, tzmax;
      compute_minmax_topography(tzmin,tzmax);
      topo_avg = 0.5*(tzmin+tzmax);
   }

   m_zmin.resize(mNumberOfGrids);
   m_global_nz.resize(mNumberOfGrids);
// Define m_zmin and m_global_nk. 
// Adjust m_global_zmin and m_global_zmax, if necessary.
   if( m_nz_base > 1 && !m_topography_exists )
   {
      // Flat 
      m_global_nz[0] = m_nz_base;
      m_zmin[0] = 0;
   }
   else if( m_nz_base > 1 && m_topography_exists )
   {
      // Curvilinear 
      int nz = static_cast<int>(1 + round((m_global_zmax-m_topo_zmax)/m_h_base));
      m_global_zmax = m_topo_zmax+(nz-1)*m_h_base;
      m_global_nz[0] = nz;
      m_zmin[0] = m_topo_zmax;

      m_global_nz[1] = static_cast<int>(1 + round((m_topo_zmax - topo_avg)/m_h_base));
      m_zmin[1] = 1e38;
   }
   else
      if( m_myrank == 0 )
	 cout << "ERROR in defineDimensionsZ, elastic domain could not be defined" << endl;
// Define local z-dimension arrays
   m_kStart.resize(mNumberOfGrids);
   m_kEnd.resize(mNumberOfGrids);
   m_kStartInt.resize(mNumberOfGrids);
   m_kEndInt.resize(mNumberOfGrids);
   for( int g = 0 ; g < mNumberOfGrids; g++ )
   {
      m_kStart[g]    = 1-m_ghost_points;
      m_kEnd[g]      = m_global_nz[g] + m_ghost_points;
      m_kStartInt[g] = 1;
      m_kEndInt[g]   = m_global_nz[g];
   }
   if (mVerbose >= 1 && m_myrank == 0)
      cout << "Extent of the computational domain xmax=" << m_global_xmax << " ymax=" << m_global_ymax << 
	 " zmin = " << m_global_zmin << " zmax=" <<  m_global_zmax << endl;
}

//-----------------------------------------------------------------------
void EW::allocateTopoArrays()
{
   if( m_topography_exists )
   {
      int ifirst = m_iStart[mNumberOfGrids-1];
      int ilast  = m_iEnd[mNumberOfGrids-1];
      int jfirst = m_jStart[mNumberOfGrids-1];
      int jlast  = m_jEnd[mNumberOfGrids-1];
// Two versions of the topography:
      mTopo.define(ifirst,ilast,jfirst,jlast,1,1); // true topography/bathymetry, read directly 
// smoothed version of true topography, with an extended number (4 instead of 2 ) of ghost points.
      m_ext_ghost_points = 2;
      mTopoGridExt.define(ifirst-m_ext_ghost_points,ilast+m_ext_ghost_points,
			  jfirst-m_ext_ghost_points,jlast+m_ext_ghost_points,1,1);
   }
}

//-----------------------------------------------------------------------
void EW::allocateArrays()
{
   for( int g=0 ; g < mNumberOfGrids ; g++ )
   {
      int ifirst = m_iStart[g];
      int ilast  = m_iEnd[g];
      int jfirst = m_jStart[g];
      int jlast  = m_jEnd[g];
      int kfirst = m_kStart[g];
      int klast  = m_kEnd[g];
    // Material data
      mMu[g].define(ifirst,ilast,jfirst,jlast,kfirst,klast);
      mRho[g].define(ifirst,ilast,jfirst,jlast,kfirst,klast);
      mLambda[g].define(ifirst,ilast,jfirst,jlast,kfirst,klast);
    // initialize the material coefficients to -1
      mMu[g].set_to_minusOne();
      mRho[g].set_to_minusOne();
      mLambda[g].set_to_minusOne();

    // Supergrid arrays
      m_sg_dc_x[g]     = newmanaged(ilast-ifirst+1);
      m_sg_dc_y[g]     = newmanaged(jlast-jfirst+1);
      m_sg_dc_z[g]     = newmanaged(klast-kfirst+1);
      m_sg_str_x[g]    = newmanaged(ilast-ifirst+1);
      m_sg_str_y[g]    = newmanaged(jlast-jfirst+1);
      m_sg_str_z[g]    = newmanaged(klast-kfirst+1);
      m_sg_corner_x[g] = newmanaged(ilast-ifirst+1);
      m_sg_corner_y[g] = newmanaged(jlast-jfirst+1);
      m_sg_corner_z[g] = newmanaged(klast-kfirst+1);
      //#pragma omp parallel for
      //      for( int k=kfirst ; k<= klast ; k++)
      //	 for( int j=jfirst ; j <= jlast ; j++ )
      //	    for( int i=ifirst ; i <= ilast ; i++ )
      //	    {
      //	       m_sg_dc_x[g][i-ifirst] =0;
      //	       m_sg_dc_y[g][j-jfirst] =0;
      //	       m_sg_dc_z[g][k-kfirst] =0;
      //	       m_sg_str_x[g][i-ifirst] =0;
      //	       m_sg_str_y[g][j-jfirst] =0;
      //	       m_sg_str_z[g][k-kfirst] =0;
      //	       m_sg_corner_x[g][i-ifirst] =0;
      //	       m_sg_corner_y[g][j-jfirst] =0;
      //	       m_sg_corner_z[g][k-kfirst] =0;
      //	    }
      if( m_topography_exists && g == mNumberOfGrids-1 )
      {
    // Grid and metric
	 mJ.define(ifirst,ilast,jfirst,jlast,kfirst,klast);
	 mX.define(ifirst,ilast,jfirst,jlast,kfirst,klast);
	 mY.define(ifirst,ilast,jfirst,jlast,kfirst,klast);
	 mZ.define(ifirst,ilast,jfirst,jlast,kfirst,klast);
	 mMetric.define(4,ifirst,ilast,jfirst,jlast,kfirst,klast);
     // Initialization, to touch memory in case OpenMP is in use
	 mJ.set_to_zero();
	 mX.set_to_zero();
	 mY.set_to_zero();
	 mZ.set_to_zero();
	 mMetric.set_to_zero();
      }
   }
}

//-----------------------------------------------------------------------
void EW::printGridSizes() const
{
   if (m_myrank == 0)
   {
      int nx, ny, nz;
      float_sw4 nTot=0.;
      printf("\nGlobal grid sizes (without ghost points)\n");
//             1234  12345679  12345679  12345679  12345679
      printf("Grid         h        Nx        Ny        Nz       Points\n");
      for (int g = 0; g < mNumberOfGrids; g++)
      {
	 nx = m_global_nx[g];
	 ny = m_global_ny[g];
	 nz = m_kEnd[g] - m_ghost_points;
	 nTot += ((long long int)nx)*ny*nz;
	 printf("%4i %9g %9i %9i %9i %12lld\n", g, mGridSize[g], nx, ny, nz, ((long long int)nx)*ny*nz);
      }
      printf("Total number of grid points (without ghost points): %g\n\n", nTot);
   }
}

//-----------------------------------------------------------------------
bool EW::parseInputFile( const string& filename )
{

   char buffer[256];
   bool foundGrid = false;
   MPI_Barrier(MPI_COMM_WORLD);

   ifstream inputFile;
   inputFile.open(filename.c_str());
   if (!inputFile.is_open())
   {
      if (m_myrank == 0)
	 cerr << endl << "ERROR: Failure opening input file: " << filename << endl;
      return false;
   }

   while (!inputFile.eof())
   {    
      inputFile.getline(buffer, 256);
      if( startswith("grid", buffer) )
      {
	 foundGrid = true;
	 processGrid(buffer);
      }
      // Need process developer before setupMPICommunication, because of array ordering m_corder
      else if(startswith("developer", buffer))
	 processDeveloper(buffer);
      else if (startswith("topography", buffer))
	 processTopography(buffer);
      else if( startswith("fileio",buffer))
	 processFileIO(buffer);
   }   
   if (!foundGrid)
      if (m_myrank == 0)
      {
	 cerr << "ERROR: No grid found in input file: " << filename << endl;
	 return false; 
      }
   defineDimensionsGXY();
   if( m_topography_exists )
   {
      allocateTopoArrays();
      if( m_topoInputStyle == EW::GaussianHill )
	 buildGaussianHillTopography(m_GaussianAmp, m_GaussianLx, m_GaussianLy, m_GaussianXc, m_GaussianYc);
   }
   defineDimensionsZ();
   setupMPICommunications();
   allocateArrays();

   if( m_topography_exists )
   {
      generate_grid();
      setup_metric();
   }

// output grid size info
   printGridSizes();

// set default boundary conditions,
   default_bcs();
   
   inputFile.clear();
   inputFile.seekg(0, ios::beg); // reset file pointer to the beginning of the input file
   while (!inputFile.eof())
   {
      inputFile.getline(buffer, 256);
      if (strlen(buffer) > 0) // empty lines produce this
      {
	 if (startswith("#", buffer) || 
	     startswith("grid", buffer) ||
             startswith("developer", buffer) ||
             startswith("topography", buffer) ||
             startswith("fileio", buffer) ||
             startswith("\n", buffer) ||
	     startswith("\r", buffer) )
	 {
	 }
	 else if(startswith("time", buffer))
	    processTime(buffer);
	 else if( startswith("source",buffer))
             processSource(buffer);
         else if( startswith("supergrid",buffer))
	    processSuperGrid(buffer);
	 else if(startswith("testpointsource", buffer))
	    processTestPointSource(buffer);
	 //	 else if(startswith("developer", buffer))
	 //	    processDeveloper(buffer);
	 else if( startswith("checkpoint",buffer))
	    processCheckPoint(buffer);
	 else if( startswith("restart",buffer))
	    processRestart(buffer);
	 else if( startswith("rec",buffer))
	    processReceiver(buffer);
	 else if( startswith("block",buffer))
	    processMaterialBlock(buffer);
         else if( startswith("dgalerkin", buffer) )
         {
             m_use_dg=true;
             if (m_myrank == 0){
                 cout << "Using DG solver" << endl;
             }
             processdGalerkin(buffer);
         }
         else if (!inputFile.eof() && m_myrank == 0)
	 {
	    cout << "*** Ignoring command: '" << buffer << "'" << endl;
	 }
      }
   }
   inputFile.close();
   if( m_myrank == 0 )
      cout << "Done reading input file " << endl;
   MPI_Barrier(MPI_COMM_WORLD);
   return true;
}

//-----------------------------------------------------------------------
void EW::setupRun()
{
// Assign values to material data arrays mRho,mMu,mLambda
   setup_materials();
// Check if any GPUs are available
   find_cuda_device( );

   m_cuobj->initialize_gpu(m_myrank);

// setup coefficients for SBP operators
   setupSBPCoeff();
// Check that f.d. operators fit inside the domains
   check_dimensions();
// Initialize IO
   create_output_directory( );
// Set up supergrid
   setup_supergrid( );
   assign_supergrid_damping_arrays();
// Copy material to GPU
   copy_material_to_device();
// B.C. data structures
   assign_local_bcs(); 
   setup_boundary_arrays();
// Time step
   computeDT( );
// Set up sources:
   for( int s=0 ; s < m_globalUniqueSources.size() ; s++)
   {
      m_globalUniqueSources[s]->set_grid_point_sources4( this, m_point_sources );
   }
  // Sorting sources on grid index will allow more efficient parallel code with multi-core
   sort_grid_point_sources();
   if( m_myrank == 0 && m_globalUniqueSources.size() > 0 )
      cout << "setup of sources done" << endl;

   if( m_cuobj->has_gpu() )
   {
      copy_point_sources_to_gpu( );
      init_point_sourcesCU( );
   }
// Setup I/O in check points
   if( m_restart_check_point != CheckPoint::nil )
      m_restart_check_point->setup_sizes();
   for( int c = 0 ; c < m_check_points.size() ; c++ )
      m_check_points[c]->setup_sizes();
   if( m_myrank == 0 && (m_restart_check_point != CheckPoint::nil || m_check_points.size() > 0) )
      cout << "setup of check point file done" << endl;
}

//-----------------------------------------------------------------------
void EW::timesteploop( vector<Sarray>& U, vector<Sarray>& Um )
{
   // input: U,Um,mMu,mLambda,mRho,

   // local arrays: F, Up, Lu, Uacc
   vector<Sarray> F, Lu, Uacc, Up;
   Sarray* dev_F;

   // Do all timing in double, time differences have to much cancellation for float.
   double time_start_solve = MPI_Wtime();
   bool saveerror = false;

   // Define local arrays
   F.resize(mNumberOfGrids);
   Lu.resize(mNumberOfGrids);
   Uacc.resize(mNumberOfGrids);
   Up.resize(mNumberOfGrids);
   U.resize(mNumberOfGrids);
   Um.resize(mNumberOfGrids);
   for( int g=0 ; g < mNumberOfGrids ; g++ )
   {
      int ifirst = m_iStart[g], ilast = m_iEnd[g];
      int jfirst = m_jStart[g], jlast = m_jEnd[g];
      int kfirst = m_kStart[g], klast = m_kEnd[g];
      F[g].define(3,ifirst,ilast,jfirst,jlast,kfirst,klast);
      Lu[g].define(3,ifirst,ilast,jfirst,jlast,kfirst,klast);
      Uacc[g].define(3,ifirst,ilast,jfirst,jlast,kfirst,klast);
      Up[g].define(3,ifirst,ilast,jfirst,jlast,kfirst,klast);
      U[g].define(3,ifirst,ilast,jfirst,jlast,kfirst,klast);
      Um[g].define(3,ifirst,ilast,jfirst,jlast,kfirst,klast);
   }

   // Set up boundary data array
   vector<float_sw4**> BCForcing;
   BCForcing.resize(mNumberOfGrids);
   for( int g = 0; g <mNumberOfGrids; g++ )
   {
      BCForcing[g] = new float_sw4*[6];
      for (int side=0; side < 6; side++)
      {
	 BCForcing[g][side]=NULL;
	 if (m_bcType[g][side] == bStressFree || m_bcType[g][side] == bDirichlet || m_bcType[g][side] == bSuperGrid)
	 {
	   BCForcing[g][side] = newmanaged(3*m_NumberOfBCPoints[g][side]);
	 }
      }
   }
   
   // Initial data, touch all memory even in
   // arrays that do not need values, in order
   // to initialize OpenMP with good memory access
   for( int g=0 ; g < mNumberOfGrids ; g++ )
   {
      U[g].set_value(0.0);
      Um[g].set_value(0.0);
      F[g].set_value(0.0);
      Up[g].set_value(0.0);
      Uacc[g].set_value(0.0);
      Lu[g].set_value(0.0);
   }

   int beginCycle = 0;
   float_sw4 t = mTstart;
   if( m_restart_check_point != CheckPoint::nil )
   {
      m_restart_check_point->read_checkpoint( t, beginCycle, Um, U );
      for(int g=0 ; g < mNumberOfGrids ; g++ )
      {
	 communicate_array( U[g], g );
	 communicate_array( Um[g], g );
      }
      cartesian_bc_forcing( t, BCForcing, m_globalUniqueSources );
      enforceBC( U, mMu, mLambda, t, BCForcing );
      cartesian_bc_forcing( t-mDt, BCForcing, m_globalUniqueSources );
      enforceBC( Um, mMu, mLambda, t-mDt, BCForcing );
   }
   beginCycle++;

   double time_measure[20];
   double time_sum[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};


   for (int ts=0; ts<m_GlobalTimeSeries.size(); ts++)
      m_GlobalTimeSeries[ts]->allocateRecordingArrays( mNumberOfTimeSteps+1, mTstart, mDt);

   if (m_myrank == 0)
   {
      cout << "Running on " << m_nprocs << " MPI tasks" << endl;
   }
   
#ifdef SW4_OPENMP
#pragma omp parallel
   {
      if( omp_get_thread_num() == 0 &&  m_myrank == 0  )
      {
	 int nth=omp_get_num_threads();
	 cout << "Using OpenMP with " << nth << " thread";
	 if( nth > 1 )
	    cout << "s";
	 cout << " per MPI task" << endl;
      }
   }
#endif
#ifdef SW4_CUDA
   cudaError_t retval = cudaMalloc( (void**)&dev_F, sizeof(Sarray)*mNumberOfGrids);
   for( int g=0 ; g < mNumberOfGrids ; g++ )
   {
      //      Lu[g].allocate_on_device(m_cuobj);
      //      Up[g].allocate_on_device(m_cuobj);
      Lu[g].copy_to_device(m_cuobj);
      Up[g].copy_to_device(m_cuobj);
      Um[g].copy_to_device(m_cuobj);
      Uacc[g].copy_to_device(m_cuobj);
      F[g].copy_to_device(m_cuobj);
      F[g].page_lock(m_cuobj);
      U[g].page_lock(m_cuobj);
      Um[g].page_lock(m_cuobj);
      Up[g].page_lock(m_cuobj);
      //      Sarray* har = F[g].create_copy_on_device(m_cuobj);
   }
   retval = cudaMemcpy( dev_F, &F[0], mNumberOfGrids*sizeof(Sarray), cudaMemcpyHostToDevice );
   //   retval = cudaMemcpy( &dev_F[g], har, sizeof(Sarray), cudaMemcpyHostToDevice );
      if( retval != cudaSuccess )
	 cout << "Error in memcpy to dev_F in timestep loop retval = " <<
	    cudaGetErrorString(retval) << endl;
      //}
#endif

// save initial data on receiver records
  vector<float_sw4> uRec;
  for (int ts=0; ts<m_GlobalTimeSeries.size(); ts++)
  {
// can't compute a 2nd order accurate time derivative at this point
// therefore, don't record anything related to velocities for the initial data
    if (m_GlobalTimeSeries[ts]->getMode() != TimeSeries::Velocity && m_GlobalTimeSeries[ts]->myPoint())
    {
      int i0 = m_GlobalTimeSeries[ts]->m_i0;
      int j0 = m_GlobalTimeSeries[ts]->m_j0;
      int k0 = m_GlobalTimeSeries[ts]->m_k0;
      int grid0 = m_GlobalTimeSeries[ts]->m_grid0;

      extractRecordData(m_GlobalTimeSeries[ts]->getMode(), i0, j0, k0, grid0, 
			uRec, Um, U); 
      m_GlobalTimeSeries[ts]->recordData(uRec);
    }
  }

   if( m_myrank == 0 )
      cout << "starting at time " << t << " at cycle " << beginCycle << endl;


   double* trdata;
   if( m_save_trace )
   {
      trdata = new double[12*(mNumberOfTimeSteps+1)];
      MPI_Barrier(m_cartesian_communicator);
   }
// Begin time stepping loop
   for( int currentTimeStep = beginCycle; currentTimeStep <= mNumberOfTimeSteps; currentTimeStep++ )
   {    
     PUSH_RANGE("TIME STEP",0);
      time_measure[0] = MPI_Wtime();
      // Predictor 
      // Need U on device for evalRHS,
      for( int g=0; g < mNumberOfGrids ; g++ )
	 U[g].copy_to_device(m_cuobj,true,0);

// all types of forcing...
      if( m_cuobj->has_gpu() )
	 ForceCU( t, dev_F, false, 1 );
      else
	 ForceOffload( t, F, m_point_sources, false );
 // Need F on device for predictor, will make this asynchronous:
      //      for( int g=0; g < mNumberOfGrids ; g++ )
      //	 F[g].copy_to_device(m_cuobj,true,1);

      if( m_checkfornan )
      {
	//std::cout<<"HEREEEEEEEEE\n";
	check_for_nan( F, 1, "F" );
	check_for_nan( U, 1, "U" );
	//std::cout<<"HEREEEEEEEEE\n";
      }
      time_measure[1] = MPI_Wtime();

// evaluate right hand side
      PUSH_RANGE("evalRHS",2);
      if( m_cuobj->has_gpu() )
	 evalRHSCU( U, mMu, mLambda, Lu, 0 ); // save Lu in composite grid 'Lu'
      else
	 evalRHS( U, mMu, mLambda, Lu ); // save Lu in composite grid 'Lu'
      POP_RANGE;
      if( m_checkfornan )
#ifdef SW4_CUDA
	 check_for_nan_GPU( Lu, 1, "Lu pred. " );
#else      
	 check_for_nan( Lu, 1, "Lu pred. " );
#endif

// take predictor step, store in Up
      m_cuobj->sync_stream( 0 );
      if( m_cuobj->has_gpu() )
	 evalPredictorCU( Up, U, Um, mRho, Lu, F, 1 );    
      else
	 evalPredictor( Up, U, Um, mRho, Lu, F );    

      for( int g=0; g < mNumberOfGrids ; g++ )
      {
	 Up[g].copy_from_device(m_cuobj,true,1);
      }
      m_cuobj->sync_stream(1);

      time_measure[2] = MPI_Wtime();

// communicate across processor boundaries
      for(int g=0 ; g < mNumberOfGrids ; g++ )
	 communicate_array( Up[g], g );

      time_measure[3] = MPI_Wtime();

// calculate boundary forcing at time t+mDt
      PUSH_RANGE("CART-BC",5);
      cartesian_bc_forcing( t+mDt, BCForcing, m_globalUniqueSources );
      POP_RANGE;
      PUSH_RANGE("enforceBC",0);
      enforceBC( Up, mMu, mLambda, t+mDt, BCForcing );
      POP_RANGE;
      PUSH_RANGE("CHK-NAN",1);
      if( m_checkfornan )
	 check_for_nan( Up, 1, "U pred. " );
      POP_RANGE;

      //      time_measure[3] = MPI_Wtime();
      time_measure[4] = MPI_Wtime();

      for( int g=0; g < mNumberOfGrids ; g++ )
	 Up[g].copy_to_device(m_cuobj,true,0);

      // Corrector
      PUSH_RANGE("ForceOffload",3);
      if( m_cuobj->has_gpu() )
	 ForceCU( t, dev_F, true, 1 );
      else
	ForceOffload( t, F, m_point_sources, true );
      POP_RANGE;
      //      for( int g=0; g < mNumberOfGrids ; g++ )
      //	 F[g].copy_to_device(m_cuobj,true,1);

      //      time_measure[4] = MPI_Wtime();
      time_measure[5] = MPI_Wtime();

      PUSH_RANGE("evalDpDmInTIme",0);
      if( m_cuobj->has_gpu() )
	 evalDpDmInTimeCU( Up, U, Um, Uacc, 0 ); // store result in Uacc
      else
	 evalDpDmInTime( Up, U, Um, Uacc ); // store result in Uacc
      POP_RANGE;
      if( m_checkfornan )
#ifdef SW4_CUDA
	 check_for_nan_GPU( Uacc, 1, "uacc " );
#else
	 check_for_nan( Uacc, 1, "uacc " );
#endif      
      PUSH_RANGE("evalRHS-2",1);
      if( m_cuobj->has_gpu() )
	 evalRHSCU( Uacc, mMu, mLambda, Lu, 0 );
      else
       	 evalRHS( Uacc, mMu, mLambda, Lu );
      POP_RANGE;
      if( m_checkfornan )
#ifdef SW4_CUDA
	 check_for_nan_GPU( Lu, 1, "L(uacc) " );
#else
	 check_for_nan( Lu, 1, "L(uacc) " );
#endif

      m_cuobj->sync_stream(0);
      if( m_cuobj->has_gpu() )
	 evalCorrectorCU( Up, mRho, Lu, F, 1 );
      else
	 evalCorrector( Up, mRho, Lu, F );
      //      time_measure[5] = MPI_Wtime();
      time_measure[6] = MPI_Wtime();

      // add in super-grid damping terms
      PUSH_RANGE("addSuperGridDamping",3);
      if ( m_use_supergrid )
      {
	 if( m_cuobj->has_gpu() )
	    addSuperGridDampingCU( Up, U, Um, mRho, 1 );
	 else
	    addSuperGridDamping( Up, U, Um, mRho );

      }
      POP_RANGE;
      for( int g=0; g < mNumberOfGrids ; g++ )
	 Up[g].copy_from_device(m_cuobj,true,1);

      m_cuobj->sync_stream(1);

      //      time_measure[6] = MPI_Wtime();
      time_measure[7] = MPI_Wtime();

// also check out EW::update_all_boundaries 
// communicate across processor boundaries
      for(int g=0 ; g < mNumberOfGrids ; g++ )
	 communicate_array( Up[g], g );

      time_measure[8] = MPI_Wtime();

// calculate boundary forcing at time t+mDt (do we really need to call this fcn again???)
      PUSH_RANGE("CART-BC2",5);
      cartesian_bc_forcing( t+mDt, BCForcing, m_globalUniqueSources );
      POP_RANGE;
      PUSH_RANGE("enforceBC2",0);
      enforceBC( Up, mMu, mLambda, t+mDt, BCForcing );
      POP_RANGE;
      if( m_checkfornan )
	 check_for_nan( Up, 1, "Up" );

// increment time
      t += mDt;

      //      time_measure[7] = MPI_Wtime();	  
      time_measure[9] = MPI_Wtime();	  

// periodically, print time stepping info to stdout
      printTime( currentTimeStep, t, currentTimeStep == mNumberOfTimeSteps ); 
// Images have to be written before the solution arrays are cycled, because both Up and Um are needed
// to compute a centered time derivative
//
      double time_chkpt, time_chkpt_tmp;
      bool wrote=false;
      //      MPI_Barrier( MPI_COMM_WORLD );
      time_chkpt=MPI_Wtime();
      for( int c=0 ; c < m_check_points.size() ; c++ )
	 if( m_check_points[c]->timeToWrite( t, currentTimeStep, mDt) )
	 {
	    m_check_points[c]->write_checkpoint( t, currentTimeStep, U, Up );
	    wrote=true;
	 }
      time_chkpt_tmp =MPI_Wtime()-time_chkpt;
      MPI_Allreduce( &time_chkpt_tmp, &time_chkpt, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
      if( m_myrank == 0 && wrote )
	 cout << "Cpu time to write check point file " << time_chkpt << " seconds " << endl;
//
// save the current solution on receiver records (time-derivative require Up and Um for a 2nd order
// approximation, so do this before cycling the arrays)
      for (int ts=0; ts<m_GlobalTimeSeries.size(); ts++)
      {
	 if (m_GlobalTimeSeries[ts]->myPoint())
	 {
	    int i0 = m_GlobalTimeSeries[ts]->m_i0;
	    int j0 = m_GlobalTimeSeries[ts]->m_j0;
	    int k0 = m_GlobalTimeSeries[ts]->m_k0;
	    int grid0 = m_GlobalTimeSeries[ts]->m_grid0;
//
// note that the solution on the new time step is in Up
// also note that all quantities related to velocities lag by one time step; they are not
// saved before the time stepping loop started
	    extractRecordData(m_GlobalTimeSeries[ts]->getMode(), i0, j0, k0, grid0, 
			      uRec, Um, Up);
	    m_GlobalTimeSeries[ts]->recordData(uRec);
	 }
      }

// // Energy evaluation, requires all three time levels present, do before cycle arrays.
//      if( m_energy_test )
//	 compute_energy( mDt, currentTimeStep == mNumberOfTimeSteps, Um, U, Up, currentTimeStep  );

// cycle the solution arrays
      cycleSolutionArrays(Um, U, Up);

      //      time_measure[8] = MPI_Wtime();	  
      time_measure[10] = MPI_Wtime();	  
// evaluate error for some test cases
//      if (m_lamb_test || m_point_source_test || m_rayleigh_wave_test )
      if ( m_point_source_test && saveerror )
      {
	 float_sw4 errInf=0, errL2=0, solInf=0; //, solL2=0;
	 exactSol( t, Up, m_globalUniqueSources ); // store exact solution in Up
//	 //	 if (m_lamb_test)
//	 //	    normOfSurfaceDifference( Up, U, errInf, errL2, solInf, solL2, a_Sources);
	 PUSH_RANGE("NormOfDiff",1);
	 normOfDifference( Up, U, errInf, errL2, solInf, m_globalUniqueSources );
	 POP_RANGE;
         if ( m_myrank == 0 )
	    cout << t << " " << errInf << " " << errL2 << " " << solInf << endl;
      }
      //      time_measure[9] = MPI_Wtime();	  	
      time_measure[11] = MPI_Wtime();	  
// // See if it is time to write a restart file
// //      if (mRestartDumpInterval > 0 &&  currentTimeStep % mRestartDumpInterval == 0)
// //        serialize(currentTimeStep, U, Um);  
      if( currentTimeStep > 1 )
      {
	 time_sum[0] += time_measure[1]-time_measure[0] + time_measure[5]-time_measure[4]; // F
	 time_sum[1] += time_measure[2]-time_measure[1] + time_measure[6]-time_measure[5]; // RHS
	 time_sum[2] += time_measure[3]-time_measure[2] + time_measure[8]-time_measure[7]; // bc comm.
	 time_sum[3] += time_measure[4]-time_measure[3] + time_measure[9]-time_measure[8]; // bc phys.
	 time_sum[4] += time_measure[7]-time_measure[6]; // super grid damping
	 time_sum[5] += time_measure[10]-time_measure[9]; // print outs
	 time_sum[6] += time_measure[11]-time_measure[10]; //  compute exact solution
	 time_sum[7] += time_measure[11]-time_measure[0]; // total measured
      }
      if( m_save_trace )
	 for( int s = 0 ; s < 12 ; s++ )
	    trdata[s+12*(currentTimeStep-beginCycle)]= time_measure[s];
      POP_RANGE;
   } // end time stepping loop
   double time_end_solve = MPI_Wtime();
   print_execution_time( time_start_solve, time_end_solve, "solver phase" );
   if( m_output_detailed_timing )
      print_execution_times( time_sum );

   if ( m_point_source_test )
   {
      float_sw4 errInf=0, errL2=0, solInf=0;//, solL2=0;
      exactSol( t, Up, m_globalUniqueSources ); // store exact solution in Up
//	 //	 if (m_lamb_test)
//	 //	    normOfSurfaceDifference( Up, U, errInf, errL2, solInf, solL2, a_Sources);
      normOfDifference( Up, U, errInf, errL2, solInf, m_globalUniqueSources );
      //      std::cout<<"VALUES"<<Up[10]<<"\n";
      if ( m_myrank == 0 )
	 cout << "Errors at time " << t << " Linf = " << errInf << " L2 = " << errL2 << " noorm of solution = " << solInf << endl;
   }
   for (int ts=0; ts<m_GlobalTimeSeries.size(); ts++)
      m_GlobalTimeSeries[ts]->writeFile();

   for( int g= 0 ; g < mNumberOfGrids ; g++ )
   {
      F[g].page_unlock(m_cuobj);
      U[g].page_unlock(m_cuobj);
      Um[g].page_unlock(m_cuobj);
      Up[g].page_unlock(m_cuobj);
   }
   m_cuobj->reset_gpu();
   if( m_save_trace )
   {
      char fname[255];
      snprintf(fname,255,"%s/trfile%04d.bin",mPath.c_str(),m_myrank);
      int fd = open(fname, O_WRONLY|O_TRUNC|O_CREAT, 0660);
      int twelve=12;
      int nsteps= mNumberOfTimeSteps-beginCycle+1;
      size_t nr=write(fd,&twelve,sizeof(int));
      nr=write(fd,&nsteps,sizeof(int));
      nr=write(fd,trdata,sizeof(double)*twelve*nsteps);
      close(fd);
   }  
}

//-----------------------------------------------------------------------
bool EW::proc_decompose_2d( int ni, int nj, int nproc, int proc_max[2] )
{
   // This routine determines a decomposition of nproc processors into
   // a 2D processor array  proc_max[0] x proc_max[1], which gives minimal 
   // communication boundary for a grid with ni x nj points.

   float_sw4 fmin = ni+nj;
   bool first  = true;
   int p1max   = ni/m_ppadding;
   int p2max   = nj/m_ppadding;
   for( int p1 = 1 ; p1 <= nproc; p1++)
      if( nproc%p1 == 0 )
      {
        int p2 = nproc/p1;
        if( p1 <= p1max && p2 <= p2max )
        {
           // int w1 = p1==1?0:1;
           // int w2 = p2==1?0:1;
           // float_sw4 f = w2*(float_sw4)(ni)/p1 + w1*(float_sw4)(nj)/p2;
// try to make each subdomain as square as possible
	  float_sw4 f = fabs((float_sw4)(ni)/p1 - (float_sw4)(nj)/p2);
           if( f < fmin || first )
           {
              fmin = f;
              proc_max[0]   = p1;
              proc_max[1]   = p2;
              first= false;
           }
        }
      }
   return !first;
}

//-----------------------------------------------------------------------
void EW::decomp1d( int nglobal, int myid, int nproc, int& s, int& e )
//
// Decompose index space 1 <= i <= nglobal into nproc blocks
// returns start and end indices for block nr. myid, 
//          where 0 <= myid <= nproc-1
//
{
   int olap    = 2*m_ppadding;
   int nlocal  = (nglobal + (nproc-1)*olap ) / nproc;
   int deficit = (nglobal + (nproc-1)*olap ) % nproc;

   if( myid < deficit )
      s = myid*(nlocal-olap) + myid+1;
   else
      s = myid*(nlocal-olap) + deficit+1;

   if (myid < deficit)
      nlocal = nlocal + 1;

   e = s + nlocal - 1;
}

//-----------------------------------------------------------------------
void EW::setupMPICommunications()
{
   if (mVerbose >= 1 && m_myrank == 0 )
      cout << "***inside setupMPICommunications***"<< endl;

// Define MPI datatypes for communication across processor boundaries
   m_send_type1.resize(2*mNumberOfGrids);
   m_send_type3.resize(2*mNumberOfGrids);
   m_send_type4.resize(2*mNumberOfGrids);

   send_type1.resize(2*mNumberOfGrids);
   send_type3.resize(2*mNumberOfGrids);
   send_type4.resize(2*mNumberOfGrids);

   bufs_type1.resize(4*mNumberOfGrids);
   bufs_type3.resize(4*mNumberOfGrids);
   bufs_type4.resize(4*mNumberOfGrids);

   //   m_send_type21.resize(2*mNumberOfGrids);
   for( int g= 0 ; g < mNumberOfGrids ; g++ )
   {
//      int ni = mU[g].m_ni, nj=mU[g].m_nj, nk=mU[g].m_nk;
      int ni = m_iEnd[g] - m_iStart[g] + 1;
      int nj = m_jEnd[g] - m_jStart[g] + 1;
      int nk = m_kEnd[g] - m_kStart[g] + 1;

      MPI_Type_vector( nj*nk, m_ppadding, ni, m_mpifloat, &m_send_type1[2*g] );
      MPI_Type_vector( nk, m_ppadding*ni, ni*nj, m_mpifloat, &m_send_type1[2*g+1] );


      send_type1[2*g]=std::make_tuple(nj*nk, m_ppadding, ni);
      send_type1[2*g+1]=std::make_tuple(nk, m_ppadding*ni, ni*nj);
      
      bufs_type1[4*g+0]=std::make_tuple(newmanagedh(nj*nk*m_ppadding),newmanagedh(nj*nk*m_ppadding));
      bufs_type1[4*g+1]=std::make_tuple(newmanagedh(nj*nk*m_ppadding),newmanagedh(nj*nk*m_ppadding));
      
      bufs_type1[4*g+2]=std::make_tuple(newmanagedh(nk*m_ppadding*ni),newmanagedh(nk*m_ppadding*ni));
      bufs_type1[4*g+3]=std::make_tuple(newmanagedh(nk*m_ppadding*ni),newmanagedh(nk*m_ppadding*ni));
      if( m_corder )
      {
	 MPI_Type_vector( 3*nj*nk, m_ppadding, ni, m_mpifloat, &m_send_type3[2*g] );
	 MPI_Type_vector( 3*nk, m_ppadding*ni, ni*nj, m_mpifloat, &m_send_type3[2*g+1] );
	 MPI_Type_vector( 4*nj*nk, m_ppadding, ni, m_mpifloat, &m_send_type4[2*g] );
	 MPI_Type_vector( 4*nk, m_ppadding*ni, ni*nj, m_mpifloat, &m_send_type4[2*g+1] );


	 send_type3[2*g]=std::make_tuple(3*nj*nk, m_ppadding, ni);
	 send_type3[2*g+1]=std::make_tuple(3*nk, m_ppadding*ni, ni*nj);
	 send_type4[2*g]=std::make_tuple(4*nj*nk, m_ppadding, ni);
	 send_type4[2*g+1]=std::make_tuple(4*nk, m_ppadding*ni, ni*nj);

	 bufs_type3[4*g+0]=std::make_tuple(newmanagedh(3*nj*nk*m_ppadding),newmanagedh(3*nj*nk*m_ppadding));
	 bufs_type3[4*g+1]=std::make_tuple(newmanagedh(3*nj*nk*m_ppadding),newmanagedh(3*nj*nk*m_ppadding));

	 bufs_type3[4*g+2]=std::make_tuple(newmanagedh(3*nk*m_ppadding*ni),newmanagedh(3*nk*m_ppadding*ni));
	 bufs_type3[4*g+3]=std::make_tuple(newmanagedh(3*nk*m_ppadding*ni),newmanagedh(3*nk*m_ppadding*ni));
	 

	 bufs_type4[4*g+0]=std::make_tuple(newmanagedh(4*nj*nk*m_ppadding),newmanagedh(4*nj*nk*m_ppadding));
	 bufs_type4[4*g+1]=std::make_tuple(newmanagedh(4*nj*nk*m_ppadding),newmanagedh(4*nj*nk*m_ppadding));

	 bufs_type4[4*g+2]=std::make_tuple(newmanagedh(4*nk*m_ppadding*ni),newmanagedh(4*nk*m_ppadding*ni));
	 bufs_type4[4*g+3]=std::make_tuple(newmanagedh(4*nk*m_ppadding*ni),newmanagedh(4*nk*m_ppadding*ni));
      }
      else
      {
	 MPI_Type_vector( nj*nk, 3*m_ppadding, 3*ni, m_mpifloat, &m_send_type3[2*g] );
	 MPI_Type_vector( nk, 3*m_ppadding*ni, 3*ni*nj, m_mpifloat, &m_send_type3[2*g+1] );
	 MPI_Type_vector( nj*nk, 4*m_ppadding, 4*ni, m_mpifloat, &m_send_type4[2*g] );
	 MPI_Type_vector( nk, 4*m_ppadding*ni, 4*ni*nj, m_mpifloat, &m_send_type4[2*g+1] );
      }

      //      MPI_Type_vector( nj*nk, 21*m_ppadding, 21*ni, MPI_DOUBLE, &m_send_type21[2*g] );
      //      MPI_Type_vector( nk, 21*m_ppadding*ni, 21*ni*nj, MPI_DOUBLE, &m_send_type21[2*g+1] );

      MPI_Type_commit( &m_send_type1[2*g] ); 
      MPI_Type_commit( &m_send_type1[2*g+1] ); 

      MPI_Type_commit( &m_send_type3[2*g] ); 
      MPI_Type_commit( &m_send_type3[2*g+1] ); 

      MPI_Type_commit( &m_send_type4[2*g] ); 
      MPI_Type_commit( &m_send_type4[2*g+1] ); 

      //      MPI_Type_commit( &m_send_type21[2*g] ); 
      //      MPI_Type_commit( &m_send_type21[2*g+1] ); 
   }
}

//-----------------------------------------------------------------------
bool EW::check_for_nan( vector<Sarray>& a_U, int verbose, string name )
{
   bool retval = false;
   for( int g=0 ; g<mNumberOfGrids; g++ )
   {
      size_t nn=a_U[g].count_nans();
      retval = retval || nn > 0;
      if( nn > 0 && verbose == 1 )
      {
	 int cnan, inan, jnan, knan;
	 a_U[g].count_nans(cnan,inan,jnan,knan);
	 cout << "grid " << g << " array " << name << " found " << nn << "  nans. First nan at " <<
	    cnan << " " << inan << " " << jnan << " " << knan << endl;
      }
   }
   return retval;
}

//-----------------------------------------------------------------------
void EW::cycleSolutionArrays(vector<Sarray> & a_Um, vector<Sarray> & a_U,
			     vector<Sarray> & a_Up ) 
{
   for (int g=0; g<mNumberOfGrids; g++)
   {
      float_sw4 *tmp = a_Um[g].c_ptr();
      a_Um[g].reference(a_U[g].c_ptr());
      a_U[g].reference(a_Up[g].c_ptr());
      a_Up[g].reference(tmp);
      if( m_cuobj->has_gpu() )
      {
	 tmp = a_Um[g].dev_ptr();
	 a_Um[g].reference_dev( a_U[g].dev_ptr());
	 a_U[g].reference_dev( a_Up[g].dev_ptr());
	 a_Up[g].reference_dev(tmp );
      }
   }
}

//-----------------------------------------------------------------------
void EW::Force(float_sw4 a_t, vector<Sarray> & a_F, vector<GridPointSource*> point_sources,
	       bool tt )
{
  for( int g =0 ; g < mNumberOfGrids ; g++ )
     a_F[g].set_to_zero();

#pragma omp parallel for
  for( int r=0 ; r<m_identsources.size()-1 ; r++ )
  {
     int s0 = m_identsources[r];
     int g = point_sources[s0]->m_grid;
     int i = point_sources[s0]->m_i0;
     int j = point_sources[s0]->m_j0;
     int k = point_sources[s0]->m_k0;
     size_t ind1 = a_F[g].index(1,i,j,k);
     //     size_t ind2 = a_F[g].index(2,i,j,k);
     //     size_t ind3 = a_F[g].index(3,i,j,k);
     size_t oc = a_F[g].m_offc;
     float_sw4* fptr =a_F[g].c_ptr();
     for( int s=m_identsources[r]; s< m_identsources[r+1] ; s++ )
  //  for( int s = 0 ; s < point_sources.size() ; s++ )
     {
	float_sw4 fxyz[3];
	if( tt )
	   point_sources[s]->getFxyztt(a_t,fxyz);
	else
	   point_sources[s]->getFxyz(a_t,fxyz);
	fptr[ind1]      += fxyz[0];
	fptr[ind1+oc]   += fxyz[1];
	fptr[ind1+2*oc] += fxyz[2];
	//	a_F[g](1,i,j,k) += fxyz[0];
	//	a_F[g](2,i,j,k) += fxyz[1];
	//	a_F[g](3,i,j,k) += fxyz[2];
	//	a_F[g](1,point_sources[s]->m_i0,point_sources[s]->m_j0,point_sources[s]->m_k0) += fxyz[0];
	//	a_F[g](2,point_sources[s]->m_i0,point_sources[s]->m_j0,point_sources[s]->m_k0) += fxyz[1];
	//	a_F[g](3,point_sources[s]->m_i0,point_sources[s]->m_j0,point_sources[s]->m_k0) += fxyz[2];
     }
  }
}
//-----------------------------------------------------------------------
void EW::ForceOffload(float_sw4 a_t, vector<Sarray> & a_F, vector<GridPointSource*> point_sources,
	       bool tt )
{
  static int firstcall=1;
  if (firstcall){
#pragma omp parallel for
    for( int r=0 ; r<m_identsources.size()-1 ; r++ )
      {
	int index=r*3;
	int s0 = m_identsources[r];
	int g = point_sources[s0]->m_grid;
	int i = point_sources[s0]->m_i0;
	int j = point_sources[s0]->m_j0;
	int k = point_sources[s0]->m_k0;
	size_t ind1 = a_F[g].index(1,i,j,k);
	//     size_t ind2 = a_F[g].index(2,i,j,k);
	//     size_t ind3 = a_F[g].index(3,i,j,k);
	size_t oc = a_F[g].m_offc;
	float_sw4* fptr =a_F[g].c_ptr();
	ForceAddress[index]=fptr+ind1;
	ForceAddress[index+1]=fptr+ind1+oc;
	ForceAddress[index+2]=fptr+ind1+2*oc;
	//PrintPointerAttributes(fptr+ind1);
	//PrintPointerAttributes(fptr+ind1+oc);
	//PrintPointerAttributes(fptr+ind1+2*oc);
      }
    firstcall=0;
  } // First call only ends
  
  for( int g =0 ; g < mNumberOfGrids ; g++ )
    a_F[g].set_to_zero();
  
#pragma omp parallel for
  for( int r=0 ; r<m_identsources.size()-1 ; r++ )
    {
      int index=r*3;
      int s0 = m_identsources[r];
      int g = point_sources[s0]->m_grid;
      int i = point_sources[s0]->m_i0;
      int j = point_sources[s0]->m_j0;
      int k = point_sources[s0]->m_k0;
      size_t ind1 = a_F[g].index(1,i,j,k);
      //     size_t ind2 = a_F[g].index(2,i,j,k);
      //     size_t ind3 = a_F[g].index(3,i,j,k);
      size_t oc = a_F[g].m_offc;
      float_sw4* fptr =a_F[g].c_ptr();
      for (int i=0;i<3;i++)
	ForceVector[index+i]=0.0;
      for( int s=m_identsources[r]; s< m_identsources[r+1] ; s++ )
	//  for( int s = 0 ; s < point_sources.size() ; s++ )
	{
	  float_sw4 fxyz[3];
	  if( tt )
	    point_sources[s]->getFxyztt(a_t,fxyz);
	  else
	    point_sources[s]->getFxyz(a_t,fxyz);
	  for (int i=0;i<3;i++)
	    ForceVector[index+i]+=fxyz[i];
	  //fptr[ind1]      += fxyz[0];
	  //fptr[ind1+oc]   += fxyz[1];
	  //fptr[ind1+2*oc] += fxyz[2];
	  //	a_F[g](1,i,j,k) += fxyz[0];
	  //	a_F[g](2,i,j,k) += fxyz[1];
	  //	a_F[g](3,i,j,k) += fxyz[2];
	  //	a_F[g](1,point_sources[s]->m_i0,point_sources[s]->m_j0,point_sources[s]->m_k0) += fxyz[0];
	  //	a_F[g](2,point_sources[s]->m_i0,point_sources[s]->m_j0,point_sources[s]->m_k0) += fxyz[1];
	  //	a_F[g](3,point_sources[s]->m_i0,point_sources[s]->m_j0,point_sources[s]->m_k0) += fxyz[2];
	}
    }
  //#pragma omp parallel for
#ifdef CUDA_CODE
  for( int r=0 ; r<m_identsources.size()-1 ; r++ ){
    int index=r*3;
      //float_sw4* fptr =a_F[g].c_ptr();
    for(int i=0;i<3;i++){
      if (!IsManaged(ForceAddress[index+i])) std::cerr<<" NO MANAGED1\n";
      if (!IsManaged(ForceVector+index+i)) std::cerr<<" NOT MANAGED2\n";
    }
  }
#endif
  // Need the lines below becuase the object is not in managed memory and the this pointer is host only
  float_sw4 *ForceVector_copy=ForceVector;
  float_sw4 **ForceAddress_copy=ForceAddress;
  PREFETCH(ForceVector);
  PREFETCH(ForceAddress);
  forall<EXEC > (0,m_identsources.size()-1,[=] RAJA_DEVICE(int r)
    {
      int index=r*3;
      //float_sw4* fptr =a_F[g].c_ptr();
      for(int i=0;i<3;i++)
	*ForceAddress_copy[index+i]+=ForceVector_copy[index+i]; // DOes this need to be an update and not an assignment
    });
	  
}

//---------------------------------------------------------------------------
void EW::evalPredictor(vector<Sarray> & a_Up, vector<Sarray> & a_U, vector<Sarray> & a_Um,
		       vector<Sarray>& a_Rho, vector<Sarray> & a_Lu, vector<Sarray> & a_F )
{
   float_sw4 dt2 = mDt*mDt;
   for( int g=0 ; g<mNumberOfGrids; g++ )
   {
     a_Up[g].prefetch();
     a_U[g].prefetch();
     a_Um[g].prefetch();
     a_Rho[g].prefetch();
     a_Lu[g].prefetch();
     a_F[g].prefetch();
#ifdef SW4_CROUTINES
      predfort( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g], 
		a_Up[g].c_ptr(), a_U[g].c_ptr(), a_Um[g].c_ptr(),
		a_Lu[g].c_ptr(), a_F[g].c_ptr(), a_Rho[g].c_ptr(), dt2 );
#else
    F77_FUNC(predfort,PREDFORT)(&m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g],
				 &m_kStart[g], &m_kEnd[g], 
				 a_Up[g].c_ptr(), a_U[g].c_ptr(), a_Um[g].c_ptr(),
				 a_Lu[g].c_ptr(), a_F[g].c_ptr(), a_Rho[g].c_ptr(), &dt2 );    
#endif
   }
}

//---------------------------------------------------------------------------
void EW::evalCorrector(vector<Sarray> & a_Up, vector<Sarray>& a_Rho,
		       vector<Sarray> & a_Lu, vector<Sarray> & a_F )
{
   float_sw4 dt4 = mDt*mDt*mDt*mDt;  
   for( int g=0 ; g<mNumberOfGrids; g++ )
   {
     a_Up[g].prefetch();
     a_Lu[g].prefetch();
     a_F[g].prefetch();
     a_Rho[g].prefetch();
#ifdef SW4_CROUTINES
      corrfort( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g], 
		a_Up[g].c_ptr(), a_Lu[g].c_ptr(), a_F[g].c_ptr(), a_Rho[g].c_ptr(), dt4 );
#else
      F77_FUNC(corrfort,CORRFORT)(&m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g],
				  &m_kStart[g], &m_kEnd[g], 
				  a_Up[g].c_ptr(), a_Lu[g].c_ptr(), a_F[g].c_ptr(),
				  a_Rho[g].c_ptr(), &dt4 );
#endif
   }
}

//---------------------------------------------------------------------------
void EW::evalDpDmInTime(vector<Sarray> & a_Up, vector<Sarray> & a_U, vector<Sarray> & a_Um,
			vector<Sarray> & a_Uacc )
{
   float_sw4 dt2i = 1./(mDt*mDt);
   for(int g=0 ; g<mNumberOfGrids; g++ )
   {
     //a_Up[g].prefetch();
     //a_U[g].prefetch();
     //a_Um[g].prefetch();
     //a_Uacc[g].prefetch();
     // Prefetched version is slower. Needs visitiing after Ray is upgraded.
#ifdef SW4_CROUTINES
      dpdmtfort( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
		 a_Up[g].c_ptr(), a_U[g].c_ptr(), a_Um[g].c_ptr(), a_Uacc[g].c_ptr(), dt2i );
#else
      F77_FUNC(dpdmtfort,DPDMTFORT)(&m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g],
				    &m_kStart[g], &m_kEnd[g],
				    a_Up[g].c_ptr(), a_U[g].c_ptr(), a_Um[g].c_ptr(),				    
				    a_Uacc[g].c_ptr(), &dt2i );
#endif
   }
   SYNC_DEVICE;
}

//-----------------------------------------------------------------------
void EW::evalRHS(vector<Sarray> & a_U, vector<Sarray>& a_Mu, vector<Sarray>& a_Lambda,
		 vector<Sarray> & a_Uacc )
{
  PUSH_RANGE("evalRHS",3);
   for(int g=0 ; g<mNumberOfCartesianGrids; g++ )
   {
     a_U[g].prefetch();
     a_Mu[g].prefetch();
     a_Lambda[g].prefetch();
     a_Uacc[g].prefetch();
     PREFETCH(m_sg_str_x[g]);
     PREFETCH(m_sg_str_y[g]);
     PREFETCH(m_sg_str_z[g]);
     PREFETCH(m_iop);
#ifdef SW4_CROUTINES
      if( m_corder )
	 rhs4sg_rev( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], 
		     m_kStart[g], m_kEnd[g], m_global_nz[g], m_onesided[g],
		     m_acof, m_bope, m_ghcof, a_Uacc[g].c_ptr(), a_U[g].c_ptr(), 
		     a_Mu[g].c_ptr(), a_Lambda[g].c_ptr(), mGridSize[g],
		     m_sg_str_x[g], m_sg_str_y[g], m_sg_str_z[g] );
      else
	 rhs4sg( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], 
		 m_kStart[g], m_kEnd[g], m_global_nz[g], m_onesided[g],
		 m_acof, m_bope, m_ghcof, a_Uacc[g].c_ptr(), a_U[g].c_ptr(), 
		 a_Mu[g].c_ptr(), a_Lambda[g].c_ptr(), mGridSize[g],
		 m_sg_str_x[g], m_sg_str_y[g], m_sg_str_z[g] );
#ifdef DEBUG_CUDA
      printf("params = %d, %d, %d, %d, %d, %d \n %f, %f, %f, %f \n %f, %f, %f, %f \n %d \n",  
           m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], 
              m_kStart[g], m_kEnd[g],
	   (a_Uacc[g].c_ptr())[1], (a_U[g].c_ptr())[1], 
              (a_Mu[g].c_ptr())[1], (a_Lambda[g].c_ptr())[1], 
	   mGridSize[g], m_sg_str_x[g][1], m_sg_str_y[g][1], m_sg_str_z[g][1],
           m_ghost_points); 
     printf("onesided[%d](4,5) = %d, %d\n", g, m_onesided[g][4], m_onesided[g][5]);
#endif
#else
      char op = '=';    // Assign Uacc := L(u)
      F77_FUNC(rhs4th3fortsgstr,RHS4TH3FORTSGSTR)( &m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g], 
						   &m_kStart[g], &m_kEnd[g], &m_global_nz[g], m_onesided[g],
						   m_acof, m_bope, m_ghcof, a_Uacc[g].c_ptr(), a_U[g].c_ptr(), 
						   a_Mu[g].c_ptr(), a_Lambda[g].c_ptr(), &mGridSize[g],
						   m_sg_str_x[g], m_sg_str_y[g], m_sg_str_z[g], &op );
#endif	
   }
   if( m_topography_exists )
   {
      int g=mNumberOfGrids-1;
      mMetric.prefetch();
      mJ.prefetch();
      char op = '=';    // Assign Uacc := L(u)
#ifdef SW4_CROUTINES
      if( m_corder )
         rhs4sgcurv_rev( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
			 a_U[g].c_ptr(), a_Mu[g].c_ptr(), a_Lambda[g].c_ptr(), mMetric.c_ptr(),
			 mJ.c_ptr(), a_Uacc[g].c_ptr(), m_onesided[g], m_acof, m_bope, m_ghcof,
			 m_sg_str_x[g], m_sg_str_y[g] );
      else
         rhs4sgcurv( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
		     a_U[g].c_ptr(), a_Mu[g].c_ptr(), a_Lambda[g].c_ptr(), mMetric.c_ptr(),
		     mJ.c_ptr(), a_Uacc[g].c_ptr(), m_onesided[g], m_acof, m_bope, m_ghcof,
		     m_sg_str_x[g], m_sg_str_y[g] );
#else      
      F77_FUNC(curvilinear4sg,CURVILINEAR4SG)(&m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g], 
					      &m_kStart[g], &m_kEnd[g], a_U[g].c_ptr(), a_Mu[g].c_ptr(),
					      a_Lambda[g].c_ptr(), mMetric.c_ptr(), mJ.c_ptr(), a_Uacc[g].c_ptr(),
					      m_onesided[g], m_acof, m_bope, m_ghcof, m_sg_str_x[g], m_sg_str_y[g],
					      &op );
#endif
   }
   POP_RANGE;
}
void EW::communicate_array( Sarray& u, int grid )
{
  //std::this_thread::sleep_for(std::chrono::milliseconds(100));
  communicate_array_async(u,grid);
  //std::this_thread::sleep_for(std::chrono::milliseconds(100));
}
//-----------------------------------------------------------------------
void EW::communicate_array_org( Sarray& u, int grid )
{
  PUSH_RANGE("COMM_ARRAY",4);
  //u.prefetch(cudaCpuDeviceId);
   REQUIRE2( u.m_nc == 1 || u.m_nc == 3 || u.m_nc == 4,
	     "Communicate array, only implemented for nc=1,3, and 4 "
	     << " nc = " << u.m_nc );
   int ie = u.m_ie, ib=u.m_ib, je=u.m_je, jb=u.m_jb, kb=u.m_kb;//,ke=u.m_ke;
   MPI_Status status;
   if( u.m_nc == 1 )
   {
      int xtag1 = 345;
      int xtag2 = 346;
      int ytag1 = 347;
      int ytag2 = 348;
      //for (int i=0;i<4;i++)std::cout<<m_neighbor[i]<<" ";
      //std::cout<<"\n";
      // X-direction communication
      MPI_Sendrecv( &u(ie-(2*m_ppadding-1),jb,kb), 1, m_send_type1[2*grid], m_neighbor[1], xtag1,
		    &u(ib,jb,kb), 1, m_send_type1[2*grid], m_neighbor[0], xtag1,
		    m_cartesian_communicator, &status );
      MPI_Sendrecv( &u(ib+m_ppadding,jb,kb), 1, m_send_type1[2*grid], m_neighbor[0], xtag2,
		    &u(ie-(m_ppadding-1),jb,kb), 1, m_send_type1[2*grid], m_neighbor[1], xtag2,
		    m_cartesian_communicator, &status );
      // Y-direction communication
      MPI_Sendrecv( &u(ib,je-(2*m_ppadding-1),kb), 1, m_send_type1[2*grid+1], m_neighbor[3], ytag1,
		    &u(ib,jb,kb), 1, m_send_type1[2*grid+1], m_neighbor[2], ytag1,
		    m_cartesian_communicator, &status );
      MPI_Sendrecv( &u(ib,jb+m_ppadding,kb), 1, m_send_type1[2*grid+1], m_neighbor[2], ytag2,
		    &u(ib,je-(m_ppadding-1),kb), 1, m_send_type1[2*grid+1], m_neighbor[3], ytag2,
		    m_cartesian_communicator, &status );
   }
   else if( u.m_nc == 3 )
   {
      int xtag1 = 345;
      int xtag2 = 346;
      int ytag1 = 347;
      int ytag2 = 348;
      // X-direction communication
      MPI_Sendrecv( &u(1,ie-(2*m_ppadding-1),jb,kb), 1, m_send_type3[2*grid], m_neighbor[1], xtag1,
		    &u(1,ib,jb,kb), 1, m_send_type3[2*grid], m_neighbor[0], xtag1,
		    m_cartesian_communicator, &status );
      MPI_Sendrecv( &u(1,ib+m_ppadding,jb,kb), 1, m_send_type3[2*grid], m_neighbor[0], xtag2,
		    &u(1,ie-(m_ppadding-1),jb,kb), 1, m_send_type3[2*grid], m_neighbor[1], xtag2,
		    m_cartesian_communicator, &status );
      // Y-direction communication
      MPI_Sendrecv( &u(1,ib,je-(2*m_ppadding-1),kb), 1, m_send_type3[2*grid+1], m_neighbor[3], ytag1,
		    &u(1,ib,jb,kb), 1, m_send_type3[2*grid+1], m_neighbor[2], ytag1,
		    m_cartesian_communicator, &status );
      MPI_Sendrecv( &u(1,ib,jb+m_ppadding,kb), 1, m_send_type3[2*grid+1], m_neighbor[2], ytag2,
		    &u(1,ib,je-(m_ppadding-1),kb), 1, m_send_type3[2*grid+1], m_neighbor[3], ytag2,
		    m_cartesian_communicator, &status );
   }
   else if( u.m_nc == 4 )
   {
      int xtag1 = 345;
      int xtag2 = 346;
      int ytag1 = 347;
      int ytag2 = 348;
      // X-direction communication
      MPI_Sendrecv( &u(1,ie-(2*m_ppadding-1),jb,kb), 1, m_send_type4[2*grid], m_neighbor[1], xtag1,
		    &u(1,ib,jb,kb), 1, m_send_type4[2*grid], m_neighbor[0], xtag1,
		    m_cartesian_communicator, &status );
      MPI_Sendrecv( &u(1,ib+m_ppadding,jb,kb), 1, m_send_type4[2*grid], m_neighbor[0], xtag2,
		    &u(1,ie-(m_ppadding-1),jb,kb), 1, m_send_type4[2*grid], m_neighbor[1], xtag2,
		    m_cartesian_communicator, &status );
      // Y-direction communication
      MPI_Sendrecv( &u(1,ib,je-(2*m_ppadding-1),kb), 1, m_send_type4[2*grid+1], m_neighbor[3], ytag1,
		    &u(1,ib,jb,kb), 1, m_send_type4[2*grid+1], m_neighbor[2], ytag1,
		    m_cartesian_communicator, &status );
      MPI_Sendrecv( &u(1,ib,jb+m_ppadding,kb), 1, m_send_type4[2*grid+1], m_neighbor[2], ytag2,
		    &u(1,ib,je-(m_ppadding-1),kb), 1, m_send_type4[2*grid+1], m_neighbor[3], ytag2,
		    m_cartesian_communicator, &status );
   }
   u.prefetch();
   POP_RANGE;
}

void EW::communicate_array_async(Sarray& u, int grid )
{
  PUSH_RANGE("COMM_ARRAY",4);
  //u.prefetch(cudaCpuDeviceId);
   REQUIRE2( u.m_nc == 1 || u.m_nc == 3 || u.m_nc == 4,
	     "Communicate array, only implemented for nc=1,3, and 4 "
	     << " nc = " << u.m_nc );
   int ie = u.m_ie, ib=u.m_ib, je=u.m_je, jb=u.m_jb, kb=u.m_kb;//,ke=u.m_ke;
   MPI_Status status;
#ifdef THREADED_MPI
   const int threaded_mpi=1;
#else
     const int threaded_mpi=0;
#endif
   if( u.m_nc == 1 )
   {
      int xtag1 = 345;
      int xtag2 = 346;
      int ytag1 = 347;
      int ytag2 = 348;
      // X-direction communication
      //for (int i=0;i<4;i++)std::cout<<m_neighbor[i]<<" ";
      //std::cout<<"\n";
#pragma omp parallel default(shared) if (threaded_mpi)
#pragma omp sections 
      {
#pragma omp section
      AMPI_Sendrecv( &u(ie-(2*m_ppadding-1),jb,kb), 1, send_type1[2*grid], m_neighbor[1], xtag1,
		    &u(ib,jb,kb), 1, send_type1[2*grid], m_neighbor[0], xtag1,
		     bufs_type1[4*grid],
		    m_cartesian_communicator, &status );
#pragma omp section
      AMPI_Sendrecv( &u(ib+m_ppadding,jb,kb), 1, send_type1[2*grid], m_neighbor[0], xtag2,
		    &u(ie-(m_ppadding-1),jb,kb), 1, send_type1[2*grid], m_neighbor[1], xtag2,
		     bufs_type1[4*grid+1],
		    m_cartesian_communicator, &status );
      // Y-direction communication
#pragma omp section
      AMPI_Sendrecv( &u(ib,je-(2*m_ppadding-1),kb), 1, send_type1[2*grid+1], m_neighbor[3], ytag1,
		     &u(ib,jb,kb), 1, send_type1[2*grid+1], m_neighbor[2], ytag1,
		     bufs_type1[4*grid+2],
		    m_cartesian_communicator, &status );    
#pragma omp section  
      AMPI_Sendrecv( &u(ib,jb+m_ppadding,kb), 1, send_type1[2*grid+1], m_neighbor[2], ytag2,
		    &u(ib,je-(m_ppadding-1),kb), 1, send_type1[2*grid+1], m_neighbor[3], ytag2,
		     bufs_type1[4*grid+3],
		    m_cartesian_communicator, &status );
   }
   }
   else if( u.m_nc == 3 )
   {

      int xtag1 = 345;
      int xtag2 = 346;
      int ytag1 = 347;
      int ytag2 = 348;
#pragma omp parallel default(shared) if (threaded_mpi)
#pragma omp sections 
      {
      // X-direction communication
#pragma omp section  
      AMPI_Sendrecv( &u(1,ie-(2*m_ppadding-1),jb,kb), 1, send_type3[2*grid], m_neighbor[1], xtag1,
		    &u(1,ib,jb,kb), 1, send_type3[2*grid], m_neighbor[0], xtag1,
		     bufs_type3[4*grid],
		    m_cartesian_communicator, &status );
#pragma omp section  
      AMPI_Sendrecv( &u(1,ib+m_ppadding,jb,kb), 1, send_type3[2*grid], m_neighbor[0], xtag2,
		    &u(1,ie-(m_ppadding-1),jb,kb), 1, send_type3[2*grid], m_neighbor[1], xtag2,
		    bufs_type3[4*grid+1],
		    m_cartesian_communicator, &status );
      // Y-direction communication
#pragma omp section  
      AMPI_Sendrecv( &u(1,ib,je-(2*m_ppadding-1),kb), 1, send_type3[2*grid+1], m_neighbor[3], ytag1,
		    &u(1,ib,jb,kb), 1, send_type3[2*grid+1], m_neighbor[2], ytag1,
		    bufs_type3[4*grid+2],
		    m_cartesian_communicator, &status );
#pragma omp section  
      AMPI_Sendrecv( &u(1,ib,jb+m_ppadding,kb), 1, send_type3[2*grid+1], m_neighbor[2], ytag2,
		    &u(1,ib,je-(m_ppadding-1),kb), 1, send_type3[2*grid+1], m_neighbor[3], ytag2,
		    bufs_type3[4*grid+3],
		    m_cartesian_communicator, &status );
      }
   }
   else if( u.m_nc == 4 )
   {
      int xtag1 = 345;
      int xtag2 = 346;
      int ytag1 = 347;
      int ytag2 = 348;
#pragma omp parallel default(shared) if (threaded_mpi)
#pragma omp sections 
      {
      // X-direction communication
#pragma omp section  
      AMPI_Sendrecv( &u(1,ie-(2*m_ppadding-1),jb,kb), 1, send_type4[2*grid], m_neighbor[1], xtag1,
		    &u(1,ib,jb,kb), 1, send_type4[2*grid], m_neighbor[0], xtag1,
		    bufs_type4[4*grid],
		    m_cartesian_communicator, &status );
#pragma omp section  
      AMPI_Sendrecv( &u(1,ib+m_ppadding,jb,kb), 1, send_type4[2*grid], m_neighbor[0], xtag2,
		    &u(1,ie-(m_ppadding-1),jb,kb), 1, send_type4[2*grid], m_neighbor[1], xtag2,
		    bufs_type4[4*grid+1],
		    m_cartesian_communicator, &status );
#pragma omp section  
      // Y-direction communication
      AMPI_Sendrecv( &u(1,ib,je-(2*m_ppadding-1),kb), 1, send_type4[2*grid+1], m_neighbor[3], ytag1,
		    &u(1,ib,jb,kb), 1, send_type4[2*grid+1], m_neighbor[2], ytag1,
		    bufs_type4[4*grid+2],
		    m_cartesian_communicator, &status );
#pragma omp section  
      AMPI_Sendrecv( &u(1,ib,jb+m_ppadding,kb), 1, send_type4[2*grid+1], m_neighbor[2], ytag2,
		    &u(1,ib,je-(m_ppadding-1),kb), 1, send_type4[2*grid+1], m_neighbor[3], ytag2,
		    bufs_type4[4*grid+3],
		    m_cartesian_communicator, &status );
      }
   }
   //u.prefetch();
   POP_RANGE;
}


void EW::communicate_array_new( Sarray& u, int grid )
{
  // The orignal communication patterns allows corner values to propagate twice during the
  // communcation. So 2nd send_recv picks put info received during the previous send_recv.
  // In this version this does not happen and it changes the results. The original way
  // can be diplicated by calling this twice, for example.
  // DO NOT USE FOR THE TIME BEING 
  
  PUSH_RANGE("COMM_ARRAY_SPLIT",4);
  //u.prefetch(cudaCpuDeviceId);
   REQUIRE2( u.m_nc == 1 || u.m_nc == 3 || u.m_nc == 4,
	     "Communicate array, only implemented for nc=1,3, and 4 "
	     << " nc = " << u.m_nc );
   int ie = u.m_ie, ib=u.m_ib, je=u.m_je, jb=u.m_jb, kb=u.m_kb;//,ke=u.m_ke;
   MPI_Status status;
   std::vector<AMPI_Ret_type> list;
   if( u.m_nc == 1 )
   {
      int xtag1 = 345;
      int xtag2 = 346;
      int ytag1 = 347;
      int ytag2 = 348;
       
      // X-direction communication
      //for (int i=0;i<4;i++)std::cout<<m_neighbor[i]<<" ";
      //std::cout<<"\n";
      list.push_back(AMPI_SendrecvSplit( &u(ie-(2*m_ppadding-1),jb,kb), 1, send_type1[2*grid], m_neighbor[1], xtag1,
       &u(ib,jb,kb), 1, send_type1[2*grid], m_neighbor[0], xtag1,
       bufs_type1[4*grid],
       m_cartesian_communicator, &status ));
      list.push_back(AMPI_SendrecvSplit( &u(ib+m_ppadding,jb,kb), 1, send_type1[2*grid], m_neighbor[0], xtag2,
       &u(ie-(m_ppadding-1),jb,kb), 1, send_type1[2*grid], m_neighbor[1], xtag2,
       bufs_type1[4*grid+1],
       m_cartesian_communicator, &status ));
      // Y-direction communication
      list.push_back(AMPI_SendrecvSplit( &u(ib,je-(2*m_ppadding-1),kb), 1, send_type1[2*grid+1], m_neighbor[3], ytag1,
       &u(ib,jb,kb), 1, send_type1[2*grid+1], m_neighbor[2], ytag1,
       bufs_type1[4*grid+2],
       m_cartesian_communicator, &status ));    
      list.push_back(AMPI_SendrecvSplit( &u(ib,jb+m_ppadding,kb), 1, send_type1[2*grid+1], m_neighbor[2], ytag2,
       &u(ib,je-(m_ppadding-1),kb), 1, send_type1[2*grid+1], m_neighbor[3], ytag2,
       bufs_type1[4*grid+3],
       m_cartesian_communicator, &status ));
      AMPI_SendrecvSync(list);
   }
   else if( u.m_nc == 3 )
   {

      int xtag1 = 345;
      int xtag2 = 346;
      int ytag1 = 347;
      int ytag2 = 348;
      
      // X-direction communication
      list.push_back(AMPI_SendrecvSplit( &u(1,ie-(2*m_ppadding-1),jb,kb), 1, send_type3[2*grid], m_neighbor[1], xtag1,
      		    &u(1,ib,jb,kb), 1, send_type3[2*grid], m_neighbor[0], xtag1,
      		     bufs_type3[4*grid],
      		    m_cartesian_communicator, &status ));
      list.push_back(AMPI_SendrecvSplit( &u(1,ib+m_ppadding,jb,kb), 1, send_type3[2*grid], m_neighbor[0], xtag2,
      		    &u(1,ie-(m_ppadding-1),jb,kb), 1, send_type3[2*grid], m_neighbor[1], xtag2,
      		    bufs_type3[4*grid+1],
      		    m_cartesian_communicator, &status ));
      // Y-direction communication
      list.push_back(AMPI_SendrecvSplit( &u(1,ib,je-(2*m_ppadding-1),kb), 1, send_type3[2*grid+1], m_neighbor[3], ytag1,
		    &u(1,ib,jb,kb), 1, send_type3[2*grid+1], m_neighbor[2], ytag1,
		    bufs_type3[4*grid+2],
		    m_cartesian_communicator, &status ));
      list.push_back(AMPI_SendrecvSplit( &u(1,ib,jb+m_ppadding,kb), 1, send_type3[2*grid+1], m_neighbor[2], ytag2,
		    &u(1,ib,je-(m_ppadding-1),kb), 1, send_type3[2*grid+1], m_neighbor[3], ytag2,
		    bufs_type3[4*grid+3],
		    m_cartesian_communicator, &status ));
      AMPI_SendrecvSync(list);
      //buffdiff(std::get<1>(bufs_typed[4*grid]),std::get<1>(bufs_type3[4*grid]),send_type3[2*grid]);
   }
   else if( u.m_nc == 4 )
   {
      int xtag1 = 345;
      int xtag2 = 346;
      int ytag1 = 347;
      int ytag2 = 348;

      // X-direction communication
      list.push_back(AMPI_SendrecvSplit( &u(1,ie-(2*m_ppadding-1),jb,kb), 1, send_type4[2*grid], m_neighbor[1], xtag1,
		    &u(1,ib,jb,kb), 1, send_type4[2*grid], m_neighbor[0], xtag1,
		    bufs_type4[4*grid],
		    m_cartesian_communicator, &status ));
      list.push_back(AMPI_SendrecvSplit( &u(1,ib+m_ppadding,jb,kb), 1, send_type4[2*grid], m_neighbor[0], xtag2,
		    &u(1,ie-(m_ppadding-1),jb,kb), 1, send_type4[2*grid], m_neighbor[1], xtag2,
		    bufs_type4[4*grid+1],
		    m_cartesian_communicator, &status ));
      // Y-direction communication
      list.push_back(AMPI_SendrecvSplit( &u(1,ib,je-(2*m_ppadding-1),kb), 1, send_type4[2*grid+1], m_neighbor[3], ytag1,
		    &u(1,ib,jb,kb), 1, send_type4[2*grid+1], m_neighbor[2], ytag1,
		    bufs_type4[4*grid+2],
		    m_cartesian_communicator, &status ));
      list.push_back(AMPI_SendrecvSplit( &u(1,ib,jb+m_ppadding,kb), 1, send_type4[2*grid+1], m_neighbor[2], ytag2,
		    &u(1,ib,je-(m_ppadding-1),kb), 1, send_type4[2*grid+1], m_neighbor[3], ytag2,
		    bufs_type4[4*grid+3],
		    m_cartesian_communicator, &status ));
      AMPI_SendrecvSync(list);
   }
   //u.prefetch();
   POP_RANGE;
}

//-----------------------------------------------------------------------
void EW::cartesian_bc_forcing( float_sw4 t, vector<float_sw4**> & a_BCForcing,
			      vector<Source*>& a_sources )
// assign the boundary forcing arrays a_BCForcing[g][side]
{
   for(int g=0 ; g<mNumberOfGrids; g++ )
   {
      if( m_point_source_test )
      {
	
	for( int side=0 ; side < 6 ; side++ ){
	  float_sw4 *tmp=a_BCForcing[g][side];
	    if( m_bcType[g][side] == bDirichlet )
	       get_exact_point_source2( a_BCForcing[g][side], t, g, *a_sources[0], &m_BndryWindow[g][6*side] );
	    else
	      forall<EXEC> (0,3*m_NumberOfBCPoints[g][side],[=] RAJA_DEVICE(int q){
		  tmp[q] = 0;});
	}
      }
      else
      {
	 // no boundary forcing
	 // we can do the same loop for all types of bc. For bParallel boundaries, numberOfBCPoints=0
	for( int side=0 ; side < 6 ; side++ ){
	  float_sw4 *tmp=a_BCForcing[g][side];
	  forall<EXEC> (0,3*m_NumberOfBCPoints[g][side],[=] RAJA_DEVICE(int q){
		  tmp[q] = 0;});
	}
      }
   }
   SYNC_DEVICE
}

//-----------------------------------------------------------------------
void EW::setup_boundary_arrays( )
{
   m_BndryWindow.resize(mNumberOfGrids);
   m_NumberOfBCPoints.resize(mNumberOfGrids);
   for (int g=0; g<mNumberOfGrids; g++ )
   {
      m_BndryWindow[g]      = new int[36];
      m_NumberOfBCPoints[g] = new int[6];
      for(int side=0; side<6 ; side++ )
      {
	 m_NumberOfBCPoints[g][side] = 0;
	 for (int qq=0; qq<6; qq+=2) // 0, 2, 4
	    m_BndryWindow[g][qq + side*6]= 999;
	 for (int qq=1; qq<6; qq+=2) // 1, 3, 5
	    m_BndryWindow[g][qq + side*6]= -999;
      }
      int wind[6];
      for(int side=0; side<6 ; side++ )
      {
	 if (m_bcType[g][side] == bStressFree || m_bcType[g][side] == bDirichlet || 
	     m_bcType[g][side] == bSuperGrid  || m_bcType[g][side] == bPeriodic)
	 {
	    // modify the window for stress free bc to only hold one plane
	    if (m_bcType[g][side] == bStressFree)
	    {
	       side_plane( g, side, wind, 1 );
// when calling side_plane with nGhost=1, you get the outermost grid plane
// for Free surface conditions, we apply the forcing on the boundary itself, i.e., just 
// inside the ghost points
// add/subtract the ghost point offset
	       if( side == 0 )
	       {
		  wind[0] += m_ghost_points;   wind[1] = wind[0];
	       }
	       else if( side == 1 )
	       {
		  wind[0] -= m_ghost_points;   wind[1] = wind[0];
	       }
	       else if( side == 2 )
	       {
		  wind[2] += m_ghost_points; wind[3] = wind[2];
	       }
	       else if( side == 3 )
	       {
		  wind[2]  -= m_ghost_points; wind[3] = wind[2];
	       }
	       else if( side == 4 )
	       {
		  wind[4] += m_ghost_points;
		  wind[5] = wind[4];
	       }
	       else
	       {
		  wind[4] -= m_ghost_points;
		  wind[5] = wind[4];
	       }
	    }
	    else // for Dirichlet, super grid, and periodic conditions, we
	       // apply the forcing directly on the ghost points
	    {
	       side_plane( g, side, wind, m_ghost_points );
	    }
	    int npts = (wind[5]-wind[4]+1)*
	       (wind[3]-wind[2]+1)*
	       (wind[1]-wind[0]+1);

	    for (int qq=0; qq<6; qq++)
	       m_BndryWindow[g][qq+side*6]=wind[qq];

	    m_NumberOfBCPoints[g][side] = npts;
	 }
      }
   }
}

//-----------------------------------------------------------------------
void EW::side_plane( int g, int side, int wind[6], int nGhost )
{
   wind[0] = m_iStart[g];
   wind[1] = m_iEnd[g];
   wind[2] = m_jStart[g];
   wind[3] = m_jEnd[g];
   wind[4] = m_kStart[g];
   wind[5] = m_kEnd[g];
   if( side == 0 )
     wind[1] = wind[0] + (nGhost-1);
   else if( side == 1 )
     wind[0] = wind[1] - (nGhost-1);
   else if( side == 2 )
     wind[3] = wind[2] + (nGhost-1);
   else if( side == 3 )
     wind[2] = wind[3] - (nGhost-1);
   else if( side == 4 )
     wind[5] = wind[4] + (nGhost-1);
   else
     wind[4] = wind[5] - (nGhost-1);
}

//-----------------------------------------------------------------------
void EW::enforceBC( vector<Sarray> & a_U, vector<Sarray>& a_Mu, vector<Sarray>& a_Lambda,
		    float_sw4 t, vector<float_sw4**> & a_BCForcing )
{
   float_sw4 om=0, ph=0, cv=0;
   for(int g=0 ; g<mNumberOfGrids; g++ )
   {
     PUSH_RANGE("enforceBC::PREFETCH",1);
     a_U[g].prefetch();
     a_Mu[g].prefetch();
     a_Lambda[g].prefetch();
     PREFETCH(m_sg_str_x[g]);
     PREFETCH(m_sg_str_y[g]);
     for(int j=0;j<6;j++)
       PREFETCH(a_BCForcing[g][j]);
     POP_RANGE;
      //      int topo=topographyExists() && g == mNumberOfGrids-1;
#ifdef SW4_CROUTINES
      if( m_corder )
	 bcfortsg_indrev( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
		m_BndryWindow[g], m_global_nx[g], m_global_ny[g], m_global_nz[g], a_U[g].c_ptr(),
		mGridSize[g], m_bcType[g], m_sbop, a_Mu[g].c_ptr(), a_Lambda[g].c_ptr(),
		t, a_BCForcing[g][0], a_BCForcing[g][1], a_BCForcing[g][2],
		a_BCForcing[g][3], a_BCForcing[g][4], a_BCForcing[g][5],
		om, ph, cv, m_sg_str_x[g], m_sg_str_y[g] );
      else
	 bcfortsg( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
		m_BndryWindow[g], m_global_nx[g], m_global_ny[g], m_global_nz[g], a_U[g].c_ptr(),
		mGridSize[g], m_bcType[g], m_sbop, a_Mu[g].c_ptr(), a_Lambda[g].c_ptr(),
		t, a_BCForcing[g][0], a_BCForcing[g][1], a_BCForcing[g][2],
		a_BCForcing[g][3], a_BCForcing[g][4], a_BCForcing[g][5],
		om, ph, cv, m_sg_str_x[g], m_sg_str_y[g] );
#else      
      F77_FUNC(bcfortsg, BCFORTSG)( &m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g], 
				    &m_kStart[g], &m_kEnd[g], m_BndryWindow[g], &m_global_nx[g],
				    &m_global_ny[g], &m_global_nz[g], a_U[g].c_ptr(),
				    &mGridSize[g], m_bcType[g], m_sbop, a_Mu[g].c_ptr(), a_Lambda[g].c_ptr(),
				    &t, a_BCForcing[g][0], a_BCForcing[g][1], a_BCForcing[g][2],
				    a_BCForcing[g][3], a_BCForcing[g][4], a_BCForcing[g][5],
				    &om, &ph, &cv, m_sg_str_x[g], m_sg_str_y[g] );
#endif
      if( m_topography_exists && g == mNumberOfGrids-1 && m_bcType[g][4] == bStressFree )
      {
	 int side = 5;
#ifdef SW4_CROUTINES
	 if( m_corder )
	    freesurfcurvisg_rev( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g], 
				 m_global_nz[g], side, a_U[g].c_ptr(), a_Mu[g].c_ptr(),
				 a_Lambda[g].c_ptr(), mMetric.c_ptr(), m_sbop,
			         a_BCForcing[g][4], m_sg_str_x[g], m_sg_str_y[g] );
	 else
	    freesurfcurvisg( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g], 
			     m_global_nz[g], side, a_U[g].c_ptr(), a_Mu[g].c_ptr(),
			     a_Lambda[g].c_ptr(), mMetric.c_ptr(), m_sbop,
			     a_BCForcing[g][4], m_sg_str_x[g], m_sg_str_y[g] );
#else
	 F77_FUNC(freesurfcurvisg,FREESURFCURVISG)( &m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g], 
						    &m_kStart[g], &m_kEnd[g], &m_global_nz[g], &side,
						    a_U[g].c_ptr(), a_Mu[g].c_ptr(), a_Lambda[g].c_ptr(),
						    mMetric.c_ptr(), m_sbop, a_BCForcing[g][4], m_sg_str_x[g], m_sg_str_y[g] );
#endif
      }

   }
   enforceCartTopo( a_U );
}

//-----------------------------------------------------------------------
void EW::enforceCartTopo( vector<Sarray>& a_U )
{
// interface between curvilinear and top Cartesian grid
   if (m_topography_exists)
   {
      int nc = 3;
      int g = mNumberOfCartesianGrids-1;
      int gc = mNumberOfGrids-1;
      int q, i, j;
// inject solution values between lower boundary of gc and upper boundary of g
      for( j = m_jStart[g] ; j <= m_jEnd[g]; j++ )
	 for( i = m_iStart[g]; i <= m_iEnd[g]; i++ )
	 {
// assign ghost points in the Cartesian grid
	    for (q = 0; q < m_ghost_points; q++) // only once when m_ghost_points==1
	    {
	       for( int c = 1; c <= nc ; c++ )
		  a_U[g](c,i,j,m_kStart[g] + q) = a_U[gc](c,i,j,m_kEnd[gc]-2*m_ghost_points + q);
	    }
// assign ghost points in the Curvilinear grid
	    for (q = 0; q <= m_ghost_points; q++) // twice when m_ghost_points==1 (overwrites solution on the common grid line)
	    {
	       for( int c = 1; c <= nc ; c++ )
		  a_U[gc](c,i,j,m_kEnd[gc]-q) = a_U[g](c,i,j,m_kStart[g]+2*m_ghost_points - q);
	    }
	 }
   }
}

//-----------------------------------------------------------------------
void EW::addSuperGridDamping(vector<Sarray> & a_Up, vector<Sarray> & a_U,
			     vector<Sarray> & a_Um, vector<Sarray> & a_Rho )
{
   for(int g=0 ; g<mNumberOfCartesianGrids; g++ )
   {
     a_Up[g].prefetch();
     a_U[g].prefetch();
     a_Um[g].prefetch();
     a_Rho[g].prefetch();

     PREFETCH(m_sg_dc_x[g]);
     PREFETCH(m_sg_dc_y[g]);
     PREFETCH(m_sg_dc_z[g]);

     PREFETCH(m_sg_str_x[g]);
     PREFETCH(m_sg_str_y[g]);
     PREFETCH(m_sg_str_z[g]);

     PREFETCH(m_sg_corner_x[g]);
     PREFETCH(m_sg_corner_y[g]);
     PREFETCH(m_sg_corner_z[g]);

      if( m_sg_damping_order == 4 )
      {
#ifdef SW4_CROUTINES	 
	 if( m_corder )
	    addsgd4fort_indrev( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
		      a_Up[g].c_ptr(), a_U[g].c_ptr(), a_Um[g].c_ptr(), a_Rho[g].c_ptr(),
		      m_sg_dc_x[g], m_sg_dc_y[g], m_sg_dc_z[g], m_sg_str_x[g], 
		      m_sg_str_y[g], m_sg_str_z[g], m_sg_corner_x[g], m_sg_corner_y[g],
		      m_sg_corner_z[g], m_supergrid_damping_coefficient );
	 else
	    addsgd4fort( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
		      a_Up[g].c_ptr(), a_U[g].c_ptr(), a_Um[g].c_ptr(), a_Rho[g].c_ptr(),
		      m_sg_dc_x[g], m_sg_dc_y[g], m_sg_dc_z[g], m_sg_str_x[g], 
		      m_sg_str_y[g], m_sg_str_z[g], m_sg_corner_x[g], m_sg_corner_y[g],
		      m_sg_corner_z[g], m_supergrid_damping_coefficient );
#else
	 F77_FUNC(addsgd4,ADDSGD4) ( &mDt, &mGridSize[g], a_Up[g].c_ptr(), a_U[g].c_ptr(), a_Um[g].c_ptr(), 
				     a_Rho[g].c_ptr(), m_sg_dc_x[g], m_sg_dc_y[g], m_sg_dc_z[g], m_sg_str_x[g], 
				     m_sg_str_y[g], m_sg_str_z[g], m_sg_corner_x[g], m_sg_corner_y[g],
				     m_sg_corner_z[g], &m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g], 
				    &m_kStart[g], &m_kEnd[g], &m_supergrid_damping_coefficient );
#endif
      }
      else if(  m_sg_damping_order == 6 )
      {
	std::cout<<" addsgd6fort_indrev needs offloading\n";
#ifdef SW4_CROUTINES	 
	 if( m_corder )
	    addsgd6fort_indrev( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
		      a_Up[g].c_ptr(), a_U[g].c_ptr(), a_Um[g].c_ptr(), a_Rho[g].c_ptr(),
		      m_sg_dc_x[g], m_sg_dc_y[g], m_sg_dc_z[g], m_sg_str_x[g], 
		      m_sg_str_y[g], m_sg_str_z[g], m_sg_corner_x[g], m_sg_corner_y[g],
		      m_sg_corner_z[g], m_supergrid_damping_coefficient );
	 else
	    addsgd6fort( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
		      a_Up[g].c_ptr(), a_U[g].c_ptr(), a_Um[g].c_ptr(), a_Rho[g].c_ptr(),
		      m_sg_dc_x[g], m_sg_dc_y[g], m_sg_dc_z[g], m_sg_str_x[g], 
		      m_sg_str_y[g], m_sg_str_z[g], m_sg_corner_x[g], m_sg_corner_y[g],
		      m_sg_corner_z[g], m_supergrid_damping_coefficient );
#else
	 F77_FUNC(addsgd6,ADDSGD6) ( &mDt, &mGridSize[g], a_Up[g].c_ptr(), a_U[g].c_ptr(), a_Um[g].c_ptr(), 
				     a_Rho[g].c_ptr(), m_sg_dc_x[g], m_sg_dc_y[g], m_sg_dc_z[g], m_sg_str_x[g], 
				     m_sg_str_y[g], m_sg_str_z[g], m_sg_corner_x[g], m_sg_corner_y[g],
				     m_sg_corner_z[g], &m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g], 
				     &m_kStart[g], &m_kEnd[g], &m_supergrid_damping_coefficient );
#endif
      }
   }
   if( m_topography_exists )
   {
      int g=mNumberOfGrids-1;
      if( m_sg_damping_order == 4 )
      {
#ifdef SW4_CROUTINES	 
	 if( m_corder )
	    addsgd4cfort_indrev( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
				 a_Up[g].c_ptr(), a_U[g].c_ptr(), a_Um[g].c_ptr(), a_Rho[g].c_ptr(),
				 m_sg_dc_x[g], m_sg_dc_y[g], m_sg_str_x[g], m_sg_str_y[g],
				 mJ.c_ptr(), m_sg_corner_x[g], m_sg_corner_y[g],
				 m_supergrid_damping_coefficient );
	 else
	    addsgd4cfort( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
			  a_Up[g].c_ptr(), a_U[g].c_ptr(), a_Um[g].c_ptr(), a_Rho[g].c_ptr(),
			  m_sg_dc_x[g], m_sg_dc_y[g], m_sg_str_x[g], m_sg_str_y[g], 
			  mJ.c_ptr(), m_sg_corner_x[g], m_sg_corner_y[g], m_supergrid_damping_coefficient );
#else
	 F77_FUNC(addsgd4c,ADDSGD4C) ( &mDt, a_Up[g].c_ptr(), a_U[g].c_ptr(), a_Um[g].c_ptr(), 
				     a_Rho[g].c_ptr(), m_sg_dc_x[g], m_sg_dc_y[g], m_sg_str_x[g], 
				       m_sg_str_y[g], mJ.c_ptr(), m_sg_corner_x[g], m_sg_corner_y[g],
				      &m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g], 
				       &m_kStart[g], &m_kEnd[g], &m_supergrid_damping_coefficient );
#endif
      }
      else if(  m_sg_damping_order == 6 )
      {
#ifdef SW4_CROUTINES	 
	 if( m_corder )
	    addsgd6cfort_indrev( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
				 a_Up[g].c_ptr(), a_U[g].c_ptr(), a_Um[g].c_ptr(), a_Rho[g].c_ptr(),
				 m_sg_dc_x[g], m_sg_dc_y[g], m_sg_str_x[g], m_sg_str_y[g],
				 mJ.c_ptr(), m_sg_corner_x[g], m_sg_corner_y[g],
				 m_supergrid_damping_coefficient );
	 else
	    addsgd6cfort( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
				 a_Up[g].c_ptr(), a_U[g].c_ptr(), a_Um[g].c_ptr(), a_Rho[g].c_ptr(),
				 m_sg_dc_x[g], m_sg_dc_y[g], m_sg_str_x[g], m_sg_str_y[g],
				 mJ.c_ptr(), m_sg_corner_x[g], m_sg_corner_y[g],
				 m_supergrid_damping_coefficient );
#else
	 F77_FUNC(addsgd6c,ADDSGD6C) ( &mDt, a_Up[g].c_ptr(), a_U[g].c_ptr(), a_Um[g].c_ptr(), 
				       a_Rho[g].c_ptr(), m_sg_dc_x[g], m_sg_dc_y[g], m_sg_str_x[g], 
				       m_sg_str_y[g], mJ.c_ptr(), m_sg_corner_x[g], m_sg_corner_y[g],
				       &m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g], 
				       &m_kStart[g], &m_kEnd[g], &m_supergrid_damping_coefficient );
#endif
      }

   }
}

//-----------------------------------------------------------------------
void EW::printTime( int cycle, float_sw4 t, bool force ) const 
{
   if (!mQuiet && m_myrank == 0 && (force || mPrintInterval == 1 ||
			(cycle % mPrintInterval) == 1 ||
				    cycle == 1) )
// string big enough for >1 million time steps
      cout << "Time step " << cycle << " t= " << t << endl;
   //      printf("Time step %7i  t = %15.7e\n", cycle, t);
}

//-----------------------------------------------------------------------
bool EW::exactSol( float_sw4 a_t, vector<Sarray> & a_U, vector<Source*>& sources )
{
  bool retval=false;
  if( m_point_source_test )
  {
     for( int g=0 ; g < mNumberOfGrids; g++ ) 
     {
	size_t npts = static_cast<size_t>(m_iEnd[g]-m_iStart[g]+1)*(m_jEnd[g]-m_jStart[g]+1)*(m_kEnd[g]-m_kStart[g]+1);
	float_sw4* utmp = newmanaged(npts*3);
	   //	get_exact_point_source( a_U[g].c_ptr(), a_t, g, *sources[0] );
	get_exact_point_source2( utmp, a_t, g, *sources[0] );
	a_U[g].assign( utmp, 0 );
	delmanaged(utmp);
     }
     retval = true;
  }
  return retval;
}

//-----------------------------------------------------------------------
// smooth wave for time dependence to test point force term with 
RAJA_HOST_DEVICE float_sw4 EW::SmoothWave(float_sw4 t, float_sw4 R, float_sw4 c)
{
  float_sw4 temp = R;
  float_sw4 c0 = 2187./8., c1 = -10935./8., c2 = 19683./8., c3 = -15309./8., c4 = 2187./4.;

  //  temp = where ( (t-R/c) > 0 && (t-R/c) < 1, (c0*pow(t-R/c,3)+c1*pow(t-R/c,4)+c2*pow(t-R/c,5)+c3*pow(t-R/c,6)+c4*pow(t-R/c,7)), 0);
  if( (t-R/c) > 0 && (t-R/c) < 1 )
     temp = (c0*pow(t-R/c,3)+c1*pow(t-R/c,4)+c2*pow(t-R/c,5)+c3*pow(t-R/c,6)+c4*pow(t-R/c,7));
  else
     temp = 0;
  return temp;
}

//-----------------------------------------------------------------------
// very smooth bump for time dependence for further testing of point force 
RAJA_HOST_DEVICE float_sw4 EW::VerySmoothBump(float_sw4 t, float_sw4 R, float_sw4 c)
{
  float_sw4 temp = R;
  float_sw4 c0 = 1024, c1 = -5120, c2 = 10240, c3 = -10240, c4 = 5120, c5 = -1024;

  //  temp = where ( (t-R/c) > 0 && (t-R/c) < 1, (c0*pow(t-R/c,5)+c1*pow(t-R/c,6)+c2*pow(t-R/c,7)+c3*pow(t-R/c,8)+c4*pow(t-R/c,9)+c5*pow(t-R/c,10)), 0);
  if( (t-R/c) > 0 && (t-R/c) < 1 )
     temp = (c0*pow(t-R/c,5)+c1*pow(t-R/c,6)+c2*pow(t-R/c,7)+c3*pow(t-R/c,8)+c4*pow(t-R/c,9)+c5*pow(t-R/c,10));
  else
     temp = 0;
  return temp;
}

//-----------------------------------------------------------------------
// C6 smooth bump for time dependence for further testing of point force 
RAJA_HOST_DEVICE float_sw4 EW::C6SmoothBump(float_sw4 t, float_sw4 R, float_sw4 c)
{
  float_sw4 retval = 0;
  if( (t-R/c) > 0 && (t-R/c) < 1 )
     retval = 51480.0*pow( (t-R/c)*(1-t+R/c), 7 );
  return retval;
}

//-----------------------------------------------------------------------
// derivative of smooth wave 
RAJA_HOST_DEVICE float_sw4 EW::d_SmoothWave_dt(float_sw4 t, float_sw4 R, float_sw4 c)
{
  float_sw4 temp = R;
  float_sw4 c0 = 2187./8., c1 = -10935./8., c2 = 19683./8., c3 = -15309./8., c4 = 2187./4.;

  //  temp = where ( (t-R/c) > 0 && (t-R/c) < 1, (3*c0*pow(t-R/c,2)+4*c1*pow(t-R/c,3)+5*c2*pow(t-R/c,4)+6*c3*pow(t-R/c,5)+7*c4*pow(t-R/c,6)), 0);
  if( (t-R/c) > 0 && (t-R/c) < 1 )
     temp = (3*c0*pow(t-R/c,2)+4*c1*pow(t-R/c,3)+5*c2*pow(t-R/c,4)+6*c3*pow(t-R/c,5)+7*c4*pow(t-R/c,6));
  else
     temp = 0;
  return temp;
}

//-----------------------------------------------------------------------
// very smooth bump for time dependence to further testing of point force 
RAJA_HOST_DEVICE float_sw4 EW::d_VerySmoothBump_dt(float_sw4 t, float_sw4 R, float_sw4 c)
{
  float_sw4 temp = R;
  float_sw4 c0 = 1024, c1 = -5120, c2 = 10240, c3 = -10240, c4 = 5120, c5 = -1024;

  //  temp = where ( (t-R/c) > 0 && (t-R/c) < 1, (5*c0*pow(t-R/c,4)+6*c1*pow(t-R/c,5)+7*c2*pow(t-R/c,6)+8*c3*pow(t-R/c,7)+9*c4*pow(t-R/c,8))+10*c5*pow(t-R/c,9), 0);
  if( (t-R/c) > 0 && (t-R/c) < 1 )
     temp = (5*c0*pow(t-R/c,4)+6*c1*pow(t-R/c,5)+7*c2*pow(t-R/c,6)+8*c3*pow(t-R/c,7)+9*c4*pow(t-R/c,8))+10*c5*pow(t-R/c,9);
  else
     temp = 0;
  return temp;
}

//-----------------------------------------------------------------------
// C6 smooth bump for time dependence to further testing of point force 
RAJA_HOST_DEVICE float_sw4 EW::d_C6SmoothBump_dt(float_sw4 t, float_sw4 R, float_sw4 c)
{
  float_sw4 retval=0;
  if( (t-R/c) > 0 && (t-R/c) < 1 )
     retval = 51480.0*7*(1-2*(t-R/c))*pow((t-R/c)*(1-t+R/c),6);
  return retval;
}

//-----------------------------------------------------------------------
// Primitive function (for T) of SmoothWave(t-T)*T
RAJA_HOST_DEVICE float_sw4 EW::SWTP(float_sw4 Lim, float_sw4 t)
{
  float_sw4 temp = Lim;

  float_sw4 c0 = 2187./8., c1 = -10935./8., c2 = 19683./8., c3 = -15309./8., c4 = 2187./4.;

  temp = (pow(t,3)*(c0 + c1*t + c2*pow(t,2) + c3*pow(t,3) + c4*pow(t,4))*pow(Lim,2))/2. - 
    (pow(t,2)*(3*c0 + 4*c1*t + 5*c2*pow(t,2) + 6*c3*pow(t,3) + 7*c4*pow(t,4))*pow(Lim,3))/3. + 
    (t*(3*c0 + 6*c1*t + 10*c2*pow(t,2) + 15*c3*pow(t,3) + 21*c4*pow(t,4))*pow(Lim,4))/4. + 
    ((-c0 - 4*c1*t - 10*c2*pow(t,2) - 20*c3*pow(t,3) - 35*c4*pow(t,4))*pow(Lim,5))/5. + 
    ((c1 + 5*c2*t + 15*c3*pow(t,2) + 35*c4*pow(t,3))*pow(Lim,6))/6. + 
    ((-c2 - 6*c3*t - 21*c4*pow(t,2))*pow(Lim,7))/7. + ((c3 + 7*c4*t)*pow(Lim,8))/8. - (c4*pow(Lim,9))/9.;

  return temp;
}

//-----------------------------------------------------------------------
// Primitive function (for T) of VerySmoothBump(t-T)*T
RAJA_HOST_DEVICE float_sw4 EW::VSBTP(float_sw4 Lim, float_sw4 t)
{
  float_sw4 temp = Lim;
  float_sw4 f = 1024., g = -5120., h = 10240., i = -10240., j = 5120., k = -1024.;

  temp = (pow(Lim,11)*(-25200*k*t-2520*j)+2310*k*pow(Lim,12)+(124740*k*pow(t,2)
							  +24948*j*t+2772*i)*pow(Lim,10)+(-369600*k*pow(t,3)-110880*j*pow(t,2)-24640*i*t-3080*h)*pow(Lim,9)+(727650*k*pow(t,4)+291060*j*pow(t,3)+97020*i*pow(t,2)+24255*h*t+3465*g)*pow(Lim,8)+(-997920*k*pow(t,5)-498960*j*pow(t,4)-221760*i*pow(t,3)-83160*h*pow(t,2)-23760*g*t-3960*f)*pow(Lim,7)+(970200*k*pow(t,6)+582120*j*pow(t,5)+323400*i*pow(t,4)+161700*h*pow(t,3)+69300*g*pow(t,2)+23100*f*t)*pow(Lim,6)+(-665280*k*pow(t,7)-465696*j*pow(t,6)-310464*i*pow(t,5)-194040*h*pow(t,4)-110880*g*pow(t,3)-55440*f*pow(t,2))*pow(Lim,5)+
	  (311850*k*pow(t,8)+249480*j*pow(t,7)+194040*i*pow(t,6)+145530*h*pow(t,5)+103950*g*pow(t,4)+69300*f*pow(t,3))*pow(Lim,4)+(-92400*
																   k*pow(t,9)-83160*j*pow(t,8)-73920*i*pow(t,7)-64680*h*pow(t,6)-55440*g*pow(t,5)-46200*f*pow(t,4))*pow(Lim,3)+(13860*k*pow(t,10)+13860*j*pow(t,9)+13860*i*pow(t,8)+13860*h*pow(t,7)+13860*g*pow(t,6)+13860*f*pow(t,5))*pow(Lim,2))/27720.0;

  return temp;
}
//-----------------------------------------------------------------------
// Primitive function (for T) of C6SmoothBump(t-T)*T
RAJA_HOST_DEVICE float_sw4 EW::C6SBTP(float_sw4 Lim, float_sw4 t)
{
  float_sw4 x = t-Lim;
  return pow(x,8)*(-3217.5*pow(x,8)+3432.0*(7+t)*pow(x,7)-25740.0*(3+t)*pow(x,6)
		   +27720.0*(5+3*t)*pow(x,5)-150150.0*(t+1)*x*x*x*x +
		   32760.0*(3+5*t)*x*x*x-36036.0*(1+3*t)*x*x+5720.0*(1+7*t)*x-6435.0*t);
}

//-----------------------------------------------------------------------
// Integral of H(t-T)*H(1-t+T)*SmoothWave(t-T)*T from R/alpha to R/beta
RAJA_HOST_DEVICE float_sw4 EW::SmoothWave_x_T_Integral(float_sw4 t, float_sw4 R, float_sw4 alpha, float_sw4 beta)
{
  float_sw4 temp = R;

  float_sw4 lowL, hiL;
  
  //  lowL = where(R / alpha > t - 1, R/alpha, t - 1); hiL = where(R / beta < t, R / beta, t);
  if( (R / alpha > t - 1 ) )
     lowL = R/alpha;
  else
     lowL = t-1;
  if( R / beta < t )
     hiL = R/beta;
  else
     hiL = t;
  
  //  temp = where (lowL < t && hiL > t - 1, SWTP(hiL, t) - SWTP(lowL, t), 0.0);
  if( lowL < t && hiL > t - 1 )
     temp = SWTP(hiL, t) - SWTP(lowL, t);
  else
     temp = 0;
  
  return temp;
}

//-----------------------------------------------------------------------
// Integral of H(t-T)*H(1-t+T)*VerySmoothBump(t-T)*T from R/alpha to R/beta
RAJA_HOST_DEVICE float_sw4 EW::VerySmoothBump_x_T_Integral(float_sw4 t, float_sw4 R, float_sw4 alpha, float_sw4 beta)
{
  float_sw4 temp = R;

  float_sw4 lowL, hiL;
  
  //  lowL = where(R / alpha > t - 1, R/alpha, t - 1); hiL = where(R / beta < t, R / beta, t);
  if( R / alpha > t - 1 )
     lowL = R/alpha;
  else
     lowL = t-1;
  if( R / beta < t )
     hiL = R/beta;
  else
     hiL = t;

  //  temp = where (lowL < t && hiL > t - 1, VSBTP(hiL, t) - VSBTP(lowL, t), 0.0);
  if( lowL < t && hiL > t - 1 )
     temp = VSBTP(hiL, t) - VSBTP(lowL, t);
  else
     temp = 0;
  return temp;
}

//-----------------------------------------------------------------------
// Integral of H(t-T)*H(1-t+T)*C6SmoothBump(t-T)*T from R/alpha to R/beta
RAJA_HOST_DEVICE float_sw4 EW::C6SmoothBump_x_T_Integral(float_sw4 t, float_sw4 R, float_sw4 alpha, float_sw4 beta)
{
  float_sw4 temp = R;

  float_sw4 lowL, hiL;
  
  //  lowL = where(R / alpha > t - 1, R/alpha, t - 1); hiL = where(R / beta < t, R / beta, t);
  if( R / alpha > t - 1 )
     lowL = R/alpha;
  else
     lowL = t-1;
  if( R / beta < t )
     hiL = R/beta;
  else
     hiL = t;

  //  temp = where (lowL < t && hiL > t - 1, VSBTP(hiL, t) - VSBTP(lowL, t), 0.0);
  if( lowL < t && hiL > t - 1 )
     temp = C6SBTP(hiL, t) - C6SBTP(lowL, t);
  else
     temp = 0;
  return temp;
}

//-----------------------------------------------------------------------
RAJA_HOST_DEVICE float_sw4 EW::Gaussian(float_sw4 t, float_sw4 R, float_sw4 c, float_sw4 f )
{
  float_sw4 temp = R;
  temp = 1 /(f* sqrt(2*M_PI))*exp(-pow(t-R/c,2) / (2*f*f));
  return temp;
}

//-----------------------------------------------------------------------
RAJA_HOST_DEVICE float_sw4 EW::d_Gaussian_dt(float_sw4 t, float_sw4 R, float_sw4 c, float_sw4 f)
{
  float_sw4 temp = R;
  temp = 1 /(f* sqrt(2*M_PI))*(-exp(-pow(t-R/c,2)/(2*f*f))*(t-R/c))/pow(f,2);
  return temp;
}

//-----------------------------------------------------------------------
RAJA_HOST_DEVICE float_sw4 EW::Gaussian_x_T_Integral(float_sw4 t, float_sw4 R, float_sw4 f, float_sw4 alpha, float_sw4 beta)
{
  float_sw4 temp = R;
  temp = -0.5*t*(erf( (t-R/beta)/(sqrt(2.0)*f))     - erf( (t-R/alpha)/(sqrt(2.0)*f)) ) -
     f/sqrt(2*M_PI)*( exp(-pow(t-R/beta,2)/(2*f*f) ) - exp( -pow(t-R/alpha,2)/(2*f*f) )  ) ;
     //  temp = 1/(f*sqrt(2*M_PI))*( f*f*(-exp(-pow(t-R/beta,2)/(2*f*f))+exp(-pow(t-R/alpha,2)/(2*f*f)) ) +
     //	     t*0.5*sqrt(M_PI*2)*f*( erf((t-R/alpha)/(sqrt(2.0)*f)) - erf((t-R/beta)/(sqrt(2.0)*f)) ) );
  //  temp = 1 /(f*sqrt(2*M_PI))*(f*( (-exp(-pow(t-R / alpha,2)/pow(f,2)) + exp(-pow(t-R / beta,2)/pow(f,2)) )*f + sqrt(M_PI)*t*(-erf((t-R / alpha) / f) + erf(R / beta / f))))/2.;
  return temp;
}

//-----------------------------------------------------------------------
void EW::get_exact_point_source( float_sw4* up, float_sw4 t, int g, Source& source,
				 int* wind )
{
   timeDep tD;
   if(!( source.getName() == "SmoothWave" || source.getName() == "VerySmoothBump" ||
	 source.getName() == "C6SmoothBump" || source.getName()== "Gaussian") )
   {
      cout << "EW::get_exact_point_source: Error, time dependency must be SmoothWave, VerySmoothBump, C6SmoothBump, or Gaussian, not "
	   << source.getName() << endl;
      return;
   }
   else if( source.getName() == "SmoothWave" )
      tD = iSmoothWave;
   else if( source.getName() == "VerySmoothBump" )
      tD = iVerySmoothBump;
   else if( source.getName() == "C6SmoothBump" )
      tD = iC6SmoothBump;
   else
      tD = iGaussian;

   //   u.set_to_zero();
   // Assume constant material, sample it in middle of domain
   int imid = (m_iStart[g]+m_iEnd[g])/2;
   int jmid = (m_jStart[g]+m_jEnd[g])/2;
   int kmid = (m_kStart[g]+m_kEnd[g])/2;
   float_sw4 rho   = mRho[g](imid,jmid,kmid);
   float_sw4 beta  =  sqrt( mMu[g](imid,jmid,kmid)/rho);
   float_sw4 alpha =  sqrt( (2*mMu[g](imid,jmid,kmid)+mLambda[g](imid,jmid,kmid))/rho);

      //   double alpha = m_point_source_test->m_cp;
      //   double beta  = m_point_source_test->m_cs;
      //   double rho   = m_point_source_test->m_rho;

   float_sw4 x0    = source.getX0();
   float_sw4 y0    = source.getY0();
   float_sw4 z0    = source.getZ0();
   float_sw4 fr=source.getFrequency();
   float_sw4 time = (t-source.getOffset()) * source.getFrequency();
   if( tD == iGaussian )
   {
      fr = 1/fr;
      time = time*fr;
   }
   bool ismomentsource = source.isMomentSource();
   float_sw4 fx, fy, fz;
   float_sw4 mxx, myy, mzz, mxy, mxz, myz, m0;

   if( !ismomentsource )
   {
      source.getForces( fx, fy, fz );
   }
   else
   {
      source.getMoments( mxx, mxy, mxz, myy, myz, mzz );
      //      m0  = source.getAmplitude();
      m0 = 1;
   }
   //   bool curvilinear = topographyExists() && g == mNumberOfGrids-1;
   bool curvilinear = false;
   //   float_sw4* up = u.c_ptr();
   float_sw4 h   = mGridSize[g];
   float_sw4 eps = 1e-3*h;
   size_t ind = 0;
   int imax, imin, jmax, jmin, kmax, kmin;
   if( wind == 0 )
   {
      imin = m_iStart[g];
      imax = m_iEnd[g];
      jmin = m_jStart[g];
      jmax = m_jEnd[g];
      kmin = m_kStart[g];
      kmax = m_kEnd[g];
   }
   else
   {
      imin = wind[0];
      imax = wind[1];
      jmin = wind[2];
      jmax = wind[3];
      kmin = wind[4];
      kmax = wind[5];
   }
   // Note: Use of ind, assumes loop is over the domain over which u is defined.
   //   for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
   //      for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
   //	 for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
   //   for( int k=kmin ; k <= kmax ; k++ )
   //   for( int j=jmin ; j <= jmax ; j++ )
   //	 for( int i=imin ; i <= imax ; i++ )
   int m_zmin_local= m_zmin[g];
   //PrintPointerAttributes(up);
   forallN<EXEC_CARTBC, int, int,int>(
				  RangeSegment(kmin,kmax+1),
				  RangeSegment(jmin,jmax+1),
				  RangeSegment(imin,imax+1),
				  [=] RAJA_DEVICE(int k, int j, int i) {
	 
            float_sw4 x,y,z;
	    size_t ind = (i-imin)+(j-jmin)*(imax-imin+1)+(k-kmin)*(jmax-jmin+1)*(imax-imin+1);
	    //	    if( curvilinear )
	    //	    {
	    //               x = mX(i,j,k);
	    //	       y = mY(i,j,k);
	    //	       z = mZ(i,j,k);
	    //	    }
	    //	    else
	    {
	       x = (i-1)*h;
	       y = (j-1)*h;
	       z = (k-1)*h + m_zmin_local;
	    }
	    if( !ismomentsource )
	    {
	       float_sw4 R = sqrt( (x - x0)*(x - x0) + (y - y0)*(y - y0) + (z - z0)*(z - z0) );
	       if( R < eps )
		  up[3*ind] = up[3*ind+1] = up[3*ind+2] = 0;
	       else
	       {
		  float_sw4 A, B;
		  if (tD == iSmoothWave)
		  {
		     A = ( 1/pow(alpha,2) * SmoothWave(time, fr*R, alpha) - 1/pow(beta,2) * SmoothWave(time, fr*R, beta) +
			   3/pow(fr*R,2) * SmoothWave_x_T_Integral(time, fr*R, alpha, beta) ) / (4*M_PI*rho*R*R*R)  ;
	  
		     B = ( 1/pow(beta,2) * SmoothWave(time, fr*R, beta) -
			   1/pow(fr*R,2) * SmoothWave_x_T_Integral(time, fr*R, alpha, beta) ) / (4*M_PI*rho*R) ;
		  }
		  else if (tD == iVerySmoothBump)
		  {
		     A = ( 1/pow(alpha,2) * VerySmoothBump(time, fr*R, alpha) - 1/pow(beta,2) * VerySmoothBump(time, fr*R, beta) +
			   3/pow(fr*R,2) * VerySmoothBump_x_T_Integral(time, fr*R, alpha, beta) ) / (4*M_PI*rho*R*R*R)  ;
		     
		     B = ( 1/pow(beta,2) * VerySmoothBump(time, fr*R, beta) -
			   1/pow(fr*R,2) * VerySmoothBump_x_T_Integral(time, fr*R, alpha, beta) ) / (4*M_PI*rho*R) ;
		  }
		  else if (tD == iC6SmoothBump)
		  {
		     A = ( 1/pow(alpha,2) * C6SmoothBump(time, fr*R, alpha) - 1/pow(beta,2) * C6SmoothBump(time, fr*R, beta) +
			   3/pow(fr*R,2) * C6SmoothBump_x_T_Integral(time, fr*R, alpha, beta) ) / (4*M_PI*rho*R*R*R)  ;
		     
		     B = ( 1/pow(beta,2) * C6SmoothBump(time, fr*R, beta) -
			   1/pow(fr*R,2) * C6SmoothBump_x_T_Integral(time, fr*R, alpha, beta) ) / (4*M_PI*rho*R) ;
		  }
                  else if( tD == iGaussian )
		  {
		     A = ( 1/pow(alpha,2) * Gaussian(time, R, alpha,fr) - 1/pow(beta,2) * Gaussian(time, R, beta,fr) +
			   3/pow(R,2) * Gaussian_x_T_Integral(time, R, fr, alpha, beta) ) / (4*M_PI*rho*R*R*R)  ;
		     
		     B = ( 1/pow(beta,2) * Gaussian(time, R, beta,fr) -
			   1/pow(R,2) * Gaussian_x_T_Integral(time, R, fr, alpha, beta) ) / (4*M_PI*rho*R) ;
		  }
		  up[3*ind]   = ( (x - x0)*(x - x0)*fx + (x - x0)*(y - y0)*fy + (x - x0)*(z - z0)*fz )*A + fx*B;
		  up[3*ind+1] = ( (y - y0)*(x - x0)*fx + (y - y0)*(y - y0)*fy + (y - y0)*(z - z0)*fz )*A + fy*B;
		  up[3*ind+2] = ( (z - z0)*(x - x0)*fx + (z - z0)*(y - y0)*fy + (z - z0)*(z - z0)*fz )*A + fz*B;
	       }
	    }
	    else 
	    {
	       up[3*ind] = up[3*ind+1] = up[3*ind+2] = 0;
	       // Here, ismomentsource == true
	       float_sw4 R = sqrt( (x - x0)*(x - x0) + (y - y0)*(y - y0) + (z - z0)*(z - z0) );
	       if( R < eps )
	       {
		  up[3*ind] = up[3*ind+1] = up[3*ind+2] = 0;
	       }
	       else
	       {
		  float_sw4 A, B, C, D, E;
		  if (tD == iSmoothWave)
		  {
		     A = SmoothWave(time, R, alpha);
		     B = SmoothWave(time, R, beta);
		     C = SmoothWave_x_T_Integral(time, R, alpha, beta);
		     D = d_SmoothWave_dt(time, R, alpha) / pow(alpha,3) / R;
		     E = d_SmoothWave_dt(time, R, beta) / pow(beta,3) / R;
		  }
		  else if (tD == iVerySmoothBump)
		  {
		     A = VerySmoothBump(time, R, alpha);
		     B = VerySmoothBump(time, R, beta);
		     C = VerySmoothBump_x_T_Integral(time, R, alpha, beta);
		     D = d_VerySmoothBump_dt(time, R, alpha) / pow(alpha,3) / R;
		     E = d_VerySmoothBump_dt(time, R, beta) / pow(beta,3) / R;
		  }
		  else if (tD == iC6SmoothBump)
		  {
		     A = C6SmoothBump(time, R, alpha);
		     B = C6SmoothBump(time, R, beta);
		     C = C6SmoothBump_x_T_Integral(time, R, alpha, beta);
		     D = d_C6SmoothBump_dt(time, R, alpha) / pow(alpha,3) / R;
		     E = d_C6SmoothBump_dt(time, R, beta) / pow(beta,3) / R;
		  }
		  else if (tD == iGaussian)
		  {
		     A = Gaussian(time, R, alpha,fr);
		     B = Gaussian(time, R, beta,fr);
		     C = Gaussian_x_T_Integral(time, R, fr,alpha, beta);
		     D = d_Gaussian_dt(time, R, alpha,fr) / pow(alpha,3) / R;
		     E = d_Gaussian_dt(time, R, beta,fr) / pow(beta,3) / R;
		  }
		  up[3*ind] += 
	// m_xx*G_xx,x
		     + m0*mxx/(4*M_PI*rho)*
		     ( 
		      + 3*(x-x0)*(x-x0)*(x-x0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      - 2*(x-x0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      + 3*(x-x0)*(x-x0) / pow(R,5) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))
	 
		      + ( 15*(x-x0)*(x-x0)*(x-x0) / pow(R,7) - 6*(x-x0) / pow(R,5) ) * C
	 
		      + (x-x0)*(x-x0) / pow(R,3)* ((x-x0)*D - (x-x0)*E)
	 
		      - 1 / pow(R,3) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))

		      - 3*(x-x0) / pow(R,5) * C

		      + (x-x0) / (pow(R,3)*pow(beta,2)) * B

		      + 1 / R * (x-x0)*E
		      );
		  up[3*ind] +=
		     // m_yy*G_xy,y
		     + m0*myy/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(y-y0)*(y-y0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      - (x-x0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (x-x0)*(y-y0) / pow(R,3)* ((y-y0)*D - (y-y0)*E)

		      + 3*(x-x0)*(y-y0) / pow(R,5) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))

		      + ( 15*(x-x0)*(y-y0)*(y-y0) / pow(R,7) - 3*(x-x0) / pow(R,5) ) * C
		      );
		  up[3*ind] +=
		     // m_zz*G_xz,z
		     + m0*mzz/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(z-z0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))

		      - (x-x0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (x-x0)*(z-z0) / pow(R,3)* ((z-z0)*D - (z-z0)*E)

		      + 3*(x-x0)*(z-z0) / pow(R,5) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))

		      + ( 15*(x-x0)*(z-z0)*(z-z0) / pow(R,7) - 3*(x-x0) / pow(R,5) ) * C
		      );
		  up[3*ind] +=
		     // m_xy*G_xy,x
		     + m0*mxy/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(x-x0)*(y-y0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))

		      - (y-y0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (x-x0)*(y-y0) / pow(R,3)* ((x-x0)*D - (x-x0)*E)

		      + 3*(x-x0)*(y-y0) / pow(R,5) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))

		      + ( 15*(x-x0)*(x-x0)*(y-y0) / pow(R,7) - 3*(y-y0) / pow(R,5) ) * C
		      );
		  up[3*ind] +=
		     // m_xy*G_xx,y
		     + m0*mxy/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(x-x0)*(y-y0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      + 3*(x-x0)*(x-x0) / pow(R,5) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))
	 
		      + 15*(x-x0)*(x-x0)*(y-y0) / pow(R,7) * C
	 
		      + (x-x0)*(x-x0) / pow(R,3)* ((y-y0)*D - (y-y0)*E)
	 
		      - 1 / pow(R,3) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))

		      - 3*(y-y0) / pow(R,5) * C

		      + (y-y0) / (pow(R,3)*pow(beta,2)) * B

		      + 1 / R * (y-y0)*E
		      );
		  up[3*ind] +=
		     // m_xz*G_xz,x
		     + m0*mxz/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(x-x0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))

		      - (z-z0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (x-x0)*(z-z0) / pow(R,3)* ((x-x0)*D - (x-x0)*E)

		      + 3*(x-x0)*(z-z0) / pow(R,5) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))

		      + ( 15*(x-x0)*(x-x0)*(z-z0) / pow(R,7) - 3*(z-z0) / pow(R,5) ) * C
		      );
		  up[3*ind] +=
		     // m_yz*G_xz,y
		     + m0*myz/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(y-y0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (x-x0)*(z-z0) / pow(R,3)* ((y-y0)*D - (y-y0)*E)

		      + 3*(x-x0)*(z-z0) / pow(R,5) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))

		      + 15*(x-x0)*(y-y0)*(z-z0) / pow(R,7) * C
		      );
		  up[3*ind] +=
		     // m_xz*G_xx,z
		     + m0*mxz/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(x-x0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      + 3*(x-x0)*(x-x0) / pow(R,5) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))
	 
		      + 15*(x-x0)*(x-x0)*(z-z0) / pow(R,7) * C
	 
		      + (x-x0)*(x-x0) / pow(R,3)* ((z-z0)*D - (z-z0)*E)
	 
		      - 1 / pow(R,3) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))

		      - 3*(z-z0) / pow(R,5) * C

		      + (z-z0) / (pow(R,3)*pow(beta,2)) * B

		      + 1 / R * (z-z0)*E
		      );
		  up[3*ind] +=
		     // m_yz*G_yx,z
		     + m0*myz/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(y-y0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (x-x0)*(y-y0) / pow(R,3)* ((z-z0)*D - (z-z0)*E)

		      + 3*(x-x0)*(y-y0) / pow(R,5) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))

		      + 15*(x-x0)*(y-y0)*(z-z0) / pow(R,7) * C
		      );
		  //------------------------------------------------------------
		  up[3*ind+1] += 
		     // m_xx*G_xy,x
		     m0*mxx/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(x-x0)*(y-y0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))

		      - (y-y0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (x-x0)*(y-y0) / pow(R,3)* ((x-x0)*D - (x-x0)*E)

		      + 3*(x-x0)*(y-y0) / pow(R,5) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))

		      + ( 15*(x-x0)*(x-x0)*(y-y0) / pow(R,7) - 3*(y-y0) / pow(R,5) ) * C
		      );
		  up[3*ind+1] += 
		     // m_yy**G_yy,y
		     + m0*myy/(4*M_PI*rho)*
		     ( 
		      + 3*(y-y0)*(y-y0)*(y-y0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      - 2*(y-y0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      + 3*(y-y0)*(y-y0) / pow(R,5) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))
	 
		      + ( 15*(y-y0)*(y-y0)*(y-y0) / pow(R,7) - 6*(y-y0) / pow(R,5) ) * C
	 
		      + (y-y0)*(y-y0) / pow(R,3)* ((y-y0)*D - (y-y0)*E)
	 
		      - 1 / pow(R,3) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))

		      - 3*(y-y0) / pow(R,5) * C

		      + (y-y0) / (pow(R,3)*pow(beta,2)) * B

		      + 1 / R * (y-y0)*E
		      );
		  up[3*ind+1] += 
		     // m_zz*G_zy,z
		     + m0*mzz/(4*M_PI*rho)*
		     (
		      + 3*(z-z0)*(z-z0)*(y-y0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))

		      - (y-y0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (z-z0)*(y-y0) / pow(R,3)* ((z-z0)*D - (z-z0)*E)

		      + 3*(z-z0)*(y-y0) / pow(R,5) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))

		      + ( 15*(z-z0)*(z-z0)*(y-y0) / pow(R,7) - 3*(y-y0) / pow(R,5) ) * C
		      );
		  up[3*ind+1] += 
		     // m_xy*G_yy,x
		     + m0*mxy/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(y-y0)*(y-y0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      + 3*(y-y0)*(y-y0) / pow(R,5) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))
	  
		      + 15*(x-x0)*(y-y0)*(y-y0) / pow(R,7) * C
	  
		      + (y-y0)*(y-y0) / pow(R,3)* ((x-x0)*D - (x-x0)*E)
	  
		      - 1 / pow(R,3) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))
	  
		      - 3*(x-x0) / pow(R,5) * C
	  
		      + (x-x0) / (pow(R,3)*pow(beta,2)) * B
	  
		      + 1 / R * (x-x0)*E
		      );
		  up[3*ind+1] += 
		     // m_xz*G_zy,x
		     + m0*mxz/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(y-y0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	  
		      + (y-y0)*(z-z0) / pow(R,3)* ((x-x0)*D - (x-x0)*E)
	  
		      + 3*(y-y0)*(z-z0) / pow(R,5) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))
	  
		      + 15*(x-x0)*(y-y0)*(z-z0) / pow(R,7) * C
		      );
		  up[3*ind+1] += 
		     // m_xy*G_xy,y
		     + m0*mxy/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(y-y0)*(y-y0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	  
		      - (x-x0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))
	  
		      + (x-x0)*(y-y0) / pow(R,3)* ((y-y0)*D - (y-y0)*E)
	  
		      + 3*(x-x0)*(y-y0) / pow(R,5) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))
	  
		      + ( 15*(x-x0)*(y-y0)*(y-y0) / pow(R,7) - 3*(x-x0) / pow(R,5) ) * C
		      );
		  up[3*ind+1] += 
		     // m_yz*G_zy,y
		     + m0*myz/(4*M_PI*rho)*
		     (
		      + 3*(z-z0)*(y-y0)*(y-y0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	  
		      - (z-z0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))
	  
		      + (z-z0)*(y-y0) / pow(R,3)* ((y-y0)*D - (y-y0)*E)
	  
		      + 3*(z-z0)*(y-y0) / pow(R,5) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))
	  
		      + ( 15*(z-z0)*(y-y0)*(y-y0) / pow(R,7) - 3*(z-z0) / pow(R,5) ) * C
		      );
		  up[3*ind+1] += 
		     // m_xz*G_xy,z
		     + m0*mxz/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(y-y0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	  
		      + (x-x0)*(y-y0) / pow(R,3)* ((z-z0)*D - (z-z0)*E)
	  
		      + 3*(x-x0)*(y-y0) / pow(R,5) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))
	  
		      + 15*(x-x0)*(y-y0)*(z-z0) / pow(R,7) * C
		      );
		  up[3*ind+1] += 
		     // m_yz*G_yy,z
		     + m0*myz/(4*M_PI*rho)*
		     (
		      + 3*(z-z0)*(y-y0)*(y-y0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      + 3*(y-y0)*(y-y0) / pow(R,5) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))
	  
		      + 15*(z-z0)*(y-y0)*(y-y0) / pow(R,7) * C
	  
		      + (y-y0)*(y-y0) / pow(R,3)* ((z-z0)*D - (z-z0)*E)
	  
		      - 1 / pow(R,3) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))
	  
		      - 3*(z-z0) / pow(R,5) * C
	  
		      + (z-z0) / (pow(R,3)*pow(beta,2)) * B
	  
		      + 1 / R * (z-z0)*E
		      );
		  //------------------------------------------------------------
		  up[3*ind+2] += 
		     // m_xx*G_zx,x
		     + m0*mxx/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(x-x0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))

		      - (z-z0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (x-x0)*(z-z0) / pow(R,3)* ((x-x0)*D - (x-x0)*E)

		      + 3*(x-x0)*(z-z0) / pow(R,5) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))

		      + ( 15*(x-x0)*(x-x0)*(z-z0) / pow(R,7) - 3*(z-z0) / pow(R,5) ) * C
		      );
		  up[3*ind+2] += 
		     // m_yy*G_zy,y
		     + m0*myy/(4*M_PI*rho)*
		     (
		      + 3*(y-y0)*(y-y0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))

		      - (z-z0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (y-y0)*(z-z0) / pow(R,3)* ((y-y0)*D - (y-y0)*E)

		      + 3*(y-y0)*(z-z0) / pow(R,5) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))

		      + ( 15*(y-y0)*(y-y0)*(z-z0) / pow(R,7) - 3*(z-z0) / pow(R,5) ) * C
		      );
		  up[3*ind+2] += 
		     // m_zz**G_zz,z
		     + m0*mzz/(4*M_PI*rho)*
		     ( 
		      + 3*(z-z0)*(z-z0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      - 2*(z-z0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      + 3*(z-z0)*(z-z0) / pow(R,5) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))
	 
		      + ( 15*(z-z0)*(z-z0)*(z-z0) / pow(R,7) - 6*(z-z0) / pow(R,5) ) * C
	 
		      + (z-z0)*(z-z0) / pow(R,3)* ((z-z0)*D - (z-z0)*E)
	 
		      - 1 / pow(R,3) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))

		      - 3*(z-z0) / pow(R,5) * C

		      + (z-z0) / (pow(R,3)*pow(beta,2)) * B

		      + 1 / R * (z-z0)*E
		      );
		  up[3*ind+2] += 
		     // m_xy*G_zy,x
		     + m0*mxy/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(y-y0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	  
		      + (y-y0)*(z-z0) / pow(R,3)* ((x-x0)*D - (x-x0)*E)
	  
		      + 3*(y-y0)*(z-z0) / pow(R,5) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))
	  
		      + 15*(x-x0)*(y-y0)*(z-z0) / pow(R,7) * C
		      );
		  up[3*ind+2] += 
		     // m_xz**G_zz,x
		     + m0*mxz/(4*M_PI*rho)*
		     ( 
		      + 3*(x-x0)*(z-z0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      + 3*(z-z0)*(z-z0) / pow(R,5) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))
	 
		      + 15*(x-x0)*(z-z0)*(z-z0) / pow(R,7) * C
	 
		      + (z-z0)*(z-z0) / pow(R,3)* ((x-x0)*D - (x-x0)*E)
	 
		      - 1 / pow(R,3) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))

		      - 3*(x-x0) / pow(R,5) * C

		      + (x-x0) / (pow(R,3)*pow(beta,2)) * B

		      + 1 / R * (x-x0)*E
		      );
		  up[3*ind+2] += 
		     // m_xy*G_xz,y
		     + m0*mxy/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(y-y0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (x-x0)*(z-z0) / pow(R,3)* ((y-y0)*D - (y-y0)*E)

		      + 3*(x-x0)*(z-z0) / pow(R,5) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))

		      + 15*(x-x0)*(y-y0)*(z-z0) / pow(R,7) * C
		      );
		  up[3*ind+2] += 
		     // m_yz*G_zz,y
		     + m0*myz/(4*M_PI*rho)*
		     ( 
		      + 3*(y-y0)*(z-z0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      + 3*(z-z0)*(z-z0) / pow(R,5) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))
	 
		      + 15*(y-y0)*(z-z0)*(z-z0) / pow(R,7) * C
	 
		      + (z-z0)*(z-z0) / pow(R,3)* ((y-y0)*D - (y-y0)*E)
	 
		      - 1 / pow(R,3) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))

		      - 3*(y-y0) / pow(R,5) * C

		      + (y-y0) / (pow(R,3)*pow(beta,2)) * B

		      + 1 / R * (y-y0)*E
		      );
		  up[3*ind+2] += 
		     // m_xz*G_xz,z
		     + m0*mxz/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(z-z0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      - (x-x0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      + (x-x0)*(z-z0) / pow(R,3)* ((z-z0)*D - (z-z0)*E)
	 
		      + 3*(x-x0)*(z-z0) / pow(R,5) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))
	 
		      + ( 15*(x-x0)*(z-z0)*(z-z0) / pow(R,7) - 3*(x-x0) / pow(R,5) ) * C
		      );
		  up[3*ind+2] += 
		     // m_yz*G_yz,z
		     + m0*myz/(4*M_PI*rho)*
		     (
		      + 3*(z-z0)*(z-z0)*(y-y0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))

		      - (y-y0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (z-z0)*(y-y0) / pow(R,3)* ((z-z0)*D - (z-z0)*E)

		      + 3*(z-z0)*(y-y0) / pow(R,5) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))

		      + ( 15*(z-z0)*(z-z0)*(y-y0) / pow(R,7) - 3*(y-y0) / pow(R,5) ) * C
		      );
		  //printf("VALUE = %lf \n",up[3*ind+2]);
	       }
	    }
	    //ind++;
				  });
#ifdef CUDA_CODE
   cudaDeviceSynchronize();
#endif
}

//-----------------------------------------------------------------------
void EW::normOfDifference( vector<Sarray> & a_Uex,  vector<Sarray> & a_U, float_sw4 &diffInf, 
                           float_sw4 &diffL2, float_sw4 &xInf, vector<Source*>& a_globalSources )
{
   float_sw4 linfLocal=0, l2Local=0, diffInfLocal=0, diffL2Local=0;
   float_sw4 xInfLocal=0, xInfGrid=0;
   float_sw4 htop = mGridSize[mNumberOfGrids-1];
   float_sw4 hbot = mGridSize[0];
   for(int g=0 ; g<mNumberOfGrids; g++ )
   {
      float_sw4 radius =-1, x0=0, y0=0, z0=0;
      float_sw4 h = mGridSize[g];
      int nsgxy = (int)(0.5+m_sg_gp_thickness*htop/h);
      int nsgz  = (int)(0.5+m_sg_gp_thickness*hbot/h);
      int imin, imax, jmin, jmax, kmin, kmax;

      // Remove supergrid layers
      if (mbcGlobalType[0] == bSuperGrid)
	 imin = max(m_iStartInt[g], nsgxy+1);
      else
	 imin = m_iStartInt[g];
  
      if (mbcGlobalType[1] == bSuperGrid)
	 imax = min(m_iEndInt[g], m_global_nx[g] - nsgxy);
      else
	 imax = m_iEndInt[g];

      if (mbcGlobalType[2] == bSuperGrid)
	 jmin = max(m_jStartInt[g], nsgxy+1);
      else
	 jmin = m_jStartInt[g];

      if (mbcGlobalType[3] == bSuperGrid)
	 jmax = min(m_jEndInt[g], m_global_ny[g] - nsgxy);
      else
	 jmax = m_jEndInt[g];

// Can not test on global type when there is more than one grid in the z-direction
// if uppermost grid has layer on top boundary, the fine grid spacing is used for the s.g. layer width
      if (m_bcType[g][4] == bSuperGrid)
	 kmin = max(m_kStartInt[g], nsgxy+1);
      else
	 kmin = m_kStartInt[g];
   // The lowermost grid has the s.g. layer width based on the spacing of the coarsest grid
      if (m_bcType[g][5] == bSuperGrid)
	 kmax = min(m_kEndInt[g], m_global_nz[g] - nsgz);
      else
	 kmax = m_kEndInt[g];
      if( m_point_source_test )
      {
	 radius = 4*h;
	 x0 = a_globalSources[0]->getX0();
	 y0 = a_globalSources[0]->getY0();
	 z0 = a_globalSources[0]->getZ0();
      }
// need to exclude parallel overlap from L2 calculation
//    int usesg = usingSupergrid();
    //    if( topographyExists() && g == mNumberOfGrids-1 )
    //       F77_FUNC(solerr3c, SOLERR3C)( &m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g], 
    //				    &m_kStart[g], &m_kEnd[g],  
    //				     a_Uex[g].c_ptr(), a_U[g].c_ptr(), mX.c_ptr(), mY.c_ptr(), mZ.c_ptr(), mJ.c_ptr(),
    //				     &linfLocal, &l2Local, &xInfGrid, &x0, &y0, &z0, &radius,
    //				     &imin, &imax, &jmin, &jmax, &kmin, &kmax,
    //				     &usesg, m_sg_str_x[g], m_sg_str_y[g] );
    //    else
#ifdef SW4_CROUTINES
      solerr3fort( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
  	           h, a_Uex[g].c_ptr(),  a_U[g].c_ptr(), linfLocal, l2Local, xInfGrid, 
		   m_zmin[g], x0, y0, z0, radius, imin, imax, jmin, jmax, kmin, kmax );
#else
      F77_FUNC(solerr3, SOLERR3)( &m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g], 
				  &m_kStart[g], &m_kEnd[g], &h, a_Uex[g].c_ptr(),
				  a_U[g].c_ptr(), &linfLocal, &l2Local, &xInfGrid, &m_zmin[g], 
				  &x0, &y0, &z0, &radius,
				  &imin, &imax, &jmin, &jmax, &kmin, &kmax );
#endif
      if (linfLocal > diffInfLocal) diffInfLocal = linfLocal;
      if (xInfGrid > xInfLocal) xInfLocal = xInfGrid;
      diffL2Local += l2Local;
   }
// communicate local results for global errors
   MPI_Allreduce( &diffInfLocal, &diffInf, 1, m_mpifloat, MPI_MAX, m_cartesian_communicator );
   MPI_Allreduce( &xInfLocal,    &xInf,    1, m_mpifloat, MPI_MAX, m_cartesian_communicator );
   MPI_Allreduce( &diffL2Local,  &diffL2,  1, m_mpifloat, MPI_SUM, m_cartesian_communicator );
   diffL2 = sqrt(diffL2);
}


//-----------------------------------------------------------------------
void EW::check_dimensions()
{
   for( int g= 0 ; g < mNumberOfGrids ; g++ )
   {
      int nz=m_kEndInt[g]-m_kStartInt[g]+1;
      int nzmin;
      if( m_onesided[g][4] && m_onesided[g][5] )
	 nzmin = 12;
      else if( m_onesided[g][4] || m_onesided[g][5] )
	 nzmin = 8;
      else
	 nzmin = 1;
      REQUIRE2( nz >= nzmin, "The number of grid points (not counting ghost pts) in the z-direction in grid " << g <<
		" must be >= " << nzmin << " current value is " << nz );
      int nx = m_iEndInt[g]-m_iStartInt[g]+1;
      REQUIRE2( nx >= 1, "No grid points left (not counting ghost pts) in the x-direction in grid " << g );
      int ny = m_jEndInt[g]-m_jStartInt[g]+1;
      REQUIRE2( ny >= 1, "No grid points left (not counting ghost pts) in the y-direction in grid " << g );
   }
}

//-----------------------------------------------------------------------
void EW::setup_supergrid( )
{
   if (mVerbose >= 3 && m_myrank == 0 )
      cout << "*** Inside setup_supergrid ***" << endl;
// check to see if there are any supergrid boundary conditions
   m_use_supergrid = false;
   for( int side=0 ; side < 6 ; side++ )
      if( mbcGlobalType[side] == bSuperGrid )
	 m_use_supergrid = true;
   if (mVerbose && m_myrank == 0 && m_use_supergrid)
      cout << "Detected at least one boundary with supergrid conditions" << endl;
   int gTop = mNumberOfCartesianGrids-1;
   int gBot = 0;
   m_supergrid_taper_z.resize(mNumberOfGrids);
   m_supergrid_taper_x.define_taper( (mbcGlobalType[0] == bSuperGrid), 0.0, (mbcGlobalType[1] == bSuperGrid), 
				     m_global_xmax, m_sg_gp_thickness*mGridSize[gTop] );
   m_supergrid_taper_y.define_taper( (mbcGlobalType[2] == bSuperGrid), 0.0, (mbcGlobalType[3] == bSuperGrid), 
				     m_global_ymax, m_sg_gp_thickness*mGridSize[gTop] );
   if( mNumberOfGrids == 1 )
      m_supergrid_taper_z[0].define_taper( !m_topography_exists && (mbcGlobalType[4] == bSuperGrid), 0.0,
					   (mbcGlobalType[5] == bSuperGrid), m_global_zmax,
					   m_sg_gp_thickness*mGridSize[gBot] );
   else
   {
      m_supergrid_taper_z[mNumberOfGrids-1].define_taper( !m_topography_exists && (mbcGlobalType[4] == bSuperGrid),
							  0.0, false, m_global_zmax,
							  m_sg_gp_thickness*mGridSize[gTop] );
      m_supergrid_taper_z[0].define_taper( false, 0.0, mbcGlobalType[5]==bSuperGrid, m_global_zmax,
					  m_sg_gp_thickness*mGridSize[gBot] );
      for( int g=1 ; g < mNumberOfGrids-1 ; g++ )
	 m_supergrid_taper_z[g].define_taper( false, 0.0, false, 0.0, m_sg_gp_thickness*mGridSize[gBot] );
   }
}

//-----------------------------------------------------------------------
void EW::assign_supergrid_damping_arrays()
{
  int i, j, k, topCartesian;
  float_sw4 x, y, z;
  
  //// resize the vectors for the pointers
  //  m_sg_dc_x.resize(mNumberOfGrids);
  //  m_sg_dc_y.resize(mNumberOfGrids);
  //  m_sg_dc_z.resize(mNumberOfGrids);

  //  m_sg_str_x.resize(mNumberOfGrids);
  //  m_sg_str_y.resize(mNumberOfGrids);
  //  m_sg_str_z.resize(mNumberOfGrids);

// new corner taper functions to reduce strength of damping near the edges and corners
//  m_sg_corner_x.resize(mNumberOfGrids);
//  m_sg_corner_y.resize(mNumberOfGrids);
//  m_sg_corner_z.resize(mNumberOfGrids);
  
// allocate storage for 1-D damping coefficients on each grid
//  for( int g=0 ; g<mNumberOfGrids; g++) 
//  {
//    m_sg_dc_x[g]  = new float_sw4[m_iEnd[g]-m_iStart[g]+1];
//    m_sg_dc_y[g]  = new float_sw4[m_jEnd[g]-m_jStart[g]+1];
//    m_sg_dc_z[g]  = new float_sw4[m_kEnd[g]-m_kStart[g]+1];

//    m_sg_str_x[g] = new float_sw4[m_iEnd[g]-m_iStart[g]+1];
//    m_sg_str_y[g] = new float_sw4[m_jEnd[g]-m_jStart[g]+1];
//    m_sg_str_z[g] = new float_sw4[m_kEnd[g]-m_kStart[g]+1];

// new corner taper functions to reduce strength of damping near the edges and corners
//    m_sg_corner_x[g] = new float_sw4[m_iEnd[g]-m_iStart[g]+1];
//    m_sg_corner_y[g] = new float_sw4[m_jEnd[g]-m_jStart[g]+1];
//    m_sg_corner_z[g] = new float_sw4[m_kEnd[g]-m_kStart[g]+1];
//  }

#define dcx(i,g) (m_sg_dc_x[g])[i-m_iStart[g]]
#define dcy(j,g) (m_sg_dc_y[g])[j-m_jStart[g]]
#define dcz(k,g) (m_sg_dc_z[g])[k-m_kStart[g]]

#define strx(i,g) (m_sg_str_x[g])[i-m_iStart[g]]
#define stry(j,g) (m_sg_str_y[g])[j-m_jStart[g]]
#define strz(k,g) (m_sg_str_z[g])[k-m_kStart[g]]

#define cornerx(i,g) (m_sg_corner_x[g])[i-m_iStart[g]]
#define cornery(j,g) (m_sg_corner_y[g])[j-m_jStart[g]]
#define cornerz(k,g) (m_sg_corner_z[g])[k-m_kStart[g]]

  //  topCartesian = mNumberOfCartesianGrids-1;
// Note: compared to WPP2, we don't need to center the damping coefficients on the half-point anymore,
// because the damping term is now 4th order: D+D-( a(x) D+D- ut(x) )

  topCartesian = mNumberOfCartesianGrids-1;
  if( m_use_supergrid )
  {
//       printf("SG: standard case!\n");
     for( int g=0 ; g<mNumberOfGrids; g++)  
     {
	for( i = m_iStart[g] ; i <= m_iEnd[g] ; i++ )
	{
	   x = (i-1)*mGridSize[g];
	   dcx(i,g)  = m_supergrid_taper_x.dampingCoeff(x);
	   strx(i,g) = m_supergrid_taper_x.stretching(x);
	   cornerx(i,g)  = m_supergrid_taper_x.cornerTaper(x);
	}
	for( j = m_jStart[g] ; j <= m_jEnd[g] ; j++ )
	{
	   y = (j-1)*mGridSize[g];
	   dcy(j,g)  = m_supergrid_taper_y.dampingCoeff(y);
	   stry(j,g) = m_supergrid_taper_y.stretching(y);
	   cornery(j,g)  = m_supergrid_taper_y.cornerTaper(y);
	}
	if (g > topCartesian || (0 < g && g < mNumberOfGrids-1)  ) // Curvilinear or refinement grid
	{
// No supergrid damping in the vertical (k-) direction on a curvilinear or refinement grid.
	   for( k = m_kStart[g] ; k <= m_kEnd[g] ; k++ )
	   {
	      dcz(k,g) = 0.;
	      strz(k,g) = 1;
	      cornerz(k,g) = 1.;
	   }
	}
	else
	{
	   for( k = m_kStart[g] ; k <= m_kEnd[g] ; k++ )
	   {
	      z = m_zmin[g] + (k-1)*mGridSize[g];
	      dcz(k,g)  = m_supergrid_taper_z[g].dampingCoeff(z);
	      strz(k,g) = m_supergrid_taper_z[g].stretching(z);
	      cornerz(k,g) = m_supergrid_taper_z[g].cornerTaper(z);
	   }
	}
     } // end for g...
  } // end if m_use_supergrid  
  else //
  {
//       printf("SG: supergrid not used!\n");
// Supergrid not used, but define arrays to simplify coding in some places.
     for( int g=0 ; g < mNumberOfGrids ; g++ )
     {
	for( i = m_iStart[g] ; i <= m_iEnd[g] ; i++ )
	{
	   dcx(i,g)  = 0;
	   strx(i,g) = 1;
	   cornerx(i,g) = 1.;
	}
	for( j = m_jStart[g] ; j <= m_jEnd[g] ; j++ )
	{
	   dcy(j,g)  = 0;
	   stry(j,g) = 1;
	   cornery(j,g) = 1.;
	}
	for( k = m_kStart[g] ; k <= m_kEnd[g] ; k++ )
	{
	   dcz(k,g)  = 0.;
	   strz(k,g) = 1;
	   cornerz(k,g) = 1.;
	}
     }
  }
  copy_supergrid_arrays_to_device();
#undef dcx
#undef dcy
#undef dcz
#undef strx
#undef stry
#undef strz
#undef cornerx
#undef cornery
#undef cornerz
}


//-----------------------------------------------------------------------
void EW::assign_local_bcs( )
{
// This routine assigns m_bcType[g][b], b=0,1,2,3, based on mbcGlobalType, taking parallel overlap boundaries into account

   int top=mNumberOfGrids-1; // index of the top grid in the arrays m_iStart, m_iEnd, etc
  
// horizontal bc's are the same for all grids
   for( int g= 0 ; g < mNumberOfGrids ; g++ )
   {
// start by copying the global bc's
      for (int b=0; b<=3; b++)
	 m_bcType[g][b] = mbcGlobalType[b];
  
      if (m_iStart[top]+m_ghost_points > 1)
      {
	 m_bcType[g][0] = bProcessor;
      }
      if (m_iEnd[top]-m_ghost_points < m_global_nx[top])
      {
	 m_bcType[g][1] = bProcessor;
      }
      if (m_jStart[top]+m_ghost_points > 1)
      {
	 m_bcType[g][2] = bProcessor;
      }
      if (m_jEnd[top]-m_ghost_points < m_global_ny[top])
      {
	 m_bcType[g][3] = bProcessor;
      }
   }
  
// vertical bc's are interpolating except at the bottom and the top, where they equal the global conditions
//   ( Only preliminary support for acoustic/elastic, not fully implemented)
   m_bcType[top][4] = mbcGlobalType[4];
   for( int g = 0 ; g < mNumberOfGrids-1 ; g++ )
  {
     if( m_is_curvilinear[g+1] && !m_is_curvilinear[g] ) // Elastic case only
	m_bcType[g][4] = bCCInterface;
     if( !m_is_curvilinear[g+1] && !m_is_curvilinear[g] ) // Two Cartesian grids, must be refinement bndry.
	m_bcType[g][4] = bRefInterface;
     if( !m_is_curvilinear[g+1] && m_is_curvilinear[g] ) // Acoustic case only
	m_bcType[g][4] = bCCInterface;
     if( m_is_curvilinear[g+1] && m_is_curvilinear[g] ) // Acoustic/Elastic interface
	m_bcType[g][4] = bAEInterface;
  }
  m_bcType[0][5] = mbcGlobalType[5];
  for( int g = 1 ; g < mNumberOfGrids ; g++ )
  {
     if( m_is_curvilinear[g] && !m_is_curvilinear[g-1] ) // Elastic case
	m_bcType[g][5] = bCCInterface;
     if( !m_is_curvilinear[g] && !m_is_curvilinear[g-1] ) // Two Cartesian grids, must be refinement bndry.
	m_bcType[g][5] = bRefInterface;
     if( !m_is_curvilinear[g] && m_is_curvilinear[g-1] ) // Acoustic case
	m_bcType[g][5] = bCCInterface;
     if( m_is_curvilinear[g] && m_is_curvilinear[g-1] ) // Acoustic/Elastic interface
	m_bcType[g][5] = bAEInterface;
  }

// Find out which boundaries need one sided approximation in mixed derivatives
  for( int g= 0 ; g < mNumberOfGrids ; g++ )
     for(int side=4 ; side < 6 ; side++ )
	m_onesided[g][side] = (m_bcType[g][side] == bStressFree) ||
	   (m_bcType[g][side] == bRefInterface) || (m_bcType[g][side] == bAEInterface); 
}

//-----------------------------------------------------------------------
void EW::create_output_directory( )
{
   if (m_myrank == 0 ) 
   {
      cout << "----------------------------------------------------" << endl
	   << " Making Output Directory: " << mPath << endl
	   << "\t\t" << endl;

      // Create directory where all these files will be written.
      int err = mkdirs(mPath);
      if (err == 0)
	cout << "... Done!" << endl
	     << "----------------------------------------------------" << endl;
      else
      {
// fatal error
	 cerr << endl << "******** Failed to create the output directory *******" << endl << endl;
	 MPI_Abort(MPI_COMM_WORLD,1);
      }

// check that we have write permission on the directory
      if (access(mPath.c_str(),W_OK)!=0)
      {
	 // fatal error
	 cerr << endl << "Error: No write permission on output directory: " << mPath << endl;
	 MPI_Abort(MPI_COMM_WORLD,1);
      }
      
   }
  // Let processor 0 finish first!
   cout.flush();  cerr.flush();
   MPI_Barrier(MPI_COMM_WORLD);

// Check that the mPath directory exists from all processes
   struct stat statBuf;
   int statErr = stat(mPath.c_str(), &statBuf);
   CHECK_INPUT(statErr == 0 && S_ISDIR(statBuf.st_mode), "Error: " << mPath << " is not a directory" << endl);
   
// check that all processes have write permission on the directory
   CHECK_INPUT(access(mPath.c_str(),W_OK)==0,
	   "Error: No write permission on output directory: " << mPath << endl);
}

//-----------------------------------------------------------------------
int EW::mkdirs(const string& path)
{
   //   string pathTemp(path.begin(), path.end()); 
   string pathTemp = path;
   //-----------------------------------------------------------------
   // Recursively call stat and then mkdir on each sub-directory in 'path'
   //-----------------------------------------------------------------
   string sep = "/";
   char * pathtemparg = new char[pathTemp.length()+1];
   strcpy(pathtemparg,pathTemp.c_str());
   char* token = strtok( pathtemparg, sep.c_str() );
//   char* token = strtok(const_cast<char*>(pathTemp.c_str()), sep.c_str());

   stringstream pathsofar;

// for checking the status:
   struct stat statBuf;
   int statErr;
   
   // If there's a leading slash, put it back on...
   if (strncmp(pathTemp.c_str(), sep.c_str(), 1) == 0) pathsofar << sep;

   while (token != NULL)
   {
      pathsofar << token << sep;

// test: check the status of the path so far...
//      cout << "Calling stat() on path: " << pathsofar.str() << endl;
      statErr = stat(pathsofar.str().c_str(), &statBuf);
      if (statErr == 0)
      {
//	cout << "stat() returned successfully." << endl;
	if ( S_ISDIR(statBuf.st_mode) )
	{
//	  cout << "stat() says: '" << pathsofar.str() << "' is a directory." << endl;
// it already exists, this is okay, let's get the next directory in the string and skip to the while statement
	  token = strtok(NULL, sep.c_str());
	  continue;
	}
	else
	{
	  cerr << "stat() says: '" << pathsofar.str() << "' is not a directory." << endl;
	// real error, let's bail...
	  delete[] pathtemparg;
	  return -1;
	}
	
      }
      else
      {
//	cerr << "stat() returned an error code." << endl;
	if (errno == EACCES)
	{
	  cerr << "Error: **Search permission is denied for one of the directories in the path prefix of " << pathsofar.str() << endl;
	  delete[] pathtemparg;
	  return -1;
	}
	else if (errno == ENOTDIR)
	{
	  cerr << "Error: **A component of the path '" <<  pathsofar.str() << "' is not a directory. " << endl;
	  delete[] pathtemparg;
	  return -1;
	}
 	else if (errno == ENOENT)
 	{
// this means that we need to call mkdir to create the directory
	  if (mVerbose >=2) 
	    cout << "Info: **stat returned ENOENT (the path does not exist, or the path " << endl
		 << "      is an empty string) " << pathsofar.str() << endl;
 	}
	else
	{
	  if (mVerbose >=2) 
	    cout << "Info: **stat returned other error code for path: " << pathsofar.str() << endl;
	}
      }

// if we got this far, then 'pathsofar' does not exists

// tmp
      if (mVerbose >=2) cout << "Calling mkdir() on path: " << pathsofar.str() << endl;
// old code for recursively making the output directory
       if (mkdir(pathsofar.str().c_str(), 
                S_IWUSR | S_IXUSR | S_IRUSR | S_IRGRP | S_IXGRP ) // why do we need group permissions?
          == -1)
      {
	if (mVerbose >=2) cout << "mkdir() returned an error code." << endl;
         // check error conditions
	if (errno == EEXIST)
	{
// can this ever happen since we called stat(), which said that the directory did not exist ???
	  if (mVerbose >=2) cout << "Info: ** The directory already exists:" << pathsofar.str() << endl;
	  
	  // it already exists, this is okay!
	  token = strtok(NULL, sep.c_str());
	  continue;
	}
	else if (errno == EACCES)
	  cerr << "Error: **Write permission is denied for the parent directory in which the new directory is to be added." << pathsofar.str() << endl;
	else if (errno == EMLINK)
	  cerr << "Error: **The parent directory has too many links (entries)." << 
	    pathsofar.str() << endl;
	else if (errno == ENOSPC)
	  cerr << "Error: **The file system doesn't have enough room to create the new directory." <<
	    pathsofar.str() << endl;
	else if (errno == EROFS)
	  cerr << "Error: **  The parent directory of the directory being created is on a read-only file system and cannot be modified." << pathsofar.str() << endl;
	else if (errno == ENOSPC)
	  cerr << "Error: ** The new directory cannot be created because the user's disk quota is exhausted." << pathsofar.str() << endl;
	// real error, let's bail...
	delete[] pathtemparg;
	return -1;
      }
      else
      {
	if (mVerbose >=2) cout << "mkdir() returned successfully." << endl;

// are there more directories to be made?
	token = strtok(NULL, sep.c_str());
      }
   }
   delete[] pathtemparg;
   return 0;
}

//-----------------------------------------------------------------------
void EW::computeDT()
{
   if (!mQuiet && mVerbose >= 1 && m_myrank == 0 )
      printf("*** computing the time step ***\n");

   float_sw4 dtloc=1.e10;
   for (int g=0; g<mNumberOfCartesianGrids; g++)
   {
      float_sw4 eigmax = -1;
      for (int k=m_kStart[g]; k<=m_kEnd[g]; k++)
	 for (int j=m_jStart[g]; j<=m_jEnd[g]; j++)
	    for (int i=m_iStart[g]; i<=m_iEnd[g]; i++)
	    {
	       float_sw4 loceig = (4*mMu[g](i,j,k) + mLambda[g](i,j,k) )/mRho[g](i,j,k);
	       eigmax = loceig > eigmax ? loceig:eigmax;
		  //	       dtGP = mCFL*mGridSize[g]/sqrt( loceig );
		  //	       dtloc = dtloc < dtGP ? dtloc : dtGP;
	    }
      float_sw4 ieigmax = 1/sqrt(eigmax);
      dtloc = dtloc < mCFL*mGridSize[g]*ieigmax ? dtloc : mCFL*mGridSize[g]*ieigmax;
   }
   if( m_topography_exists )
   {
#define SQR(x) (x)*(x)
// Curvilinear grid
      float_sw4 dtCurv;
      int g = mNumberOfGrids-1;
      float_sw4  la, mu, la2mu;
      int N=3, LDZ=1, INFO;
      char JOBZ='N', UPLO='L';
      float_sw4 eigmax = -1;
      // always use double precision version of lapack routine, for simplicity
      double Amat[6], W[3], Z[1], WORK[9];
// do consider ghost points (especially the ghost line above the topography might be important)
      for (int k=m_kStart[g]; k<=m_kEnd[g]; k++)
	 for (int j=m_jStart[g]; j<=m_jEnd[g]; j++)
	    for (int i=m_iStart[g]; i<=m_iEnd[g]; i++)
	    {
	       la = mLambda[g](i,j,k);
	       mu = mMu[g](i,j,k);
	       //	       for( int a = 0 ; a < m_number_mechanisms ; a++ )
	       //	       {
	       //		  la += mLambdaVE[g][a](i,j,k);
	       //		  mu += mMuVE[g][a](i,j,k);
	       //	       }
	       la2mu = la + 2.*mu;
	       float_sw4 jinv = 1/mJ(i,j,k);
// A11
	       Amat[0] = -4*(SQR(mMetric(1,i,j,k))*la2mu + SQR(mMetric(1,i,j,k))*mu + 
			 SQR(mMetric(2,i,j,k))*la2mu + SQR(mMetric(3,i,j,k))*mu + SQR(mMetric(4,i,j,k))*mu)*jinv;
// A21 = A12
	       Amat[1] = -4.*mMetric(2,i,j,k)*mMetric(3,i,j,k)*(mu+la)*jinv;
// A31 = A13	   
	       Amat[2] = -4.*mMetric(2,i,j,k)*mMetric(4,i,j,k)*(mu+la)*jinv;
// A22	   
	       Amat[3] = -4.*(SQR(mMetric(1,i,j,k))*mu + SQR(mMetric(1,i,j,k))*la2mu +
		        + SQR(mMetric(2,i,j,k))*mu + SQR(mMetric(3,i,j,k))*la2mu + SQR(mMetric(4,i,j,k))*mu)*jinv;
// A32 = A23
	       Amat[4] = -4.*mMetric(3,i,j,k)*mMetric(4,i,j,k)*(mu+la)*jinv;
// A33
	       Amat[5] = -4.*(SQR(mMetric(1,i,j,k))*mu + SQR(mMetric(1,i,j,k))*mu
			+ SQR(mMetric(2,i,j,k))*mu + SQR(mMetric(3,i,j,k))*mu + SQR(mMetric(4,i,j,k))*la2mu)*jinv;
// calculate eigenvalues of symmetric matrix
#ifndef SW4_CUDA
	       F77_FUNC(dspev,DSPEV)(JOBZ, UPLO, N, Amat, W, Z, LDZ, WORK, INFO);
#endif
	       if (INFO != 0)
	       {
		  printf("ERROR: computeDT: dspev returned INFO = %i for grid point (%i, %i, %i)\n", INFO, i, j, k);
		  printf("lambda = %e, mu = %e\n", la, mu);
		  printf("Jacobian = %15.7g \n",mJ(i,j,k));
		  printf("Matrix = \n");
		  printf(" %15.7g  %15.7g %15.7g \n",Amat[0],Amat[1],Amat[2]);
		  printf(" %15.7g  %15.7g %15.7g \n",Amat[1],Amat[3],Amat[4]);
		  printf(" %15.7g  %15.7g %15.7g \n",Amat[2],Amat[4],Amat[5]);
		  MPI_Abort(MPI_COMM_WORLD, 1);
	       }
// eigenvalues in ascending order: W[0] < W[1] < W[2]
	       if (W[0] >= 0.)
	       {
		  printf("ERROR: computeDT: determining eigenvalue is non-negative; W[0] = %e at curvilinear grid point (%i, %i, %i)\n", W[0], i, j, k);
		  MPI_Abort(MPI_COMM_WORLD, 1);
	       }
	       float_sw4 loceig = (-W[0])/(4.*mRho[g](i,j,k));
	       eigmax = loceig > eigmax ? loceig:eigmax;
	    }
      float_sw4 ieigmax = 1/sqrt(eigmax);
      dtCurv = mCFL*ieigmax;
      dtloc = dtloc<dtCurv ? dtloc: dtCurv;
#undef SQR
   } // end if topographyExists()
   mDt = dtloc;
// compute the global minima
   MPI_Allreduce( &dtloc, &mDt, 1, m_mpifloat, MPI_MIN, m_cartesian_communicator);
   if (!mQuiet && mVerbose >= 1 && m_myrank == 0 )
      cout << " CFL= " << mCFL << " prel. time step=" << mDt << endl;

   if( mTimeIsSet )
   {
// constrain the dt based on the goal time
      mNumberOfTimeSteps = static_cast<int> ((mTmax - mTstart) / mDt + 0.5); 
      mNumberOfTimeSteps = (mNumberOfTimeSteps==0)? 1: mNumberOfTimeSteps;
// the resulting mDt could be slightly too large, because the numberOfTimeSteps is rounded to the nearest int
      mDt = (mTmax - mTstart) / mNumberOfTimeSteps;
   }
}

//-----------------------------------------------------------------------
void EW::computeNearestGridPoint(int & a_i, 
                                   int & a_j, 
                                   int & a_k,
                                   int & a_g, // grid on which indices are located
                                   float_sw4 a_x, 
                                   float_sw4 a_y, 
                                   float_sw4 a_z)
{
  bool breakLoop = false;
  
  for (int g = 0; g < mNumberOfGrids; g++)
    {
      if (a_z > m_zmin[g] || g == mNumberOfGrids-1) // We can not trust zmin for the curvilinear grid, since it doesn't mean anything
        {
          a_i = (int)floor(a_x/mGridSize[g])+1;
          if (a_x-((a_i-0.5)*mGridSize[g]) > 0.) (a_i)++;
          
          a_j = (int)floor(a_y/mGridSize[g])+1;
          if (a_y-((a_j-0.5)*mGridSize[g]) > 0.) (a_j)++;
          
          a_k = (int)floor((a_z-m_zmin[g])/mGridSize[g])+1;  //Note: this component will be garbage for g=curvilinear grid
          if (a_z-(m_zmin[g]+(a_k-0.5)*mGridSize[g]) > 0.)   (a_k)++;
          
          a_g = g                                        ;
          
          breakLoop = true;
        }
      else if (a_z == m_zmin[g]) // testing for equality between doubles is kind of pointless...
        {
           // Point is located on top surface if g=finest grid, else the location is on
	   // a grid/grid interface, and point is flagged as located on the finer (upper) grid.
          if (g == mNumberOfGrids-1)
            {
              a_i = (int)floor(a_x/mGridSize[g])+1;
              if (a_x-((a_i-0.5)*mGridSize[g]) > 0.) (a_i)++;
              
              a_j = (int)floor(a_y/mGridSize[g])+1;
              if (a_y-((a_j-0.5)*mGridSize[g]) > 0.) (a_j)++;
              
              a_k = 1;
              
              a_g = g;
            }
          else
            {
              a_i = (int)floor(a_x/mGridSize[g+1])+1;
              if (a_x-((a_i-0.5)*mGridSize[g+1]) > 0.) (a_i)++;
              
              a_j = (int)floor(a_y/mGridSize[g+1])+1;
              if (a_y-((a_j-0.5)*mGridSize[g+1]) > 0.) (a_j)++;
              
              a_k = (int)floor((a_z-m_zmin[g+1])/mGridSize[g+1])+1; // Here, I know I am on a grid line
              
              a_g = g+1                                    ;
            }
          breakLoop = true;
        }
      
      if (breakLoop)
        {
              break;
        } 
    }
  
//  if (m_topography_exists && (a_g == mNumberOfGrids-1)) // The curvilinear grid will always be the one with the highest number. 
//    {
// tmp
//      printf("EW/computeNearestGridPt: You are in the curvilinear part of the grid, but we do compute the gridpt index using only the Cartesian grid\n");
//    }

// if z > zmax in grid 0 because the coordinate has not yet been corrected for topography, we simply set a_k to m_kEnd
  if (m_topography_exists && a_z >= m_global_zmax)
  {
    a_k = m_kEnd[0];
    a_g = 0;
  }

  if (!m_topography_exists || (m_topography_exists && a_g < mNumberOfCartesianGrids))
    {
      VERIFY2(a_i >= 1-m_ghost_points && a_i <= m_global_nx[a_g]+m_ghost_points,
              "Grid Error: i (" << a_i << ") is out of bounds: ( " << 1 << "," 
              << m_global_nx[a_g] << ")" << " x,y,z = " << a_x << " " << a_y << " " << a_z);
      VERIFY2(a_j >= 1-m_ghost_points && a_j <= m_global_ny[a_g]+m_ghost_points,
              "Grid Error: j (" << a_j << ") is out of bounds: ( " << 1 << ","
              << m_global_ny[a_g] << ")" << " x,y,z = " << a_x << " " << a_y << " " << a_z);
      VERIFY2(a_k >= m_kStart[a_g] && a_k <= m_kEnd[a_g],
              "Grid Error: k (" << a_k << ") is out of bounds: ( " << 1 << "," 
              << m_kEnd[a_g]-m_ghost_points << ")" << " x,y,z = " << a_x << " " << a_y << " " << a_z);
    }
}

//-----------------------------------------------------------------------
bool EW::interior_point_in_proc(int a_i, int a_j, int a_g)
{
// NOT TAKING PARALLEL GHOST POINTS INTO ACCOUNT!
// Determine if grid point with index (a_i, a_j) on grid a_g is an interior grid point on this processor 

   bool retval = false;
   if (a_g >=0 && a_g < mNumberOfGrids){
     retval = (a_i >= m_iStartInt[a_g]) && (a_i <= m_iEndInt[a_g]) &&   
       (a_j >= m_jStartInt[a_g]) && (a_j <= m_jEndInt[a_g]);
   }
   return retval; 
}

//-----------------------------------------------------------------------
bool EW::point_in_proc(int a_i, int a_j, int a_g)
{
// TAKING PARALLEL GHOST POINTS INTO ACCOUNT!
// Determine if grid point with index (a_i, a_j) on grid a_g is a grid point on this processor 

   bool retval = false;
   if (a_g >=0 && a_g < mNumberOfGrids){
     retval = (a_i >= m_iStart[a_g] && a_i <= m_iEnd[a_g] &&   
               a_j >= m_jStart[a_g] && a_j <= m_jEnd[a_g] );
   }
   return retval; 
}

//-----------------------------------------------------------------------
bool EW::point_in_proc_ext(int a_i, int a_j, int a_g)
{
// TAKING PARALLEL GHOST POINTS+EXTRA GHOST POINTS INTO ACCOUNT!
// Determine if grid point with index (a_i, a_j) on grid a_g is a grid point on this processor 

   bool retval = false;
   if (a_g >=0 && a_g < mNumberOfGrids){
     retval = (a_i >= m_iStart[a_g]-m_ext_ghost_points && a_i <= m_iEnd[a_g]+m_ext_ghost_points &&   
               a_j >= m_jStart[a_g]-m_ext_ghost_points && a_j <= m_jEnd[a_g]+m_ext_ghost_points );
   }
   return retval; 
}

//-----------------------------------------------------------------------
bool EW::is_onesided( int g, int side ) const
{
   return m_onesided[g][side] == 1;
}

//-----------------------------------------------------------------------
void EW::print_execution_time( double t1, double t2, string msg )
{
//   if( !mQuiet && proc_zero() )
   if( m_myrank == 0 )
   {
      double s = t2 - t1;
      int h = static_cast<int>(s/3600.0);
      s = s - h*3600;
      int m = static_cast<int>(s/60.0);
      s = s - m*60;
      cout << "   Execution time, " << msg << " ";
      if( h > 1 )
	 cout << h << " hours ";
      else if( h > 0 )
	 cout << h << " hour  ";

      if( m > 1 )
	 cout << m << " minutes ";
      else if( m > 0 )
	 cout << m << " minute  ";

      if( s > 0 )
	 cout << s << " seconds " ;
      cout << endl;
   }
}

//-----------------------------------------------------------------------
void EW::print_execution_times( double times[8] )
{
   double* time_sums =new double[8*m_nprocs];
   MPI_Gather( times, 8, MPI_DOUBLE, time_sums, 8, MPI_DOUBLE, 0, MPI_COMM_WORLD );
   bool printavgs = true;
   if( m_myrank == 0 )
   {
      double avgs[8]={0,0,0,0,0,0,0,0};
      for( int p= 0 ; p < m_nprocs ; p++ )
	 for( int c=0 ; c < 8 ; c++ )
	    avgs[c] += time_sums[8*p+c];
      for( int c=0 ; c < 8 ; c++ )
	 avgs[c] /= m_nprocs;
      
      cout << "\n----------------------------------------" << endl;
      cout << "          Execution time summary " << endl;
//      cout << "Processor  Total      BC total   Step   Image&Time series  Comm.ref   Comm.bndry BC impose  "
      if( printavgs )
      {
	 cout << "  Total      BC comm    BC phys    Scheme     Supergrid  Forcing "
	   <<endl;
	 cout.setf(ios::left);
	 cout.precision(5);
	 cout.width(11);
	 cout << avgs[7];
	 cout.width(11);
	 cout << avgs[2];
	 cout.width(11);
	 cout << avgs[3];
	 cout.width(11);
	 cout << avgs[1];
	 cout.width(11);
	 cout << avgs[4];
	 cout.width(11);
	 cout << avgs[0];
	 cout.width(11);
      }
      else
      {
	 cout << "Processor  Total      BC comm    BC phys    Scheme     Supergrid  Forcing "
	      <<endl;
	 cout.setf(ios::left);
	 cout.precision(5);
	 for( int p= 0 ; p < m_nprocs ; p++ )
	 {
	    cout.width(11);
	    cout << p;
	    cout.width(11);
	    cout << time_sums[8*p+7];
	    cout.width(11);
	    cout << time_sums[8*p+2];
	    cout.width(11);
	    cout << time_sums[8*p+3];
	    cout.width(11);
	    cout << time_sums[8*p+1];
	    cout.width(11);
	    cout << time_sums[8*p+4];
	    cout.width(11);
	    cout << time_sums[8*p];
	    cout.width(11);
	 //	 cout << time_sums[7*p+4];
	 //	 cout.width(11);
	 //	 cout << time_sums[7*p+5];
	 //	 cout.width(11);
	 //	 cout << time_sums[7*p+6];
	    cout << endl;
	 }
      }
      //
      // << "|" << time_sums[p*7+3] << "|\t" << time_sums[p*7+1] << "|\t" << time_sums[p*7]
      //	      << "|\t " << time_sums[7*p+2] << "|\t" << time_sums[p*7+4] << "|\t" << time_sums[p*7+5]
      //	      << "|\t" << time_sums[7*p+6]<<endl;
      cout << "Clock tick is " << MPI_Wtick() << " seconds" << endl;
      //      cout << "MPI_Wtime is ";
      //      int flag;
      //      bool wtime_is_global;
      // MPI_Comm_get_attr( MPI_COMM_WORLD, MPI_WTIME_IS_GLOBAL, &wtime_is_global, &flag );
      // if( wtime_is_global )
      // 	 cout << "global";
      // else
      // 	 cout << "local";
      // cout << endl;
      cout << "----------------------------------------\n" << endl;
      cout.setf(ios::right);
      cout.precision(6);

      // Save timings to file 
      string fname = mPath+"timings.bin";
      int fd=open( fname.c_str(), O_TRUNC|O_CREAT|O_WRONLY, 0660 );
      if( fd == -1 )
	 cout << "Error opening " << fname.c_str()  << " for writing execution times" << endl;
      size_t nr=write(fd,&m_nprocs,sizeof(int));
      if( nr != sizeof(int) )
	 cout << "Error writing nprocs on " << fname.c_str() << " nr = " << nr << " bytes" << endl;
      nr = write(fd, time_sums, 7*m_nprocs*sizeof(double));
      if( nr != 7*m_nprocs*sizeof(double) )
	 cout << "Error writing time_sums on " << fname.c_str() << " nr = " << nr << " bytes" << endl;
      close(fd);
   }
   delete[] time_sums;
}

//-----------------------------------------------------------------------
bool EW::check_for_match_on_cpu_gpu( vector<Sarray>& a_U, int verbose, string name )
{

   bool retval=false;

   if( m_cuobj->has_gpu() )
   {
      retval = false;
      for( int g=0 ; g<mNumberOfGrids; g++ )
      {
         size_t nn=a_U[g].check_match_cpu_gpu( m_cuobj, name );
         retval = retval || nn > 0;
         if( nn > 0 && verbose == 1 )
         {
            int cnan, inan, jnan, knan;
            a_U[g].check_match_cpu_gpu( m_cuobj, cnan, inan, jnan, knan, name );
            cout << "grid " << g << " array " << name << " found " << nn << "  dismatch. First dismatch at " <<
                    cnan << " " << inan << " " << jnan << " " << knan << endl;
         }
      }
   }
   return retval;
}

//-----------------------------------------------------------------------
void EW::setup_materials()
{
   // Point source test sets material directly in processTestPointSource
   if( !m_point_source_test )
   {
      // Undefined q-factors, attenutation not yet implemented
      vector<Sarray> Qs(mNumberOfGrids), Qp(mNumberOfGrids);
      for( int b=0 ; b < m_mtrlblocks.size() ; b++ )
	 m_mtrlblocks[b]->set_material_properties( mRho, mMu, mLambda, Qs, Qp );
      // Here mMu contains cs, and mLambda contains cp
      int g = mNumberOfGrids-1;
      extrapolateInZ( g, mRho[g],    true, false ); 
      extrapolateInZ( g, mLambda[g], true, false ); 
      extrapolateInZ( g, mMu[g],     true, false );
      g = 0;
      extrapolateInZ( g, mRho[g],    false, true ); 
      extrapolateInZ( g, mLambda[g], false, true ); 
      extrapolateInZ( g, mMu[g],     false, true );
      extrapolateInXY( mRho );
      extrapolateInXY( mMu );
      extrapolateInXY( mLambda );
      // Convert mMu to mu, and mLambda to lambda
      convert_material_to_mulambda( );
   }
}

//-----------------------------------------------------------------------
void EW::convert_material_to_mulambda( )
{
  for( int g = 0 ; g < mNumberOfGrids; g++)
  {
      // On input, we have stored cs in MU, cp in Lambda
      // use mu = rho*cs*cs and lambda = rho*cp*cp  - 2*mu
      for( int k = m_kStart[g] ; k <= m_kEnd[g]; k++ )
      {
          for( int j = m_jStart[g] ; j <= m_jEnd[g]; j++ )
          {
              for( int i = m_iStart[g] ; i <= m_iEnd[g] ; i++ )
              {
                  mMu[g](i,j,k)     = mRho[g](i,j,k)*mMu[g](i,j,k)*mMu[g](i,j,k);                  
                  mLambda[g](i,j,k) = mRho[g](i,j,k)*mLambda[g](i,j,k)*mLambda[g](i,j,k)-2*mMu[g](i,j,k);
              }
          }
      }
  }
}

//-----------------------------------------------------------------------
void EW::extrapolateInXY( vector<Sarray>& field )
{
   for( int g= 0; g < mNumberOfGrids ; g++ )
   {
      if( m_iStartInt[g] == 1 )
         for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
	    for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
	       for( int i=m_iStart[g] ; i < 1 ; i++ )
	       {
		  if( field[g](i,j,k) == -1 )
		     field[g](i,j,k) = field[g](1,j,k);
	       }
      if( m_iEndInt[g] == m_global_nx[g] )
         for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
	    for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
	       for( int i=m_iEndInt[g]+1 ; i <= m_iEnd[g] ; i++ )
	       {
		  if( field[g](i,j,k) == -1 )
		     field[g](i,j,k) = field[g](m_iEndInt[g],j,k);
	       }
      if( m_jStartInt[g] == 1 )
         for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
	    for( int j=m_jStart[g] ; j < 1 ; j++ )
	       for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
	       {
		  if( field[g](i,j,k) == -1 )
		     field[g](i,j,k) = field[g](i,1,k);
	       }
      if( m_jEndInt[g] == m_global_ny[g] )
         for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
	    for( int j=m_jEndInt[g]+1 ; j <= m_jEnd[g] ; j++ )
	       for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
	       {
		  if( field[g](i,j,k) == -1 )
		     field[g](i,j,k) = field[g](i,m_jEndInt[g],k);
	       }
// corners not necessary to treat explicitly???
      
   }
}

//-----------------------------------------------------------------------
void EW::extrapolateInZ( int g, Sarray& field, bool lowk, bool highk )
{
   if( lowk )
      for( int k=m_kStart[g] ; k < 1 ; k++ )
	 for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
	    for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
	       if( field(i,j,k) == -1 )
		  field(i,j,k) = field(i,j,1);
   if( highk )
      for( int k=m_kEndInt[g]+1 ; k <= m_kEnd[g] ; k++ )
	 for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
	    for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
	       if( field(i,j,k) == -1 )
		  field(i,j,k) = field(i,j,m_kEndInt[g]);
}

//-----------------------------------------------------------------------
void EW::getGlobalBoundingBox(float_sw4 bbox[6])
{
  bbox[0] = 0.;
  bbox[1] = m_global_xmax;
  bbox[2] = 0.;
  bbox[3] = m_global_ymax;
  bbox[4] = m_global_zmin;
  bbox[5] = m_global_zmax;
}

//-----------------------------------------------------------------------
bool EW::getDepth( float_sw4 x, float_sw4 y, float_sw4 z, float_sw4 & depth )
{
   if( !m_topography_exists )
   {
      depth = z;
      return true;
   }
   else
   {
      float_sw4 ztopo=0;
      if( find_topo_zcoord_owner(x,y,ztopo) )
      {
	 depth = z-ztopo;
	 return true;
      }
      else
	 return false;
   }
}

//-----------------------------------------------------------------------
bool EW::interpolate_topography( float_sw4 q, float_sw4 r, float_sw4 & Z0, bool smoothed)
{
// Interpolate the smoothed or raw topography

// Assume that (q,r) are indices in the curvilinear grid.

// if (q,r) is on this processor (need a 2x2 interval in (i,j)-index space:
// Return true and assign Z0 corresponding to (q,r)

// Returns false if 
// 1) (q,r) is outside the global parameter domain (expanded by ghost points)
// 2) (q,r) is not on this processor

// NOTE:
// The parameters are normalized such that 1 <= q <= Nx is the full domain (without ghost points),
//  1 <= r <= Ny.

// 0. No topography, easy case:
   if( !topographyExists() )
   {
      Z0 = 0;
      return true;
   }
// 1. Check that the point is inside the domain
   int g = mNumberOfGrids-1;
   float_sw4 h = mGridSize[g];
   float_sw4 qMin = (float_sw4) (1- m_ghost_points);
   float_sw4 qMax = (float_sw4) (m_global_nx[g] + m_ghost_points);
   float_sw4 rMin = (float_sw4) (1- m_ghost_points);
   float_sw4 rMax = (float_sw4) (m_global_ny[g] + m_ghost_points);
   if (!(q >= qMin && q <= qMax && r >= rMin && r <= rMax))
   {
      Z0 = 0;
      return false;
   }
// 2. Compute elevation at (q,r)
   float_sw4 tau; // holds the elevation at (q,r). Recall that elevation=-z
   if (m_analytical_topo)
   {
      float_sw4 X0   = (q-1.0)*h;
      float_sw4 Y0   = (r-1.0)*h;
      float_sw4 igx2 = 1.0/(m_GaussianLx*m_GaussianLx);
      float_sw4 igy2 = 1.0/(m_GaussianLy*m_GaussianLy);
      tau = m_GaussianAmp*exp(-(X0-m_GaussianXc)*(X0-m_GaussianXc)*igx2
			      -(Y0-m_GaussianYc)*(Y0-m_GaussianYc)*igy2 );
   }
   else
   {
// 3.Compute nearest grid point
      int iNear = static_cast<int>(round(q));
      int jNear = static_cast<int>(round(r));
      bool smackOnTop = (fabs(iNear-q) < 1.e-9 && fabs(jNear-r)) < 1.e-9;
      if (smackOnTop && point_in_proc(iNear,jNear,g)) 
      {
	 // 3a. (q,r) coincides with a grid point. Get elevation at that point.
	 if (smoothed)
	    tau = mTopoGridExt(iNear,jNear,1);
	 else
	    tau = mTopo(iNear,jNear,1);
      }
      else
      {
	 // 3b. (q,r) not at a grid point. Interpolate to get the elevation.
	 // Nearest lower grid point:
	 int i = static_cast<int>(floor(q));
	 int j = static_cast<int>(floor(r));
	 if( point_in_proc_ext(i-3,j-3,g) && point_in_proc_ext(i+4,j+4,g) )
	 {
	    float_sw4 a6cofi[8], a6cofj[8];
	    gettopowgh( q-i, a6cofi );
	    gettopowgh( r-j, a6cofj );
	    tau = 0;
	    for( int l=j-3 ; l <= j+4 ; l++ )
	       for( int k=i-3 ; k <= i+4 ; k++ )
		  tau += a6cofi[k-i+3]*a6cofj[l-j+3]*mTopoGridExt(k,l,1);
	 }
	 else
	 {
	    Z0 = 0;
	    return false;
	 }
      }
   }
   Z0 = -tau;
   return true;
}

//-----------------------------------------------------------------------
void EW::gettopowgh( float_sw4 ai, float_sw4 wgh[8] ) const
{
   float_sw4 pol = ai*ai*ai*ai*ai*ai*ai*(-251+135*ai+25*ai*ai-
                                      33*ai*ai*ai+6*ai*ai*ai*ai)/720;
   wgh[0] = -1.0/60*ai + 1.0/180*ai*ai + 1.0/48*ai*ai*ai + 23.0/144*ai*ai*ai*ai 
      - (17.0*ai + 223.0)*ai*ai*ai*ai*ai/720 - pol;
   wgh[1] = 3.0/20*ai -3.0/40*ai*ai -1.0/6*ai*ai*ai - 13.0/12*ai*ai*ai*ai + 
      97.0/45*ai*ai*ai*ai*ai + 1.0/6*ai*ai*ai*ai*ai*ai + 7*pol;
   wgh[2] = -0.75*ai +0.75*ai*ai+(13.0+155*ai)*ai*ai*ai/48 -103.0/16*ai*ai*ai*ai*ai
      - 121.0/240*ai*ai*ai*ai*ai*ai - 21*pol;
   wgh[3] = 1 - 49.0/36*ai*ai - 49.0/9*ai*ai*ai*ai+385.0/36*ai*ai*ai*ai*ai +
      61.0/72*ai*ai*ai*ai*ai*ai + 35*pol;
   wgh[4] = 0.75*ai + 0.75*ai*ai - 13.0/48*ai*ai*ai + 89.0/16*ai*ai*ai*ai - 
         1537.0/144*ai*ai*ai*ai*ai - 41.0/48*ai*ai*ai*ai*ai*ai - 35*pol;
   wgh[5] = -3.0/20*ai - 3.0/40*ai*ai + 1.0/6*ai*ai*ai - 41.0/12*ai*ai*ai*ai
      + 6.4*ai*ai*ai*ai*ai + 31.0/60*ai*ai*ai*ai*ai*ai + 21*pol;
   wgh[6] = 1.0/60*ai + 1.0/180*ai*ai - 1.0/48*ai*ai*ai + 167.0/144*ai*ai*ai*ai -
      1537.0/720*ai*ai*ai*ai*ai- 25.0/144*ai*ai*ai*ai*ai*ai - 7*pol;
   wgh[7] = -1.0/6*ai*ai*ai*ai + 11.0/36*ai*ai*ai*ai*ai + 1.0/40*ai*ai*ai*ai*ai*ai + pol;
}

//-----------------------------------------------------------------------
void EW::grid_mapping( float_sw4 q, float_sw4 r, float_sw4 s, float_sw4& x,
		       float_sw4& y, float_sw4& z )
{
   int g=mNumberOfGrids-1;
   float_sw4 h=mGridSize[g];
   x = (q-1)*h;
   y = (r-1)*h;
   float_sw4 ztopo;
   if( interpolate_topography(q,r,ztopo,true) )
   {
      int nz = m_global_nz[g];
      float_sw4 izb = 1.0/(m_zetaBreak*(nz-1));
      float_sw4 sa  = (s-1)*izb;
      float_sw4 omsm = (1-sa);
      for( int l=2 ; l <= m_grid_interpolation_order ; l++ )
	 omsm *= (1-sa);
      if( sa >= 1 )
	 z = m_topo_zmax - (nz-s)*h;
      else
	 z = m_topo_zmax - (nz-s)*h - omsm*(m_topo_zmax-(nz-1)*h-ztopo);
   }
   else
      z = -1e38;
   //   double zloc = z;
   //   MPI_Allreduce( &zloc, &z, 1, MPI_DOUBLE, MPI_MAX, m_cartesian_communicator );
}

//-----------------------------------------------------------------------
bool EW::invert_grid_mapping( int g, float_sw4 x, float_sw4 y, float_sw4 z, 
			      float_sw4& q, float_sw4& r, float_sw4& s )
{
// Translates (x,y,z) to grid indices on grid g.
// Successful only if (x,y) is in my processor, will return false if
// the point is outside the processor.
//
   bool success=true;
   q = x/mGridSize[g]+1;
   r = y/mGridSize[g]+1;
   if( g < mNumberOfCartesianGrids )
      s = (z-m_zmin[g])/mGridSize[g]+1;
   else
   {
      // Grid is curvilinear
      // Maximum number of iterations, and error tolerance
      // for Newton iterations
      int    maxit = 10;
      float_sw4 tol = 1e-9;
      float_sw4 zTopo;
      if( interpolate_topography(q, r, zTopo, true ) )      
      {
	 int    nz  = m_global_nz[g];
	 float_sw4 h   = mGridSize[g];
	 float_sw4 izb = 1.0/m_zetaBreak;
	    // Elastic region top grid, sun is s normalized to 0 < sun < 1
	 float_sw4 sun = 1-(m_topo_zmax-z)/((nz-1)*h);
	 if( sun >= m_zetaBreak )
	  // In uniform part of grid
	    s = (nz-1)*sun+1;
	 else
	 {
	  // Non-uniform, solve for s by Newton iteration
	    int it = 0;
	    float_sw4 numerr=tol+1;
	    while( numerr > tol && it < maxit )
	    {
	       float_sw4 omsm = (1-izb*sun);
	       for( int l=2 ; l <= m_grid_interpolation_order-1 ; l++ )
		  omsm *= (1-izb*sun);
	       float_sw4 fp =  h*(nz-1) + izb*m_grid_interpolation_order*omsm*(m_topo_zmax - (nz-1)*h - zTopo);
	       omsm *= (1-izb*sun);
	       float_sw4 f  = m_topo_zmax - (nz-1)*h*(1-sun) - omsm*(m_topo_zmax-(nz-1)*h-zTopo)-z;
	       float_sw4 ds= f/fp;
	       numerr = fabs(ds);
	       sun = sun - ds;
	       it++;
	    }
	    s = (nz-1)*sun+1;
	    if( numerr >= tol )
	    {
	       cout << "EW::invert_grid_mapping: WARNING no convergence err=" << numerr << " tol = " << tol << endl;
	       s = -1e38;
	       success = false;
	    }
	 }
      }
      else
      {
      // point not in processor, could not evaluate topography
	 s = -1e38;
	 success = false;
      }
   }
// Communicate the result to all processors
//   double sloc = s;
//   MPI_Allreduce( &sloc, &s, 1, MPI_DOUBLE, MPI_MAX, m_cartesian_communicator );
//   success = s != -1e38;
   return success;
}

//-----------------------------------------------------------------------
void EW::computeGeographicCoord(float_sw4 x, float_sw4 y, float_sw4 & longitude, float_sw4 & latitude)
{
  // conversion factor between degrees and radians
   //   if( m_geoproj == 0 )
   {
      float_sw4 deg2rad = M_PI/180.0;
      float_sw4 phi = mGeoAz * deg2rad;
      latitude = mLatOrigin + 
	 (x*cos(phi) - y*sin(phi))/mMetersPerDegree;
      if (mConstMetersPerLongitude)
      {
	 longitude = mLonOrigin + 
	    (x*sin(phi) + y*cos(phi))/(mMetersPerLongitude);
      }
      else
      {
	 longitude = mLonOrigin + 
	    (x*sin(phi) + y*cos(phi))/(mMetersPerDegree*cos(latitude*deg2rad));
      }
   }
   //   else
      //      m_geoproj->computeGeographicCoord( x, y, longitude, latitude );
}

//-----------------------------------------------------------------------
void EW::computeCartesianCoord(float_sw4 &x, float_sw4 &y, float_sw4 lon, float_sw4 lat)
{
  // -----------------------------------------------------------------
  // Compute the cartesian coordinate given the geographic coordinate
  // -----------------------------------------------------------------
   //   if( m_geoproj == 0 )
   //  // compute x and y
   {
      float_sw4 deg2rad = M_PI/180.0;
      float_sw4 phi = mGeoAz * deg2rad;
     //     x = mMetersPerDegree*(cos(phi)*(lat-mLatOrigin) + cos(lat*deg2rad)*(lon-mLonOrigin)*sin(phi));
     //     y = mMetersPerDegree*(-sin(phi)*(lat-mLatOrigin) + cos(lat*deg2rad)*(lon-mLonOrigin)*cos(phi));
      if (mConstMetersPerLongitude)
      {
	 x = mMetersPerDegree*cos(phi)*(lat-mLatOrigin)    + mMetersPerLongitude*(lon-mLonOrigin)*sin(phi);
	 y = mMetersPerDegree*(-sin(phi))*(lat-mLatOrigin) + mMetersPerLongitude*(lon-mLonOrigin)*cos(phi);
      }
      else
      {
	 x = mMetersPerDegree*(cos(phi)*(lat-mLatOrigin) + cos(lat*deg2rad)*(lon-mLonOrigin)*sin(phi));
	 y = mMetersPerDegree*(-sin(phi)*(lat-mLatOrigin) + cos(lat*deg2rad)*(lon-mLonOrigin)*cos(phi));
      }
   }
   //   else
      //      m_geoproj->computeCartesianCoord(x,y,lon,lat);
}

//-----------------------------------------------------------------------
void EW::get_utc( int utc[7] ) const
{
   for( int c=0 ; c < 7 ; c++ )
      utc[c] = m_utc0[c];
}

//-----------------------------------------------------------------------
void EW::extractRecordData(TimeSeries::receiverMode mode, int i0, int j0, int k0, int g0, 
			   vector<float_sw4> &uRec, vector<Sarray> &Um2, vector<Sarray> &U)
{
  if (mode == TimeSeries::Displacement)
  {
    uRec.resize(3);
    uRec[0] = U[g0](1, i0, j0, k0);
    uRec[1] = U[g0](2, i0, j0, k0);
    uRec[2] = U[g0](3, i0, j0, k0);
  }
  else if (mode == TimeSeries::Velocity)
  {
    uRec.resize(3);
    uRec[0] = (U[g0](1, i0, j0, k0) - Um2[g0](1, i0, j0, k0))/(2*mDt);
    uRec[1] = (U[g0](2, i0, j0, k0) - Um2[g0](2, i0, j0, k0))/(2*mDt);
    uRec[2] = (U[g0](3, i0, j0, k0) - Um2[g0](3, i0, j0, k0))/(2*mDt);
  }
  else if(mode == TimeSeries::Div)
  {
    uRec.resize(1);
    if (g0 < mNumberOfCartesianGrids) // must be a Cartesian grid
    {
//      int i=m_i0, j=m_j0, k=m_k0, g=m_grid0;
      float_sw4 factor = 1.0/(2*mGridSize[g0]);
      uRec[0] = ((U[g0](1,i0+1, j0, k0) - U[g0](1,i0-1, j0, k0)+
		  U[g0](2,i0, j0+1, k0) - U[g0](2,i0, j0-1, k0)+
		  U[g0](3,i0, j0, k0+1) - U[g0](3,i0, j0, k0-1))*factor);
    }
    else // must be curvilinear
    {
//      int i=m_i0, j=m_j0, k=m_k0, g=m_grid0;
       float_sw4 factor = 0.5/sqrt(mJ(i0,j0,k0));
       uRec[0] = ( ( mMetric(1,i0,j0,k0)*(U[g0](1,i0+1,j0,k0) - U[g0](1,i0-1,j0,k0))+
		     mMetric(1,i0,j0,k0)*(U[g0](2,i0,j0+1,k0) - U[g0](2,i0,j0-1,k0))+
		     mMetric(2,i0,j0,k0)*(U[g0](1,i0,j0,k0+1) - U[g0](1,i0,j0,k0-1))+
		     mMetric(3,i0,j0,k0)*(U[g0](2,i0,j0,k0+1) - U[g0](2,i0,j0,k0-1))+
		     mMetric(4,i0,j0,k0)*(U[g0](3,i0,j0,k0+1) - U[g0](3,i0,j0,k0-1))  )*factor);
    }
  } // end div
  else if(mode == TimeSeries::Curl)
  {
    uRec.resize(3);
    if (g0 < mNumberOfCartesianGrids) // must be a Cartesian grid
    {
//       int i=m_i0, j=m_j0, k=m_k0, g=m_grid0;
      float_sw4 factor = 1.0/(2*mGridSize[g0]);
      float_sw4 duydx = (U[g0](2,i0+1,j0,k0) - U[g0](2,i0-1,j0,k0))*factor;
      float_sw4 duzdx = (U[g0](3,i0+1,j0,k0) - U[g0](3,i0-1,j0,k0))*factor;
      float_sw4 duxdy = (U[g0](1,i0,j0+1,k0) - U[g0](1,i0,j0-1,k0))*factor;
      float_sw4 duzdy = (U[g0](3,i0,j0+1,k0) - U[g0](3,i0,j0-1,k0))*factor;
      float_sw4 duxdz = (U[g0](1,i0,j0,k0+1) - U[g0](1,i0,j0,k0-1))*factor;
      float_sw4 duydz = (U[g0](2,i0,j0,k0+1) - U[g0](2,i0,j0,k0-1))*factor;
//       if( m_xycomponent )
//       {
      uRec[0] = ( duzdy-duydz );
      uRec[1] = ( duxdz-duzdx );
      uRec[2] = ( duydx-duxdy );
//       }
//       else
//       {
// 	 float_sw4 uns = m_thynrm*(duzdy-duydz)-m_thxnrm*(duxdz-duzdx);
// 	 float_sw4 uew = m_salpha*(duzdy-duydz)+m_calpha*(duxdz-duzdx);
// 	 mRecordedUX.push_back( uew );
// 	 mRecordedUY.push_back( uns );
// 	 mRecordedUZ.push_back( -(duydx-duxdy) );
//       }
    }
    else // must be curvilinear
    {
//       int i=m_i0, j=m_j0, k=m_k0, g=m_grid0;
      float_sw4 factor = 0.5/sqrt(mJ(i0,j0,k0));
      float_sw4 duxdq = (U[g0](1,i0+1,j0,k0) - U[g0](1,i0-1,j0,k0));
      float_sw4 duydq = (U[g0](2,i0+1,j0,k0) - U[g0](2,i0-1,j0,k0));
      float_sw4 duzdq = (U[g0](3,i0+1,j0,k0) - U[g0](3,i0-1,j0,k0));
      float_sw4 duxdr = (U[g0](1,i0,j0+1,k0) - U[g0](1,i0,j0-1,k0));
      float_sw4 duydr = (U[g0](2,i0,j0+1,k0) - U[g0](2,i0,j0-1,k0));
      float_sw4 duzdr = (U[g0](3,i0,j0+1,k0) - U[g0](3,i0,j0-1,k0));
      float_sw4 duxds = (U[g0](1,i0,j0,k0+1) - U[g0](1,i0,j0,k0-1));
      float_sw4 duyds = (U[g0](2,i0,j0,k0+1) - U[g0](2,i0,j0,k0-1));
      float_sw4 duzds = (U[g0](3,i0,j0,k0+1) - U[g0](3,i0,j0,k0-1));
      float_sw4 duzdy = mMetric(1,i0,j0,k0)*duzdr+mMetric(3,i0,j0,k0)*duzds;
      float_sw4 duydz = mMetric(4,i0,j0,k0)*duyds;
      float_sw4 duxdz = mMetric(4,i0,j0,k0)*duxds;
      float_sw4 duzdx = mMetric(1,i0,j0,k0)*duzdq+mMetric(2,i0,j0,k0)*duzds;
      float_sw4 duydx = mMetric(1,i0,j0,k0)*duydq+mMetric(2,i0,j0,k0)*duyds;
      float_sw4 duxdy = mMetric(1,i0,j0,k0)*duxdr+mMetric(3,i0,j0,k0)*duxds;
//       if( m_xycomponent )
//       {
      uRec[0] = (duzdy-duydz)*factor;
      uRec[1] = (duxdz-duzdx)*factor;
      uRec[2] = (duydx-duxdy)*factor;
//       }
//       else
//       {
// 	 float_sw4 uns = m_thynrm*(duzdy-duydz)-m_thxnrm*(duxdz-duzdx);
// 	 float_sw4 uew = m_salpha*(duzdy-duydz)+m_calpha*(duxdz-duzdx);
// 	 mRecordedUX.push_back( uew*factor );
// 	 mRecordedUY.push_back( uns*factor );
// 	 mRecordedUZ.push_back( -(duydx-duxdy)*factor );
//       }
    }
  } // end Curl
  else if(mode == TimeSeries::Strains )
  {
     uRec.resize(6);
    if (g0 < mNumberOfCartesianGrids) // must be a Cartesian grid
    {
//       int i=m_i0, j=m_j0, k=m_k0, g=m_grid0;
      float_sw4 factor = 1.0/(2*mGridSize[g0]);
      float_sw4 duydx = (U[g0](2,i0+1,j0,k0) - U[g0](2,i0-1,j0,k0))*factor;
      float_sw4 duzdx = (U[g0](3,i0+1,j0,k0) - U[g0](3,i0-1,j0,k0))*factor;
      float_sw4 duxdy = (U[g0](1,i0,j0+1,k0) - U[g0](1,i0,j0-1,k0))*factor;
      float_sw4 duzdy = (U[g0](3,i0,j0+1,k0) - U[g0](3,i0,j0-1,k0))*factor;
      float_sw4 duxdz = (U[g0](1,i0,j0,k0+1) - U[g0](1,i0,j0,k0-1))*factor;
      float_sw4 duydz = (U[g0](2,i0,j0,k0+1) - U[g0](2,i0,j0,k0-1))*factor;
      float_sw4 duxdx = (U[g0](1,i0+1,j0,k0) - U[g0](1,i0-1,j0,k0))*factor;
      float_sw4 duydy = (U[g0](2,i0,j0+1,k0) - U[g0](2,i0,j0-1,k0))*factor;
      float_sw4 duzdz = (U[g0](3,i0,j0,k0+1) - U[g0](3,i0,j0,k0-1))*factor;
      uRec[0] = ( duxdx );
      uRec[1] = ( duydy );
      uRec[2] = ( duzdz );
      uRec[3] = ( 0.5*(duydx+duxdy) );
      uRec[4] = ( 0.5*(duzdx+duxdz) );
      uRec[5] = ( 0.5*(duydz+duzdy) );
   }
    else // must be curvilinear
   {
//       int i=m_i0, j=m_j0, k0=m_k00, g0=m_grid0;
      float_sw4 factor = 0.5/sqrt(mJ(i0,j0,k0));
      float_sw4 duxdq = (U[g0](1,i0+1,j0,k0) - U[g0](1,i0-1,j0,k0));
      float_sw4 duydq = (U[g0](2,i0+1,j0,k0) - U[g0](2,i0-1,j0,k0));
      float_sw4 duzdq = (U[g0](3,i0+1,j0,k0) - U[g0](3,i0-1,j0,k0));
      float_sw4 duxdr = (U[g0](1,i0,j0+1,k0) - U[g0](1,i0,j0-1,k0));
      float_sw4 duydr = (U[g0](2,i0,j0+1,k0) - U[g0](2,i0,j0-1,k0));
      float_sw4 duzdr = (U[g0](3,i0,j0+1,k0) - U[g0](3,i0,j0-1,k0));
      float_sw4 duxds = (U[g0](1,i0,j0,k0+1) - U[g0](1,i0,j0,k0-1));
      float_sw4 duyds = (U[g0](2,i0,j0,k0+1) - U[g0](2,i0,j0,k0-1));
      float_sw4 duzds = (U[g0](3,i0,j0,k0+1) - U[g0](3,i0,j0,k0-1));
      float_sw4 duzdy = (mMetric(1,i0,j0,k0)*duzdr+mMetric(3,i0,j0,k0)*duzds)*factor;
      float_sw4 duydz = (mMetric(4,i0,j0,k0)*duyds)*factor;
      float_sw4 duxdz = (mMetric(4,i0,j0,k0)*duxds)*factor;
      float_sw4 duzdx = (mMetric(1,i0,j0,k0)*duzdq+mMetric(2,i0,j0,k0)*duzds)*factor;
      float_sw4 duydx = (mMetric(1,i0,j0,k0)*duydq+mMetric(2,i0,j0,k0)*duyds)*factor;
      float_sw4 duxdy = (mMetric(1,i0,j0,k0)*duxdr+mMetric(3,i0,j0,k0)*duxds)*factor;
      float_sw4 duxdx = ( mMetric(1,i0,j0,k0)*(U[g0](1,i0+1,j0,k0) - U[g0](1,i0-1,j0,k0))+
		       mMetric(2,i0,j0,k0)*(U[g0](1,i0,j0,k0+1) - U[g0](1,i0,j0,k0-1)) )*factor;
      float_sw4 duydy = ( mMetric(1,i0,j0,k0)*(U[g0](2,i0,j0+1,k0) - U[g0](2,i0,j0-1,k0))+
		       mMetric(3,i0,j0,k0)*(U[g0](2,i0,j0,k0+1) - U[g0](2,i0,j0,k0-1)) )*factor;
      float_sw4 duzdz = ( mMetric(4,i0,j0,k0)*(U[g0](3,i0,j0,k0+1) - U[g0](3,i0,j0,k0-1)) )*factor;
      uRec[0] = ( duxdx );
      uRec[1] = ( duydy );
      uRec[2] = ( duzdz );
      uRec[3] = ( 0.5*(duydx+duxdy) );
      uRec[4] = ( 0.5*(duzdx+duxdz) );
      uRec[5] = ( 0.5*(duydz+duzdy) );
   }
  } // end Strains
  else if(mode == TimeSeries::DisplacementGradient )
  {
     uRec.resize(9);
     if (g0 < mNumberOfCartesianGrids) // must be a Cartesian grid
     {
//       int i=m_i0, j=m_j0, k=m_k0, g=m_grid0;
	float_sw4 factor = 1.0/(2*mGridSize[g0]);
	float_sw4 duydx = (U[g0](2,i0+1,j0,k0) - U[g0](2,i0-1,j0,k0))*factor;
	float_sw4 duzdx = (U[g0](3,i0+1,j0,k0) - U[g0](3,i0-1,j0,k0))*factor;
	float_sw4 duxdy = (U[g0](1,i0,j0+1,k0) - U[g0](1,i0,j0-1,k0))*factor;
	float_sw4 duzdy = (U[g0](3,i0,j0+1,k0) - U[g0](3,i0,j0-1,k0))*factor;
	float_sw4 duxdz = (U[g0](1,i0,j0,k0+1) - U[g0](1,i0,j0,k0-1))*factor;
	float_sw4 duydz = (U[g0](2,i0,j0,k0+1) - U[g0](2,i0,j0,k0-1))*factor;
	float_sw4 duxdx = (U[g0](1,i0+1,j0,k0) - U[g0](1,i0-1,j0,k0))*factor;
	float_sw4 duydy = (U[g0](2,i0,j0+1,k0) - U[g0](2,i0,j0-1,k0))*factor;
	float_sw4 duzdz = (U[g0](3,i0,j0,k0+1) - U[g0](3,i0,j0,k0-1))*factor;
	uRec[0] =  duxdx;
	uRec[1] =  duxdy;
	uRec[2] =  duxdz;
	uRec[3] =  duydx;
	uRec[4] =  duydy;
	uRec[5] =  duydz;
	uRec[6] =  duzdx;
	uRec[7] =  duzdy;
	uRec[8] =  duzdz;
     }
     else // must be curvilinear
     {
//       int i=m_i0, j=m_j0, k0=m_k00, g0=m_grid0;
	float_sw4 factor = 0.5/sqrt(mJ(i0,j0,k0));
	float_sw4 duxdq = (U[g0](1,i0+1,j0,k0) - U[g0](1,i0-1,j0,k0));
	float_sw4 duydq = (U[g0](2,i0+1,j0,k0) - U[g0](2,i0-1,j0,k0));
	float_sw4 duzdq = (U[g0](3,i0+1,j0,k0) - U[g0](3,i0-1,j0,k0));
	float_sw4 duxdr = (U[g0](1,i0,j0+1,k0) - U[g0](1,i0,j0-1,k0));
	float_sw4 duydr = (U[g0](2,i0,j0+1,k0) - U[g0](2,i0,j0-1,k0));
	float_sw4 duzdr = (U[g0](3,i0,j0+1,k0) - U[g0](3,i0,j0-1,k0));
	float_sw4 duxds = (U[g0](1,i0,j0,k0+1) - U[g0](1,i0,j0,k0-1));
	float_sw4 duyds = (U[g0](2,i0,j0,k0+1) - U[g0](2,i0,j0,k0-1));
	float_sw4 duzds = (U[g0](3,i0,j0,k0+1) - U[g0](3,i0,j0,k0-1));
	float_sw4 duzdy = (mMetric(1,i0,j0,k0)*duzdr+mMetric(3,i0,j0,k0)*duzds)*factor;
	float_sw4 duydz = (mMetric(4,i0,j0,k0)*duyds)*factor;
	float_sw4 duxdz = (mMetric(4,i0,j0,k0)*duxds)*factor;
	float_sw4 duzdx = (mMetric(1,i0,j0,k0)*duzdq+mMetric(2,i0,j0,k0)*duzds)*factor;
	float_sw4 duydx = (mMetric(1,i0,j0,k0)*duydq+mMetric(2,i0,j0,k0)*duyds)*factor;
	float_sw4 duxdy = (mMetric(1,i0,j0,k0)*duxdr+mMetric(3,i0,j0,k0)*duxds)*factor;
	float_sw4 duxdx = ( mMetric(1,i0,j0,k0)*(U[g0](1,i0+1,j0,k0) - U[g0](1,i0-1,j0,k0))+
		       mMetric(2,i0,j0,k0)*(U[g0](1,i0,j0,k0+1) - U[g0](1,i0,j0,k0-1)) )*factor;
	float_sw4 duydy = ( mMetric(1,i0,j0,k0)*(U[g0](2,i0,j0+1,k0) - U[g0](2,i0,j0-1,k0))+
		       mMetric(3,i0,j0,k0)*(U[g0](2,i0,j0,k0+1) - U[g0](2,i0,j0,k0-1)) )*factor;
	float_sw4 duzdz = ( mMetric(4,i0,j0,k0)*(U[g0](3,i0,j0,k0+1) - U[g0](3,i0,j0,k0-1)) )*factor;
	uRec[0] =  duxdx;
	uRec[1] =  duxdy;
	uRec[2] =  duxdz;
	uRec[3] =  duydx;
	uRec[4] =  duydy;
	uRec[5] =  duydz;
	uRec[6] =  duzdx;
	uRec[7] =  duzdy;
	uRec[8] =  duzdz;
     }

  } // end DisplacementGradient
  return;
}

//-----------------------------------------------------------------------
void EW::default_bcs( )
{
   for( int side=0 ; side < 6 ; side++ )
      mbcGlobalType[side] = bSuperGrid;
   mbcGlobalType[4] = bStressFree; // low-z is normally free surface
}

//-----------------------------------------------------------------------
void EW::buildGaussianHillTopography(float_sw4 amp, float_sw4 Lx, float_sw4 Ly, float_sw4 x0, float_sw4 y0)
{
   if (mVerbose >= 1 && (m_myrank == 0 ) )
    cout << "***inside buildGaussianHillTopography***"<< endl;

#define SQR(x) (x)*(x)
  int topLevel = mNumberOfGrids-1;

  float_sw4 x, y;

// copy data
  m_analytical_topo = true;
//  m_analytical_topo = false;
  m_GaussianAmp = amp;
  m_GaussianLx = Lx;
  m_GaussianLy = Ly;
  m_GaussianXc = x0;
  m_GaussianYc = y0;

  for (int i = m_iStart[topLevel]; i <= m_iEnd[topLevel]; ++i)
    for (int j = m_jStart[topLevel]; j <= m_jEnd[topLevel]; ++j)
    {
      x = (i-1)*mGridSize[topLevel];
      y = (j-1)*mGridSize[topLevel];
// positive topography  is up (negative z)
      mTopo(i,j,1) = m_GaussianAmp*exp(-SQR((x-m_GaussianXc)/m_GaussianLx) 
			               -SQR((y-m_GaussianYc)/m_GaussianLy));
    }

  for (int i = mTopoGridExt.m_ib ; i <= mTopoGridExt.m_ie ; ++i)
     for (int j = mTopoGridExt.m_jb ; j <= mTopoGridExt.m_je; ++j)
     {
	x = (i-1)*mGridSize[topLevel];
	y = (j-1)*mGridSize[topLevel];
// positive topography  is up (negative z)
	mTopoGridExt(i,j,1) = m_GaussianAmp*exp(-SQR((x-m_GaussianXc)/m_GaussianLx) 
						-SQR((y-m_GaussianYc)/m_GaussianLy)); 
     }
#undef SQR
}

//-----------------------------------------------------------------------
void EW::compute_minmax_topography( float_sw4& topo_zmin, float_sw4& topo_zmax )
{
   if( m_topography_exists )
   {
      int g = mNumberOfGrids-1;

      int i=m_iStart[g], j=m_jEnd[g];
// The z-coordinate points downwards, so positive topography (above sea level)
// gets negative z-values
      float_sw4 zMinLocal, zMaxLocal;
      zMaxLocal = zMinLocal = -mTopoGridExt(i,j,1);
// tmp
      int i_min_loc=i, i_max_loc=i;
      int j_min_loc=j, j_max_loc=j;
// end tmp
      int imin = mTopoGridExt.m_ib;
      int imax = mTopoGridExt.m_ie;
      int jmin = mTopoGridExt.m_jb;
      int jmax = mTopoGridExt.m_je;
      for (i= imin ; i<=imax ; i++)
	 for (j=jmin; j<=jmax ; j++)
	 {
	    if (-mTopoGridExt(i,j,1) > zMaxLocal)
	    {
	       zMaxLocal = -mTopoGridExt(i,j,1);
	       i_max_loc = i;
	       j_max_loc = j;
	    }
	    if (-mTopoGridExt(i,j,1) < zMinLocal)
	    {
	       zMinLocal = -mTopoGridExt(i,j,1);
	       i_min_loc = i;
	       j_min_loc = j;
	    }
	 }
      MPI_Allreduce( &zMinLocal, &topo_zmin, 1, m_mpifloat, MPI_MIN, m_cartesian_communicator);
      MPI_Allreduce( &zMaxLocal, &topo_zmax, 1, m_mpifloat, MPI_MAX, m_cartesian_communicator);
   }
   else
   {
      topo_zmin = topo_zmax = 0;
   }
}

//-----------------------------------------------------------------------
void EW::generate_grid()
{
   // Generate grid on domain: topography <= z <= zmax, 
   // The 2D grid on z=zmax, is given by ifirst <= i <= ilast, jfirst <= j <= jlast
   // spacing h. 
  if (!m_topography_exists ) return;
  
//  m_grid_interpolation_order = a_order;

  if (mVerbose >= 1 && (m_myrank==0) )
     cout << "***inside generate_grid***"<< endl;

// get the size from the top Cartesian grid
  int g = mNumberOfCartesianGrids-1;
  int ifirst = m_iStart[g];
  int ilast  = m_iEnd[g];
  int jfirst = m_jStart[g];
  int jlast  = m_jEnd[g];

  float_sw4 h = mGridSize[g]; // grid size must agree with top cartesian grid
  float_sw4 zMaxCart = m_zmin[g]; // bottom z-level for curvilinear grid

  int i, j;
  int gTop = mNumberOfGrids-1;
  int Nz = m_kEnd[gTop] - m_ghost_points;

  if(mVerbose > 4 &&  (m_myrank == 0 ) )
  {
    printf("generate_grid: Number of grid points in curvilinear grid = %i, kStart = %i, kEnd = %i\n", 
	Nz, m_kStart[gTop], m_kEnd[gTop]);
  }

// generate the grid by calling the curvilinear mapping function
  float_sw4 X0, Y0, Z0;
  int k;
  for (k=m_kStart[gTop]; k<=m_kEnd[gTop]; k++)
    for (j=m_jStart[gTop]; j<=m_jEnd[gTop]; j++)
      for (i=m_iStart[gTop]; i<=m_iEnd[gTop]; i++)
      {
	grid_mapping((float_sw4) i, (float_sw4) j, (float_sw4) k, X0, Y0, Z0);
	mX(i,j,k) = X0;
	mY(i,j,k) = Y0;
	mZ(i,j,k) = Z0;
      }
  communicate_array( mZ, gTop );

// calculate min and max((mZ(i,j,k)-mZ(i,j,k-1))/h) for k=Nz
  k = Nz;
  float_sw4 hRatio;
  float_sw4 mZmin = 1.0e9, mZmax=0;
  for (j=m_jStart[gTop]; j<=m_jEnd[gTop]; j++)
    for (i=m_iStart[gTop]; i<=m_iEnd[gTop]; i++)
    {
       hRatio = (mZ(i,j,k)-mZ(i,j,k-1))/mGridSize[gTop];
	if (hRatio < mZmin) mZmin = hRatio;
	if (hRatio > mZmax) mZmax = hRatio;
    }
  float_sw4 zMinGlobal, zMaxGlobal;
  MPI_Allreduce( &mZmin, &zMinGlobal, 1, m_mpifloat, MPI_MIN, m_cartesian_communicator);
  MPI_Allreduce( &mZmax, &zMaxGlobal, 1, m_mpifloat, MPI_MAX, m_cartesian_communicator);
   printf("Curvilinear/Cartesian interface (k=Nz-1): Min grid size ratio - 1 = %e, max ratio z - 1 = %e, top grid # = %i Local = %e %e Type = %s\n", 
	  zMinGlobal-1., zMaxGlobal-1., gTop,mZmin,mZmax,mZ.status());
  if(mVerbose > 3 &&  (m_myrank == 0) )
  {
    printf("Curvilinear/Cartesian interface (k=Nz-1): Min grid size ratio - 1 = %e, max ratio z - 1 = %e, top grid # = %i\n", 
           zMinGlobal-1., zMaxGlobal-1., gTop);
  }
}

//---------------------------------------------------------
void EW::setup_metric()
{
   if (!m_topography_exists ) return;
   if (mVerbose >= 1 && (m_myrank == 0))
      cout << "***inside setup_metric***"<< endl;
   int g=mNumberOfGrids-1;
   int Bx=m_iStart[g];
   int By=m_jStart[g];
   int Bz=m_kStart[g];
   int Nx=m_iEnd[g];
   int Ny=m_jEnd[g];
   int Nz=m_kEnd[g];

   if( m_analytical_topo && m_use_analytical_metric )
   {
      // Gaussian hill topography, analytical expressions for metric derivatives.
      int nxg = m_global_nx[g];
      int nyg = m_global_ny[g];
      int nzg = m_global_nz[g];
      float_sw4 h= mGridSize[g];   
      float_sw4 zmax = m_zmin[g-1] - (nzg-1)*h*(1-m_zetaBreak);
#ifdef SW4_CROUTINES
      if( m_corder )
	 metricexgh_rev( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
			 m_global_nz[g], mX.c_ptr(), mY.c_ptr(), mZ.c_ptr(), mMetric.c_ptr(), mJ.c_ptr(),
			 m_grid_interpolation_order, m_zetaBreak, zmax, m_GaussianAmp, m_GaussianXc, 
			 m_GaussianYc, m_GaussianLx, m_GaussianLy );
      else
	 metricexgh( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
		     m_global_nz[g], mX.c_ptr(), mY.c_ptr(), mZ.c_ptr(), mMetric.c_ptr(), mJ.c_ptr(),
		     m_grid_interpolation_order, m_zetaBreak, zmax, m_GaussianAmp, m_GaussianXc, 
		     m_GaussianYc, m_GaussianLx, m_GaussianLy );
#else
      F77_FUNC(metricexgh,METRICEXGH)( &Bx, &Nx, &By, &Ny, &Bz, &Nz, &nxg, &nyg, &nzg, mX.c_ptr(), mY.c_ptr(), mZ.c_ptr(),
				       mMetric.c_ptr(), mJ.c_ptr(), &m_grid_interpolation_order, &m_zetaBreak, &zmax, 
				       &m_GaussianAmp, &m_GaussianXc, &m_GaussianYc, &m_GaussianLx, &m_GaussianLy );
#endif      
   }
   else
   {
     int ierr=0;
#ifdef SW4_CROUTINES
     if( m_corder )
	ierr = metric_rev( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
		    mX.c_ptr(), mY.c_ptr(), mZ.c_ptr(), mMetric.c_ptr(), mJ.c_ptr() );
     else
	ierr = metric( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
		mX.c_ptr(), mY.c_ptr(), mZ.c_ptr(), mMetric.c_ptr(), mJ.c_ptr() );
#else
     F77_FUNC(metric,METRIC)( &Bx, &Nx, &By, &Ny, &Bz, &Nz, mX.c_ptr(), mY.c_ptr(), mZ.c_ptr(),
			      mMetric.c_ptr(), mJ.c_ptr(), &ierr );
#endif
     CHECK_INPUT(ierr==0, "Problems calculating the metric coefficients");
   }
   communicate_array( mMetric, mNumberOfGrids-1 );
   communicate_array( mJ, mNumberOfGrids-1 );

   //   if( m_analytical_topo && !m_use_analytical_metric && mVerbose > 3 )
   //      // Test metric derivatives if available
   //      metric_derivatives_test( );

   float_sw4 minJ, maxJ;
#ifdef SW4_CROUTINES
   gridinfo( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
	     mMetric.c_ptr(), mJ.c_ptr(), minJ, maxJ );
#else   
   F77_FUNC(gridinfo,GRIDINFO)(&Bx, &Nx, &By, &Ny, &Bz, &Nz, mMetric.c_ptr(), mJ.c_ptr(), &minJ, &maxJ );
#endif
   float_sw4 minJglobal, maxJglobal;
   MPI_Allreduce( &minJ, &minJglobal, 1, m_mpifloat, MPI_MIN, m_cartesian_communicator);
   MPI_Allreduce( &maxJ, &maxJglobal, 1, m_mpifloat, MPI_MAX, m_cartesian_communicator);
   if (mVerbose>3 && (m_myrank == 0))
      printf("*** Jacobian of metric: minJ = %e maxJ = %e\n", minJglobal, maxJglobal);
   //// just save the results for now... do the sanity check later
   //   m_minJacobian= minJglobal;
   //   m_maxJacobian= maxJglobal;
}

//-----------------------------------------------------------------------
bool EW::find_topo_zcoord_owner( float_sw4 X, float_sw4 Y, float_sw4& Ztopo )
{
   bool success = true;
   if ( m_topography_exists )
   {
      float_sw4 h = mGridSize[mNumberOfGrids-1];
      float_sw4 q, r;
      q = X/h + 1.0;
      r = Y/h + 1.0;
// evaluate elevation of topography on the grid
      if (!interpolate_topography(q, r, Ztopo, true))
      {
	 cerr << "Unable to evaluate topography at" << " X= " << X << " Y= " << Y << endl;
	 cerr << "Setting topography to ZERO" << endl;
	 Ztopo = 0;
	 success = false;
      }
   }
   else
   {
      Ztopo = 0; // no topography
   }
   return success;
}

//-----------------------------------------------------------------------
bool EW::find_topo_zcoord_all( float_sw4 X, float_sw4 Y, float_sw4& Ztopo )
{
   bool success = true;
   if (m_topography_exists )
   {
      float_sw4 h = mGridSize[mNumberOfGrids-1];
      float_sw4 q, r;
      q = X/h + 1.0;
      r = Y/h + 1.0;
      float_sw4 Ztopoloc;
// evaluate elevation of topography on the grid
      if (!interpolate_topography(q, r, Ztopoloc, true))
      {
	 Ztopoloc = -1e38;
      }
      MPI_Allreduce( &Ztopoloc, &Ztopo, 1, m_mpifloat, MPI_MAX, m_cartesian_communicator );
      success = Ztopo > -1e38;
   }
   else
   {
      Ztopo = 0; // no topography
      success = true;
   }
   return success;
}

//-----------------------------------------------------------------------
bool less_than( GridPointSource* ptsrc1, GridPointSource* ptsrc2 )
{
   return ptsrc1->m_key < ptsrc2->m_key;
}

//-----------------------------------------------------------------------
void EW::sort_grid_point_sources()
{
   size_t* gptr = new size_t[mNumberOfGrids];
   gptr[0] = 0;
   for(int g=0 ; g < mNumberOfGrids-1 ; g++ )
   {
      gptr[g+1] = gptr[g] + static_cast<size_t>((m_iEnd[g]-m_iStart[g]+1))*
	 (m_jEnd[g]-m_jStart[g]+1)*(m_kEnd[g]-m_kStart[g]+1);
   }
   size_t* ni   = new size_t[mNumberOfGrids];
   size_t* nij  = new size_t[mNumberOfGrids];
   for(int g=0 ; g < mNumberOfGrids ; g++ )
   {
      ni[g] = (m_iEnd[g]-m_iStart[g]+1);
      nij[g] = ni[g]*(m_jEnd[g]-m_jStart[g]+1);
   }
   for( int s=0 ; s < m_point_sources.size() ; s++ )
   {
      int g = m_point_sources[s]->m_grid;
      size_t key = gptr[g] + (m_point_sources[s]->m_i0-m_iStart[g]) +
	 ni[g]*(m_point_sources[s]->m_j0-m_jStart[g]) +
	 nij[g]*(m_point_sources[s]->m_k0-m_kStart[g]);
      m_point_sources[s]->set_sort_key(key);
   }
   delete[] gptr;
   delete[] ni;
   delete[] nij;

   std::sort(m_point_sources.begin(), m_point_sources.end(), less_than );
   // set up array detecting sources belonging to idential points
   m_identsources.resize(1);
   m_identsources[0] = 0;
   int k = 0;
   while( m_identsources[k] < m_point_sources.size() )
   {
      int m = m_identsources[k];
      size_t key = m_point_sources[m]->m_key;
      while( m+1 < m_point_sources.size() && m_point_sources[m+1]->m_key == key )
	 m++;
      m_identsources.push_back(m+1);
      k++;
   }

   // Test   
   int nrsrc =m_point_sources.size();
   int nrunique = m_identsources.size()-1;

#ifdef CUDA_CODE
   ForceVector = newmanaged(nrunique*3);
   cudaMallocManaged(&ForceAddress, 3*nrunique*sizeof(float_sw4*));
   map[ForceAddress]= 3*nrunique*sizeof(float_sw4*);
   cudaDeviceSynchronize();
#else
   ForceVector = new double[nrunique*3];
   ForceAddress= new float_sw4*[nrunique*3];
#endif
   int nrsrctot, nruniquetot;
   MPI_Reduce( &nrsrc, &nrsrctot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
   MPI_Reduce( &nrunique, &nruniquetot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
   if( m_myrank == 0 )
   {
      cout << "number of grid point  sources = " << nrsrctot << endl;
      cout << "number of unique g.p. sources = " << nruniquetot << endl;
   }

   //   for( int s=0 ; s<m_identsources.size()-1 ; s++ )
   //      for( int i=m_identsources[s]; i< m_identsources[s+1] ; i++ )
   //	 std::cout << "src= " << i << " key=" << m_point_sources[i]->m_key <<
   //	    "grid= " << m_point_sources[i]->m_grid << " (i,j,k) = " <<
   //	    m_point_sources[i]->m_i0 << " " << 
   //	    m_point_sources[i]->m_j0 << " " << 
   //	    m_point_sources[i]->m_k0 << std::endl;
}

//-----------------------------------------------------------------------
void EW::copy_point_sources_to_gpu()
{
#ifdef SW4_CUDA
   cudaError_t retcode=cudaMalloc( (void**)&dev_point_sources, sizeof(GridPointSource*)*m_point_sources.size() );
   if( cudaSuccess != retcode )
      cout << "Error EW::copy_point_sources_to_gpu, cudaMalloc, retcode = " << cudaGetErrorString(retcode) << endl;
   GridPointSource** hdev_src = new GridPointSource*[m_point_sources.size()];
   for( int s=0 ; s < m_point_sources.size() ; s++ )
   {
      GridPointSource* dev_local;
      retcode=cudaMalloc( (void**)&dev_local, sizeof(GridPointSource) );
      if( cudaSuccess != retcode )
	 cout << "Error EW::copy_point_sources_to_gpu, cudaMalloc, 2, retcode = " << cudaGetErrorString(retcode) << endl;
      retcode = cudaMemcpy( dev_local, m_point_sources[s], sizeof(GridPointSource),
			    cudaMemcpyHostToDevice );
      if( cudaSuccess != retcode )
	 cout << "Error EW::copy_point_sources_to_gpu, cudaMemcpy, retcode = " <<
	    cudaGetErrorString(retcode) << endl;
      hdev_src[s] = dev_local;
   }
   retcode = cudaMemcpy( dev_point_sources, hdev_src, sizeof(GridPointSource*)*m_point_sources.size(), cudaMemcpyHostToDevice);
   if( cudaSuccess != retcode )
      cout << "Error EW::copy_point_sources_to_gpu, cudaMemcpy, 2, retcode = " <<
	 cudaGetErrorString(retcode) << endl;
   retcode = cudaMalloc( (void**)&dev_identsources, sizeof(int)*m_identsources.size() );
   if( cudaSuccess != retcode )
      cout << "Error EW::copy_point_sources_to_gpu, cudaMalloc, 3, retcode = " << cudaGetErrorString(retcode) << endl;
   retcode = cudaMemcpy( dev_identsources, &m_identsources[0], sizeof(int)*m_identsources.size(), cudaMemcpyHostToDevice );
   if( cudaSuccess != retcode )
      cout << "Error EW::copy_point_sources_to_gpu, cudaMemcpy, 3, retcode = " <<
	 cudaGetErrorString(retcode) << endl;
#endif
}
float_sw4* EW::newmanaged(size_t len){
   void *ptr;
   //std::cout<<"******Using overloaded new**********\n";
   //cudaMallocManaged(&ptr, len*sizeof(float_sw4),cudaMemAttachHost);
#ifdef CUDA_CODE
   cudaMallocManaged(&ptr, len*sizeof(float_sw4),cudaMemAttachGlobal);
   map[ptr]=len*sizeof(float_sw4);
   prefetched[ptr]=false;
   cudaDeviceSynchronize();
#else
   ptr=malloc(len*sizeof(float_sw4));
#endif
   return (float_sw4*)ptr;
}
float_sw4* EW::newmanagedh(size_t len){
   void *ptr;
   //std::cout<<"******Using overloaded new**********\n";
   //cudaMallocManaged(&ptr, len*sizeof(float_sw4),cudaMemAttachHost);
#ifdef CUDA_CODE
   cudaMalloc(&ptr, len*sizeof(float_sw4));
   //cudaMallocManaged(&ptr, len*sizeof(float_sw4),cudaMemAttachGlobal);
   map[ptr]=len*sizeof(float_sw4);
   //prefetched[ptr]=false;
   cudaDeviceSynchronize();
#else
   ptr=malloc(len*sizeof(float_sw4));
#endif
   return (float_sw4*)ptr;
}
void EW::delmanaged(float_sw4* &dptr){
#ifdef CUDA_CODE
  if (map.find(dptr)!=map.end())
    map.erase(dptr);
  else{
    std::cerr<<"ERROR EW::delamanged trying to delete non-managed pointer\n";
    abort();
  }
  cudaDeviceSynchronize();
  cudaFree((void*)dptr);
#else
  free((void*)dptr);
#endif
  dptr=NULL;
}
int EW::prefetch(void *ptr,int device){
#ifdef CUDA_CODE
  if (ptr==NULL) return 0;
  if (prefetched[ptr]) return 0;
  prefetched[ptr]=true;
  PUSH_RANGE("EW::PREFETCH",2);

  if (map.find(ptr)!=map.end()){
    cudaMemPrefetchAsync(ptr, 
			 map[ptr],
			 device,
			 0);

    SYNC_DEVICE
      POP_RANGE;
      return 0;
  }
  else{
    //std::cout<<"BAD PREFETCH\n";
    return 1;
  }
#endif
}
int EW::prefetchforced(void *ptr,int device){
#ifdef CUDA_CODE
  if (ptr==NULL) return 0;
  prefetched[ptr]=true;
  PUSH_RANGE("EW::PREFETCHFORCED",2);

  if (map.find(ptr)!=map.end()){
    cudaMemPrefetchAsync(ptr, 
			 map[ptr],
			 device,
			 0);

    SYNC_DEVICE
      POP_RANGE;
      return 0;
  }
  else{
    //std::cout<<"BAD PREFETCH\n";
    return 1;
  }
#endif
}

void EW::AMPI_Sendrecv(float_sw4* a, int scount, std::tuple<int,int,int> &sendt, int sendto, int stag,
		       float_sw4* b, int rcount, std::tuple<int,int,int> &recvt, int recvfrom, int rtag,
		       std::tuple<float_sw4*,float_sw4*> &buf,
		       MPI_Comm comm, MPI_Status *status){
  PUSH_RANGE_PAYLOAD("AMPI_SENDRECV",stag%6,stag);
  MPI_Request send_req,recv_req;
  
  int recv_count=std::get<0>(recvt)*std::get<1>(recvt);
  int send_count=std::get<0>(sendt)*std::get<1>(sendt);
  
  // std::cout<<"send_count "<<send_count<<" recv_count "<<recv_count<<"\n";
  if (recvfrom!=MPI_PROC_NULL)
    if (MPI_Irecv(std::get<1>(buf),recv_count,MPI_DOUBLE,recvfrom,rtag,comm,&recv_req)!=MPI_SUCCESS) std::cerr<<"MPI_Irecv failed in EW::AMPI_Sendrecv\n";

  if (sendto!=MPI_PROC_NULL){
    getbuffer_device(a,std::get<0>(buf),sendt);
    //std::cout<<"SENDING :: "<<sendto<<" ";
    //for(int i=0;i<10;i++) std::cout<<std::get<0>(buf)[i]<<" ";
    //std::cout<<"\n";
    if (MPI_Isend(std::get<0>(buf),send_count,MPI_DOUBLE,sendto,stag,comm,&send_req)!=MPI_SUCCESS) std::cerr<<"MPI_Isend failed in EW::AMPI_Sendrecv\n";
  }
  MPI_Status send_status,recv_status;

  if (recvfrom!=MPI_PROC_NULL){
    if (MPI_Wait(&recv_req,&recv_status)!=MPI_SUCCESS) std::cerr<<"MPI_WAIT RECV FAILED IN AMPI_SENDrecv\n";
    putbuffer_device(b,std::get<1>(buf),recvt);
    //std::cout<<"RECEIVING :: "<<recvfrom<<" ";
    //for(int i=0;i<10;i++) std::cout<<std::get<1>(buf)[i]<<" ";
    //std::cout<<"\n";
  }
  
  if (sendto!=MPI_PROC_NULL)
    if (MPI_Wait(&send_req,&send_status)!=MPI_SUCCESS) std::cerr<<"MPI_WAIT SEND FAILED IN AMPI_SENDrecv\n";

  POP_RANGE;
}


AMPI_Ret_type
EW::AMPI_SendrecvSplit(float_sw4* a, int scount, std::tuple<int,int,int> &sendt, int sendto, int stag,
		       float_sw4* b, int rcount, std::tuple<int,int,int> &recvt, int recvfrom, int rtag,
		       std::tuple<float_sw4*,float_sw4*> &buf,
		       MPI_Comm comm, MPI_Status *status){
  PUSH_RANGE_PAYLOAD("AMPI_SR_SPLIT",stag%6,stag);
  MPI_Request *send_req,*recv_req;
  send_req = new MPI_Request;
  recv_req = new MPI_Request;

  //*send_req = MPI_REQUEST_NULL;
  //*recv_req = MPI_REQUEST_NULL;

  int recv_count=std::get<0>(recvt)*std::get<1>(recvt);
  int send_count=std::get<0>(sendt)*std::get<1>(sendt);
  
  // std::cout<<"send_count "<<send_count<<" recv_count "<<recv_count<<"\n";
  if (recvfrom!=MPI_PROC_NULL){
    if (MPI_Irecv(std::get<1>(buf),recv_count,MPI_DOUBLE,recvfrom,rtag,comm,recv_req)!=MPI_SUCCESS) std::cerr<<"MPI_Irecv failed in EW::AMPI_SendrecvSplit\n";
  } else *recv_req = MPI_REQUEST_NULL;

  if (sendto!=MPI_PROC_NULL){
    getbuffer_device(a,std::get<0>(buf),sendt);
    //std::cout<<"SENDING :: "<<sendto<<" ";
    //for(int i=0;i<10;i++) std::cout<<std::get<0>(buf)[i]<<" ";
    //std::cout<<"\n";
    if (MPI_Isend(std::get<0>(buf),send_count,MPI_DOUBLE,sendto,stag,comm,send_req)!=MPI_SUCCESS) std::cerr<<"MPI_Isend failed in EW::AMPI_SendrecvSplit\n";
  } else *send_req = MPI_REQUEST_NULL;

  // if (recvfrom!=MPI_PROC_NULL){
  //   if (MPI_Wait(recv_req,&recv_status)!=MPI_SUCCESS) std::cerr<<"MPI_WAIT RECV FAILED IN AMPI_SENDrecv\n";
  //   putbuffer_device(b,std::get<1>(buf),recvt);
  // }
  
  //std::tuple<MPI_Status,double *, double *, std::tuple<int,int,int> &,MPI_Status> retval;
  auto retval = std::make_tuple(recv_req,b,std::get<1>(buf),recvt,send_req);
  POP_RANGE;
  return retval;
}

void EW::AMPI_SendrecvSync(std::vector<AMPI_Ret_type> &list){
  //std::cout<<"AMPI_SYNC.."<<MPI_REQUEST_NULL;
  PUSH_RANGE("AMPI_SYNC",4);
  MPI_Request *waitlist;
  int indx;
  waitlist = new MPI_Request[list.size()];
  //for (int i=0;i<list.size();i++) std::cout<<i<<"["<<*std::get<0>(list[i])<<","<<std::get<1>(list[i])<<","<<std::get<2>(list[i])<<","<<*std::get<4>(list[i])<<"]  ";
  //std::cout<<"EOL\n";
  const bool ordered=false;
  if (ordered){
    for (int i=0;i<list.size();i++){
      if (*std::get<0>(list[i])!=MPI_REQUEST_NULL){
	//std::cout<<i<<"["<<*std::get<0>(list[i])<<","<<std::get<1>(list[i])<<","<<std::get<2>(list[i])<<","<<*std::get<4>(list[i])<<"]  ";
	//std::cout<<"EOL\n";
	MPI_Status status;
	if (MPI_Wait(std::get<0>(list[i]),&status)!=MPI_SUCCESS) std::cout<<"ERROR:: MPI_Wait failed on line :"<<__LINE__<<" in "<<__FILE__<<"\n";
	putbuffer_device(std::get<1>(list[i]),std::get<2>(list[i]),std::get<3>(list[i]));
      }
    }
  } else {
    for (int i=0;i<list.size();i++) waitlist[i]=*(std::get<0>(list[i]));
    for (int i=0;i<list.size();i++){
      MPI_Status status;
      if (MPI_Waitany(list.size(),waitlist,&indx,&status)!=MPI_SUCCESS) std::cout<<"ERROR:: MPI failed on line :"<<__LINE__<<" in "<<__FILE__<<"\n";
      if (indx>=0)
	putbuffer_device(std::get<1>(list[indx]),std::get<2>(list[indx]),std::get<3>(list[indx]));
    }
  }
  for (int i=0;i<list.size();i++) waitlist[i]=*(std::get<4>(list[i]));
  MPI_Status *stats = new MPI_Status[list.size()];
  if (MPI_Waitall(list.size(),waitlist,stats)!=MPI_SUCCESS) std::cerr<<"ERROR in MPI_Waitany in"<<__LINE__<<"\n";
  delete [] stats;
  delete [] waitlist;
  for (int i=0;i<list.size();i++){
    //MPI_Request_free(std::get<0>(list[i]));
    //MPI_Request_free(std::get<4>(list[i]));
    delete std::get<0>(list[i]);
    delete std::get<4>(list[i]);
  }
  //std::cout<<"Done \n";
  POP_RANGE;
}
void EW::getbuffer(float_sw4 *data, float_sw4* buf, std::tuple<int,int,int> &mtype ){

  int count=std::get<0>(mtype);
  int bl = std::get<1>(mtype);
  int stride = std::get<2>(mtype);
  for (int i=0;i<count;i++){
    for (int k=0;k<bl;k++) buf[k+i*bl]=data[i*stride+k];
  }
}
void EW::getbuffer_device(float_sw4 *data, float_sw4* buf, std::tuple<int,int,int> &mtype ){

  int count=std::get<0>(mtype);
  int bl = std::get<1>(mtype);
  int stride = std::get<2>(mtype);
  //std::cout<<"getbuffer_device...";
  PUSH_RANGE_PAYLOAD("GET_BUFFER",1,count*bl);
  //PREFETCHFORCED(buf);
  forall<EXEC > (0,count,[=] RAJA_DEVICE(int i){
    for (int k=0;k<bl;k++) buf[k+i*bl]=data[i*stride+k];
  });

  SYNC_DEVICE;
  POP_RANGE;
  //std::cout<<"Done\n";
}

void EW::putbuffer(float_sw4 *data, float_sw4* buf, std::tuple<int,int,int> &mtype ){
  int count=std::get<0>(mtype);
  int bl = std::get<1>(mtype);
  int stride = std::get<2>(mtype);
  for (int i=0;i<count;i++){
    for (int k=0;k<bl;k++) data[i*stride+k]=buf[k+i*bl];
  }
  
}
void EW:: buffdiff(float_sw4* buf1, float_sw4*buf2,std::tuple<int,int,int> &mtype){
  int count=std::get<0>(mtype);
  int bl = std::get<1>(mtype);
  int stride = std::get<2>(mtype);
  int lc=0;
  for(int i=0;i<count;i++) for(int k=0;k<bl;k++) if (buf1[k+i*bl]!=buf2[k+i*bl]) {
      std:cout<<"DIFF "<<i<<" "<<k<<" "<<k+i*bl<<" "<<buf1[k+i*bl]<< " "<<buf2[k+i*bl]<<"\n";
	lc++;
	if (lc==2) return;
      }
}
void EW::putbuffer_device(float_sw4 *data, float_sw4* buf, std::tuple<int,int,int> &mtype ){
  int count=std::get<0>(mtype);
  int bl = std::get<1>(mtype);
  int stride = std::get<2>(mtype);
  //std::cout<<"putbuffer_device...";
  PUSH_RANGE_PAYLOAD("PUT_BUFFER",2,count*bl);
  //PREFETCHFORCED(buf);
  forall<EXEC > (0,count,[=] RAJA_DEVICE(int i){
      for (int k=0;k<bl;k++) data[i*stride+k]=buf[k+i*bl];
    });

  SYNC_DEVICE;
  //std::cout<<"Done\n";
  POP_RANGE;
}
//-----------------------------------------------------------------------
void EW::get_exact_point_source2( float_sw4* up, float_sw4 t, int g, Source& source,
				 int* wind )
{
   timeDep tD;
   if(!( source.getName() == "SmoothWave" || source.getName() == "VerySmoothBump" ||
	 source.getName() == "C6SmoothBump" || source.getName()== "Gaussian") )
   {
      cout << "EW::get_exact_point_source: Error, time dependency must be SmoothWave, VerySmoothBump, C6SmoothBump, or Gaussian, not "
	   << source.getName() << endl;
      return;
   }
   else if( source.getName() == "SmoothWave" )
      tD = iSmoothWave;
   else if( source.getName() == "VerySmoothBump" )
      tD = iVerySmoothBump;
   else if( source.getName() == "C6SmoothBump" )
      tD = iC6SmoothBump;
   else
      tD = iGaussian;

   //   u.set_to_zero();
   // Assume constant material, sample it in middle of domain
   int imid = (m_iStart[g]+m_iEnd[g])/2;
   int jmid = (m_jStart[g]+m_jEnd[g])/2;
   int kmid = (m_kStart[g]+m_kEnd[g])/2;
   float_sw4 rho   = mRho[g](imid,jmid,kmid);
   float_sw4 beta  =  sqrt( mMu[g](imid,jmid,kmid)/rho);
   float_sw4 alpha =  sqrt( (2*mMu[g](imid,jmid,kmid)+mLambda[g](imid,jmid,kmid))/rho);

      //   double alpha = m_point_source_test->m_cp;
      //   double beta  = m_point_source_test->m_cs;
      //   double rho   = m_point_source_test->m_rho;

   float_sw4 x0    = source.getX0();
   float_sw4 y0    = source.getY0();
   float_sw4 z0    = source.getZ0();
   float_sw4 fr=source.getFrequency();
   float_sw4 time = (t-source.getOffset()) * source.getFrequency();
   if( tD == iGaussian )
   {
      fr = 1/fr;
      time = time*fr;
   }
   bool ismomentsource = source.isMomentSource();
   float_sw4 fx, fy, fz;
   float_sw4 mxx, myy, mzz, mxy, mxz, myz, m0;

   if( !ismomentsource )
   {
      source.getForces( fx, fy, fz );
   }
   else
   {
      source.getMoments( mxx, mxy, mxz, myy, myz, mzz );
      //      m0  = source.getAmplitude();
      m0 = 1;
   }
   //   bool curvilinear = topographyExists() && g == mNumberOfGrids-1;
   bool curvilinear = false;
   //   float_sw4* up = u.c_ptr();
   float_sw4 h   = mGridSize[g];
   float_sw4 eps = 1e-3*h;
   //size_t ind = 0;
   int imax, imin, jmax, jmin, kmax, kmin;
   if( wind == 0 )
   {
      imin = m_iStart[g];
      imax = m_iEnd[g];
      jmin = m_jStart[g];
      jmax = m_jEnd[g];
      kmin = m_kStart[g];
      kmax = m_kEnd[g];
   }
   else
   {
      imin = wind[0];
      imax = wind[1];
      jmin = wind[2];
      jmax = wind[3];
      kmin = wind[4];
      kmax = wind[5];
   }
   // Note: Use of ind, assumes loop is over the domain over which u is defined.
   //   for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
   //      for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
   //	 for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
   //   for( int k=kmin ; k <= kmax ; k++ )
   //   for( int j=jmin ; j <= jmax ; j++ )
   //	 for( int i=imin ; i <= imax ; i++ )
   int m_zmin_local= m_zmin[g];
   //PrintPointerAttributes(up);
   if( !ismomentsource )
     {
   forallN<EXEC_CARTBC, int, int,int>(
				  RangeSegment(kmin,kmax+1),
				  RangeSegment(jmin,jmax+1),
				  RangeSegment(imin,imax+1),
				  [=] RAJA_DEVICE(int k, int j, int i) {
	 
            float_sw4 x,y,z;
	    size_t ind = (i-imin)+(j-jmin)*(imax-imin+1)+(k-kmin)*(jmax-jmin+1)*(imax-imin+1);
	    //	    if( curvilinear )
	    //	    {
	    //               x = mX(i,j,k);
	    //	       y = mY(i,j,k);
	    //	       z = mZ(i,j,k);
	    //	    }
	    //	    else
	    {
	       x = (i-1)*h;
	       y = (j-1)*h;
	       z = (k-1)*h + m_zmin_local;
	    }

	       float_sw4 R = sqrt( (x - x0)*(x - x0) + (y - y0)*(y - y0) + (z - z0)*(z - z0) );
	       if( R < eps )
		  up[3*ind] = up[3*ind+1] = up[3*ind+2] = 0;
	       else
	       {
		  float_sw4 A, B;
		  if (tD == iSmoothWave)
		  {
		     A = ( 1/pow(alpha,2) * SmoothWave(time, fr*R, alpha) - 1/pow(beta,2) * SmoothWave(time, fr*R, beta) +
			   3/pow(fr*R,2) * SmoothWave_x_T_Integral(time, fr*R, alpha, beta) ) / (4*M_PI*rho*R*R*R)  ;
	  
		     B = ( 1/pow(beta,2) * SmoothWave(time, fr*R, beta) -
			   1/pow(fr*R,2) * SmoothWave_x_T_Integral(time, fr*R, alpha, beta) ) / (4*M_PI*rho*R) ;
		  }
		  else if (tD == iVerySmoothBump)
		  {
		     A = ( 1/pow(alpha,2) * VerySmoothBump(time, fr*R, alpha) - 1/pow(beta,2) * VerySmoothBump(time, fr*R, beta) +
			   3/pow(fr*R,2) * VerySmoothBump_x_T_Integral(time, fr*R, alpha, beta) ) / (4*M_PI*rho*R*R*R)  ;
		     
		     B = ( 1/pow(beta,2) * VerySmoothBump(time, fr*R, beta) -
			   1/pow(fr*R,2) * VerySmoothBump_x_T_Integral(time, fr*R, alpha, beta) ) / (4*M_PI*rho*R) ;
		  }
		  else if (tD == iC6SmoothBump)
		  {
		     A = ( 1/pow(alpha,2) * C6SmoothBump(time, fr*R, alpha) - 1/pow(beta,2) * C6SmoothBump(time, fr*R, beta) +
			   3/pow(fr*R,2) * C6SmoothBump_x_T_Integral(time, fr*R, alpha, beta) ) / (4*M_PI*rho*R*R*R)  ;
		     
		     B = ( 1/pow(beta,2) * C6SmoothBump(time, fr*R, beta) -
			   1/pow(fr*R,2) * C6SmoothBump_x_T_Integral(time, fr*R, alpha, beta) ) / (4*M_PI*rho*R) ;
		  }
                  else if( tD == iGaussian )
		  {
		     A = ( 1/pow(alpha,2) * Gaussian(time, R, alpha,fr) - 1/pow(beta,2) * Gaussian(time, R, beta,fr) +
			   3/pow(R,2) * Gaussian_x_T_Integral(time, R, fr, alpha, beta) ) / (4*M_PI*rho*R*R*R)  ;
		     
		     B = ( 1/pow(beta,2) * Gaussian(time, R, beta,fr) -
			   1/pow(R,2) * Gaussian_x_T_Integral(time, R, fr, alpha, beta) ) / (4*M_PI*rho*R) ;
		  }
		  up[3*ind]   = ( (x - x0)*(x - x0)*fx + (x - x0)*(y - y0)*fy + (x - x0)*(z - z0)*fz )*A + fx*B;
		  up[3*ind+1] = ( (y - y0)*(x - x0)*fx + (y - y0)*(y - y0)*fy + (y - y0)*(z - z0)*fz )*A + fy*B;
		  up[3*ind+2] = ( (z - z0)*(x - x0)*fx + (z - z0)*(y - y0)*fy + (z - z0)*(z - z0)*fz )*A + fz*B;
	       }
				  });}
	    else 
	      {
		forallN<EXEC_CARTBC, int, int,int>(
						   RangeSegment(kmin,kmax+1),
						   RangeSegment(jmin,jmax+1),
						   RangeSegment(imin,imax+1),
						   [=] RAJA_DEVICE(int k, int j, int i) {
						     size_t ind = (i-imin)+(j-jmin)*(imax-imin+1)+(k-kmin)*(jmax-jmin+1)*(imax-imin+1);
						     up[3*ind] = up[3*ind+1] = up[3*ind+2] = 0;
						     // Here, ismomentsource == true
						     float_sw4 x,y,z;
						     {
						     x = (i-1)*h;
	       y = (j-1)*h;
	       z = (k-1)*h + m_zmin_local;
	    }
	       float_sw4 R = sqrt( (x - x0)*(x - x0) + (y - y0)*(y - y0) + (z - z0)*(z - z0) );
	       if( R < eps )
	       {
	       	  up[3*ind] = up[3*ind+1] = up[3*ind+2] = 0;
	       }
	       else
	       {
		  float_sw4 A, B, C, D, E;
		  if (tD == iSmoothWave)
		  {
		     A = SmoothWave(time, R, alpha);
		     B = SmoothWave(time, R, beta);
		     C = SmoothWave_x_T_Integral(time, R, alpha, beta);
		     D = d_SmoothWave_dt(time, R, alpha) / pow(alpha,3) / R;
		     E = d_SmoothWave_dt(time, R, beta) / pow(beta,3) / R;
		  }
		  else if (tD == iVerySmoothBump)
		  {
		     A = VerySmoothBump(time, R, alpha);
		     B = VerySmoothBump(time, R, beta);
		     C = VerySmoothBump_x_T_Integral(time, R, alpha, beta);
		     D = d_VerySmoothBump_dt(time, R, alpha) / pow(alpha,3) / R;
		     E = d_VerySmoothBump_dt(time, R, beta) / pow(beta,3) / R;
		  }
		  else if (tD == iC6SmoothBump)
		  {
		     A = C6SmoothBump(time, R, alpha);
		     B = C6SmoothBump(time, R, beta);
		     C = C6SmoothBump_x_T_Integral(time, R, alpha, beta);
		     D = d_C6SmoothBump_dt(time, R, alpha) / pow(alpha,3) / R;
		     E = d_C6SmoothBump_dt(time, R, beta) / pow(beta,3) / R;
		  }
		  else if (tD == iGaussian)
		  {
		     A = Gaussian(time, R, alpha,fr);
		     B = Gaussian(time, R, beta,fr);
		     C = Gaussian_x_T_Integral(time, R, fr,alpha, beta);
		     D = d_Gaussian_dt(time, R, alpha,fr) / pow(alpha,3) / R;
		     E = d_Gaussian_dt(time, R, beta,fr) / pow(beta,3) / R;
		  }
		  up[3*ind] += 
	// m_xx*G_xx,x
		     + m0*mxx/(4*M_PI*rho)*
		     ( 
		      + 3*(x-x0)*(x-x0)*(x-x0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      - 2*(x-x0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      + 3*(x-x0)*(x-x0) / pow(R,5) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))
	 
		      + ( 15*(x-x0)*(x-x0)*(x-x0) / pow(R,7) - 6*(x-x0) / pow(R,5) ) * C
	 
		      + (x-x0)*(x-x0) / pow(R,3)* ((x-x0)*D - (x-x0)*E)
	 
		      - 1 / pow(R,3) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))

		      - 3*(x-x0) / pow(R,5) * C

		      + (x-x0) / (pow(R,3)*pow(beta,2)) * B

		      + 1 / R * (x-x0)*E
		      );
		  up[3*ind] +=
		     // m_yy*G_xy,y
		     + m0*myy/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(y-y0)*(y-y0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      - (x-x0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (x-x0)*(y-y0) / pow(R,3)* ((y-y0)*D - (y-y0)*E)

		      + 3*(x-x0)*(y-y0) / pow(R,5) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))

		      + ( 15*(x-x0)*(y-y0)*(y-y0) / pow(R,7) - 3*(x-x0) / pow(R,5) ) * C
		      );
		  up[3*ind] +=
		     // m_zz*G_xz,z
		     + m0*mzz/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(z-z0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))

		      - (x-x0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (x-x0)*(z-z0) / pow(R,3)* ((z-z0)*D - (z-z0)*E)

		      + 3*(x-x0)*(z-z0) / pow(R,5) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))

		      + ( 15*(x-x0)*(z-z0)*(z-z0) / pow(R,7) - 3*(x-x0) / pow(R,5) ) * C
		      );
		  up[3*ind] +=
		     // m_xy*G_xy,x
		     + m0*mxy/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(x-x0)*(y-y0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))

		      - (y-y0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (x-x0)*(y-y0) / pow(R,3)* ((x-x0)*D - (x-x0)*E)

		      + 3*(x-x0)*(y-y0) / pow(R,5) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))

		      + ( 15*(x-x0)*(x-x0)*(y-y0) / pow(R,7) - 3*(y-y0) / pow(R,5) ) * C
		      );
		  up[3*ind] +=
		     // m_xy*G_xx,y
		     + m0*mxy/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(x-x0)*(y-y0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      + 3*(x-x0)*(x-x0) / pow(R,5) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))
	 
		      + 15*(x-x0)*(x-x0)*(y-y0) / pow(R,7) * C
	 
		      + (x-x0)*(x-x0) / pow(R,3)* ((y-y0)*D - (y-y0)*E)
	 
		      - 1 / pow(R,3) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))

		      - 3*(y-y0) / pow(R,5) * C

		      + (y-y0) / (pow(R,3)*pow(beta,2)) * B

		      + 1 / R * (y-y0)*E
		      );
		  up[3*ind] +=
		     // m_xz*G_xz,x
		     + m0*mxz/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(x-x0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))

		      - (z-z0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (x-x0)*(z-z0) / pow(R,3)* ((x-x0)*D - (x-x0)*E)

		      + 3*(x-x0)*(z-z0) / pow(R,5) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))

		      + ( 15*(x-x0)*(x-x0)*(z-z0) / pow(R,7) - 3*(z-z0) / pow(R,5) ) * C
		      );
		  up[3*ind] +=
		     // m_yz*G_xz,y
		     + m0*myz/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(y-y0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (x-x0)*(z-z0) / pow(R,3)* ((y-y0)*D - (y-y0)*E)

		      + 3*(x-x0)*(z-z0) / pow(R,5) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))

		      + 15*(x-x0)*(y-y0)*(z-z0) / pow(R,7) * C
		      );
		  up[3*ind] +=
		     // m_xz*G_xx,z
		     + m0*mxz/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(x-x0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      + 3*(x-x0)*(x-x0) / pow(R,5) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))
	 
		      + 15*(x-x0)*(x-x0)*(z-z0) / pow(R,7) * C
	 
		      + (x-x0)*(x-x0) / pow(R,3)* ((z-z0)*D - (z-z0)*E)
	 
		      - 1 / pow(R,3) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))

		      - 3*(z-z0) / pow(R,5) * C

		      + (z-z0) / (pow(R,3)*pow(beta,2)) * B

		      + 1 / R * (z-z0)*E
		      );
		  up[3*ind] +=
		     // m_yz*G_yx,z
		     + m0*myz/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(y-y0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (x-x0)*(y-y0) / pow(R,3)* ((z-z0)*D - (z-z0)*E)

		      + 3*(x-x0)*(y-y0) / pow(R,5) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))

		      + 15*(x-x0)*(y-y0)*(z-z0) / pow(R,7) * C
		      );
		  //------------------------------------------------------------
		  up[3*ind+1] += 
		     // m_xx*G_xy,x
		     m0*mxx/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(x-x0)*(y-y0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))

		      - (y-y0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (x-x0)*(y-y0) / pow(R,3)* ((x-x0)*D - (x-x0)*E)

		      + 3*(x-x0)*(y-y0) / pow(R,5) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))

		      + ( 15*(x-x0)*(x-x0)*(y-y0) / pow(R,7) - 3*(y-y0) / pow(R,5) ) * C
		      );
		  up[3*ind+1] += 
		     // m_yy**G_yy,y
		     + m0*myy/(4*M_PI*rho)*
		     ( 
		      + 3*(y-y0)*(y-y0)*(y-y0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      - 2*(y-y0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      + 3*(y-y0)*(y-y0) / pow(R,5) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))
	 
		      + ( 15*(y-y0)*(y-y0)*(y-y0) / pow(R,7) - 6*(y-y0) / pow(R,5) ) * C
	 
		      + (y-y0)*(y-y0) / pow(R,3)* ((y-y0)*D - (y-y0)*E)
	 
		      - 1 / pow(R,3) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))

		      - 3*(y-y0) / pow(R,5) * C

		      + (y-y0) / (pow(R,3)*pow(beta,2)) * B

		      + 1 / R * (y-y0)*E
		      );
		  up[3*ind+1] += 
		     // m_zz*G_zy,z
		     + m0*mzz/(4*M_PI*rho)*
		     (
		      + 3*(z-z0)*(z-z0)*(y-y0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))

		      - (y-y0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (z-z0)*(y-y0) / pow(R,3)* ((z-z0)*D - (z-z0)*E)

		      + 3*(z-z0)*(y-y0) / pow(R,5) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))

		      + ( 15*(z-z0)*(z-z0)*(y-y0) / pow(R,7) - 3*(y-y0) / pow(R,5) ) * C
		      );
		  up[3*ind+1] += 
		     // m_xy*G_yy,x
		     + m0*mxy/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(y-y0)*(y-y0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      + 3*(y-y0)*(y-y0) / pow(R,5) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))
	  
		      + 15*(x-x0)*(y-y0)*(y-y0) / pow(R,7) * C
	  
		      + (y-y0)*(y-y0) / pow(R,3)* ((x-x0)*D - (x-x0)*E)
	  
		      - 1 / pow(R,3) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))
	  
		      - 3*(x-x0) / pow(R,5) * C
	  
		      + (x-x0) / (pow(R,3)*pow(beta,2)) * B
	  
		      + 1 / R * (x-x0)*E
		      );
		  up[3*ind+1] += 
		     // m_xz*G_zy,x
		     + m0*mxz/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(y-y0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	  
		      + (y-y0)*(z-z0) / pow(R,3)* ((x-x0)*D - (x-x0)*E)
	  
		      + 3*(y-y0)*(z-z0) / pow(R,5) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))
	  
		      + 15*(x-x0)*(y-y0)*(z-z0) / pow(R,7) * C
		      );
		  up[3*ind+1] += 
		     // m_xy*G_xy,y
		     + m0*mxy/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(y-y0)*(y-y0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	  
		      - (x-x0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))
	  
		      + (x-x0)*(y-y0) / pow(R,3)* ((y-y0)*D - (y-y0)*E)
	  
		      + 3*(x-x0)*(y-y0) / pow(R,5) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))
	  
		      + ( 15*(x-x0)*(y-y0)*(y-y0) / pow(R,7) - 3*(x-x0) / pow(R,5) ) * C
		      );
		  up[3*ind+1] += 
		     // m_yz*G_zy,y
		     + m0*myz/(4*M_PI*rho)*
		     (
		      + 3*(z-z0)*(y-y0)*(y-y0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	  
		      - (z-z0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))
	  
		      + (z-z0)*(y-y0) / pow(R,3)* ((y-y0)*D - (y-y0)*E)
	  
		      + 3*(z-z0)*(y-y0) / pow(R,5) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))
	  
		      + ( 15*(z-z0)*(y-y0)*(y-y0) / pow(R,7) - 3*(z-z0) / pow(R,5) ) * C
		      );
		  up[3*ind+1] += 
		     // m_xz*G_xy,z
		     + m0*mxz/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(y-y0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	  
		      + (x-x0)*(y-y0) / pow(R,3)* ((z-z0)*D - (z-z0)*E)
	  
		      + 3*(x-x0)*(y-y0) / pow(R,5) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))
	  
		      + 15*(x-x0)*(y-y0)*(z-z0) / pow(R,7) * C
		      );
		  up[3*ind+1] += 
		     // m_yz*G_yy,z
		     + m0*myz/(4*M_PI*rho)*
		     (
		      + 3*(z-z0)*(y-y0)*(y-y0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      + 3*(y-y0)*(y-y0) / pow(R,5) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))
	  
		      + 15*(z-z0)*(y-y0)*(y-y0) / pow(R,7) * C
	  
		      + (y-y0)*(y-y0) / pow(R,3)* ((z-z0)*D - (z-z0)*E)
	  
		      - 1 / pow(R,3) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))
	  
		      - 3*(z-z0) / pow(R,5) * C
	  
		      + (z-z0) / (pow(R,3)*pow(beta,2)) * B
	  
		      + 1 / R * (z-z0)*E
		      );
		  //------------------------------------------------------------
		  up[3*ind+2] += 
		     // m_xx*G_zx,x
		     + m0*mxx/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(x-x0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))

		      - (z-z0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (x-x0)*(z-z0) / pow(R,3)* ((x-x0)*D - (x-x0)*E)

		      + 3*(x-x0)*(z-z0) / pow(R,5) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))

		      + ( 15*(x-x0)*(x-x0)*(z-z0) / pow(R,7) - 3*(z-z0) / pow(R,5) ) * C
		      );
		  up[3*ind+2] += 
		     // m_yy*G_zy,y
		     + m0*myy/(4*M_PI*rho)*
		     (
		      + 3*(y-y0)*(y-y0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))

		      - (z-z0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (y-y0)*(z-z0) / pow(R,3)* ((y-y0)*D - (y-y0)*E)

		      + 3*(y-y0)*(z-z0) / pow(R,5) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))

		      + ( 15*(y-y0)*(y-y0)*(z-z0) / pow(R,7) - 3*(z-z0) / pow(R,5) ) * C
		      );
		  up[3*ind+2] += 
		     // m_zz**G_zz,z
		     + m0*mzz/(4*M_PI*rho)*
		     ( 
		      + 3*(z-z0)*(z-z0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      - 2*(z-z0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      + 3*(z-z0)*(z-z0) / pow(R,5) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))
	 
		      + ( 15*(z-z0)*(z-z0)*(z-z0) / pow(R,7) - 6*(z-z0) / pow(R,5) ) * C
	 
		      + (z-z0)*(z-z0) / pow(R,3)* ((z-z0)*D - (z-z0)*E)
	 
		      - 1 / pow(R,3) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))

		      - 3*(z-z0) / pow(R,5) * C

		      + (z-z0) / (pow(R,3)*pow(beta,2)) * B

		      + 1 / R * (z-z0)*E
		      );
		  up[3*ind+2] += 
		     // m_xy*G_zy,x
		     + m0*mxy/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(y-y0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	  
		      + (y-y0)*(z-z0) / pow(R,3)* ((x-x0)*D - (x-x0)*E)
	  
		      + 3*(y-y0)*(z-z0) / pow(R,5) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))
	  
		      + 15*(x-x0)*(y-y0)*(z-z0) / pow(R,7) * C
		      );
		  up[3*ind+2] += 
		     // m_xz**G_zz,x
		     + m0*mxz/(4*M_PI*rho)*
		     ( 
		      + 3*(x-x0)*(z-z0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      + 3*(z-z0)*(z-z0) / pow(R,5) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))
	 
		      + 15*(x-x0)*(z-z0)*(z-z0) / pow(R,7) * C
	 
		      + (z-z0)*(z-z0) / pow(R,3)* ((x-x0)*D - (x-x0)*E)
	 
		      - 1 / pow(R,3) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))

		      - 3*(x-x0) / pow(R,5) * C

		      + (x-x0) / (pow(R,3)*pow(beta,2)) * B

		      + 1 / R * (x-x0)*E
		      );
		  up[3*ind+2] += 
		     // m_xy*G_xz,y
		     + m0*mxy/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(y-y0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (x-x0)*(z-z0) / pow(R,3)* ((y-y0)*D - (y-y0)*E)

		      + 3*(x-x0)*(z-z0) / pow(R,5) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))

		      + 15*(x-x0)*(y-y0)*(z-z0) / pow(R,7) * C
		      );
		  up[3*ind+2] += 
		     // m_yz*G_zz,y
		     + m0*myz/(4*M_PI*rho)*
		     ( 
		      + 3*(y-y0)*(z-z0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      + 3*(z-z0)*(z-z0) / pow(R,5) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))
	 
		      + 15*(y-y0)*(z-z0)*(z-z0) / pow(R,7) * C
	 
		      + (z-z0)*(z-z0) / pow(R,3)* ((y-y0)*D - (y-y0)*E)
	 
		      - 1 / pow(R,3) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))

		      - 3*(y-y0) / pow(R,5) * C

		      + (y-y0) / (pow(R,3)*pow(beta,2)) * B

		      + 1 / R * (y-y0)*E
		      );
		  up[3*ind+2] += 
		     // m_xz*G_xz,z
		     + m0*mxz/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(z-z0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      - (x-x0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      + (x-x0)*(z-z0) / pow(R,3)* ((z-z0)*D - (z-z0)*E)
	 
		      + 3*(x-x0)*(z-z0) / pow(R,5) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))
	 
		      + ( 15*(x-x0)*(z-z0)*(z-z0) / pow(R,7) - 3*(x-x0) / pow(R,5) ) * C
		      );
		  up[3*ind+2] += 
		     // m_yz*G_yz,z
		     + m0*myz/(4*M_PI*rho)*
		     (
		      + 3*(z-z0)*(z-z0)*(y-y0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))

		      - (y-y0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (z-z0)*(y-y0) / pow(R,3)* ((z-z0)*D - (z-z0)*E)

		      + 3*(z-z0)*(y-y0) / pow(R,5) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))

		      + ( 15*(z-z0)*(z-z0)*(y-y0) / pow(R,7) - 3*(y-y0) / pow(R,5) ) * C
		      );
		  //printf("VALUE = %lf \n",up[3*ind+2]);
	       }
				  });}
	    //ind++;
				  
#ifdef CUDA_CODE
   cudaDeviceSynchronize();
#endif
}
