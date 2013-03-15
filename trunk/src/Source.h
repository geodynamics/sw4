//-*-c++-*-
#ifndef EW_SOURCE_H
#define EW_SOURCE_H

#include <iostream>
#include <vector>
#include <string>
#include "TimeDep.h"

class EW;

class GridPointSource;

class Filter;

class Source
{
   friend std::ostream& operator<<(std::ostream& output, const Source& s);
public:
  
  Source(EW * a_ew, double frequency, double t0,
	 double x0, double y0, double z0,
	 double Mxx,
	 double Mxy,
	 double Mxz,
	 double Myy,
	 double Myz,
	 double Mzz,
	 timeDep tDep,
	 const char *name,int ncyc=1,
	 double* pars=NULL, int npars=0, int* ipars=NULL, int nipars=0 );

  Source(EW * a_ew, double frequency, double t0,
         double x0, double y0, double z0,
         double Fx,
         double Fy,
         double Fz,
         timeDep tDep,
         const char *name,int ncyc=1,
	 double* pars=NULL, int npars=0, int* ipars=NULL, int nipars=0 );

 ~Source();

  bool gridPointSet() const;

  int m_i0, m_j0, m_k0;
  int m_grid;

  double getX0() const;
  double getY0() const;
  double getZ0() const;
  bool ignore() const;

  // Amplitude
  double getAmplitude() const;
   //  void setAmplitude(double amp);
  
  // Offset in time
  double getOffset() const;

  // Frequency
  double getFrequency() const;
  timeDep getTfunc() const {return mTimeDependence;}
  void setMaxFrequency(double max_freq);

  // Type of source
  bool isMomentSource() const;

  // discretize a time function at each time step and change the time function to be "Discrete()"
   //  void discretizeTimeFuncAndFilter(double tStart, double dt, int nSteps, double fc);

  double dt_to_resolve( int ppw ) const;
  int ppw_to_resolve( double dt ) const;

  const std::string& getName() const { return mName; };
  void correct_Z_level( );
  void limit_frequency( int ppw, double minvsoh );
  double compute_t0_increase( double t0_min ) const;
  void adjust_t0( double dt0 );

  void set_grid_point_sources4( EW *a_EW, std::vector<GridPointSource*>& point_sources ) const;

  void set_z_is_relative_to_topography( bool tf ) { m_zRelativeToTopography = tf; };
  void exact_testmoments( int kx[3], int ky[3], int kz[3], double momexact[3] );
  void getForces( double& fx, double& fy, double& fz ) const;
  void getMoments( double& mxx, double& mxy, double& mxz, double& myy, double& myz, double& mzz ) const;
  void setMoments( double mxx, double mxy, double mxz, double myy, double myz, double mzz );
  void printPointer(){std::cout << "Source pointer = "  << mPar << std::endl;}
  void perturb( double h, int comp );
  void set_derivative( int der );
  void set_noderivative( );
  void set_dirderivative( double dir[11] );
  Source* copy( std::string a_name );
  void set_parameters( double x[11] );
  void setFrequency( double freq );
  void get_parameters( double x[11] ) const;
  void filter_timefunc( Filter* fi, double tstart, double dt, int nsteps );

 private:
  Source();

  void compute_metric_at_source( EW* a_EW, double q, double r, double s, int ic, int jc, int kc,
				 int g, double& zq, double& zr, double& zs ) const;
  int spline_interpolation( );
  void getsourcewgh(double ai, double wgh[6], double dwghda[6], double ddwghda[6] ) const;
  void getsourcedwgh(double ai, double wgh[6], double dwghda[6], double ddwghda[6] ) const;
  void getsourcewghlow(double ai, double wgh[6], double dwghda[6], double ddwghda[6] ) const;
  void getsourcedwghlow(double ai, double wgh[6], double dwghda[6], double ddwghda[6] ) const;
  void getmetwgh( double alph, double wgh[8] ) const;
  void getmetdwgh( double alph, double wgh[8] ) const;
  void getmetwgh7( double ai, double wgh[7] ) const;
  void getmetdwgh7( double ai, double wgh[7] ) const;

  double find_min_exponent() const;
  std::string mName;
  std::vector<double> mForces;
  bool mIsMomentSource;
  double mFreq, mT0;

  bool mGridPointSet;
  bool m_zRelativeToTopography;
  double mX0,mY0,mZ0;
  double* mPar;
  int* mIpar;
  int mNpar, mNipar;
  int mNcyc;
  int m_derivative;  
  timeDep mTimeDependence;
  double m_dir[11];
  bool m_is_filtered;
};

#endif
