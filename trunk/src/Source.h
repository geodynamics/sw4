//-*-c++-*-
#ifndef EW_SOURCE_H
#define EW_SOURCE_H

#include <iostream>
#include <vector>
#include <string>
#include "TimeDep.h"

class EW;

class GridPointSource;


class Source
{
   friend std::ostream& operator<<(std::ostream& output, const Source& s);
public:
  
  Source(EW * a_ew, double amplitude, double frequency, double t0,
	 double x0, double y0, double z0,
	 double Mxx,
	 double Mxy,
	 double Mxz,
	 double Myy,
	 double Myz,
	 double Mzz,
	 timeDep tDep,
	 char *name,int ncyc=1);

  Source(EW * a_ew, double amplitude, double frequency, double t0,
         double x0, double y0, double z0,
         double Fx,
         double Fy,
         double Fz,
         timeDep tDep,
         const char *name,int ncyc=1);

 ~Source();

  bool gridPointSet() const;

  int m_i0, m_j0, m_k0; //C.B: renaming mN,mM,mL
  int m_grid;

  double getX0() const;
  double getY0() const;
  double getZ0() const;
  bool ignore() const;

  // Amplitude
  double getAmplitude() const;
  void setAmplitude(double amp);
  
  // Offset in time
  double getOffset() const;

  // Frequency
  double getFrequency() const;
  timeDep getTfunc() const {return mTimeDependence;}
  void setMaxFrequency(double max_freq);

  // Type of source
  bool isMomentSource() const;

  // discretize a time function at each time step and change the time function to be "Discrete()"
  void discretizeTimeFuncAndFilter(double tStart, double dt, int nSteps, double fc);

  double dt_to_resolve( int ppw ) const;
  int ppw_to_resolve( double dt ) const;

  const std::string& getName() const { return mName; };
  void correct_Z_level( );
  void limit_frequency( int ppw, double minvsoh );
  double compute_t0_increase( double t0_min );
  void adjust_t0( double dt0 );

  void set_grid_point_sources4( EW *a_EW, std::vector<GridPointSource*>& point_sources );
  void set_grid_point_sources4b( EW *a_EW, std::vector<GridPointSource*>& point_sources );
  void set_grid_point_sources( EW *a_EW, std::vector<GridPointSource*>& point_sources );

  void distribute_source_xyplane( EW *a_EW, std::vector<GridPointSource*>& point_sources, 
				  int g, int k, double wghz );
  void distribute_source_xyplane_mom( EW *a_EW, std::vector<GridPointSource*>& point_sources,
				      int g, int k, double wghz, double dwghz );

  void set_z_is_relative_to_topography( bool tf ) { m_zRelativeToTopography = tf; };
  void exact_testmoments( int kx[3], int ky[3], int kz[3], double momexact[3] );
  void getForces( double& fx, double& fy, double& fz ) const;
  void getMoments( double& mxx, double& myy, double& mzz, double& mxy, double& mxz, double& myz ) const;
   void printPointer(){std::cout << "Source pointer = "  << mPar << std::endl;}

 private:
  Source();
  void getsourcewgh7(double ai, double wgh[7] );
  void getsourcedwgh7(double ai, double wgh[7] );
  void getsourcewgh(double ai, double wgh[6] );
  void getsourcedwgh(double ai, double wgh[6] );
  double dist_d_dx_dirac(double x);
  double distributedhat(double x);
  double hat(double x);
  double d_qubic_dx(double x);
  double d_quadric_dx(double x);
  double d_hat_dx(double x);
  int sgn(double arg);
  
//  void initializeTimeFunction();
  double find_min_exponent();
  std::string mName;
  std::vector<double> mForces;
  double mAmp;
  bool mIsMomentSource;
  double mFreq, mT0;

  bool mGridPointSet;
  bool m_zRelativeToTopography;
  double mX0,mY0,mZ0;
  double* mPar;
  int mNcyc;
  bool mIgnore;

  timeDep mTimeDependence;
//  EW * mEW;
};

#endif
