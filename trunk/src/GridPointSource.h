//-*-c++-*-
#ifndef GRID_POINT_SOURCE_H
#define GRID_POINT_SOURCE_H

#include <iostream>
#include <vector>
#include <string>
#include "TimeDep.h"
#include "Source.h"
#include "Sarray.h"

class GridPointSource
{
   friend std::ostream& operator<<(std::ostream& output, const GridPointSource& s);
public:

  GridPointSource(double amplitude, double frequency, double t0,
		  int i0, int j0, int k0, int g,
		  double Fx, double Fy, double Fz,
		  timeDep tDep, int ncyc, double* jacobian=NULL );

 ~GridPointSource();

  int m_i0,m_j0,m_k0; // grid point index
  int m_grid;

  void getFxyz( double t, double* fxyz ) const;
  void getFxyztt( double t, double* fxyz ) const;
  void getFxyz_notime( double* fxyz ) const;

  // evaluate time fcn: RENAME to evalTimeFunc
  double getTimeFunc(double t) const;
  void limitFrequency(double max_freq);
  // discretize a time function at each time step and change the time function to be "Discrete()"
  void discretizeTimeFuncAndFilter(double tStart, double dt, int nSteps, double fc);
   void add_to_gradient( std::vector<Sarray>& kappa, std::vector<Sarray> & eta,
			 double t, double dt, double gradient[11], std::vector<double> & h );

 private:

  GridPointSource();

  void initializeTimeFunction();
  double mForces[3];
  double mAmp;
  double mFreq, mT0;

  timeDep mTimeDependence;
  double (*mTimeFunc)(double f, double t,double* par);
  double (*mTimeFunc_t)(double f, double t,double* par);
  double (*mTimeFunc_tt)(double f, double t,double* par);
  double (*mTimeFunc_ttt)(double f, double t,double* par);
  double (*mTimeFunc_om)(double f, double t,double* par);
  double (*mTimeFunc_omtt)(double f, double t,double* par);

  double* mPar;
  int mNcyc;
  double m_min_exponent;
  bool m_jacobian_known;
  double m_jacobian[27];
};

#endif
