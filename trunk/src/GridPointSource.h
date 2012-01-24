//-*-c++-*-
#ifndef GRID_POINT_SOURCE_H
#define GRID_POINT_SOURCE_H

#include <iostream>
#include <vector>
#include <string>
#include "TimeDep.h"
#include "Source.h"

class EW;

class GridPointSource
{
   friend std::ostream& operator<<(std::ostream& output, const GridPointSource& s);
public:

  GridPointSource(double amplitude, double frequency, double t0,
		  int N, int M, int L,int g,
		  double Fx, double Fy, double Fz,
		  timeDep tDep, int ncyc );

 ~GridPointSource();

  int m_i0,m_j0,m_k0; // grid point index
  int m_grid;

  void getFxyz( double t, double* fxyz ) const;

  void limitFrequency(double max_freq);

  // evaluate time fcn: RENAME to evalTimeFunc
  double getTimeFunc(double t) const;

  // discretize a time function at each time step and change the time function to be "Discrete()"
  void discretizeTimeFuncAndFilter(double tStart, double dt, int nSteps, double fc);

 private:

  GridPointSource();

  void initializeTimeFunction();
  double mForces[3];
  double mAmp;
  double mFreq, mT0;

  timeDep mTimeDependence;
  double (*mTimeFunc)(double f, double t,double* par);
  double* mPar;
  int mNcyc;
  double m_min_exponent;
};

#endif
