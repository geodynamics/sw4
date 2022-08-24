//-*-c++-*-
//  SW4 LICENSE
// # ----------------------------------------------------------------------
// # SW4 - Seismic Waves, 4th order
// # ----------------------------------------------------------------------
// # Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
// # Produced at the Lawrence Livermore National Laboratory. 
// # 
// # Written by:
// # N. Anders Petersson (petersson1@llnl.gov)
// # Bjorn Sjogreen      (sjogreen2@llnl.gov)
// # 
// # LLNL-CODE-643337 
// # 
// # All rights reserved. 
// # 
// # This file is part of SW4, Version: 1.0
// # 
// # Please also read LICENCE.txt, which contains "Our Notice and GNU General Public License"
// # 
// # This program is free software; you can redistribute it and/or modify
// # it under the terms of the GNU General Public License (as published by
// # the Free Software Foundation) version 2, dated June 1991. 
// # 
// # This program is distributed in the hope that it will be useful, but
// # WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
// # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
// # conditions of the GNU General Public License for more details. 
// # 
// # You should have received a copy of the GNU General Public License
// # along with this program; if not, write to the Free Software
// # Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA 
#ifndef GRID_POINT_SOURCE_H
#define GRID_POINT_SOURCE_H

#include <iostream>
#include <vector>
#include <string>
#include "sw4.h"
#include "TimeDep.h"
#include "Source.h"
#include "Sarray.h"
#include "Filter.h"

class GridPointSource
{
   friend std::ostream& operator<<(std::ostream& output, const GridPointSource& s);
public:

  GridPointSource(float_sw4 frequency, float_sw4 t0,
		  int i0, int j0, int k0, int g,
		  float_sw4 Fx, float_sw4 Fy, float_sw4 Fz,
		  timeDep tDep, int ncyc, float_sw4* pars, int npar, int* ipars, int nipar,
		  float_sw4* jacobian=NULL, float_sw4* dddp=NULL, float_sw4* hess1=NULL,
		  float_sw4* hess2=NULL, float_sw4* hess3=NULL );

 ~GridPointSource();

  int m_i0,m_j0,m_k0; // grid point index
  int m_grid;
  size_t m_key; // For sorting sources

  void getFxyz( float_sw4 t, float_sw4* fxyz ) const;
  void getFxyztt( float_sw4 t, float_sw4* fxyz ) const;
  void getFxyz_notime( float_sw4* fxyz ) const;

  // evaluate time fcn: RENAME to evalTimeFunc
  float_sw4 getTimeFunc(float_sw4 t) const;
  float_sw4 evalTimeFunc_t(float_sw4 t) const;
  float_sw4 evalTimeFunc_tt(float_sw4 t) const;
  float_sw4 evalTimeFunc_ttt(float_sw4 t) const;
  float_sw4 evalTimeFunc_tttt(float_sw4 t) const;

  void limitFrequency(float_sw4 max_freq);

  void add_to_gradient( std::vector<Sarray>& kappa, std::vector<Sarray> & eta,
			float_sw4 t, float_sw4 dt, float_sw4 gradient[11], std::vector<float_sw4> & h,
			vector<Sarray>& Jac, bool topography_exists );
  void add_to_hessian( std::vector<Sarray> & kappa, std::vector<Sarray> & eta,
		       float_sw4 t, float_sw4 dt, float_sw4 hessian[121], std::vector<float_sw4> & h );
  void set_derivative( int der, const float_sw4 dir[11] );
  void set_noderivative( );
  void print_info() const;
  void set_sort_key( size_t key );

   //// discretize a time function at each time step and change the time function to be "Discrete()"
   //  void discretizeTimeFuncAndFilter(float_sw4 tStart, float_sw4 dt, int nSteps, Filter *filter_ptr);

// make public for simplicity
  float_sw4 mFreq, mT0;

 private:

  GridPointSource();

  void initializeTimeFunction();
  float_sw4 mForces[3];
   //  float_sw4 mAmp;

  timeDep mTimeDependence;
   float_sw4 (*mTimeFunc)(float_sw4 f, float_sw4 t,float_sw4* par, int npar, int* ipar, int nipar );
  float_sw4 (*mTimeFunc_t)(float_sw4 f, float_sw4 t,float_sw4* par, int npar, int* ipar, int nipar );
  float_sw4 (*mTimeFunc_tt)(float_sw4 f, float_sw4 t,float_sw4* par, int npar, int* ipar, int nipar );
  float_sw4 (*mTimeFunc_ttt)(float_sw4 f, float_sw4 t,float_sw4* par, int npar, int* ipar, int nipar );
  float_sw4 (*mTimeFunc_om)(float_sw4 f, float_sw4 t,float_sw4* par, int npar, int* ipar, int nipar );
  float_sw4 (*mTimeFunc_omtt)(float_sw4 f, float_sw4 t,float_sw4* par, int npar, int* ipar, int nipar );
  float_sw4 (*mTimeFunc_tttt)(float_sw4 f, float_sw4 t,float_sw4* par, int npar, int* ipar, int nipar );
  float_sw4 (*mTimeFunc_tttom)(float_sw4 f, float_sw4 t,float_sw4* par, int npar, int* ipar, int nipar );
  float_sw4 (*mTimeFunc_ttomom)(float_sw4 f, float_sw4 t,float_sw4* par, int npar, int* ipar, int nipar );
  float_sw4 (*mTimeFunc_tom)(float_sw4 f, float_sw4 t,float_sw4* par, int npar, int* ipar, int nipar );
   float_sw4 (*mTimeFunc_omom)(float_sw4 f, float_sw4 t,float_sw4* par, int npar, int* ipar, int nipar );

  float_sw4* mPar;
  int* mIpar; 
  int  mNpar, mNipar;

  int mNcyc;
   //  float_sw4 m_min_exponent;
  int m_derivative;
  bool m_jacobian_known;
  float_sw4 m_jacobian[27];
  bool m_hessian_known;
  float_sw4 m_hesspos1[9], m_hesspos2[9], m_hesspos3[9], m_dddp[9]; 
  float_sw4 m_dir[11];
};

#endif
