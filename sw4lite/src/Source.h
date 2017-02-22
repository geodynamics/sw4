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
#ifndef EW_SOURCE_H
#define EW_SOURCE_H

#include "sw4.h"

#include <iostream>
#include <vector>
#include <string>


class EW;

class GridPointSource;

//class Filter;

class Source
{
   friend std::ostream& operator<<(std::ostream& output, const Source& s);
public:
  
  Source(EW * a_ew, float_sw4 frequency, float_sw4 t0,
	 float_sw4 x0, float_sw4 y0, float_sw4 z0,
	 float_sw4 Mxx,
	 float_sw4 Mxy,
	 float_sw4 Mxz,
	 float_sw4 Myy,
	 float_sw4 Myz,
	 float_sw4 Mzz,
	 timeDep tDep,
	 const char *name,
	 bool topodepth, 
	 int ncyc=1,
	 float_sw4* pars=NULL, int npars=0, int* ipars=NULL, int nipars=0, bool correctForMu=false );

  Source(EW * a_ew, float_sw4 frequency, float_sw4 t0,
         float_sw4 x0, float_sw4 y0, float_sw4 z0,
         float_sw4 Fx,
         float_sw4 Fy,
         float_sw4 Fz,
         timeDep tDep,
         const char *name,
	 bool topodepth, 
	 int ncyc=1,
	 float_sw4* pars=NULL, int npars=0, int* ipars=NULL, int nipars=0, bool correctForMu=false );

 ~Source();

  int m_i0, m_j0, m_k0;
  int m_grid;

  float_sw4 getX0() const;
  float_sw4 getY0() const;
  float_sw4 getZ0() const;
  float_sw4 getDepth() const;
  bool ignore() const {return mIgnore;}
  bool myPoint(){ return m_myPoint; }

  // Amplitude
  float_sw4 getAmplitude() const;
   //  void setAmplitude(float_sw4 amp);
  
  // Offset in time
  float_sw4 getOffset() const;

  // Frequency
  float_sw4 getFrequency() const;
  timeDep getTfunc() const {return mTimeDependence;}
  void setMaxFrequency(float_sw4 max_freq);

  // Type of source
  bool isMomentSource() const;

  float_sw4 dt_to_resolve( int ppw ) const;
  int ppw_to_resolve( float_sw4 dt ) const;

  const std::string& getName() const { return mName; };
  void limit_frequency( int ppw, float_sw4 minvsoh );
  float_sw4 compute_t0_increase( float_sw4 t0_min ) const;
  void adjust_t0( float_sw4 dt0 );

  void set_grid_point_sources4( EW *a_EW, std::vector<GridPointSource*>& point_sources );

  void exact_testmoments( int kx[3], int ky[3], int kz[3], float_sw4 momexact[3] );
  void getForces( float_sw4& fx, float_sw4& fy, float_sw4& fz ) const;
  void getMoments( float_sw4& mxx, float_sw4& mxy, float_sw4& mxz, float_sw4& myy, float_sw4& myz, float_sw4& mzz ) const;
  void setMoments( float_sw4 mxx, float_sw4 mxy, float_sw4 mxz, float_sw4 myy, float_sw4 myz, float_sw4 mzz );
  void printPointer(){std::cout << "Source pointer = "  << mPar << std::endl;}
  void perturb( float_sw4 h, int comp );
  void set_derivative( int der );
  void set_noderivative( );
  void set_dirderivative( float_sw4 dir[11] );
  Source* copy( std::string a_name );
  void set_parameters( float_sw4 x[11] );
  void setFrequency( float_sw4 freq );
  void get_parameters( float_sw4 x[11] ) const;
   //  void filter_timefunc( Filter* fi, float_sw4 tstart, float_sw4 dt, int nsteps );
  bool get_CorrectForMu(){return mShearModulusFactor;};
  void set_CorrectForMu(bool smf){mShearModulusFactor=smf;};
  void copy_pars_to_device();
 private:
  Source();
  void adjust_zcoord( EW* a_ew );
  void compute_mapped_coordinates( EW *a_ew );
   //  void correct_Z_level( EW *a_ew );
  void compute_metric_at_source( EW* a_EW, float_sw4 q, float_sw4 r, float_sw4 s, int ic,
				 int jc, int kc, int g, float_sw4& zq, float_sw4& zr,
				 float_sw4& zs, float_sw4& zqq, float_sw4& zqr, float_sw4& zqs,
				 float_sw4& zrr, float_sw4& zrs, float_sw4& zss ) const;
   //  int spline_interpolation( );
  void getsourcewgh(float_sw4 ai, float_sw4 wgh[6], float_sw4 dwghda[6], float_sw4 ddwghda[6] ) const;
  void getsourcedwgh(float_sw4 ai, float_sw4 wgh[6], float_sw4 dwghda[6], float_sw4 ddwghda[6] ) const;
  void getsourcewghlow(float_sw4 ai, float_sw4 wgh[6], float_sw4 dwghda[6], float_sw4 ddwghda[6] ) const;
  void getsourcedwghlow(float_sw4 ai, float_sw4 wgh[6], float_sw4 dwghda[6], float_sw4 ddwghda[6] ) const;
  void getmetwgh( float_sw4 alph, float_sw4 wgh[8], float_sw4 dwgh[8], float_sw4 ddwgh[8], float_sw4 dddwgh[8] ) const;
  void getmetdwgh( float_sw4 alph, float_sw4 wgh[8] ) const;
  void getmetwgh7( float_sw4 ai, float_sw4 wgh[7] ) const;
  void getmetdwgh7( float_sw4 ai, float_sw4 wgh[7] ) const;

  void getsourcewghNM2sm6(  float_sw4 ci,  float_sw4 wghk[6] ) const;
  void getsourcedwghNM2sm6(  float_sw4 ci,  float_sw4 dwghk[6] ) const;
  void getsourcewghNM1sm6(  float_sw4 ci,  float_sw4 wghk[6] ) const;
  void getsourcedwghNM1sm6(  float_sw4 ci,  float_sw4 dwghk[6] ) const;
  void getsourcewghNsm6(  float_sw4 ci,  float_sw4 wghk[6] ) const;
  void getsourcedwghNsm6(  float_sw4 ci,  float_sw4 dwghk[6] ) const;
  void getsourcewghP1sm6(  float_sw4 ci,  float_sw4 wghk[6] ) const;
  void getsourcedwghP1sm6(  float_sw4 ci,  float_sw4 dwghk[6] ) const;

  float_sw4 find_min_exponent() const;
  std::string mName;
  std::vector<float_sw4> mForces;
  bool mIsMomentSource;
  float_sw4 mFreq, mT0;

  bool m_myPoint;
  bool m_zRelativeToTopography;
  float_sw4 mX0,mY0,mZ0;
  float_sw4 mQ0,mR0,mS0;
  float_sw4* mPar;
  float_sw4* mdevPar;
  int* mIpar;
  int* mdevIpar;
  int mNpar, mNipar;
  int mNcyc;
  int m_derivative;  
  timeDep mTimeDependence;
  float_sw4 m_dir[11];
  bool m_is_filtered;

  float_sw4 m_zTopo;
  bool mIgnore;
  bool mShearModulusFactor;

};

#endif
