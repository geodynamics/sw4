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
#ifndef EW_QSPLINE_H
#define EW_QSPLINE_H

#include <cmath>
#include <iostream>
using namespace std;

class Qspline
{
   int m_npts;
   double* m_polcof;
   double m_tmin, m_dt, m_dti;
public:
   Qspline( int npts, double* fun, double tmin, double dt, int bclow=1, int bchigh=1,
	    double s1=0, double t1=0, double sn=0, double tn=0 );
   void Qsplineold( int npts, double* fun, double tmin, double dt );
   ~Qspline();
   double* get_polycof_ptr() { return m_polcof;}
// AP: evalf and evaldd are not used by GridPointSource::getTimeFunc. 
// Instead the functions Discrete() and Discrete_tt(), defined in time_functions.C, are called
//   void evalf( double t, double& f );
   //   void evald( double t, double& f, double& f1, double& f2 );
//   void evaldd( double t, double& f, double& f1, double& f2, double& f3, double& f4 );
};

#endif
