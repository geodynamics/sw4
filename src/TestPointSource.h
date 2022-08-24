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
#ifndef TEST_POINT_SOURCE_H
#define TEST_POINT_SOURCE_H

#include "Sarray.h"

class Source;

class TestPointSource
{
public:
   double m_rho, m_cs, m_cp, m_lambda, m_mu;
   Source* m_source_ptr;

   TestPointSource( double rho, double cs, double cp, Source* source_ptr=NULL ) : 
      m_rho(rho),m_cs(cs),m_cp(cp),m_source_ptr(source_ptr)
   {
      m_mu = m_cs*m_cs*m_rho;
      m_lambda = m_cp*m_cp*m_rho-2*m_mu;
   }
   void ubnd( Sarray& u, Sarray& x, Sarray&  y, Sarray&  z, float_sw4 t, 
              float_sw4 h, int nghost, int sides[6] );
   void set_source( Source* source_ptr ){m_source_ptr=source_ptr;}

private:
   TestPointSource(const TestPointSource&);
   TestPointSource& operator=(const TestPointSource&);

   float_sw4 SmoothWave(float_sw4 t, float_sw4 R, float_sw4 c);
   float_sw4 VerySmoothBump(float_sw4 t, float_sw4 R, float_sw4 c);
   float_sw4 C6SmoothBump(float_sw4 t, float_sw4 R, float_sw4 c);
   float_sw4 d_SmoothWave_dt(float_sw4 t, float_sw4 R, float_sw4 c);
   float_sw4 d_VerySmoothBump_dt(float_sw4 t, float_sw4 R, float_sw4 c);
   float_sw4 d_C6SmoothBump_dt(float_sw4 t, float_sw4 R, float_sw4 c);
   float_sw4 SWTP(float_sw4 Lim, float_sw4 t);
   float_sw4 VSBTP(float_sw4 Lim, float_sw4 t);
   float_sw4 C6SBTP(float_sw4 Lim, float_sw4 t);
   float_sw4 SmoothWave_x_T_Integral(float_sw4 t, float_sw4 R, float_sw4 alpha, float_sw4 beta);
   float_sw4 VerySmoothBump_x_T_Integral(float_sw4 t, float_sw4 R, float_sw4 alpha, float_sw4 beta);
   float_sw4 C6SmoothBump_x_T_Integral(float_sw4 t, float_sw4 R, float_sw4 alpha, float_sw4 beta);
   float_sw4 Gaussian(float_sw4 t, float_sw4 R, float_sw4 c, float_sw4 f );
   float_sw4 d_Gaussian_dt(float_sw4 t, float_sw4 R, float_sw4 c, float_sw4 f);
   float_sw4 Gaussian_x_T_Integral(float_sw4 t, float_sw4 R, float_sw4 f, float_sw4 alpha, float_sw4 beta);
};


#endif
