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
#include "SuperGrid.h"
#include "Require.h"
#include <cstdio>
using namespace std;

SuperGrid::SuperGrid()
{
  m_left = false;
  m_right = false;
  m_x0=0.;
  m_x1=1.;
  m_width=0.1;
  m_const_width=0.;
  m_epsL = 1e-4;
  m_tw_omega = 1.2;
}

void SuperGrid::print_parameters() const
{
   printf("SuperGrid parameters left=%i, right=%i, x0=%e, x1=%e, width=%e, transition=%e epsL=%e\n", m_left, m_right, m_x0, m_x1, m_width, m_trans_width,m_epsL);
}

void SuperGrid::define_taper(bool left, double leftStart, bool right, double rightEnd, double width)
{
  m_left = left;
  m_x0 = leftStart;
  m_right = right;
  m_x1 = rightEnd;
  m_width = width;
// always use the full width for the transition, making m_const_width=0
// experimenting ...
  m_trans_width = 0.5*width;
//  m_trans_width = 1.0*width;
//  m_trans_width = transWidth;
  m_const_width = m_width - m_trans_width;
  
// sanity checks
  if (m_left || m_right)
  {
     double dlen = m_x1-m_x0;
     CHECK_INPUT(m_width > 0., "The supergrid taper width must be positive, not = " << m_width);
     CHECK_INPUT(m_width < dlen, "The supergrid taper width must be smaller than the domain, not = " << m_width);
     CHECK_INPUT(m_trans_width > 0., "The supergrid taper transition width must be positive, not = " << m_trans_width);
     CHECK_INPUT(m_const_width >= 0., "The supergrid const_width = width - trans_width must be non-negative, not = " << m_const_width);
  }
  
  if (m_left && m_right)
  {
    if (m_x0+m_width > m_x1-m_width)
    {
      print_parameters();
      CHECK_INPUT(false, "The supergrid taper functions at the left and right must be separated. Here x0+width = " << m_x0+m_width << 
		  " and x1-width = " << m_x1-m_width);
    }
    
  }
  else if( m_left )
  {
    if (m_x0+m_width > m_x1 )
    {
      print_parameters();
      CHECK_INPUT(false, "The supergrid taper functions at the left must be smaller than the domain. Here x0+width = " << m_x0+m_width << 
		  " and x1 = " << m_x1);
    }
  }    
  else if( m_right )
  {
    if (m_x0 > m_x1-m_width )
    {
      print_parameters();
      CHECK_INPUT(false, "The supergrid taper functions at the right must be smaller than the domain. Here x0 = " << m_x0 << 
		  " and x1-width = " << m_x1-m_width );
    }
  }
}

double SuperGrid::dampingCoeff(double x) const
{
  double phi = stretching(x);
// should be equivalent to PsiAux/phi
//  double f=(1-phi)/phi/(1-m_epsL);
  double f = PsiDamp(x)/phi;
// replaced PsiAux by PsiDamp, which goes to one faster
  
  return f;
}

double SuperGrid::stretching( double x ) const
{ // this function satisfies 0 < epsL <= f <= 1
  return 1-(1-m_epsL)*PsiAux(x); // PsiAux(x) = psi(x) in our papers
}

double SuperGrid::PsiAux(double x) const
{ // PsiAux = psi in our papers
// this function is zero for m_x0+m_width <= x <= m_x1-m_width
// and one for x=m_x0 and x=m_x1
  double f=0.;
  if (m_left && x < m_x0+m_width)
// the following makes the damping transition in 0 <= x <= m_width
    f=Psi0( (m_x0 + m_width - x)/m_width); 
  else if (m_right && x > m_x1-m_width)
// the following makes the damping transition in m_x1-m_width < x < m_x1 
    f=Psi0( (x - (m_x1-m_width) )/m_width);
  return f;
}

double SuperGrid::PsiDamp(double x) const
{ // PsiAux = psi in our papers
// this function is zero for m_x0+m_width <= x <= m_x1-m_width
// and one for x=m_x0 and x=m_x1
  double f=0.;
  if (m_left && x < m_x0+m_width)
// the following makes the damping transition in 0 < const_width <= x <= const_width+trans_width = m_width
// constant damping in 0 <= x <= const_width
    f=Psi0( (m_x0 + m_width - x)/m_trans_width); 
  else if (m_right && x > m_x1-m_width)
// the following makes the damping transition in m_x1-m_width < x < m_x1 - const_width < m_x1
// constant damping in m_x1 - const_width <= x <= m_x1
    f=Psi0( (x - (m_x1-m_width) )/m_trans_width);
  return f;
}

double SuperGrid::linTaper(double x) const
{ 
// this function is zero for m_x0+m_width <= x <= m_x1-m_width
// and one for x=m_x0 and x=m_x1
  double f=0.;
  if (m_left && x < m_x0+m_width)
//  linear taper from 0 to 1
    f= (m_x0 + m_width - x)/m_width; 
  else if (m_right && x > m_x1-m_width)
// linear taper from 0 to 1
    f= (x - (m_x1-m_width) )/m_width;
  return f;
}


// used for damping coefficient
double SuperGrid::Psi0(double xi) const
{
   double f;
   if (xi<=0.)
      f = 0;
   else if (xi>=1.)
      f = 1.0;
   else
//    f=xi*xi*xi*(10 - 15*xi + 6*xi*xi);
//    f = fmin + (1.-fmin)*xi*xi*xi*(10 - 15*xi + 6*xi*xi);
// C4 function
//    f = fmin + (1.-fmin)* xi*xi*xi*xi*xi*( 
//      126 - 420*xi + 540*xi*xi - 315*xi*xi*xi + 70*xi*xi*xi*xi );
// Skewed C4 fcn (p3)
//    f = xi*xi*xi*xi*xi*(-14.0 + 70.0*xi - 90.0*xi*xi + 35.0*xi*xi*xi);
// C5 function (currently the default stretching function)  (p1)
     f =  xi*xi*xi*xi*xi*xi*(
       462-1980*xi+3465*xi*xi-3080*xi*xi*xi+1386*xi*xi*xi*xi-252*xi*xi*xi*xi*xi);
// one-sided C5 fcn (p2)
//     f =  xi*xi*xi*xi*xi*xi*(84.0 - 216.0*xi + 189.0*xi*xi - 56.0*xi*xi*xi);
   
   return f;
}

double SuperGrid::cornerTaper( double x ) const
{ // this function is 1 in the interior and tapers linearly to 1/2 in the SG layers
  const double cmin=0.33;
  return 1.0 - (1.0-cmin)*linTaper(x);
}

double SuperGrid::tw_stretching( double x ) const
{
   return 1 + 0.5*sin(m_tw_omega*x);
}
  
void SuperGrid::set_twilight( double omega )
{
   m_tw_omega = omega;
}
