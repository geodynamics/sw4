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
#include <cmath>
#include <iostream>
#include "Require.h"
#include "sw4.h"
float_sw4 VerySmoothBump(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = - 1024*pow(t*freq,10) + 5120*pow(t*freq,9) - 10240*pow(t*freq,8) + 10240*pow(t*freq,7) - 5120*pow(t*freq,6) + 1024*pow(t*freq,5);
  return tmp;
}

float_sw4 VerySmoothBump_t(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = freq*( - 1024*10*pow(t*freq,9) + 5120*9*pow(t*freq,8) - 10240*8*pow(t*freq,7) + 10240*7*pow(t*freq,6) - 5120*6*pow(t*freq,5) + 1024*5*pow(t*freq,4));
  return tmp;
}

float_sw4 VerySmoothBump_om(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = t*( - 1024*10*pow(t*freq,9) + 5120*9*pow(t*freq,8) - 10240*8*pow(t*freq,7) + 10240*7*pow(t*freq,6) - 5120*6*pow(t*freq,5) + 1024*5*pow(t*freq,4));
  return tmp;
}

float_sw4 VerySmoothBump_tt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = freq*freq*( - 1024*90*pow(t*freq,8) + 5120*72*pow(t*freq,7) - 10240*56*pow(t*freq,6) + 10240*42*pow(t*freq,5) - 5120*30*pow(t*freq,4) + 1024*20*pow(t*freq,3) );
  return tmp;
}

float_sw4 VerySmoothBump_tom(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = 5120*pow(t*freq,4)*(-20*pow(t*freq,5) + 81*pow(t*freq,4) -
			       128*pow(t*freq,3) + 98*pow(t*freq,2) - 36*t*freq + 5 );
  return tmp;
}

float_sw4 VerySmoothBump_omom(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = t*t*( - 1024*90*pow(t*freq,8) + 5120*72*pow(t*freq,7) - 10240*56*pow(t*freq,6) + 10240*42*pow(t*freq,5) - 5120*30*pow(t*freq,4) + 1024*20*pow(t*freq,3) );
  return tmp;
}

float_sw4 VerySmoothBump_ttt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = freq*freq*freq*( - 1024*90*8*pow(t*freq,7) + 5120*72*7*pow(t*freq,6) - 10240*56*6*pow(t*freq,5) + 10240*42*5*pow(t*freq,4) - 5120*30*4*pow(t*freq,3) + 1024*20*3*pow(t*freq,2) );
  return tmp;
}

float_sw4 VerySmoothBump_omtt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = freq*freq*t*( - 1024*90*8*pow(t*freq,7) + 5120*72*7*pow(t*freq,6) - 10240*56*6*pow(t*freq,5) + 10240*42*5*pow(t*freq,4) - 5120*30*4*pow(t*freq,3) + 1024*20*3*pow(t*freq,2) )
+2*freq*( - 1024*90*pow(t*freq,8) + 5120*72*pow(t*freq,7) - 10240*56*pow(t*freq,6) + 10240*42*pow(t*freq,5) - 5120*30*pow(t*freq,4) + 1024*20*pow(t*freq,3) );
  return tmp;
}

float_sw4 VerySmoothBump_tttt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = pow(freq,4)*122880*
	(-42*pow(t*freq,6)+126*pow(t*freq,5)-140*pow(t*freq,4)+70*pow(t*freq,3)-15*pow(t*freq,2)+t*freq);
  return tmp;
}

float_sw4 VerySmoothBump_tttom(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = freq*freq*61440*
	(-120*pow(t*freq,7)+378*pow(t*freq,6)-448*pow(t*freq,5)+245*pow(t*freq,4)-60*pow(t*freq,3)+5*t*freq*t*freq);
  return tmp;
}

float_sw4 VerySmoothBump_ttomom(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = freq*freq*t*t*122880*
	(-42*pow(t*freq,6)+126*pow(t*freq,5)-140*pow(t*freq,4)+70*pow(t*freq,3)-15*pow(t*freq,2)+t*freq) +
20480*freq*freq*freq*t*t*t*(-153*pow(t*freq,5)+540*pow(t*freq,4)-728*pow(t*freq,3)+462*pow(t*freq,2)-135*t*freq+14);
  return tmp;
}


float_sw4 RickerWavelet(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
    return (2*factor - 1)*exp(-factor);
  else
    return 0;
}

float_sw4 RickerWavelet_t(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
     return pow(M_PI*freq,2)*t*( 6 - 4*factor )*exp(-factor);
  else
    return 0;
}

float_sw4 RickerWavelet_om(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
     return M_PI*M_PI*freq*t*t*( 6 - 4*factor )*exp(-factor);
  else
    return 0;
}

float_sw4 RickerWavelet_tt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
     return M_PI*M_PI*freq*freq*( 6-24*factor+8*factor*factor)*exp(-factor);
  else
    return 0;
}

float_sw4 RickerWavelet_ttt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
     return pow(M_PI*freq,4)*t*( -60+80*factor-16*factor*factor)*exp(-factor);
  else
    return 0;
}

float_sw4 RickerWavelet_omtt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
     return M_PI*M_PI*freq*(12-108*factor+96*factor*factor-16*factor*factor*factor)*exp(-factor);
 
  else
    return 0;
}

float_sw4 RickerInt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
    return -t*exp(-factor);
  else
    return 0;
}

float_sw4 RickerInt_t(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
     return (2*factor-1)*exp(-factor);
  else
    return 0;
}

float_sw4 RickerInt_om(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
     return 2*t*t*t*freq*M_PI*M_PI*exp(-factor);
  else
    return 0;
}

float_sw4 RickerInt_tt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
     return M_PI*M_PI*freq*freq*t*(6-4*factor)*exp(-factor);
  else
    return 0;
}

float_sw4 RickerInt_ttt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
     return M_PI*M_PI*freq*freq*(6-24*factor+8*factor*factor)*exp(-factor);
  else
    return 0;
}

float_sw4 RickerInt_omtt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
     return t*M_PI*M_PI*freq*(12-28*factor+8*factor*factor)*exp(-factor);
  else
    return 0;
}

float_sw4 Gaussian(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
    return freq / sqrt(2*M_PI)*exp(-factor);
  else
    return 0;
}

float_sw4 Gaussian_t(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
    return -freq*freq*freq*t / sqrt(2*M_PI)*exp(-factor);
  else
    return 0;
}

float_sw4 Gaussian_om(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
     return (1-2*factor)/ sqrt(2*M_PI)*exp(-factor);
  else
    return 0;
}

float_sw4 Gaussian_tt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
     return freq / sqrt(2*M_PI)* freq*freq*(2*factor-1)*exp(-factor);
  else
    return 0;
}

float_sw4 Gaussian_tom(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
     return freq*freq*t / sqrt(2*M_PI)*(-3 + 2*factor)*exp(-factor);
  else
    return 0;
}

float_sw4 Gaussian_omom(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
     return freq*t*t / sqrt(2*M_PI)*(-3 + 2*factor)*exp(-factor);
  else
    return 0;
}

float_sw4 Gaussian_ttt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
     return freq*freq*freq*freq*freq*t / sqrt(2*M_PI)*(3-2*factor)*exp(-factor);
  else
    return 0;
}

float_sw4 Gaussian_omtt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
     return freq*freq*(12*factor-3-4*factor*factor)/sqrt(2*M_PI)*exp(-factor);
  else
    return 0;
}

float_sw4 Gaussian_tttt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
     return freq*freq*freq*freq*freq / sqrt(2*M_PI)*(3-12*factor + 4*factor*factor)*exp(-factor);
  else
    return 0;
}

float_sw4 Gaussian_tttom(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
     return freq*freq*freq*freq*t / sqrt(2*M_PI)*(15-20*factor + 4*factor*factor)*exp(-factor);
  else
    return 0;
}

float_sw4 Gaussian_ttomom(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
     return freq / sqrt(2*M_PI)*(-6+54*factor-48*factor*factor+8*factor*factor*factor)*exp(-factor);
  else
    return 0;
}

float_sw4 Erf( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  return 0.5*(1+erf( freq*t/sqrt(2.0)) );
}

float_sw4 Erf_t(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
    return freq / sqrt(2*M_PI)*exp(-factor);
  else
    return 0;
}

float_sw4 Erf_om(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
    return t / sqrt(2*M_PI)*exp(-factor);
  else
    return 0;
}

float_sw4 Erf_tt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
    return -freq / sqrt(2*M_PI)* freq*freq*t*exp(-factor);
  else
    return 0;
}

float_sw4 Erf_ttt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
     return freq / sqrt(2*M_PI)* freq*freq*(2*factor-1)*exp(-factor);
  else
    return 0;
}

float_sw4 Erf_omtt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
     return t / sqrt(2*M_PI)* freq*freq*(2*factor-3)*exp(-factor);
  else
    return 0;
}

float_sw4 Ramp(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 1.0;
  else
    tmp = 0.5*(1 - cos(M_PI*t*freq));
  
  return tmp;
}

float_sw4 Ramp_t(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = 0.5*M_PI*freq*sin(M_PI*t*freq);
  
  return tmp;
}

float_sw4 Ramp_om(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = 0.5*M_PI*t*sin(M_PI*t*freq);
  
  return tmp;
}

float_sw4 Ramp_tt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = 0.5*M_PI*M_PI*freq*freq*cos(M_PI*t*freq);
  
  return tmp;
}

float_sw4 Ramp_ttt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = -0.5*M_PI*M_PI*M_PI*freq*freq*freq*sin(M_PI*t*freq);
  
  return tmp;
}

float_sw4 Ramp_omtt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = M_PI*M_PI*freq*(cos(M_PI*t*freq)-0.5*M_PI*t*freq*sin(M_PI*t*freq));
  
  return tmp;
}

float_sw4 Triangle(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = 2*freq*8./pow(M_PI,2)*(sin(M_PI*(t*freq)) - sin(3*M_PI*(t*freq))/9 + sin(5*M_PI*(t*freq))/25 - sin(7*M_PI*(t*freq))/49);

  return tmp; 
}

float_sw4 Triangle_t(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = 2*freq*8./pow(M_PI,2)*M_PI*freq*(cos(M_PI*(t*freq)) - cos(3*M_PI*(t*freq))/3 + cos(5*M_PI*(t*freq))/5 - cos(7*M_PI*(t*freq))/7);

  return tmp; 
}

float_sw4 Triangle_om(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = 2*8./pow(M_PI,2)*(M_PI*freq*t*(cos(M_PI*(t*freq)) - cos(3*M_PI*(t*freq))/3 + cos(5*M_PI*(t*freq))/5 - cos(7*M_PI*(t*freq))/7) + (sin(M_PI*(t*freq)) - sin(3*M_PI*(t*freq))/9 + sin(5*M_PI*(t*freq))/25 - sin(7*M_PI*(t*freq))/49) );

  return tmp; 
}

float_sw4 Triangle_tt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = 2*freq*8./pow(M_PI,2)*(-M_PI*M_PI*freq*freq)*
	( sin(M_PI*(t*freq)) - sin(3*M_PI*(t*freq)) + 
          sin(5*M_PI*(t*freq)) - sin(7*M_PI*(t*freq)) );

  return tmp; 
}

float_sw4 Triangle_ttt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = 2*freq*8./pow(M_PI,2)*(-M_PI*M_PI*M_PI*freq*freq*freq)*
	( cos(M_PI*(t*freq)) - 3*cos(3*M_PI*(t*freq)) + 
          5*cos(5*M_PI*(t*freq)) - 7*cos(7*M_PI*(t*freq)) );

  return tmp; 
}

float_sw4 Triangle_omtt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = 2*freq*freq*M_PI*M_PI*8/pow(M_PI,2)*( 
	(-3)*( sin(M_PI*(t*freq)) - sin(3*M_PI*(t*freq)) + 
	       sin(5*M_PI*(t*freq)) - sin(7*M_PI*(t*freq)) ) 
	-freq*t*M_PI*(cos(M_PI*(t*freq)) - 3*cos(3*M_PI*(t*freq)) + 
		 5*cos(5*M_PI*(t*freq)) - 7*cos(7*M_PI*(t*freq)) ));
  return tmp; 
}

float_sw4 Sawtooth(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = 8./pow(M_PI,2)*(sin(M_PI*(2*t*freq)) - sin(3*M_PI*(2*t*freq))/9 + sin(5*M_PI*(2*t*freq))/25 - sin(7*M_PI*(2*t*freq))/49);

  return tmp; 
}

float_sw4 Sawtooth_t(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = 8./pow(M_PI,2)*(2*M_PI*freq)*(cos(M_PI*(2*t*freq)) - cos(3*M_PI*(2*t*freq))/3 + cos(5*M_PI*(2*t*freq))/5 - cos(7*M_PI*(2*t*freq))/7);

  return tmp; 
}

float_sw4 Sawtooth_om(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = 8./pow(M_PI,2)*(2*M_PI*t)*(cos(M_PI*(2*t*freq)) - cos(3*M_PI*(2*t*freq))/3 + cos(5*M_PI*(2*t*freq))/5 - cos(7*M_PI*(2*t*freq))/7);

  return tmp; 
}

float_sw4 Sawtooth_tt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = 8./pow(M_PI,2)*(-M_PI*M_PI*2*2*freq*freq)*
              (sin(M_PI*(2*t*freq)) - sin(3*M_PI*(2*t*freq)) +
	       sin(5*M_PI*(2*t*freq)) - sin(7*M_PI*(2*t*freq)));

  return tmp; 
}

float_sw4 Sawtooth_ttt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = -64.*(M_PI*freq*freq*freq)*
              (cos(M_PI*(2*t*freq)) - 3*cos(3*M_PI*(2*t*freq)) +
	       5*cos(5*M_PI*(2*t*freq)) - 7*cos(7*M_PI*(2*t*freq)));

  return tmp; 
}

float_sw4 Sawtooth_omtt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = -64*freq*(sin(M_PI*(2*t*freq)) - sin(3*M_PI*(2*t*freq)) +
		     sin(5*M_PI*(2*t*freq)) - sin(7*M_PI*(2*t*freq))) 
          -64*M_PI*freq*freq*t*
               (cos(M_PI*(2*t*freq)) - 3*cos(3*M_PI*(2*t*freq)) +
		5*cos(5*M_PI*(2*t*freq)) - 7*cos(7*M_PI*(2*t*freq)));

  return tmp; 
}

float_sw4 SmoothWave(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 c0 = 2187./8., c1 = -10935./8., c2 = 19683./8., c3 = -15309./8., c4 = 2187./4.;
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = (c0*pow(t*freq,3)+c1*pow(t*freq,4)+c2*pow(t*freq,5)+c3*pow(t*freq,6)+c4*pow(t*freq,7));
  
  return tmp;
}

float_sw4 SmoothWave_t(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 c0 = 2187./8., c1 = -10935./8., c2 = 19683./8., c3 = -15309./8., c4 = 2187./4.;
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = freq*(c0*3*pow(t*freq,2)+c1*4*pow(t*freq,3)+c2*5*pow(t*freq,4)+c3*6*pow(t*freq,5)+c4*7*pow(t*freq,6));
  return tmp;
}

float_sw4 SmoothWave_om(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 c0 = 2187./8., c1 = -10935./8., c2 = 19683./8., c3 = -15309./8., c4 = 2187./4.;
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = t*(c0*3*pow(t*freq,2)+c1*4*pow(t*freq,3)+c2*5*pow(t*freq,4)+c3*6*pow(t*freq,5)+c4*7*pow(t*freq,6));
  return tmp;
}

float_sw4 SmoothWave_tt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 c0 = 2187./8., c1 = -10935./8., c2 = 19683./8., c3 = -15309./8., c4 = 2187./4.;
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = freq*freq*(c0*6*t*freq+c1*12*pow(t*freq,2)+c2*20*pow(t*freq,3)+c3*30*pow(t*freq,4)+c4*42*pow(t*freq,5));
  
  return tmp;
}

float_sw4 SmoothWave_ttt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 c0 = 2187./8., c1 = -10935./8., c2 = 19683./8., c3 = -15309./8., c4 = 2187./4.;
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = freq*freq*freq*(c0*6+c1*24*t*freq+c2*60*pow(t*freq,2)+c3*120*pow(t*freq,3)+c4*210*pow(t*freq,4));
  return tmp;
}

float_sw4 SmoothWave_omtt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 c0 = 2187./8., c1 = -10935./8., c2 = 19683./8., c3 = -15309./8., c4 = 2187./4.;
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = 2*freq*(c0*6*t*freq+c1*12*pow(t*freq,2)+c2*20*pow(t*freq,3)+c3*30*pow(t*freq,4)+c4*42*pow(t*freq,5)) +
         freq*freq*t*(c0*6+c1*24*t*freq+c2*60*pow(t*freq,2)+c3*120*pow(t*freq,3)+c4*210*pow(t*freq,4));
  return tmp;
}

float_sw4 Brune( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  const float_sw4 tf = t*freq;
  if( tf < 0 )
    return 0;
  else
    {
      if( -tf > par[0] )
	return 1-exp(-tf)*(1+tf);
      else
	return 1;
    }
}

float_sw4 Brune_t( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  const float_sw4 tf = t*freq;
  if( tf < 0 )
    return 0;
  else
    {
      if( -tf > par[0] )
	 return tf*freq*exp(-tf);
      else
	return 0;
    }
}

float_sw4 Brune_om( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  const float_sw4 tf = t*freq;
  if( tf < 0 )
    return 0;
  else
    {
      if( -tf > par[0] )
	 return tf*t*exp(-tf);
      else
	return 0;
    }
}

float_sw4 Brune_tt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  const float_sw4 tf = t*freq;
  if( tf < 0 )
    return 0;
  else
    {
      if( -tf > par[0] )
	 return freq*freq*(1-tf)*exp(-tf);
      else
	return 0;
    }
}

float_sw4 Brune_ttt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  const float_sw4 tf = t*freq;
  if( tf < 0 )
    return 0;
  else
    {
      if( tf < -par[0] )
	 return (tf-2)*freq*freq*freq*exp(-tf);
      else
	return 0;
    }
}

float_sw4 Brune_omtt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  const float_sw4 tf = t*freq;
  if( tf < 0 )
    return 0;
  else
    {
      if( tf < -par[0] )
	 return freq*(2-4*tf+tf*tf)*exp(-tf);
      else
	return 0;
    }
}

float_sw4 DBrune( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  const float_sw4 tf = t*freq;
  if( tf < 0 )
    return 0;
  else
    {
      if( tf < -par[0] )
	return tf*freq*exp(-tf);
      else
	return 0;
    }
}

float_sw4 DBrune_t( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  const float_sw4 tf = t*freq;
  if( tf < 0 )
    return 0;
  else
    {
      if( -tf > par[0] )
	 return freq*freq*(1-tf)*exp(-tf);
      else
	return 0;
    }
}

float_sw4 DBrune_om( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  const float_sw4 tf = t*freq;
  if( tf < 0 )
    return 0;
  else
    {
      if( -tf > par[0] )
	 return tf*(2-tf)*exp(-tf);
      else
	return 0;
    }
}

float_sw4 DBrune_tt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  const float_sw4 tf = t*freq;
  if( tf < 0 )
    return 0;
  else
    {
      if( tf < -par[0] )
	 return (tf-2)*freq*freq*freq*exp(-tf);
      else
	return 0;
    }
}

float_sw4 DBrune_ttt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  const float_sw4 tf = t*freq;
  if( tf < 0 )
    return 0;
  else
    {
      if( tf < -par[0] )
	 return (3-tf)*freq*freq*freq*freq*exp(-tf);
      else
	return 0;
    }
}

float_sw4 DBrune_omtt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  const float_sw4 tf = t*freq;
  if( tf < 0 )
    return 0;
  else
    {
      if( tf < -par[0] )
	 return (6*tf-6-tf*tf)*freq*freq*exp(-tf);
      else
	return 0;
    }
}

float_sw4 BruneSmoothed( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  const float_sw4 tf = t*freq;
  const float_sw4 h  = 2.31;
  const float_sw4 hi = 1/h;
  if( tf < 0 )
    return 0;
  else if( tf < h )
    return 1-exp(-tf)*(1 + tf + 0.5*tf*tf - 1.5*hi*tf*tf*tf + 
		       1.5*hi*hi*tf*tf*tf*tf -0.5*hi*hi*hi*tf*tf*tf*tf*tf);
  else
    {
      if( -tf > par[0] )
	return 1-exp(-tf)*(1+tf);
      else
	return 1;
    }
}

float_sw4 BruneSmoothed_t( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  const float_sw4 tf = t*freq;
  const float_sw4 h  = 2.31;
  const float_sw4 hi = 1/h;
  if( tf < 0 )
    return 0;
  else if( tf < h )
  {
     const float_sw4 c3 = - 1.5*hi;
     const float_sw4 c4 = 1.5*hi*hi;
     const float_sw4 c5 = -0.5*hi*hi*hi;
     return exp(-tf)*freq*((0.5-3*c3)*tf*tf+(c3-4*c4)*tf*tf*tf+(c4-5*c5)*tf*tf*tf*tf+c5*tf*tf*tf*tf*tf);
  }
  else
    {
      if( -tf > par[0] )
	 return tf*freq*exp(-tf);
      else
	return 0;
    }
}

float_sw4 BruneSmoothed_om( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  const float_sw4 tf = t*freq;
  const float_sw4 h  = 2.31;
  const float_sw4 hi = 1/h;
  if( tf < 0 )
    return 0;
  else if( tf < h )
  {
     const float_sw4 c3 = - 1.5*hi;
     const float_sw4 c4 = 1.5*hi*hi;
     const float_sw4 c5 = -0.5*hi*hi*hi;
     return exp(-tf)*t*((0.5-3*c3)*tf*tf+(c3-4*c4)*tf*tf*tf+(c4-5*c5)*tf*tf*tf*tf+c5*tf*tf*tf*tf*tf);
  }
  else
    {
      if( -tf > par[0] )
	 return tf*t*exp(-tf);
      else
	return 0;
    }
}

float_sw4 BruneSmoothed_tt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  const float_sw4 tf = t*freq;
  const float_sw4 h  = 2.31;
  const float_sw4 hi = 1/h;
  if( tf < 0 )
    return 0;
  else if( tf < h )
  {
     const float_sw4 c3 = - 1.5*hi;
     const float_sw4 c4 = 1.5*hi*hi;
     const float_sw4 c5 = -0.5*hi*hi*hi;
     return exp(-tf)*( freq*freq*( (1-6*c3)*tf+(-0.5+6*c3-12*c4)*tf*tf+(-c3+8*c4-20*c5)*tf*tf*tf+
				   (-c4+10*c5)*tf*tf*tf*tf -c5*tf*tf*tf*tf*tf));
				
  }
  else
    {
      if( -tf > par[0] )
	return freq*freq*(1-tf)*exp(-tf);
      else
	return 0;
    }
}

float_sw4 BruneSmoothed_ttt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  const float_sw4 tf = t*freq;
  const float_sw4 h  = 2.31;
  const float_sw4 hi = 1/h;
  if( tf < 0 )
    return 0;
  else if( tf < h )
  {
     const float_sw4 c3 = - 1.5*hi;
     const float_sw4 c4 = 1.5*hi*hi;
     const float_sw4 c5 = -0.5*hi*hi*hi;
     return exp(-tf)*freq*freq*freq*( (1-6*c3) + (18*c3-2-24*c4)*tf+(0.5-9*c3+36*c4-60*c5)*tf*tf+(c3-12*c4+60*c5)*tf*tf*tf+
					(c4-15*c5)*tf*tf*tf*tf +c5*tf*tf*tf*tf*tf);
				
  }
  else
    {
      if( -tf > par[0] )
	 return (tf-2)*freq*freq*freq*exp(-tf);
      else
	return 0;
    }
}

float_sw4 BruneSmoothed_omtt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  const float_sw4 tf = t*freq;
  const float_sw4 h  = 2.31;
  const float_sw4 hi = 1/h;
  if( tf < 0 )
    return 0;
  else if( tf < h )
  {
     const float_sw4 c3 = - 1.5*hi;
     const float_sw4 c4 = 1.5*hi*hi;
     const float_sw4 c5 = -0.5*hi*hi*hi;
     return exp(-tf)*freq*( tf*( (1-6*c3) + (12*c3-1-24*c4)*tf+(-3*c3+24*c4-60*c5)*tf*tf+
				 (-4*c4+40*c5)*tf*tf*tf -5*c5*tf*tf*tf*tf ) + 
			    (2-tf)*((1-6*c3)*tf+(-0.5+6*c3-12*c4)*tf*tf+(-c3+8*c4-20*c5)*tf*tf*tf+
				    (-c4+10*c5)*tf*tf*tf*tf -c5*tf*tf*tf*tf*tf) );
  }
  else
    {
      if( -tf > par[0] )
	 return freq*(2-4*tf+tf*tf)*exp(-tf);
      else
	return 0;
    }
}


float_sw4 GaussianWindow( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 incyc2 = 1/(par[1]*par[1]);
  const float_sw4 tf = t*freq;
  if( -0.5*tf*tf*incyc2  > par[0] )
    return sin(tf)*exp(-0.5*tf*tf*incyc2);
  else
    return 0;
}

float_sw4 GaussianWindow_t( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 incyc2 = 1/(par[1]*par[1]);
  const float_sw4 tf = t*freq;
  if( -0.5*tf*tf*incyc2  > par[0] )
     return (freq*cos(tf)-freq*tf*incyc2*sin(tf))*exp(-0.5*tf*tf*incyc2 );
  else
    return 0;
}

float_sw4 GaussianWindow_om( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 incyc2 = 1/(par[1]*par[1]);
  const float_sw4 tf = t*freq;
  if( -0.5*tf*tf*incyc2  > par[0] )
     return (t*cos(tf)-t*tf*incyc2*sin(tf))*exp(-0.5*tf*tf*incyc2 );
  else
    return 0;
}

float_sw4 GaussianWindow_tt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 incyc2 = 1/(par[1]*par[1]);
  const float_sw4 tf = t*freq;
  if( -0.5*tf*tf*incyc2  > par[0] )
     return ( (-freq*freq-freq*freq*incyc2+freq*freq*tf*tf*incyc2*incyc2)*sin(tf)-
	      tf*2*freq*freq*incyc2*cos(tf) )*exp(-0.5*tf*tf*incyc2);
  else
    return 0;
}

float_sw4 GaussianWindow_ttt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 incyc2 = 1/(par[1]*par[1]);
  const float_sw4 tf = t*freq;
  if( -0.5*tf*tf*incyc2  > par[0] )
     return ( freq*freq*freq*(3*tf*incyc2*(1+incyc2)-tf*tf*tf*incyc2*incyc2*incyc2)*sin(tf) +
	      freq*freq*freq*( 3*tf*tf*incyc2*incyc2-3*incyc2-1)*cos(tf))*exp(-0.5*tf*tf*incyc2);
  else
    return 0;
}

float_sw4 GaussianWindow_omtt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 incyc2 = 1/(par[1]*par[1]);
  const float_sw4 tf = t*freq;
  if( -0.5*tf*tf*incyc2  > par[0] )
     return ( freq*(-2-2*incyc2 + 3*incyc2*tf*tf +5*tf*tf*incyc2*incyc2 -tf*tf*tf*tf*incyc2*incyc2*incyc2)*sin(tf) +
	      t*freq*freq*( 3*tf*tf*incyc2*incyc2-7*incyc2-1)*cos(tf) )*exp(-0.5*tf*tf*incyc2);
  else
    return 0;
}

float_sw4 Liu( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
   float_sw4 tau = 2*M_PI/freq;
   float_sw4 tau1 = 0.13*tau;
   float_sw4 tau2 = tau-tau1;
   if( t < 0 )
      return 0;
   else if( t >= tau )
      return 1;
   else
   {
      float_sw4 ipi = 1.0/M_PI;
      float_sw4 cn = 1.0/(1.4*tau1+1.2*tau1*ipi + 0.3*tau2);
      if( t <= tau1 )
	 return cn*(0.7*t-0.7*tau1*ipi*sin(M_PI*t/tau1)-1.2*tau1*ipi*(cos(0.5*M_PI*t/tau1)-1));
      else if( t <= 2*tau1 )
	 return cn*(1.0*t-0.3*tau1+1.2*tau1*ipi - 0.7*tau1*ipi*sin(M_PI*t/tau1)+0.3*tau2*ipi*sin(M_PI*(t-tau1)/tau2));
      else if( t <= tau )
	 return cn*(0.3*t+1.1*tau1+1.2*tau1*ipi+0.3*tau2*ipi*sin(M_PI*(t-tau1)/tau2));
   }
   return 0.; // should never get here, but keeps compiler happy
}

float_sw4 Liu_t( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
   float_sw4 tau = 2*M_PI/freq;
   float_sw4 tau1 = 0.13*tau;
   float_sw4 tau2 = tau-tau1;
   if( t < 0 )
      return 0;
   else if( t >= tau )
      return 0;
   else
   {
      float_sw4 ipi = 1.0/M_PI;
      float_sw4 cn = 1.0/(1.4*tau1+1.2*tau1*ipi + 0.3*tau2);
      if( t <= tau1 )
	 return cn*(0.7-0.7*cos(M_PI*t/tau1)+0.6*sin(0.5*M_PI*t/tau1));
      else if( t <= 2*tau1 )
	 return cn*(1-0.7*cos(M_PI*t/tau1)+0.3*cos(M_PI*(t-tau1)/tau2));
      else if( t <= tau )
	 return cn*(0.3+0.3*cos(M_PI*(t-tau1)/tau2));
   }
   return 0.; // should never get here, but keeps compiler happy
}

float_sw4 Liu_om( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
   float_sw4 tau = 2*M_PI/freq;
   float_sw4 tau1 = 0.13*tau;
   float_sw4 tau2 = tau-tau1;
   if( t < 0 )
      return 0;
   else if( t >= tau )
      return 0;
   else
   {
      float_sw4 ipi = 1.0/M_PI;
      float_sw4 cn = t*1.0/(1.4*tau1+1.2*tau1*ipi + 0.3*tau2)/freq;
      if( t <= tau1 )
	 return cn*(0.7-0.7*cos(M_PI*t/tau1)+0.6*sin(0.5*M_PI*t/tau1));
      else if( t <= 2*tau1 )
	 return cn*(1-0.7*cos(M_PI*t/tau1)+0.3*cos(M_PI*(t-tau1)/tau2));
      else if( t <= tau )
	 return cn*(0.3+0.3*cos(M_PI*(t-tau1)/tau2));
   }
   return 0.; // should never get here, but keeps compiler happy
}

float_sw4 Liu_tt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
   float_sw4 tau = 2*M_PI/freq;
   float_sw4 tau1 = 0.13*tau;
   float_sw4 tau2 = tau-tau1;
   if( t < 0 )
      return 0;
   else if( t >= tau )
      return 0;
   else
   {
      float_sw4 ipi = 1.0/M_PI;
      float_sw4 cn = 1.0/(1.4*tau1+1.2*tau1*ipi + 0.3*tau2);
      if( t <= tau1 )
	 return cn*(0.7*M_PI*sin(M_PI*t/tau1)+0.3*M_PI*cos(0.5*M_PI*t/tau1))/tau1;
      else if( t <= 2*tau1 )
	 return cn*(0.7*M_PI*sin(M_PI*t/tau1)/tau1-0.3*M_PI*sin(M_PI*(t-tau1)/tau2)/tau2);
      else if( t <= tau )
	 return cn*(-0.3*M_PI*sin(M_PI*(t-tau1)/tau2))/tau2;
   }
   return 0.; // should never get here, but keeps compiler happy
}

float_sw4 Liu_ttt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
   float_sw4 tau = 2*M_PI/freq;
   float_sw4 tau1 = 0.13*tau;
   float_sw4 tau2 = tau-tau1;
   if( t < 0 )
      return 0;
   else if( t >= tau )
      return 0;
   else
   {
      float_sw4 ipi = 1.0/M_PI;
      float_sw4 cn = 1.0/(1.4*tau1+1.2*tau1*ipi + 0.3*tau2);
      if( t <= tau1 )
	 return cn*(0.7*M_PI*M_PI*cos(M_PI/tau1*t)-0.15*M_PI*M_PI*sin(0.5*M_PI/tau1*t))/(tau1*tau1);
      else if( t <= 2*tau1 )
	 return cn*(0.7*M_PI*M_PI*cos(M_PI*t/tau1)/(tau1*tau1)-0.3*M_PI*M_PI*cos(M_PI*(t-tau1)/tau2)/(tau2*tau2));
      else if( t <= tau )
	 return cn*(-0.3*M_PI*M_PI*cos(M_PI*(t-tau1)/tau2))/(tau2*tau2);
   }
   return 0.; // should never get here, but keeps compiler happy
}

float_sw4 Liu_omtt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
   float_sw4 tau = 2*M_PI/freq;
   float_sw4 tau1 = 0.13*tau;
   float_sw4 tau2 = tau-tau1;
   if( t < 0 )
      return 0;
   else if( t >= tau )
      return 0;
   else
   {
      float_sw4 ipi = 1.0/M_PI;
      float_sw4 cn = 1.0/(1.4*tau1+1.2*tau1*ipi + 0.3*tau2)/freq;
      if( t <= tau1 )
	 return cn*(2*(0.7*M_PI*sin(M_PI/tau1*t)+0.3*M_PI*cos(0.5*M_PI/tau1*t)) + 
        (0.7*M_PI*M_PI*t/tau1*cos(M_PI/tau1*t)-0.15*M_PI*M_PI*t/tau1*sin(0.5*M_PI/tau1*t)) )/(tau1);
      else if( t <= 2*tau1 )
	 return cn*(2*(0.7*M_PI*sin(M_PI*t/tau1)/tau1-0.3*M_PI*sin(M_PI*(t-tau1)/tau2)/tau2) + 
		    t*(0.7*M_PI*M_PI*cos(M_PI*t/tau1)/(tau1*tau1)-0.3*M_PI*M_PI*cos(M_PI*(t-tau1)/tau2)/(tau2*tau2)));
      else if( t <= tau )
	 return -0.3*M_PI*cn*( 2*sin(M_PI*(t-tau1)/tau2)/tau2 + M_PI*t*cos(M_PI*(t-tau1)/tau2)/(tau2*tau2) );
   }
   return 0.; // should never get here, but keeps compiler happy
}

float_sw4 NullFunc( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
// this function should never be called
  CHECK_INPUT(false,"The NullFunc time function was called!");
  return 0.;  
}

float_sw4 Dirac( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
   const float_sw4 c1=2.59765625;
   const float_sw4 c2=10.0625;
   const float_sw4 c3=22.875;
   const float_sw4 c4=26.6666666666667;
   const float_sw4 c5=11.6666666666667;
   const float_sw4 o12=0.0833333333333333;
   const float_sw4 o6=0.166666666666667;
   const float_sw4 a1=0.41015625;
   const float_sw4 a2=2.140625;
   const float_sw4 a3=3.4609375;
   //   delta(0)
   // freq holds 1/dt
   // Stencil from -s to p
   // k0 is center of pulse on grid given by t + k*dt

   float_sw4 kc = -t*freq;
   // stencil point of t in [-2,..,2] interval
   int k0 = (int)floor(kc+0.5);
   //   std::cout << "t="<< t << " kc=" << kc << " k0= " << k0 << std::endl;
   if( k0 < -2 || k0 > 2 )
      return 0;
   else
   {
      float_sw4 alpha =(-t*freq-k0);
      float_sw4 alpha2=alpha*alpha;
      float_sw4 pol=alpha2*alpha2*( c1 - c2*alpha2 + c3*alpha2*alpha2 -
	    c4*alpha2*alpha2*alpha2+c5*alpha2*alpha2*alpha2*alpha2);
      float_sw4 wgh;
      if( k0 == 2 )
         wgh = o12*alpha*(1-alpha2)-a1*alpha2+pol;
      else if( k0 == 1 )
         wgh = o6*alpha*(-4+alpha2)+a2*alpha2-4*pol;
      else if( k0 == 0 )
         wgh = 1-a3*alpha2 + 6*pol;
      else if( k0 == -1 )
         wgh = o6*alpha*(4-alpha2)+a2*alpha2-4*pol;
      else if( k0 == -2 )
         wgh = o12*alpha*(-1+alpha2)-a1*alpha2+pol;
      //      std::cout << "wgh = " << wgh << std::endl;
      return freq*wgh;
   }
}

float_sw4 Dirac_t( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
   const float_sw4 c1=10.390625;
   const float_sw4 c2= 60.375;
   const float_sw4 c3=183.0;
   const float_sw4 c4=266.666666666667;
   const float_sw4 c5=140;
   const float_sw4 o12=0.0833333333333333;
   const float_sw4 o6=0.166666666666667;
   const float_sw4 a1=0.8203125;
   const float_sw4 a2=4.281250;
   const float_sw4 a3=6.921875;
   //   delta'(0)
   // freq holds 1/dt
   // Stencil from -s to p
   // k0 is center of pulse on grid given by t + k*dt
   float_sw4 kc = -t*freq;
   // stencil point of t in [-2,..,2] interval
   int k0 = (int)floor(kc+0.5);
   if( k0 < -2 || k0 > 2 )
      return 0;
   else
   {
      float_sw4 wgh;
      float_sw4 alpha =(-t*freq-k0);
      float_sw4 alpha2=alpha*alpha;
      float_sw4 polp=alpha2*alpha*( c1 - c2*alpha2 + c3*alpha2*alpha2 -
	    c4*alpha2*alpha2*alpha2+c5*alpha2*alpha2*alpha2*alpha2);
      if( k0 == 2 )
         wgh = o12*(1-3*alpha2)-a1*alpha + polp;
      else if( k0 == 1 )
         wgh = o6*(-4+3*alpha2)+a2*alpha-4*polp;
      else if( k0 == 0 )
         wgh = -a3*alpha + 6*polp;
      else if( k0 == -1 )
         wgh = o6*(4-3*alpha2)+a2*alpha-4*polp;
      else if( k0 == -2 )
         wgh = o12*(-1+3*alpha2)-a1*alpha+polp;
      return freq*freq*wgh;
   }
}

float_sw4 Dirac_tt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
   const float_sw4 c1=31.171875;
   const float_sw4 c2=301.875;
   const float_sw4 a1=0.8203125;
   const float_sw4 a2=4.281250;
   const float_sw4 a3=6.921875;
   //   delta(0)
   // freq holds 1/dt
   // Stencil from -s to p
   // k0 is center of pulse on grid given by t + k*dt
   float_sw4 kc = -t*freq;
   // stencil point of t in [-2,..,2] interval
   int k0 = (int)floor(kc+0.5);
   if( k0 < -2 || k0 > 2 )
      return 0;
   else
   {
      float_sw4 wgh;
      float_sw4 alpha =(-t*freq-k0);
      float_sw4 alpha2=alpha*alpha;
      float_sw4 polpp=alpha2*( c1 - c2*alpha2 + 1281.0*alpha2*alpha2 -
	    2400.0*alpha2*alpha2*alpha2+1540*alpha2*alpha2*alpha2*alpha2);
      if( k0 == 2 )
         wgh = -0.5*alpha-a1 + polpp;
      else if( k0 == 1 )
         wgh = alpha + a2 - 4*polpp;
      else if( k0 == 0 )
         wgh = -a3 + 6*polpp;
      else if( k0 == -1 )
         wgh = -alpha+a2-4*polpp;
      else if( k0 == -2 )
         wgh = 0.5*alpha-a1+polpp;
      return freq*freq*freq*wgh;
   }
}
float_sw4 Dirac_ttt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
   const float_sw4 c1=62.34375;
   const float_sw4 c2=1207.5;
   //   delta(0)
   // freq holds 1/dt
   // Stencil from -s to p
   // k0 is center of pulse on grid given by t + k*dt
   float_sw4 kc = -t*freq;
   // stencil point of t in [-2,..,2] interval
   int k0 = (int)floor(kc+0.5);
   if( k0 < -2 || k0 > 2 )
      return 0;
   else
   {
      float_sw4 wgh;
      float_sw4 alpha =(-t*freq-k0);
      float_sw4 alpha2=alpha*alpha;
      float_sw4 polppp=alpha*( c1 - c2*alpha2 + 7686.0*alpha2*alpha2 -
	    19200.0*alpha2*alpha2*alpha2+15400*alpha2*alpha2*alpha2*alpha2);
      if( k0 == 2 )
         wgh = -0.5 + polppp;
      else if( k0 == 1 )
         wgh = 1.0 - 4*polppp;
      else if( k0 == 0 )
         wgh =  6*polppp;
      else if( k0 == -1 )
         wgh = -1.0-4*polppp;
      else if( k0 == -2 )
         wgh = 0.5+polppp;
      return freq*freq*freq*freq*wgh;
   }
}
float_sw4 Dirac_tttt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
   const float_sw4 c1=62.34375;
   const float_sw4 c2=3622.5;
   //   delta(0)
   // freq holds 1/dt
   // Stencil from -s to p
   // k0 is center of pulse on grid given by t + k*dt
   float_sw4 kc = -t*freq;
   // stencil point of t in [-2,..,2] interval
   int k0 = (int)floor(kc+0.5);
   if( k0 < -2 || k0 > 2 )
      return 0;
   else
   {
      float_sw4 wgh;
      float_sw4 alpha =(-t*freq-k0);
      float_sw4 alpha2=alpha*alpha;
      float_sw4 polpppp=( c1 - c2*alpha2 + 38430.0*alpha2*alpha2 -
	    134400.0*alpha2*alpha2*alpha2+138600.0*alpha2*alpha2*alpha2*alpha2);
      if( k0 == 2 )
         wgh = polpppp;
      else if( k0 == 1 )
         wgh = -4*polpppp;
      else if( k0 == 0 )
         wgh =  6*polpppp;
      else if( k0 == -1 )
         wgh = -4*polpppp;
      else if( k0 == -2 )
         wgh = polpppp;
      return freq*freq*freq*freq*freq*wgh;
   }
}
float_sw4 Dirac_om( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
   // This source have no omega dependence
   return 0;
}
float_sw4 Dirac_tom( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
   // This source have no omega dependence
   return 0;
}
float_sw4 Dirac_omom( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
   // This source have no omega dependence
   return 0;
}
float_sw4 Dirac_omtt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
   // This source have no omega dependence
   return 0;
}
float_sw4 Dirac_tttom( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
   // This source have no omega dependence
   return 0;
}
float_sw4 Dirac_ttomom( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
   // This source have no omega dependence
   return 0;
}

float_sw4 Discrete( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
// freq holds 1/dt
   float_sw4 tstart = par[0];
   int npts = ipar[0];

   int k = static_cast<int>(floor((t-tstart)*freq));

   if( k < 0 )
   {
      k = 0;
      t = tstart;
   }
   if( k > npts-2 )
   {
      k = npts-2;
      t = tstart+(npts-1)/freq;
   }

   float_sw4 arg=(t-tstart)*freq-k; // (t-(tstart+k*dt))/dt
//std::cout <<  "t= " << t << " npts " << npts << " k= " << k << "arg = " << arg <<  std::endl;
   return par[6*k+1] + par[2+6*k]*arg + par[3+6*k]*arg*arg + par[4+6*k]*arg*arg*arg +
       par[5+6*k]*arg*arg*arg*arg + par[6+6*k]*arg*arg*arg*arg*arg; 
}

float_sw4 Discrete_t( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
// freq holds 1/dt
   float_sw4 tstart = par[0];
   int npts = ipar[0];
   int k = static_cast<int>(floor((t-tstart)*freq));
   if( k < 0 )
   {
      k = 0;
      t = tstart;
   }
   if( k > npts-2 )
   {
      k = npts-2;
      t = tstart+(npts-1)/freq;
   }
   float_sw4 arg=(t-tstart)*freq-k; // (t-(tstart+k*dt))/dt
   return (par[2+6*k] + 2*par[3+6*k]*arg + 3*par[4+6*k]*arg*arg + 4*par[5+6*k]*arg*arg*arg+
	5*par[6+6*k]*arg*arg*arg*arg)*freq;
}

float_sw4 Discrete_tt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
   float_sw4 tstart = par[0];
   int npts = ipar[0];
   int k = static_cast<int>(floor((t-tstart)*freq));
   if( k < 0 )
   {
      k = 0;
      t = tstart;
   }
   if( k > npts-2 )
   {
      k = npts-2;
      t = tstart+(npts-1)/freq;
   }
   float_sw4 arg=(t-tstart)*freq-k; // (t-(tstart+k*dt))/dt
   return (2*par[3+6*k] + 6*par[4+6*k]*arg + 12*par[5+6*k]*arg*arg + 20*par[6+6*k]*arg*arg*arg)*freq*freq;
}

float_sw4 Discrete_ttt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
   float_sw4 tstart = par[0];
   int npts = ipar[0];
   int k = static_cast<int>(floor((t-tstart)*freq));
   if( k < 0 )
   {
      k = 0;
      t = tstart;
   }
   if( k > npts-2 )
   {
      k = npts-2;
      t = tstart+(npts-1)/freq;
   }
   float_sw4 arg=(t-tstart)*freq-k; // (t-(tstart+k*dt))/dt
   return (6*par[4+6*k] + 24*par[5+6*k]*arg + 60*par[6+6*k]*arg*arg)*freq*freq*freq;
}

float_sw4 Discrete_tttt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
   float_sw4 tstart = par[0];
   int npts = ipar[0];
   int k = static_cast<int>(floor((t-tstart)*freq));
   if( k < 0 )
   {
      t = tstart;
      k = 0;
   }
   if( k > npts-2 )
   {
      k = npts-2;
      t = tstart+(npts-1)/freq;
   }
   float_sw4 arg=(t-tstart)*freq-k; // (t-(tstart+k*dt))/dt
   return (24*par[5+6*k] + 120*par[6+6*k]*arg)*freq*freq*freq*freq;
}

float_sw4 Discrete_om( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
   // This source have no omega dependence
   return 0;
}
float_sw4 Discrete_tom( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
   // This source have no omega dependence
   return 0;
}
float_sw4 Discrete_omom( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
   // This source have no omega dependence
   return 0;
}
float_sw4 Discrete_omtt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
   // This source have no omega dependence
   return 0;
}
float_sw4 Discrete_tttom( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
   // This source have no omega dependence
   return 0;
}
float_sw4 Discrete_ttomom( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
   // This source have no omega dependence
   return 0;
}

float_sw4 C6SmoothBump(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = 51480*pow(t*freq*(1-t*freq),7);
  //    tmp = 16384*pow(t*freq*(1-t*freq),7);
  return tmp;
}

float_sw4 C6SmoothBump_t(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     //     tmp = 16384*freq*7*(1-2*t*freq)*pow(t*freq*(1-t*freq),6);
     tmp = 51480*freq*7*(1-2*t*freq)*pow(t*freq*(1-t*freq),6);
  return tmp;
}

float_sw4 C6SmoothBump_om(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     //     tmp = 16384*t*7*(1-2*t*freq)*pow(t*freq*(1-t*freq),6);
     tmp = 51480*t*7*(1-2*t*freq)*pow(t*freq*(1-t*freq),6);
  return tmp;
}

float_sw4 C6SmoothBump_tt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     //     tmp = 16384*freq*freq*7*( 6*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),5) - 2*pow(t*freq*(1-t*freq),6));
     tmp = 51480*freq*freq*7*( 6*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),5) - 2*pow(t*freq*(1-t*freq),6));
  return tmp;
}

float_sw4 C6SmoothBump_tom(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     //     tmp = 16384*(7*(1-2*t*freq)*pow(t*freq*(1-t*freq),6)+
     //	          t*freq*7*( 6*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),5) - 2*pow(t*freq*(1-t*freq),6) ));
     tmp = 51480*(7*(1-2*t*freq)*pow(t*freq*(1-t*freq),6)+
	          t*freq*7*( 6*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),5) - 2*pow(t*freq*(1-t*freq),6) ));
  return tmp;
}

float_sw4 C6SmoothBump_omom(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     //     tmp = 16384*t*t*7*( 6*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),5) - 2*pow(t*freq*(1-t*freq),6));
     tmp = 51480*t*t*7*( 6*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),5) - 2*pow(t*freq*(1-t*freq),6));
  return tmp;
}

float_sw4 C6SmoothBump_ttt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     //     tmp = 16384*freq*freq*freq*42*(1-2*t*freq)*( -6*pow(t*freq*(1-t*freq),5)+
     //					    5*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),4));
     tmp = 51480*freq*freq*freq*42*(1-2*t*freq)*( -6*pow(t*freq*(1-t*freq),5)+
					    5*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),4));
  return tmp;
}

float_sw4 C6SmoothBump_omtt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     //     tmp = 16384*(2*freq*(7*( 6*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),5) - 2*pow(t*freq*(1-t*freq),6)) ) 
     //            + freq*freq*t*(42*(1-2*t*freq)*( -6*pow(t*freq*(1-t*freq),5)+
     //					     5*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),4)) ) );
     tmp = 51480*(2*freq*(7*( 6*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),5) - 2*pow(t*freq*(1-t*freq),6)) ) 
            + freq*freq*t*(42*(1-2*t*freq)*( -6*pow(t*freq*(1-t*freq),5)+
					     5*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),4)) ) );
  return tmp;
}

float_sw4 C6SmoothBump_tttt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     //     tmp = 16384*freq*freq*freq*freq*(12*42*(pow(t*freq*(1-t*freq),5)-5*(1-2*freq*t)*(1-2*freq*t)*pow(t*freq*(1-t*freq),4))+840*pow((1-2*t*freq)*t*freq*(1-t*freq),4) );
     tmp = 51480*freq*freq*freq*freq*(12*42*(pow(t*freq*(1-t*freq),5)-5*(1-2*freq*t)*(1-2*freq*t)*pow(t*freq*(1-t*freq),4))+840*pow((1-2*t*freq)*t*freq*(1-t*freq),4) );
  return tmp;
}

float_sw4 C6SmoothBump_tttom(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     //     tmp = 16384*freq*freq*(3*(42*(1-2*t*freq)*( -6*pow(t*freq*(1-t*freq),5)+
     //		 5*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),4))) +
     //	    freq*t*(12*42*(pow(t*freq*(1-t*freq),5)-5*(1-2*freq*t)*(1-2*freq*t)*pow(t*freq*(1-t*freq),4))+
     //                    840*pow((1-2*t*freq)*t*freq*(1-t*freq),4) ) );
     tmp = 51480*freq*freq*(3*(42*(1-2*t*freq)*( -6*pow(t*freq*(1-t*freq),5)+
		 5*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),4))) +
	    freq*t*(12*42*(pow(t*freq*(1-t*freq),5)-5*(1-2*freq*t)*(1-2*freq*t)*pow(t*freq*(1-t*freq),4))+
                    840*pow((1-2*t*freq)*t*freq*(1-t*freq),4) ) );
  return tmp;
}

float_sw4 C6SmoothBump_ttomom(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar )
{
  float_sw4 tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     //     tmp = 16384*( 2*(7*( 6*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),5) - 2*pow(t*freq*(1-t*freq),6)))+
     //		   4*freq*t*( 42*(1-2*t*freq)*( -6*pow(t*freq*(1-t*freq),5)+
     //				5*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),4))) +
     //	   freq*freq*t*t*( (12*42*(pow(t*freq*(1-t*freq),5)-
     //	   5*(1-2*freq*t)*(1-2*freq*t)*pow(t*freq*(1-t*freq),4))+840*pow((1-2*t*freq)*t*freq*(1-t*freq),4) )));
     tmp = 51480*( 2*(7*( 6*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),5) - 2*pow(t*freq*(1-t*freq),6)))+
		   4*freq*t*( 42*(1-2*t*freq)*( -6*pow(t*freq*(1-t*freq),5)+
				5*(1-2*t*freq)*(1-2*t*freq)*pow(t*freq*(1-t*freq),4))) +
	   freq*freq*t*t*( (12*42*(pow(t*freq*(1-t*freq),5)-
	   5*(1-2*freq*t)*(1-2*freq*t)*pow(t*freq*(1-t*freq),4))+840*pow((1-2*t*freq)*t*freq*(1-t*freq),4) )));
  return tmp;
}
