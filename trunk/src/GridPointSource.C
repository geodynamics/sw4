#include "mpi.h"
#include "GridPointSource.h"
#include "Require.h"

#include <fenv.h>
#include <cmath>

#include  "EW.h"

using namespace std;

double VerySmoothBump(double freq, double t, double* par )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = - 1024*pow(t*freq,10) + 5120*pow(t*freq,9) - 10240*pow(t*freq,8) + 10240*pow(t*freq,7) - 5120*pow(t*freq,6) + 1024*pow(t*freq,5);
  return tmp;
}

double VerySmoothBump_t(double freq, double t, double* par )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = freq*( - 1024*10*pow(t*freq,9) + 5120*9*pow(t*freq,8) - 10240*8*pow(t*freq,7) + 10240*7*pow(t*freq,6) - 5120*6*pow(t*freq,5) + 1024*5*pow(t*freq,4));
  return tmp;
}

double VerySmoothBump_om(double freq, double t, double* par )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = t*( - 1024*10*pow(t*freq,9) + 5120*9*pow(t*freq,8) - 10240*8*pow(t*freq,7) + 10240*7*pow(t*freq,6) - 5120*6*pow(t*freq,5) + 1024*5*pow(t*freq,4));
  return tmp;
}

double VerySmoothBump_tt(double freq, double t, double* par )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = freq*freq*( - 1024*90*pow(t*freq,8) + 5120*72*pow(t*freq,7) - 10240*56*pow(t*freq,6) + 10240*42*pow(t*freq,5) - 5120*30*pow(t*freq,4) + 1024*20*pow(t*freq,3) );
  return tmp;
}

double VerySmoothBump_tom(double freq, double t, double* par )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = 5120*pow(t*freq,4)*(-20*pow(t*freq,5) + 81*pow(t*freq,4) -
			       128*pow(t*freq,3) + 98*pow(t*freq,2) - 36*t*freq + 5 );
  return tmp;
}

double VerySmoothBump_omom(double freq, double t, double* par )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = t*t*( - 1024*90*pow(t*freq,8) + 5120*72*pow(t*freq,7) - 10240*56*pow(t*freq,6) + 10240*42*pow(t*freq,5) - 5120*30*pow(t*freq,4) + 1024*20*pow(t*freq,3) );
  return tmp;
}

double VerySmoothBump_ttt(double freq, double t, double* par )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = freq*freq*freq*( - 1024*90*8*pow(t*freq,7) + 5120*72*7*pow(t*freq,6) - 10240*56*6*pow(t*freq,5) + 10240*42*5*pow(t*freq,4) - 5120*30*4*pow(t*freq,3) + 1024*20*3*pow(t*freq,2) );
  return tmp;
}

double VerySmoothBump_omtt(double freq, double t, double* par )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = freq*freq*t*( - 1024*90*8*pow(t*freq,7) + 5120*72*7*pow(t*freq,6) - 10240*56*6*pow(t*freq,5) + 10240*42*5*pow(t*freq,4) - 5120*30*4*pow(t*freq,3) + 1024*20*3*pow(t*freq,2) )
+2*freq*( - 1024*90*pow(t*freq,8) + 5120*72*pow(t*freq,7) - 10240*56*pow(t*freq,6) + 10240*42*pow(t*freq,5) - 5120*30*pow(t*freq,4) + 1024*20*pow(t*freq,3) );
  return tmp;
}

double VerySmoothBump_tttt(double freq, double t, double* par )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = pow(freq,4)*122880*
	(-42*pow(t*freq,6)+126*pow(t*freq,5)-140*pow(t*freq,4)+70*pow(t*freq,3)-15*pow(t*freq,2)+t*freq);
  return tmp;
}

double VerySmoothBump_tttom(double freq, double t, double* par )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = freq*freq*61440*
	(-120*pow(t*freq,7)+378*pow(t*freq,6)-448*pow(t*freq,5)+245*pow(t*freq,4)-60*pow(t*freq,3)+5*t*freq*t*freq);
  return tmp;
}

double VerySmoothBump_ttomom(double freq, double t, double* par )
{
  double tmp;
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


double RickerWavelet(double freq, double t, double* par )
{
  double factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
    return (2*factor - 1)*exp(-factor);
  else
    return 0;
}

double RickerWavelet_t(double freq, double t, double* par )
{
  double factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
     return pow(M_PI*freq,2)*t*( 6 - 4*factor )*exp(-factor);
  else
    return 0;
}

double RickerWavelet_om(double freq, double t, double* par )
{
  double factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
     return M_PI*M_PI*freq*t*t*( 6 - 4*factor )*exp(-factor);
  else
    return 0;
}

double RickerWavelet_tt(double freq, double t, double* par )
{
  double factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
     return M_PI*M_PI*freq*freq*( 6-24*factor+8*factor*factor)*exp(-factor);
  else
    return 0;
}

double RickerWavelet_ttt(double freq, double t, double* par )
{
  double factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
     return pow(M_PI*freq,4)*t*( -60+80*factor-16*factor*factor)*exp(-factor);
  else
    return 0;
}

double RickerWavelet_omtt(double freq, double t, double* par )
{
  double factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
     return M_PI*M_PI*freq*(12-108*factor+96*factor*factor-16*factor*factor*factor)*exp(-factor);
 
  else
    return 0;
}

double RickerInt(double freq, double t, double* par )
{
  double factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
    return -t*exp(-factor);
  else
    return 0;
}

double RickerInt_t(double freq, double t, double* par )
{
  double factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
     return (2*factor-1)*exp(-factor);
  else
    return 0;
}

double RickerInt_om(double freq, double t, double* par )
{
  double factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
     return 2*t*t*t*freq*M_PI*M_PI*exp(-factor);
  else
    return 0;
}

double RickerInt_tt(double freq, double t, double* par )
{
  double factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
     return M_PI*M_PI*freq*freq*t*(6-4*factor)*exp(-factor);
  else
    return 0;
}

double RickerInt_ttt(double freq, double t, double* par )
{
  double factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
     return M_PI*M_PI*freq*freq*(6-24*factor+8*factor*factor)*exp(-factor);
  else
    return 0;
}

double RickerInt_omtt(double freq, double t, double* par )
{
  double factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
     return t*M_PI*M_PI*freq*(12-28*factor+8*factor*factor)*exp(-factor);
  else
    return 0;
}

double Gaussian(double freq, double t, double* par )
{
  double factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
    return freq / sqrt(2*M_PI)*exp(-factor);
  else
    return 0;
}

double Gaussian_t(double freq, double t, double* par )
{
  double factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
    return -freq*freq*freq*t / sqrt(2*M_PI)*exp(-factor);
  else
    return 0;
}

double Gaussian_om(double freq, double t, double* par )
{
  double factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
     return (1-2*factor)/ sqrt(2*M_PI)*exp(-factor);
  else
    return 0;
}

double Gaussian_tt(double freq, double t, double* par )
{
  double factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
     return freq / sqrt(2*M_PI)* freq*freq*(2*factor-1)*exp(-factor);
  else
    return 0;
}

double Gaussian_tom(double freq, double t, double* par )
{
  double factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
     return freq*freq*t / sqrt(2*M_PI)*(-3 + 2*factor)*exp(-factor);
  else
    return 0;
}

double Gaussian_omom(double freq, double t, double* par )
{
  double factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
     return freq*t*t / sqrt(2*M_PI)*(-3 + 2*factor)*exp(-factor);
  else
    return 0;
}

double Gaussian_ttt(double freq, double t, double* par )
{
  double factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
     return freq*freq*freq*freq*freq*t / sqrt(2*M_PI)*(3-2*factor)*exp(-factor);
  else
    return 0;
}

double Gaussian_omtt(double freq, double t, double* par )
{
  double factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
     return freq*freq*(12*factor-3-4*factor*factor)/sqrt(2*M_PI)*exp(-factor);
  else
    return 0;
}

double Gaussian_tttt(double freq, double t, double* par )
{
  double factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
     return freq*freq*freq*freq*freq / sqrt(2*M_PI)*(3-12*factor + 4*factor*factor)*exp(-factor);
  else
    return 0;
}

double Gaussian_tttom(double freq, double t, double* par )
{
  double factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
     return freq*freq*freq*freq*t / sqrt(2*M_PI)*(15-20*factor + 4*factor*factor)*exp(-factor);
  else
    return 0;
}

double Gaussian_ttomom(double freq, double t, double* par )
{
  double factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
     return freq / sqrt(2*M_PI)*(-6+54*factor-48*factor*factor+8*factor*factor*factor)*exp(-factor);
  else
    return 0;
}

double Erf( double freq, double t, double* par )
{
  return 0.5*(1+erf( freq*t/sqrt(2.0)) );
}

double Erf_t(double freq, double t, double* par )
{
  double factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
    return freq / sqrt(2*M_PI)*exp(-factor);
  else
    return 0;
}

double Erf_om(double freq, double t, double* par )
{
  double factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
    return t / sqrt(2*M_PI)*exp(-factor);
  else
    return 0;
}

double Erf_tt( double freq, double t, double* par )
{
  double factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
    return -freq / sqrt(2*M_PI)* freq*freq*t*exp(-factor);
  else
    return 0;
}

double Erf_ttt(double freq, double t, double* par )
{
  double factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
     return freq / sqrt(2*M_PI)* freq*freq*(2*factor-1)*exp(-factor);
  else
    return 0;
}

double Erf_omtt(double freq, double t, double* par )
{
  double factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
     return t / sqrt(2*M_PI)* freq*freq*(2*factor-3)*exp(-factor);
  else
    return 0;
}

double Ramp(double freq, double t, double* par )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 1.0;
  else
    tmp = 0.5*(1 - cos(M_PI*t*freq));
  
  return tmp;
}

double Ramp_t(double freq, double t, double* par )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = 0.5*M_PI*freq*sin(M_PI*t*freq);
  
  return tmp;
}

double Ramp_om(double freq, double t, double* par )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = 0.5*M_PI*t*sin(M_PI*t*freq);
  
  return tmp;
}

double Ramp_tt(double freq, double t, double* par )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = 0.5*M_PI*M_PI*freq*freq*cos(M_PI*t*freq);
  
  return tmp;
}

double Ramp_ttt(double freq, double t, double* par )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = -0.5*M_PI*M_PI*M_PI*freq*freq*freq*sin(M_PI*t*freq);
  
  return tmp;
}

double Ramp_omtt(double freq, double t, double* par )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = M_PI*M_PI*freq*(cos(M_PI*t*freq)-0.5*M_PI*t*freq*sin(M_PI*t*freq));
  
  return tmp;
}

double Triangle(double freq, double t, double* par )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = 2*freq*8./pow(M_PI,2)*(sin(M_PI*(t*freq)) - sin(3*M_PI*(t*freq))/9 + sin(5*M_PI*(t*freq))/25 - sin(7*M_PI*(t*freq))/49);

  return tmp; 
}

double Triangle_t(double freq, double t, double* par )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = 2*freq*8./pow(M_PI,2)*M_PI*freq*(cos(M_PI*(t*freq)) - cos(3*M_PI*(t*freq))/3 + cos(5*M_PI*(t*freq))/5 - cos(7*M_PI*(t*freq))/7);

  return tmp; 
}

double Triangle_om(double freq, double t, double* par )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = 2*8./pow(M_PI,2)*(M_PI*freq*t*(cos(M_PI*(t*freq)) - cos(3*M_PI*(t*freq))/3 + cos(5*M_PI*(t*freq))/5 - cos(7*M_PI*(t*freq))/7) + (sin(M_PI*(t*freq)) - sin(3*M_PI*(t*freq))/9 + sin(5*M_PI*(t*freq))/25 - sin(7*M_PI*(t*freq))/49) );

  return tmp; 
}

double Triangle_tt(double freq, double t, double* par )
{
  double tmp;
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

double Triangle_ttt(double freq, double t, double* par )
{
  double tmp;
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

double Triangle_omtt(double freq, double t, double* par )
{
  double tmp;
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

double Sawtooth(double freq, double t, double* par )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = 8./pow(M_PI,2)*(sin(M_PI*(2*t*freq)) - sin(3*M_PI*(2*t*freq))/9 + sin(5*M_PI*(2*t*freq))/25 - sin(7*M_PI*(2*t*freq))/49);

  return tmp; 
}

double Sawtooth_t(double freq, double t, double* par )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = 8./pow(M_PI,2)*(2*M_PI*freq)*(cos(M_PI*(2*t*freq)) - cos(3*M_PI*(2*t*freq))/3 + cos(5*M_PI*(2*t*freq))/5 - cos(7*M_PI*(2*t*freq))/7);

  return tmp; 
}

double Sawtooth_om(double freq, double t, double* par )
{
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = 8./pow(M_PI,2)*(2*M_PI*t)*(cos(M_PI*(2*t*freq)) - cos(3*M_PI*(2*t*freq))/3 + cos(5*M_PI*(2*t*freq))/5 - cos(7*M_PI*(2*t*freq))/7);

  return tmp; 
}

double Sawtooth_tt(double freq, double t, double* par )
{
  double tmp;
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

double Sawtooth_ttt(double freq, double t, double* par )
{
  double tmp;
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

double Sawtooth_omtt(double freq, double t, double* par )
{
  double tmp;
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

double SmoothWave(double freq, double t, double* par )
{
  double c0 = 2187./8., c1 = -10935./8., c2 = 19683./8., c3 = -15309./8., c4 = 2187./4.;
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = (c0*pow(t*freq,3)+c1*pow(t*freq,4)+c2*pow(t*freq,5)+c3*pow(t*freq,6)+c4*pow(t*freq,7));
  
  return tmp;
}

double SmoothWave_t(double freq, double t, double* par )
{
  double c0 = 2187./8., c1 = -10935./8., c2 = 19683./8., c3 = -15309./8., c4 = 2187./4.;
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = freq*(c0*3*pow(t*freq,2)+c1*4*pow(t*freq,3)+c2*5*pow(t*freq,4)+c3*6*pow(t*freq,5)+c4*7*pow(t*freq,6));
  return tmp;
}

double SmoothWave_om(double freq, double t, double* par )
{
  double c0 = 2187./8., c1 = -10935./8., c2 = 19683./8., c3 = -15309./8., c4 = 2187./4.;
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
    tmp = t*(c0*3*pow(t*freq,2)+c1*4*pow(t*freq,3)+c2*5*pow(t*freq,4)+c3*6*pow(t*freq,5)+c4*7*pow(t*freq,6));
  return tmp;
}

double SmoothWave_tt(double freq, double t, double* par )
{
  double c0 = 2187./8., c1 = -10935./8., c2 = 19683./8., c3 = -15309./8., c4 = 2187./4.;
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = freq*freq*(c0*6*t*freq+c1*12*pow(t*freq,2)+c2*20*pow(t*freq,3)+c3*30*pow(t*freq,4)+c4*42*pow(t*freq,5));
  
  return tmp;
}

double SmoothWave_ttt(double freq, double t, double* par )
{
  double c0 = 2187./8., c1 = -10935./8., c2 = 19683./8., c3 = -15309./8., c4 = 2187./4.;
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = freq*freq*freq*(c0*6+c1*24*t*freq+c2*60*pow(t*freq,2)+c3*120*pow(t*freq,3)+c4*210*pow(t*freq,4));
  return tmp;
}

double SmoothWave_omtt(double freq, double t, double* par )
{
  double c0 = 2187./8., c1 = -10935./8., c2 = 19683./8., c3 = -15309./8., c4 = 2187./4.;
  double tmp;
  if (t*freq < 0)
    tmp = 0.0;
  else if (t*freq > 1)
    tmp = 0.0;
  else
     tmp = 2*freq*(c0*6*t*freq+c1*12*pow(t*freq,2)+c2*20*pow(t*freq,3)+c3*30*pow(t*freq,4)+c4*42*pow(t*freq,5)) +
         freq*freq*t*(c0*6+c1*24*t*freq+c2*60*pow(t*freq,2)+c3*120*pow(t*freq,3)+c4*210*pow(t*freq,4));
  return tmp;
}

double Brune( double freq, double t, double* par )
{
  const double tf = t*freq;
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

double Brune_t( double freq, double t, double* par )
{
  const double tf = t*freq;
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

double Brune_om( double freq, double t, double* par )
{
  const double tf = t*freq;
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

double Brune_tt( double freq, double t, double* par )
{
  const double tf = t*freq;
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

double Brune_ttt( double freq, double t, double* par )
{
  const double tf = t*freq;
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

double Brune_omtt( double freq, double t, double* par )
{
  const double tf = t*freq;
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

double DBrune( double freq, double t, double* par )
{
  const double tf = t*freq;
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

double DBrune_t( double freq, double t, double* par )
{
  const double tf = t*freq;
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

double DBrune_om( double freq, double t, double* par )
{
  const double tf = t*freq;
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

double DBrune_tt( double freq, double t, double* par )
{
  const double tf = t*freq;
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

double DBrune_ttt( double freq, double t, double* par )
{
  const double tf = t*freq;
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

double DBrune_omtt( double freq, double t, double* par )
{
  const double tf = t*freq;
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

double BruneSmoothed( double freq, double t, double* par )
{
  const double tf = t*freq;
  const double h  = 2.31;
  const double hi = 1/h;
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

double BruneSmoothed_t( double freq, double t, double* par )
{
  const double tf = t*freq;
  const double h  = 2.31;
  const double hi = 1/h;
  if( tf < 0 )
    return 0;
  else if( tf < h )
  {
     const double c3 = - 1.5*hi;
     const double c4 = 1.5*hi*hi;
     const double c5 = -0.5*hi*hi*hi;
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

double BruneSmoothed_om( double freq, double t, double* par )
{
  const double tf = t*freq;
  const double h  = 2.31;
  const double hi = 1/h;
  if( tf < 0 )
    return 0;
  else if( tf < h )
  {
     const double c3 = - 1.5*hi;
     const double c4 = 1.5*hi*hi;
     const double c5 = -0.5*hi*hi*hi;
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

double BruneSmoothed_tt( double freq, double t, double* par )
{
  const double tf = t*freq;
  const double h  = 2.31;
  const double hi = 1/h;
  if( tf < 0 )
    return 0;
  else if( tf < h )
  {
     const double c3 = - 1.5*hi;
     const double c4 = 1.5*hi*hi;
     const double c5 = -0.5*hi*hi*hi;
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

double BruneSmoothed_ttt( double freq, double t, double* par )
{
  const double tf = t*freq;
  const double h  = 2.31;
  const double hi = 1/h;
  if( tf < 0 )
    return 0;
  else if( tf < h )
  {
     const double c3 = - 1.5*hi;
     const double c4 = 1.5*hi*hi;
     const double c5 = -0.5*hi*hi*hi;
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

double BruneSmoothed_omtt( double freq, double t, double* par )
{
  const double tf = t*freq;
  const double h  = 2.31;
  const double hi = 1/h;
  if( tf < 0 )
    return 0;
  else if( tf < h )
  {
     const double c3 = - 1.5*hi;
     const double c4 = 1.5*hi*hi;
     const double c5 = -0.5*hi*hi*hi;
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


double GaussianWindow( double freq, double t, double* par )
{
  double incyc2 = 1/(par[1]*par[1]);
  const double tf = t*freq;
  if( -0.5*tf*tf*incyc2  > par[0] )
    return sin(tf)*exp(-0.5*tf*tf*incyc2);
  else
    return 0;
}

double GaussianWindow_t( double freq, double t, double* par )
{
  double incyc2 = 1/(par[1]*par[1]);
  const double tf = t*freq;
  if( -0.5*tf*tf*incyc2  > par[0] )
     return (freq*cos(tf)-freq*tf*incyc2*sin(tf))*exp(-0.5*tf*tf*incyc2 );
  else
    return 0;
}

double GaussianWindow_om( double freq, double t, double* par )
{
  double incyc2 = 1/(par[1]*par[1]);
  const double tf = t*freq;
  if( -0.5*tf*tf*incyc2  > par[0] )
     return (t*cos(tf)-t*tf*incyc2*sin(tf))*exp(-0.5*tf*tf*incyc2 );
  else
    return 0;
}

double GaussianWindow_tt( double freq, double t, double* par )
{
  double incyc2 = 1/(par[1]*par[1]);
  const double tf = t*freq;
  if( -0.5*tf*tf*incyc2  > par[0] )
     return ( (-freq*freq-freq*freq*incyc2+freq*freq*tf*tf*incyc2*incyc2)*sin(tf)-
	      tf*2*freq*freq*incyc2*cos(tf) )*exp(-0.5*tf*tf*incyc2);
  else
    return 0;
}

double GaussianWindow_ttt( double freq, double t, double* par )
{
  double incyc2 = 1/(par[1]*par[1]);
  const double tf = t*freq;
  if( -0.5*tf*tf*incyc2  > par[0] )
     return ( freq*freq*freq*(3*tf*incyc2*(1+incyc2)-tf*tf*tf*incyc2*incyc2*incyc2)*sin(tf) +
	      freq*freq*freq*( 3*tf*tf*incyc2*incyc2-3*incyc2-1)*cos(tf))*exp(-0.5*tf*tf*incyc2);
  else
    return 0;
}

double GaussianWindow_omtt( double freq, double t, double* par )
{
  double incyc2 = 1/(par[1]*par[1]);
  const double tf = t*freq;
  if( -0.5*tf*tf*incyc2  > par[0] )
     return ( freq*(-2-2*incyc2 + 3*incyc2*tf*tf +5*tf*tf*incyc2*incyc2 -tf*tf*tf*tf*incyc2*incyc2*incyc2)*sin(tf) +
	      t*freq*freq*( 3*tf*tf*incyc2*incyc2-7*incyc2-1)*cos(tf) )*exp(-0.5*tf*tf*incyc2);
  else
    return 0;
}

double Liu( double freq, double t, double* par )
{
   double tau = 2*M_PI/freq;
   double tau1 = 0.13*tau;
   double tau2 = tau-tau1;
   if( t < 0 )
      return 0;
   else if( t >= tau )
      return 1;
   else
   {
      double ipi = 1.0/M_PI;
      double cn = 1.0/(1.4*tau1+1.2*tau1*ipi + 0.3*tau2);
      if( t <= tau1 )
	 return cn*(0.7*t-0.7*tau1*ipi*sin(M_PI*t/tau1)-1.2*tau1*ipi*(cos(0.5*M_PI*t/tau1)-1));
      else if( t <= 2*tau1 )
	 return cn*(1.0*t-0.3*tau1+1.2*tau1*ipi - 0.7*tau1*ipi*sin(M_PI*t/tau1)+0.3*tau2*ipi*sin(M_PI*(t-tau1)/tau2));
      else if( t <= tau )
	 return cn*(0.3*t+1.1*tau1+1.2*tau1*ipi+0.3*tau2*ipi*sin(M_PI*(t-tau1)/tau2));
   }
   return 0.; // should never get here, but keeps compiler happy
}

double Liu_t( double freq, double t, double* par )
{
   double tau = 2*M_PI/freq;
   double tau1 = 0.13*tau;
   double tau2 = tau-tau1;
   if( t < 0 )
      return 0;
   else if( t >= tau )
      return 0;
   else
   {
      double ipi = 1.0/M_PI;
      double cn = 1.0/(1.4*tau1+1.2*tau1*ipi + 0.3*tau2);
      if( t <= tau1 )
	 return cn*(0.7-0.7*cos(M_PI*t/tau1)+0.6*sin(0.5*M_PI*t/tau1));
      else if( t <= 2*tau1 )
	 return cn*(1-0.7*cos(M_PI*t/tau1)+0.3*cos(M_PI*(t-tau1)/tau2));
      else if( t <= tau )
	 return cn*(0.3+0.3*cos(M_PI*(t-tau1)/tau2));
   }
   return 0.; // should never get here, but keeps compiler happy
}

double Liu_om( double freq, double t, double* par )
{
   double tau = 2*M_PI/freq;
   double tau1 = 0.13*tau;
   double tau2 = tau-tau1;
   if( t < 0 )
      return 0;
   else if( t >= tau )
      return 0;
   else
   {
      double ipi = 1.0/M_PI;
      double cn = t*1.0/(1.4*tau1+1.2*tau1*ipi + 0.3*tau2)/freq;
      if( t <= tau1 )
	 return cn*(0.7-0.7*cos(M_PI*t/tau1)+0.6*sin(0.5*M_PI*t/tau1));
      else if( t <= 2*tau1 )
	 return cn*(1-0.7*cos(M_PI*t/tau1)+0.3*cos(M_PI*(t-tau1)/tau2));
      else if( t <= tau )
	 return cn*(0.3+0.3*cos(M_PI*(t-tau1)/tau2));
   }
   return 0.; // should never get here, but keeps compiler happy
}

double Liu_tt( double freq, double t, double* par )
{
   double tau = 2*M_PI/freq;
   double tau1 = 0.13*tau;
   double tau2 = tau-tau1;
   if( t < 0 )
      return 0;
   else if( t >= tau )
      return 0;
   else
   {
      double ipi = 1.0/M_PI;
      double cn = 1.0/(1.4*tau1+1.2*tau1*ipi + 0.3*tau2);
      if( t <= tau1 )
	 return cn*(0.7*M_PI*sin(M_PI*t/tau1)+0.3*M_PI*cos(0.5*M_PI*t/tau1))/tau1;
      else if( t <= 2*tau1 )
	 return cn*(0.7*M_PI*sin(M_PI*t/tau1)/tau1-0.3*M_PI*sin(M_PI*(t-tau1)/tau2)/tau2);
      else if( t <= tau )
	 return cn*(-0.3*M_PI*sin(M_PI*(t-tau1)/tau2))/tau2;
   }
   return 0.; // should never get here, but keeps compiler happy
}

double Liu_ttt( double freq, double t, double* par )
{
   double tau = 2*M_PI/freq;
   double tau1 = 0.13*tau;
   double tau2 = tau-tau1;
   if( t < 0 )
      return 0;
   else if( t >= tau )
      return 0;
   else
   {
      double ipi = 1.0/M_PI;
      double cn = 1.0/(1.4*tau1+1.2*tau1*ipi + 0.3*tau2);
      if( t <= tau1 )
	 return cn*(0.7*M_PI*M_PI*cos(M_PI/tau1*t)-0.15*M_PI*M_PI*sin(0.5*M_PI/tau1*t))/(tau1*tau1);
      else if( t <= 2*tau1 )
	 return cn*(0.7*M_PI*M_PI*cos(M_PI*t/tau1)/(tau1*tau1)-0.3*M_PI*M_PI*cos(M_PI*(t-tau1)/tau2)/(tau2*tau2));
      else if( t <= tau )
	 return cn*(-0.3*M_PI*M_PI*cos(M_PI*(t-tau1)/tau2))/(tau2*tau2);
   }
   return 0.; // should never get here, but keeps compiler happy
}

double Liu_omtt( double freq, double t, double* par )
{
   double tau = 2*M_PI/freq;
   double tau1 = 0.13*tau;
   double tau2 = tau-tau1;
   if( t < 0 )
      return 0;
   else if( t >= tau )
      return 0;
   else
   {
      double ipi = 1.0/M_PI;
      double cn = 1.0/(1.4*tau1+1.2*tau1*ipi + 0.3*tau2)/freq;
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

double NullFunc( double freq, double t, double* par )
{
// this function should never be called
  CHECK_INPUT(false,"The NullFunc time function was called!");
  return 0.;  
}

double Discrete( double freq, double t, double* par )
{
// note that 
// t holds EW::mTime - tStart
// freq holds 1/dt
// number of time steps:
  int nSteps = (int) par[0];
  
  int k = (int) (t*freq + 0.5);
  if (k>=0 && k<=nSteps)
    return par[k+1]; // offset by 1 since par[0] holds the number of time steps
  else
    return 0;
  
}

double Discrete_tt( double freq, double t, double* par )
{
// note that 
// t holds EW::mTime - tStart
// freq holds 1/dt
// number of time steps:
  int nSteps = (int) par[0];

  int kt = (int) (t*freq + 0.5);
  int k = kt+1; // offset by 1 since we save nSteps in par[0]
  
//
// centered differences in time
//
  double d2t=0;
  if (kt>0 && kt<nSteps)
  {
// freq=1/dt
    d2t = (par[k-1] - 2*par[k] + par[k+1])*freq*freq; // 2nd order for now
  }
  else if (kt<=0)
  {
    k = 2;
    d2t = (par[k-1] - 2*par[k] + par[k+1])*freq*freq; // 1st order for now
  }
  else
  {
    k = nSteps;
    d2t = (par[k-1] - 2*par[k] + par[k+1])*freq*freq; // 1st order for now
  }

  return d2t;
}

void butter2(double *b, double *a, double Wn)
{
  double fr = 2/Wn;
  double omegac = tan(M_PI/fr);
  double c = 1. + sqrt(2.)*omegac + omegac*omegac;

  b[0] = omegac*omegac / c;
  b[1] = 2*b[0];
  b[2] = b[0];
  a[0] = 1.;
  a[1] = 2.*(omegac*omegac - 1.)/c;
  a[2] = (1. - sqrt(2.)*omegac + omegac*omegac)/c;
}

// this routine will be replaced by a general direct form II algorithm
// (in the Filter class)
void filtfilt(double *b, double *a, double *u, int N)
{
// note this routine overwrites the input signal 'u'
  double op;
  int i;
  
// forwards 
  double x1=u[0];
  double x2=u[0];
  double y1=u[0];
  double y2=u[0];
  for (i=0; i<N; i++)
  {
    op = b[0]*u[i] + b[1]*x1 + b[2]*x2 - (a[1]*y1 + a[2]*y2);
    y2=y1;
    y1=op;
    x2=x1;
    x1=u[i];
    u[i]=op;
  }

// backwards
  x1=u[N-1];
  x2=u[N-1];
  y1=u[N-1];
  y2=u[N-1];
  for (i=N-1; i>=0; i--)
  {
    op = b[0]*u[i] + b[1]*x1 + b[2]*x2 - (a[1]*y1 + a[2]*y2);
    y2=y1;
    y1=op;
    x2=x1;
    x1=u[i];
    u[i]=op;
  }
}


GridPointSource::GridPointSource(double amplitude, double frequency, double t0,
				 int N, int M, int L, int G,
				 double Fx, double Fy, double Fz,
				 timeDep tDep,
				 int ncyc, double* jacobian,
				 double* dddp, double* hess1, double* hess2, double* hess3 ):
  mAmp(amplitude),
  mFreq(frequency),
  mT0(t0),
  m_i0(N), m_j0(M), m_k0(L),m_grid(G),
  mTimeDependence(tDep),
  mNcyc(ncyc),
  m_min_exponent(-700)
{
  mPar = new double[2];
  mForces[0] = Fx;
  mForces[1] = Fy;
  mForces[2] = Fz;
  initializeTimeFunction();
  m_derivative = -1;

  m_jacobian_known = jacobian != NULL;
  if( jacobian != NULL )
     for( int m=0 ; m < 27 ; m++ )
	m_jacobian[m] = jacobian[m];

  m_hessian_known = false;
  if( dddp != NULL && hess1 != NULL && hess2 != NULL && hess3 != NULL )
  {
     for( int m=0 ; m < 9 ; m++ )
     {
	m_dddp[m] = dddp[m];
        m_hesspos1[m] = hess1[m];
        m_hesspos2[m] = hess2[m];
        m_hesspos3[m] = hess3[m];
     }
     m_hessian_known = true;
  }
}

GridPointSource::~GridPointSource()
{
  //   if( mPar != NULL )
  delete[] mPar;
}

void GridPointSource::initializeTimeFunction()
{
  mPar[0] = m_min_exponent;
  switch(mTimeDependence)
    {
    case iRicker :
      mTimeFunc = RickerWavelet;
      mTimeFunc_t = RickerWavelet_t;
      mTimeFunc_tt = RickerWavelet_tt;
      mTimeFunc_ttt = RickerWavelet_ttt;
      mTimeFunc_om = RickerWavelet_om;
      mTimeFunc_omtt = RickerWavelet_omtt;
      break;
    case iGaussian :
      mTimeFunc   = Gaussian;
      mTimeFunc_t = Gaussian_t;
      mTimeFunc_tt = Gaussian_tt;
      mTimeFunc_ttt = Gaussian_ttt;
      mTimeFunc_om = Gaussian_om;
      mTimeFunc_omtt = Gaussian_omtt;
      break;
    case iRamp :
      mTimeFunc = Ramp;
      mTimeFunc_t = Ramp_t;
      mTimeFunc_tt = Ramp_tt;
      mTimeFunc_ttt = Ramp_ttt;
      mTimeFunc_om = Ramp_om;
      mTimeFunc_omtt = Ramp_omtt;
      break;
    case iTriangle :
      mTimeFunc = Triangle;
      mTimeFunc_t = Triangle_t;
      mTimeFunc_tt = Triangle_tt;
      mTimeFunc_ttt = Triangle_ttt;
      mTimeFunc_om = Triangle_om;
      mTimeFunc_omtt = Triangle_omtt;
      break;
    case iSawtooth :
      mTimeFunc = Sawtooth;
      mTimeFunc_t = Sawtooth_t;
      mTimeFunc_tt = Sawtooth_tt;
      mTimeFunc_ttt = Sawtooth_ttt;
      mTimeFunc_om = Sawtooth_om;
      mTimeFunc_omtt = Sawtooth_omtt;
      break;
    case iSmoothWave :
      mTimeFunc = SmoothWave;
      mTimeFunc_t = SmoothWave_t;
      mTimeFunc_tt = SmoothWave_tt;
      mTimeFunc_ttt = SmoothWave_ttt;
      mTimeFunc_om = SmoothWave_om;
      mTimeFunc_omtt = SmoothWave_omtt;
      break;
    case iErf :
      mTimeFunc = Erf;
      mTimeFunc_t = Erf_t;
      mTimeFunc_tt = Erf_tt;
      mTimeFunc_ttt = Erf_ttt;
      mTimeFunc_om = Erf_om;
      mTimeFunc_omtt = Erf_omtt;
      break;
    case iVerySmoothBump :
      mTimeFunc = VerySmoothBump;
      mTimeFunc_t = VerySmoothBump_t;
      mTimeFunc_tt = VerySmoothBump_tt;
      mTimeFunc_ttt = VerySmoothBump_ttt;
      mTimeFunc_om = VerySmoothBump_om;
      mTimeFunc_omtt = VerySmoothBump_omtt;
      break;
    case iRickerInt :
      mTimeFunc = RickerInt;
      mTimeFunc_t = RickerInt_t;
      mTimeFunc_tt = RickerInt_tt;
      mTimeFunc_ttt = RickerInt_ttt;
      mTimeFunc_om = RickerInt_om;
      mTimeFunc_omtt = RickerInt_omtt;
      break;
    case iBrune :
      mTimeFunc = Brune;
      mTimeFunc_t = Brune_t;
      mTimeFunc_tt = Brune_tt;
      mTimeFunc_ttt = Brune_ttt;
      mTimeFunc_om = Brune_om;
      mTimeFunc_omtt = Brune_omtt;
      break;
    case iBruneSmoothed :
      mTimeFunc = BruneSmoothed;
      mTimeFunc_t = BruneSmoothed_t;
      mTimeFunc_tt = BruneSmoothed_tt;
      mTimeFunc_ttt = BruneSmoothed_ttt;
      mTimeFunc_om = BruneSmoothed_om;
      mTimeFunc_omtt = BruneSmoothed_omtt;
      break;
    case iDBrune :
      mTimeFunc = DBrune;
      mTimeFunc_t = DBrune_t;
      mTimeFunc_tt = DBrune_tt;
      mTimeFunc_ttt = DBrune_ttt;
      mTimeFunc_om = DBrune_om;
      mTimeFunc_omtt = DBrune_omtt;
      break;
    case iGaussianWindow :
      //      if( mPar == NULL )
      //          mPar = new double[1];
      mPar[1] = mNcyc;
      mTimeFunc = GaussianWindow;
      mTimeFunc_t = GaussianWindow_t;
      mTimeFunc_tt = GaussianWindow_tt;
      mTimeFunc_ttt = GaussianWindow_ttt;
      mTimeFunc_om = GaussianWindow_om;
      mTimeFunc_omtt = GaussianWindow_omtt;
      break;
    case iLiu :
       mTimeFunc = Liu;
       mTimeFunc_t = Liu_t;
       mTimeFunc_tt = Liu_tt;
       mTimeFunc_ttt = Liu_ttt;
       mTimeFunc_om = Liu_om;
       mTimeFunc_omtt = Liu_omtt;
       break;
    default :
      std::cout << "incorrect argument to GridPointSource constructor : default RickerWavelet used " << std::endl;
      mTimeFunc = RickerWavelet;
      mTimeFunc_t = RickerWavelet_t;
      mTimeFunc_tt = RickerWavelet_tt;
      mTimeFunc_ttt = RickerWavelet_ttt;
      mTimeFunc_om = RickerWavelet_om;
      mTimeFunc_omtt = RickerWavelet_omtt;
    }
  // Treat fourth derivatives in special 'switch', because not (yet?) implemented for all time functions
  switch( mTimeDependence )
  {
  case iVerySmoothBump :
     mTimeFunc_tttt = VerySmoothBump_tttt;
     mTimeFunc_tttom = VerySmoothBump_tttom;
     mTimeFunc_ttomom = VerySmoothBump_ttomom;
     mTimeFunc_tom = VerySmoothBump_tom;
     mTimeFunc_omom = VerySmoothBump_omom;
     break;
  case iGaussian :
     mTimeFunc_tttt = Gaussian_tttt;
     mTimeFunc_tttom = Gaussian_tttom;
     mTimeFunc_ttomom = Gaussian_ttomom;
     mTimeFunc_tom = Gaussian_tom;
     mTimeFunc_omom = Gaussian_omom;
     break;
  default: 
// tmp
// std::cout << "High derivatives not implemented for time fuction:" << mTimeDependence <<
//   " default Gaussian used for tttt, ttt-omega derivatives, etc " << std::endl;
     mTimeFunc_tttt = Gaussian_tttt;
     mTimeFunc_tttom = Gaussian_tttom;
     mTimeFunc_ttomom = Gaussian_ttomom;
     mTimeFunc_tom = Gaussian_tom;
     mTimeFunc_omom = Gaussian_omom;
  }
}

void GridPointSource::getFxyz( double t, double* fxyz ) const
{
  double afun = mAmp*mTimeFunc(mFreq,t-mT0,mPar);
  if( m_derivative==-1)
  {
     fxyz[0] = mForces[0]*afun;
     fxyz[1] = mForces[1]*afun;
     fxyz[2] = mForces[2]*afun;
  }
  else if( m_derivative >= 0 && m_derivative <= 8 )
  {
     fxyz[0] = m_jacobian[m_derivative*3]*afun;
     fxyz[1] = m_jacobian[m_derivative*3+1]*afun;
     fxyz[2] = m_jacobian[m_derivative*3+2]*afun;
  }
  else if( m_derivative == 9 )
  {
     afun = -mAmp*mTimeFunc_t(mFreq,t-mT0,mPar);
     fxyz[0] = mForces[0]*afun;
     fxyz[1] = mForces[1]*afun;
     fxyz[2] = mForces[2]*afun;
  }
  else if( m_derivative == 10 )
  {
     afun = mAmp*mTimeFunc_om(mFreq,t-mT0,mPar);
     fxyz[0] = mForces[0]*afun;
     fxyz[1] = mForces[1]*afun;
     fxyz[2] = mForces[2]*afun;
  }
  else if( m_derivative == 11 )
  {
     fxyz[0]= fxyz[1]=fxyz[2]=0;
     int i;
     for( i=0 ; i < 9; i++ )
     {
        fxyz[0] += afun*m_jacobian[i*3]*m_dir[i];
        fxyz[1] += afun*m_jacobian[i*3+1]*m_dir[i];
        fxyz[2] += afun*m_jacobian[i*3+2]*m_dir[i];
     }
     i = 9;
     afun = -mAmp*mTimeFunc_t(mFreq,t-mT0,mPar);
     fxyz[0] += afun*mForces[0]*m_dir[i];
     fxyz[1] += afun*mForces[1]*m_dir[i];
     fxyz[2] += afun*mForces[2]*m_dir[i];
     i = 10;
     afun =  mAmp*mTimeFunc_om(mFreq,t-mT0,mPar);
     fxyz[0] += afun*mForces[0]*m_dir[i];
     fxyz[1] += afun*mForces[1]*m_dir[i];
     fxyz[2] += afun*mForces[2]*m_dir[i];
  }
}

void GridPointSource::getFxyz_notime( double* fxyz ) const
{
// For source spatial discretization testing
  fxyz[0] = mForces[0]*mAmp;
  fxyz[1] = mForces[1]*mAmp;
  fxyz[2] = mForces[2]*mAmp;
}

void GridPointSource::getFxyztt( double t, double* fxyz ) const
{
  double afun = mAmp*mTimeFunc_tt(mFreq,t-mT0,mPar);
  if( m_derivative==-1)
  {
     fxyz[0] = mForces[0]*afun;
     fxyz[1] = mForces[1]*afun;
     fxyz[2] = mForces[2]*afun;
  }
  else if( m_derivative >= 0 && m_derivative <= 8 )
  {
     fxyz[0] = m_jacobian[m_derivative*3]*afun;
     fxyz[1] = m_jacobian[m_derivative*3+1]*afun;
     fxyz[2] = m_jacobian[m_derivative*3+2]*afun;
  }
  else if( m_derivative == 9 )
  {
     afun = -mAmp*mTimeFunc_ttt(mFreq,t-mT0,mPar);
     fxyz[0] = mForces[0]*afun;
     fxyz[1] = mForces[1]*afun;
     fxyz[2] = mForces[2]*afun;
  }
  else if( m_derivative == 10 )
  {
     afun = mAmp*mTimeFunc_omtt(mFreq,t-mT0,mPar);
     fxyz[0] = mForces[0]*afun;
     fxyz[1] = mForces[1]*afun;
     fxyz[2] = mForces[2]*afun;
  }
  else if( m_derivative == 11 )
  {
     fxyz[0]= fxyz[1]=fxyz[2]=0;
     int i;
     for( i=0 ; i < 9; i++ )
     {
        fxyz[0] += afun*m_jacobian[i*3]*m_dir[i];
        fxyz[1] += afun*m_jacobian[i*3+1]*m_dir[i];
        fxyz[2] += afun*m_jacobian[i*3+2]*m_dir[i];
     }
     i = 9;
     afun = -mAmp*mTimeFunc_ttt(mFreq,t-mT0,mPar);
     fxyz[0] += afun*mForces[0]*m_dir[i];
     fxyz[1] += afun*mForces[1]*m_dir[i];
     fxyz[2] += afun*mForces[2]*m_dir[i];
     i = 10;
     afun =  mAmp*mTimeFunc_omtt(mFreq,t-mT0,mPar);
     fxyz[0] += afun*mForces[0]*m_dir[i];
     fxyz[1] += afun*mForces[1]*m_dir[i];
     fxyz[2] += afun*mForces[2]*m_dir[i];
  }
}

//-----------------------------------------------------------------------
void GridPointSource::set_derivative( int der, double dir[11] )
{
   if( der >= 0 && der <= 11 )
      m_derivative = der;
   for( int i=0 ; i < 11 ; i++ )
      m_dir[i] = dir[i];
}

//-----------------------------------------------------------------------
void GridPointSource::set_noderivative( )
{
   m_derivative = -1;
}

//-----------------------------------------------------------------------
void GridPointSource::limitFrequency(double max_freq)
{
  if (mFreq > max_freq)
    mFreq=max_freq;
}

//-----------------------------------------------------------------------
double GridPointSource::getTimeFunc(double t) const
{
  return mTimeFunc(mFreq, t - mT0, mPar);
}

//-----------------------------------------------------------------------
void GridPointSource::discretizeTimeFuncAndFilter(double tStart, double dt, int nSteps, Filter *filter_ptr)
{
// allocate new parameter array of size(nSteps+1)
  double *new_par= new double[nSteps+2];

// save the number of time steps 
  new_par[0] = (double) nSteps+0.5; // add the 0.5 here so we don't have to worry 
                                    // about round-off when casting it back to an int
  
// tmp
  // int ipar0 = (int) new_par[0];
  // printf("discretizeTimeFunc: test: (int) new_par[0] = %i, nSteps=%i\n", ipar0, nSteps);

// evaluate current time function at t_k = tStart + k*dt for k=0,1,2,...,nSteps
  double t=tStart;
  for (int k=0; k<=nSteps; k++)
  {
    new_par[k+1] = getTimeFunc(tStart+k*dt);
  }
  
// free current parameter array
  delete[] mPar;
  
// point old parameter pointer to new parameter array
  mPar = new_par;

// save 1/dt in mFreq
  mFreq = 1./dt; // correct dimension and also convenient for evaluation purposes: see function Discrete()

// save tStart in mT0
  mT0 = tStart;
  
// update the fcn pointer and type
  mTimeFunc = Discrete;
  mTimeFunc_tt = Discrete_tt;
// the following are not (yet) defined
  mTimeFunc_t = NullFunc;   // not defined
  mTimeFunc_ttt = NullFunc; // not defined
  mTimeFunc_om = NullFunc;  // not defined
  mTimeFunc_omtt = NullFunc;// not defined
       
  mTimeDependence = iDiscrete;

//
// the following will be replaced by a Filter object that defines more general low-pass and band-pass filters
//
//  double fc = filter_ptr->get_corner_freq2();
// calculate coefficients in the butterworthfilter
//  double a[3], b[3];
//  butter2(b,a,2*dt*fc);
// tmp
//   printf("Butterworth coeff for dt=%.18e, fc=%e\n", dt, fc);
//   printf("a[0]=%e, a[1]=%e, a[2]=%e\n", a[0], a[1], a[2]);
//   printf("b[0]=%e, b[1]=%e, b[2]=%e\n", b[0], b[1], b[2]);
  
// now filter (forwards + backwards) the time series with a lowpass butterworth filter of order 2
//  filtfilt(b, a, &mPar[1], nSteps+1);

// this call corresponds to passes=2
  filter_ptr->zerophase(nSteps+1, &mPar[1], &mPar[1]);
  
}

//-----------------------------------------------------------------------
ostream& operator<<( ostream& output, const GridPointSource& s )
{
   output << "GridPointSource at (i,j,k) = " << s.m_i0 << "," << s.m_j0 << "," << s.m_k0 << 
     " in grid no " << s.m_grid << endl;
   output << "   Strength " << s.mAmp;
   output << " Fx Fy Fz = " << s.mForces[0] << " " << s.mForces[1] << " " << s.mForces[2] << endl;

   return output;
}

//-----------------------------------------------------------------------
void GridPointSource::add_to_gradient( std::vector<Sarray> & kappa, std::vector<Sarray> & eta,
				       double t, double dt, double gradient[11], std::vector<double> & h )
{
   if( m_jacobian_known )
   {
      double dt2o12 = dt*dt/12.0;
      double g0= mTimeFunc( mFreq, t-mT0, mPar );
      double g = g0 + dt2o12*mTimeFunc_tt( mFreq, t-mT0, mPar);


      // save some work by accessing array elements only once:
      double kap1 = kappa[m_grid](1,m_i0,m_j0,m_k0);
      double kap2 = kappa[m_grid](2,m_i0,m_j0,m_k0);
      double kap3 = kappa[m_grid](3,m_i0,m_j0,m_k0);
      double eta1 = eta[m_grid](1,m_i0,m_j0,m_k0);
      double eta2 = eta[m_grid](2,m_i0,m_j0,m_k0);
      double eta3 = eta[m_grid](3,m_i0,m_j0,m_k0);
      double h3   = h[m_grid]*h[m_grid]*h[m_grid];
      //      double h3 = 1.0;

      // derivative wrt. position (m=0,1,2) and moment tensor components (m=3,..,8)
      for( int m= 0 ; m < 9 ; m++ )
      {
	 gradient[m] -= g*mAmp*(  kap1*m_jacobian[3*m] + kap2*m_jacobian[3*m+1] +
				  kap3*m_jacobian[3*m+2]  )*h3;
	 gradient[m] -= dt2o12*g0*mAmp*( eta1*m_jacobian[3*m] + eta2*m_jacobian[3*m+1] +
					 eta3*m_jacobian[3*m+2]  )*h3;
      }

      // derivative wrt. (t0, freq)
      double dgt0 = -mTimeFunc_t(mFreq,t-mT0,mPar);
      double dgom =  mTimeFunc_om(mFreq,t-mT0,mPar);
      gradient[9]   -= dt2o12*dgt0*mAmp*( eta1*mForces[0] + eta2*mForces[1] + eta3*mForces[2])*h3;
      gradient[10]  -= dt2o12*dgom*mAmp*( eta1*mForces[0] + eta2*mForces[1] + eta3*mForces[2])*h3;

      dgt0 = dgt0 - dt2o12*mTimeFunc_ttt(mFreq,t-mT0,mPar);
      dgom = dgom + dt2o12*mTimeFunc_omtt(mFreq,t-mT0,mPar);
      gradient[9]  -= dgt0*mAmp*( kap1*mForces[0] + kap2*mForces[1] + kap3*mForces[2])*h3;
      gradient[10] -= dgom*mAmp*( kap1*mForces[0] + kap2*mForces[1] + kap3*mForces[2])*h3;
   }
}

//-----------------------------------------------------------------------
void GridPointSource::add_to_hessian( std::vector<Sarray> & kappa, std::vector<Sarray> & eta,
				      double t, double dt, double hessian[121],
				      std::vector<double> & h )
// Add upper part of symmetric matrix
{
   if( m_hessian_known && m_jacobian_known )
   {
      double dt2o12 = dt*dt/12.0;
      double g0= mTimeFunc( mFreq, t-mT0, mPar );
      double g = g0 + dt2o12*mTimeFunc_tt( mFreq, t-mT0, mPar);

      // save some work by accessing array elements only once:
      double kap1 = kappa[m_grid](1,m_i0,m_j0,m_k0);
      double kap2 = kappa[m_grid](2,m_i0,m_j0,m_k0);
      double kap3 = kappa[m_grid](3,m_i0,m_j0,m_k0);
      double eta1 = eta[m_grid](1,m_i0,m_j0,m_k0);
      double eta2 = eta[m_grid](2,m_i0,m_j0,m_k0);
      double eta3 = eta[m_grid](3,m_i0,m_j0,m_k0);
      double h3   = h[m_grid]*h[m_grid]*h[m_grid];

      double c1 = mAmp*g*h3;
      double c2 = mAmp*g0*dt2o12*h3;

      // (pos,pos)
      for( int m=0 ; m < 3 ; m++ )
	 for( int j=m ; j < 3 ; j++ )
	 {
	    hessian[m+11*j] -= (kap1*c1+eta1*c2)*m_hesspos1[m+3*j]+(kap2*c1+eta2*c2)*m_hesspos2[m+3*j]+
	       (kap3*c1+eta3*c2)*m_hesspos3[m+3*j];
	 }
      // (pos,mij)
      for( int m = 0; m < 3 ; m++ )
      {
	 int j=3;
	 hessian[m+11*j] -= (c1*kap1+eta1*c2)*m_dddp[m];
	 j = 4;
	 hessian[m+11*j] -= (c1*kap1+eta1*c2)*m_dddp[m+3] + (c1*kap2+eta2*c2)*m_dddp[m];
	 j = 5;
	 hessian[m+11*j] -= (c1*kap1+eta1*c2)*m_dddp[m+6] + (c1*kap3+eta3*c2)*m_dddp[m];      
	 j = 6;
	 hessian[m+11*j] -= (c1*kap2+eta2*c2)*m_dddp[m+3];
	 j = 7;
	 hessian[m+11*j] -= (c1*kap2+eta2*c2)*m_dddp[m+6] + (c1*kap3+eta3*c2)*m_dddp[m+3];
	 j = 8;
	 hessian[m+11*j] -= (c1*kap3+eta3*c2)*m_dddp[m+6];
      }
      // (pos,t0)
      double dgt0 = -mTimeFunc_t(mFreq,t-mT0,mPar);
      double dgom =  mTimeFunc_om(mFreq,t-mT0,mPar);

      double c2t0  = dgt0*mAmp*dt2o12;
      double c2om0 = dgom*mAmp*dt2o12;

      dgt0 = dgt0 - dt2o12*mTimeFunc_ttt(mFreq,t-mT0,mPar);
      dgom = dgom + dt2o12*mTimeFunc_omtt(mFreq,t-mT0,mPar);
      double c1t0  = dgt0*mAmp;
      double c1om0 = dgom*mAmp;

      //define c1t0,c2t0, c1om0, c2om0
      for( int m=0 ; m < 9 ; m++ )
      {
	 int j = 9;
	 hessian[m+11*j] -= ((c1t0*kap1+c2t0*eta1)*m_jacobian[3*m] + (c1t0*kap2+c2t0*eta2)*m_jacobian[3*m+1] +
			     (c1t0*kap3+c2t0*eta3)*m_jacobian[3*m+2])*h3;
	 j = 10;
	 hessian[m+11*j] -= ((c1om0*kap1+c2om0*eta1)*m_jacobian[3*m] + (c1om0*kap2+c2om0*eta2)*m_jacobian[3*m+1] +
			     (c1om0*kap3+c2om0*eta3)*m_jacobian[3*m+2])*h3;
      }
      double cmfact0 = mAmp*(kap1*mForces[0]+kap2*mForces[1]+kap3*mForces[2]);
      double cmfact  = mAmp*((kap1*mForces[0]+kap2*mForces[1]+kap3*mForces[2])+
		      dt2o12*(eta1*mForces[0]+eta2*mForces[1]+eta3*mForces[2]));
      // Second derivatives of time function
      double d2gdt02   =  mTimeFunc_tt(mFreq,t-mT0,mPar);
      double d2gdt0dom = -mTimeFunc_tom(mFreq,t-mT0,mPar);
      double d2gdomdom =  mTimeFunc_omom(mFreq,t-mT0,mPar);

      double dgdttt0t0  =  dt2o12*mTimeFunc_tttt(mFreq,t-mT0,mPar);
      double dgdttt0om  = -dt2o12*mTimeFunc_tttom(mFreq,t-mT0,mPar);
      double dgdttomom  =  dt2o12*mTimeFunc_ttomom(mFreq,t-mT0,mPar);

      int m = 9;
      int j = 9;
      hessian[m+11*j] -= (cmfact*d2gdt02 + cmfact0*dgdttt0t0)*h3;
      j = 10;
      hessian[m+11*j] -= (cmfact*d2gdt0dom + cmfact0*dgdttt0om)*h3;
      m = 10;
      j = 10;
      hessian[m+11*j] -= (cmfact*d2gdomdom + cmfact0*dgdttomom)*h3;
   }
}

//-----------------------------------------------------------------------
void GridPointSource::print_info()
{
   cout << "-----------------------------------------------------------------------"<<endl;
   cout << " position " << m_i0 << " " << m_j0 << " " << m_k0 << endl;
   cout << "Forces = " << mForces[0] << " " << mForces[1] << " " << mForces[2] << endl;
   cout << " jac = \n";
   for( int c=0 ; c < 3 ; c++ )
   {
       for( int m=0 ; m < 9 ; m++ )
	  cout << m_jacobian[c+3*m] << " " ;
       cout << endl;
   }
   cout << "-----------------------------------------------------------------------"<<endl;
}
