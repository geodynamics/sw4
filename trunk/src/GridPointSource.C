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

double RickerWavelet(double freq, double t, double* par )
{
  double factor = pow(M_PI*freq*t,2);
  if( -factor > par[0] )
    return (2*factor - 1)*exp(-factor);
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

double Gaussian(double freq, double t, double* par )
{
  double factor=pow(t*freq,2) / 2;
  if( -factor > par[0] )
    return freq / sqrt(2*M_PI)*exp(-factor);
  else
    return 0;
}

double Erf( double freq, double t, double* par )
{
  return 0.5*(1+erf( freq*t/sqrt(2.0)) );
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


double GaussianWindow( double freq, double t, double* par )
{
  double incyc2 = 1/(par[1]*par[1]);
  const double tf = t*freq;
  if( -0.5*tf*tf*incyc2  > par[0] )
    return sin(tf)*exp(-0.5*tf*tf*incyc2);
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

double Discrete( double freq, double t, double* par )
{
// note that 
// t holds EW::mTime - tStart
// freq holds 1/dt
  int k = (int) (t*freq + 0.5);
  return par[k];
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
				 int N, int M, int L,int G,
				 double Fx, double Fy, double Fz,
				 timeDep tDep,
				 int ncyc ):
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
      break;
    case iGaussian :
      mTimeFunc = Gaussian;
      break;
    case iRamp :
      mTimeFunc = Ramp;
      break;
    case iTriangle :
      mTimeFunc = Triangle;
      break;
    case iSawtooth :
      mTimeFunc = Sawtooth;
      break;
    case iSmoothWave :
      mTimeFunc = SmoothWave;
      break;
    case iErf :
      mTimeFunc = Erf;
      break;
    case iVerySmoothBump :
      mTimeFunc = VerySmoothBump;
      break;
    case iRickerInt :
      mTimeFunc = RickerInt;
      break;
    case iBrune :
      mTimeFunc = Brune;
      break;
    case iBruneSmoothed :
      mTimeFunc = BruneSmoothed;
      break;
    case iDBrune :
      mTimeFunc = DBrune;
      break;
    case iGaussianWindow :
      //      if( mPar == NULL )
      //          mPar = new double[1];
      mPar[1] = mNcyc;
      mTimeFunc = GaussianWindow;
      break;
    case iLiu :
       mTimeFunc = Liu;
       break;
    default :
      std::cout << "erroneous argument to GridPointSource constructor : default RickerWavelet used " << std::endl;
      mTimeFunc = RickerWavelet;
    }
}

void GridPointSource::getFxyz( double t, double* fxyz ) const
{
  double afun = mAmp*mTimeFunc(mFreq,t-mT0,mPar);
  fxyz[0] = mForces[0]*afun;
  fxyz[1] = mForces[1]*afun;
  fxyz[2] = mForces[2]*afun;
}

void GridPointSource::limitFrequency(double max_freq)
{
  if (mFreq > max_freq)
    mFreq=max_freq;
}

double GridPointSource::getTimeFunc(double t) const
{
  return mTimeFunc(mFreq, t - mT0, mPar);
}

void GridPointSource::discretizeTimeFuncAndFilter(double tStart, double dt, int nSteps, double fc)
{
// allocate new parameter array of size(nSteps)
  double *new_par= new double[nSteps];

// evaluate current time function at t_k = tStart + k*dt for k=0,1,2,...,nSteps-1
  double t=tStart;
  for (int k=0; k< nSteps; k++)
  {
    new_par[k] = getTimeFunc(t);
    t += dt;
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
  mTimeDependence = iDiscrete;

// calculate coefficients in the butterworthfilter
  double a[3], b[3];
  butter2(b,a,2*dt*fc);
// tmp
//   printf("Butterworth coeff for dt=%.18e, fc=%e\n", dt, fc);
//   printf("a[0]=%e, a[1]=%e, a[2]=%e\n", a[0], a[1], a[2]);
//   printf("b[0]=%e, b[1]=%e, b[2]=%e\n", b[0], b[1], b[2]);
  
// now filter (forwards + backwards) the time series with a lowpass butterworth filter of order 2
  filtfilt(b, a, mPar, nSteps);
}

//---- print operator --------------
//-----------------------------------------------------------------------
ostream& operator<<( ostream& output, const GridPointSource& s )
{
   output << "GridPointSource at (i,j,k) = " << s.m_i0 << "," << s.m_j0 << "," << s.m_k0 << 
     " in grid no " << s.m_grid << endl;
   output << "   Strength " << s.mAmp;
   output << " Fx Fy Fz = " << s.mForces[0] << " " << s.mForces[1] << " " << s.mForces[2] << endl;

   return output;
}
