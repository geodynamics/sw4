#include "mpi.h"
#include "GridPointSource.h"
#include "Source.h"
#include "Require.h"

#include <fenv.h>
#include <cmath>

#include  "EW.h"
#include "Filter.h"
#include "Qspline.h"

#include "time_functions.h"


using namespace std;

#define SQR(x) ((x)*(x))

//-----------------------------------------------------------------------
// Constructor, 
//
//    ncyc is only used in the 'GaussianWindow' time function
//
//    pars is only used in the 'Discrete' time function
//        when pars[1],..pars[npts] should contain the discrete function on a uniform
//        grid with spacing dt=1/freq, and pars[0] is the first time, thus the grid is
//            t_k = pars[0] + dt*k, k=0,1,..,npts-1
//    ipar should have size 1, with ipar[0] containing npts.
//
//    When the source time function is not 'Discrete', the input pars and ipars will
//    not be used.
//
Source::Source(EW *a_ew, 
	       double frequency, 
	       double t0,
	       double x0, 
	       double y0, 
	       double z0,
	       double Mxx,
	       double Mxy,
	       double Mxz,
	       double Myy,
	       double Myz,
	       double Mzz,
	       timeDep tDep,
	       const char *name,
	       int ncyc, double* pars, int npar, int* ipars, int nipar ):
  mIsMomentSource(true),
  mFreq(frequency),
  mT0(t0),
  m_zRelativeToTopography(false),
  mX0(x0), mY0(y0), mZ0(z0),
  mGridPointSet(false),
  mTimeDependence(tDep),
  mNcyc(ncyc),
  m_derivative(-1),
  m_is_filtered(false)
{
   mForces.resize(6);
   mForces[0] = Mxx;
   mForces[1] = Mxy;
   mForces[2] = Mxz;
   mForces[3] = Myy;
   mForces[4] = Myz;
   mForces[5] = Mzz;
   mName = name;

   a_ew->computeNearestGridPoint(m_i0,m_j0,m_k0,m_grid,mX0,mY0,mZ0);

   mNpar = npar;
   if( mNpar > 0 )
   {
      mPar   = new double[mNpar];
      for( int i= 0 ; i < mNpar ; i++ )
	 mPar[i] = pars[i];
   }
   else
   {
      mNpar = 2;
      mPar = new double[2];
   }
   mNipar = nipar;
   if( mNipar > 0 )
   {
      mIpar = new int[mNipar];
      for( int i= 0 ; i < mNipar ; i++ )
         mIpar[i] = ipars[i];
   }
   else
   {
      mNipar = 1;
      mIpar  = new int[1];
   }

   if( mTimeDependence == iDiscrete )
      spline_interpolation();
   else
   {
      mPar[0] = find_min_exponent();
      mPar[1] = mNcyc;
   }
}

//-----------------------------------------------------------------------
Source::Source(EW *a_ew, double frequency, double t0,
	       double x0, double y0, double z0,
	       double Fx,
	       double Fy,
	       double Fz,
	       timeDep tDep,
	       const char *name, int ncyc, double* pars, int npar, int* ipars, int nipar ):
  mIsMomentSource(false),
  mFreq(frequency),
  mT0(t0),
  m_zRelativeToTopography(false),
  mX0(x0), mY0(y0), mZ0(z0),
  mGridPointSet(false),
  mTimeDependence(tDep),
  mNcyc(ncyc),
  m_derivative(-1),
  m_is_filtered(false)
{
  mForces.resize(3);
  mForces[0] = Fx;
  mForces[1] = Fy;
  mForces[2] = Fz;
  mName = name;

  a_ew->computeNearestGridPoint(m_i0,m_j0,m_k0,m_grid,mX0,mY0,mZ0);

  mNpar = npar;
  if( mNpar > 0 )
  {
     mPar   = new double[mNpar];
     for( int i= 0 ; i < mNpar ; i++ )
	mPar[i] = pars[i];
  }
  else
  {
     mNpar = 2;
     mPar = new double[2];
  }

  mNipar = nipar;
  if( mNipar > 0 )
  {
     mIpar = new int[mNipar];
     for( int i= 0 ; i < mNipar ; i++ )
        mIpar[i] = ipars[i];
  }
  else
  {
     mNipar = 1;
     mIpar  = new int[1];
  }
  if( mTimeDependence == iDiscrete )
     spline_interpolation();
  else
  {
     mPar[0] = find_min_exponent();
     mPar[1] = mNcyc;
  }
}

//-----------------------------------------------------------------------
Source::Source()
{

}

//-----------------------------------------------------------------------
Source::~Source()
{
   delete[] mPar;
   delete[] mIpar;
}

//-----------------------------------------------------------------------
double Source::getX0() const
{
  return mX0;
}

//-----------------------------------------------------------------------
double Source::getY0() const
{
  return mY0;
}

//-----------------------------------------------------------------------
double Source::getZ0() const
{
  return mZ0;
}

//-----------------------------------------------------------------------
double Source::getOffset() const
{
  return mT0;
}

//-----------------------------------------------------------------------
double Source::getFrequency() const
{
  return mFreq;
}

//-----------------------------------------------------------------------
void Source::setMaxFrequency(double max_freq)
{
  if (mFreq > max_freq)
    mFreq=max_freq;
}

//-----------------------------------------------------------------------
bool Source::isMomentSource() const
{
  return mIsMomentSource;
}

//-----------------------------------------------------------------------
void Source::getForces( double& fx, double& fy, double& fz ) const
{
   if( !mIsMomentSource )
   {
      fx = mForces[0];
      fy = mForces[1];
      fz = mForces[2];
   }
   else
      fx = fy = fz = 0;
}

//-----------------------------------------------------------------------
void Source::getMoments( double& mxx, double& mxy, double& mxz, double& myy, double& myz, double& mzz ) const
{
   if( mIsMomentSource )
   {
      mxx = mForces[0];
      mxy = mForces[1];
      mxz = mForces[2];
      myy = mForces[3];
      myz = mForces[4];
      mzz = mForces[5];
   }
   else
      mxx = mxy = mxz = myy = myz = mzz = 0;
}

//-----------------------------------------------------------------------
void Source::setMoments( double mxx, double mxy, double mxz, double myy, double myz, double mzz )
{
   if( mIsMomentSource )
   {
      
      mForces[0] = mxx;
      mForces[1] = mxy;
      mForces[2] = mxz;
      mForces[3] = myy;
      mForces[4] = myz;
      mForces[5] = mzz;
   }
   else
   {
      mForces[0] = mxx;
      mForces[1] = myy;
      mForces[2] = mzz;
   }
}

//-----------------------------------------------------------------------
double Source::getAmplitude() const
{
  double amplitude=0;
  if (mIsMomentSource)
  {
    double msqr=0;
    for (int q=0; q<6; q++)
      msqr += SQR(mForces[q]);
    //    amplitude = mAmp*sqrt(msqr/2.);
    msqr += SQR(mForces[1])+SQR(mForces[2])+SQR(mForces[4]);
    amplitude = sqrt(0.5*msqr);
  }
  else
  {
    double fsqr=0;
    for (int q=0; q<3; q++)
      fsqr += SQR(mForces[q]);
    //    amplitude = mAmp*sqrt(fsqr);
    amplitude = sqrt(fsqr);
  }
  return amplitude;
}

//-----------------------------------------------------------------------
ostream& operator<<( ostream& output, const Source& s ) 
{
  output << s.mName << (s.isMomentSource()? " moment":" force") << " source term" << endl;
   output << "   Location (X,Y,Z) = " << s.mX0 << "," << s.mY0 << "," << s.mZ0 << " in grid no " << s.m_grid << endl;
   output << "   Strength " << s.getAmplitude();
   output << "   t0 = " << s.mT0 << " freq = " << s.mFreq << endl;
   if( s.mIsMomentSource )
   {
      output << " Mxx Mxy Myy Mxz Myz Mzz = " << s.mForces[0] << " " << s.mForces[1] << " " << s.mForces[3] <<
	 " " << s.mForces[2] << " " << s.mForces[4] << " " << s.mForces[5] << endl;
   }
   else
   {
      output << " Fx Fy Fz = " << s.mForces[0] << " " << s.mForces[1] << " " << s.mForces[2] << endl;
   }
   return output;
}

//-----------------------------------------------------------------------
void Source::limit_frequency( int ppw, double minvsoh )
{
   double freqlim = minvsoh/(ppw);

   if( mTimeDependence == iBrune     || mTimeDependence == iBruneSmoothed || mTimeDependence == iDBrune ||
       mTimeDependence == iGaussian  || mTimeDependence == iErf || 
       mTimeDependence == iVerySmoothBump || mTimeDependence == iSmoothWave || 
       mTimeDependence == iLiu || mTimeDependence == iC6SmoothBump )
   {
      if( mFreq > 2*M_PI*freqlim )
	 mFreq = 2*M_PI*freqlim;
   }
   else
   {
      if( mFreq > freqlim )
	 mFreq = freqlim;
   }      
}

//-----------------------------------------------------------------------
double Source::compute_t0_increase(double t0_min) const
{
// Gaussian, GaussianInt=Erf, Ricker and RickerInt are all centered around mT0
  if( mTimeDependence == iGaussian  || mTimeDependence == iErf )
    return t0_min + 6.0/mFreq-mT0; // translating these by at least 6*sigma = 6/freq
  else if( mTimeDependence == iRicker  || mTimeDependence == iRickerInt ) 
    return t0_min + 1.9/mFreq-mT0; // 1.9 ?
  else
    return t0_min - mT0; // the rest of the time functions are zero for t<mT0
}

//-----------------------------------------------------------------------
void Source::adjust_t0( double dt0 )
{
   if( dt0 > 0 && !m_is_filtered )
      mT0 += dt0;
}

//-----------------------------------------------------------------------
double Source::dt_to_resolve( int ppw ) const
{
  double dt_resolved = 0;
  if( mTimeDependence == iBrune || mTimeDependence == iBruneSmoothed ||  mTimeDependence == iDBrune)
    {
      const double t95 = 4.744/mFreq;
      dt_resolved = t95/ppw;
    }
  else
    {

    }
  return dt_resolved;
}

//-----------------------------------------------------------------------
int Source::ppw_to_resolve( double dt ) const
{
  int ppw = 1;
  if( mTimeDependence == iBrune || mTimeDependence == iBruneSmoothed ||  mTimeDependence == iDBrune)
    {
      const double t95 = 4.744/mFreq;
      ppw = static_cast<int>(t95/dt);
    }
  else
    {

    }
  return ppw;
}


//-----------------------------------------------------------------------
void Source::set_derivative( int der )
{
   if( der >= 0 && der <= 10 )
      m_derivative = der;
}

//-----------------------------------------------------------------------
void Source::set_noderivative( )
{
   m_derivative = -1;
}

//-----------------------------------------------------------------------
void Source::set_dirderivative( double dir[11] )
{
   for( int i=0 ; i < 11 ; i++ )
      m_dir[i] = dir[i];
   m_derivative = 11;
}

//-----------------------------------------------------------------------
void Source::set_parameters( double x[11] )
{
   if( mIsMomentSource )
   {
      mX0 = x[0];
      mY0 = x[1];
      mZ0 = x[2];
      mForces[0] = x[3];
      mForces[1] = x[4];
      mForces[2] = x[5];
      mForces[3] = x[6];
      mForces[4] = x[7];
      mForces[5] = x[8];
      mT0 = x[9];
      mFreq = x[10];
   }
   else
      cout << "Error in Source::set_parameters(), " <<
             "function only implemented for moment sources" << endl;
}

//-----------------------------------------------------------------------
void Source::get_parameters( double x[11] ) const
{
   if( mIsMomentSource )
   {
      x[0] = mX0;
      x[1] = mY0;
      x[2] = mZ0;
      x[3] = mForces[0];
      x[4] = mForces[1];
      x[5] = mForces[2];
      x[6] = mForces[3];
      x[7] = mForces[4];
      x[8] = mForces[5];
      x[9] = mT0;
      x[10]= mFreq;
   }
   else
      cout << "Error in Source::get_parameters(), " <<
             "function only implemented for moment sources" << endl;
}

//-----------------------------------------------------------------------
void Source::setFrequency( double freq )
{
   mFreq=freq;
}

//-----------------------------------------------------------------------
void Source::correct_Z_level( )
{
// not yet functional
//    int i,j,k,g;
//    mEW->computeNearestGridPoint( i, j, k, g, mX0, mY0, mZ0 );
//    double q, r, s, minDepth=0; // this is a hack...
//    double lat, lon;
      
//    if( g == mEW->mNumberOfGrids-1 && mEW->topographyExists() )
//    {
// // Curvilinear
//      bool canBeInverted = mEW->invert_curvilinear_grid_mapping( mX0, mY0, mZ0, q, r, s );
//      if (canBeInverted) // the source location lives on this processor and needs to be corrected
//      {
//        double zMinTilde, zMin, zMax, zSource;
// // evaluate smoothed && raw topography
//        if (mEW->interpolate_topography(q, r, zMinTilde, true) && mEW->interpolate_topography(q, r, zMin, false))
//        {
// // if the topodepth source specification was used, we need to add the z-level of the raw topography to mZ0
// 	 if (m_zRelativeToTopography)
// 	 {
// 	   mZ0 += zMin;
// 	 }
	 
// 	 if (mZ0 > zMin + minDepth) // should be + minDepth because z is down
// 	 {
// 	   zMax = mEW->m_zmin[mEW->mNumberOfCartesianGrids-1];
// 	   zSource = zMinTilde + (mZ0 - zMin) * (zMax - zMinTilde)/(zMax - zMin);
// 	   if (mEW->getVerbosity()>=1)
// 	     printf("Correcting source z-level from %e (relative=%i) to %e. Raw topo z-level = %e, smoothed topo z-level = %e\n", 
// 		  mZ0, m_zRelativeToTopography, zSource, zMin, zMinTilde);
// 	   if (! mEW->invert_curvilinear_grid_mapping( mX0, mY0, zSource, q, r, s ))
// 	   {
// 	     mEW->computeGeographicCoord(mX0, mY0, lon, lat);
// 	     cerr << "The corrected source location zSource = " << zSource << " could not be inverted" << endl
// 		  << " Raw topography z = " << zMin << " Smoothed topography z = " << zMinTilde << endl
// 		  << " Source at lat = " << lat << " lon = " << lon << " mZ0 = " << mZ0 << " zSource = " << zSource << endl;
	     
// 	     mIgnore = true;
// 	   }
// 	   else if (s<0.)
// 	   {
// 	     mEW->computeGeographicCoord(mX0, mY0, lon, lat);
// 	     cerr << "Inverting the corrected source location zSource = " << zSource << " gives s = " << s << endl 
// 		  << " Raw topography z = " << zMin << " Smoothed topography z = " << zMinTilde << endl
// 		  << " Source at lat = " << lat << " lon = " << lon << " mZ0 = " << mZ0 << " zSource = " << zSource << endl;
// 	     mIgnore = true;
// 	   }
// 	   mZ0 = zSource;
// 	 }
// 	 else
// 	 {
// 	   mIgnore = true;
// 	   mEW->computeGeographicCoord(mX0, mY0, lon, lat);
// 	   cerr << "Ignoring source at lat = " << lat << " lon = " << lon << " z-level = " << mZ0 
// 		<< " Too close or above raw topography z = " << zMin << endl;
// 	 }
//        }
//        else
//        {
// 	 mIgnore = true;
// 	 mEW->computeGeographicCoord(mX0, mY0, lon, lat);
// 	 cerr << "Ignoring source at lat = " << lat << " lon = " << lon << " z-level = " << mZ0 << endl
// 	      << " Failed to evaluate smoothed and/or raw topography at q= " << q << ", r= " << r << endl;
//        }
//      }
// // ignoring sources not on this proc destroys the calculation of the moment magnitude, which is done from proc 0
// //      else // if the location can not be inverted, the source is not on this processor
// //      {
// //        mIgnore = true;
// //      }
     
//    } // if in top grid which is curvilinear
   
}


//-----------------------------------------------------------------------
void Source::getsourcewgh( double ai, double wgh[6], double dwghda[6], double ddwghda[6] ) const
{
   // Moments k=0,1,2,3,4 exact, two cont. derivatives wrt. position
   double p5 = ai*ai*ai*ai*ai*(5.0/3-7.0/24*ai -17/12.0*ai*ai+1.125*ai*ai*ai-0.25*ai*ai*ai*ai);
   wgh[0] = 1.0/24*(2*ai-ai*ai-2*ai*ai*ai-19*ai*ai*ai*ai) + p5;
   wgh[1] = 1.0/6*(-4*ai+4*ai*ai+ai*ai*ai)+4*ai*ai*ai*ai -5*p5;
   wgh[2] = 1-1.25*ai*ai-97.0/12*ai*ai*ai*ai + 10*p5;
   wgh[3] = 1.0/6*( 4*ai+4*ai*ai-ai*ai*ai+49*ai*ai*ai*ai)-10*p5;
   wgh[4] = 1.0/24*(-2*ai-ai*ai+2*ai*ai*ai)-4.125*ai*ai*ai*ai+5*p5;
   wgh[5] = 5.0/6*ai*ai*ai*ai - p5;

   // Derivatives of wgh wrt. ai:
   p5 = 5*ai*ai*ai*ai*(5.0/3-7.0/24*ai -17/12.0*ai*ai+1.125*ai*ai*ai-0.25*ai*ai*ai*ai) +
      ai*ai*ai*ai*ai*(-7.0/24 -17/6.0*ai+3*1.125*ai*ai-ai*ai*ai); 
   dwghda[0] = 1.0/24*(2-2*ai-6*ai*ai-19*4*ai*ai*ai) + p5;
   dwghda[1] = 1.0/6*(-4+8*ai+3*ai*ai)+ 16*ai*ai*ai - 5*p5;
   dwghda[2] = -2.5*ai-97.0/3*ai*ai*ai + 10*p5;
   dwghda[3] = 1.0/6*(4+8*ai-3*ai*ai+49*4*ai*ai*ai) - 10*p5;
   dwghda[4] = 1.0/24*(-2-2*ai+6*ai*ai)-4.125*4*ai*ai*ai + 5*p5;
   dwghda[5] = 20.0/6*ai*ai*ai - p5;

   // Second derivatives of wgh wrt. ai:
   p5 = ai*ai*ai*(100.0/3-8.75*ai-59.5*ai*ai+63*ai*ai*ai-18*ai*ai*ai*ai);

   ddwghda[0] = -1.0/12- 0.5*ai-9.5*ai*ai + p5;
   ddwghda[1] =  4.0/3 + ai     + 48*ai*ai - 5*p5;
   ddwghda[2] =  -2.5           - 97*ai*ai + 10*p5;
   ddwghda[3] =   4.0/3 - ai      + 98*ai*ai- 10*p5;
   ddwghda[4] =  -1.0/12 + 0.5*ai -49.5*ai*ai + 5*p5;
   ddwghda[5] =                    10*ai*ai - p5;
      
}

//-----------------------------------------------------------------------
void Source::getsourcedwgh( double ai, double wgh[6], double dwghda[6], double ddwghda[6] ) const
{
   // Moments k=0,1,2,3,4 exact, two cont. derivatives wrt. position
   double p5 = ai*ai*ai*ai*(-25.0/12-0.75*ai + 59.0/12*ai*ai - 4*ai*ai*ai + ai*ai*ai*ai);
   wgh[0] = 1.0/12*(-1+ai+3*ai*ai+8*ai*ai*ai) + p5;
   wgh[1] = 2.0/3*(1-2*ai) - 0.5*ai*ai-3.5*ai*ai*ai -5*p5;
   wgh[2] = 2.5*ai + 22.0/3*ai*ai*ai + 10*p5;
   wgh[3] = 2.0/3*(-1-2*ai)+0.5*ai*ai-23.0/3*ai*ai*ai-10*p5;
   wgh[4] = (1+ai)/12-0.25*ai*ai+4*ai*ai*ai + 5*p5;
   wgh[5] = -5.0/6*ai*ai*ai - p5;

   // Derivatives of wgh wrt. ai:
   p5 = 4*ai*ai*ai*(-25.0/12-0.75*ai + 59.0/12*ai*ai - 4*ai*ai*ai + ai*ai*ai*ai) +
      ai*ai*ai*ai*(-0.75 + 59.0/6*ai - 12*ai*ai + 4*ai*ai*ai);
   dwghda[0] = 1.0/12*(1+6*ai+24*ai*ai) + p5;
   dwghda[1] = 2.0/3*(-2) - ai-3*3.5*ai*ai -5*p5;
   dwghda[2] = 2.5 + 22.0*ai*ai + 10*p5;
   dwghda[3] = 2.0/3*(-2)+ai-23.0*ai*ai-10*p5;
   dwghda[4] = 1.0/12-0.5*ai+12*ai*ai + 5*p5;
   dwghda[5] = -5.0/2*ai*ai - p5;

   // Second derivatives of wgh wrt. ai:
   p5 = ai*ai*(-25-15*ai+147.5*ai*ai-168*ai*ai*ai+56*ai*ai*ai*ai);

   ddwghda[0] =  0.5 + 4*ai + p5;
   ddwghda[1] =  -1  -21*ai -5*p5;
   ddwghda[2] =       44*ai + 10*p5;
   ddwghda[3] =  1   -46*ai            -10*p5;
   ddwghda[4] =  -0.5 + 24*ai + 5*p5;
   ddwghda[5] =        -5*ai    - p5;

}

//-----------------------------------------------------------------------
void Source::getsourcewghlow( double ai, double wgh[6], double dwghda[6], double ddwghda[6] ) const
{
   // Lower component stencil, to use at lower boundaries
   // Moments k=0,1,2,3,4 exact, two cont. derivatives wrt. position

   wgh[0] = (2*ai-ai*ai-2*ai*ai*ai+ai*ai*ai*ai)/24;
   wgh[1] = (-4*ai+4*ai*ai+ai*ai*ai-ai*ai*ai*ai)/6;
   wgh[2] = 1-1.25*ai*ai+0.25*ai*ai*ai*ai;
   wgh[3] = (4*ai+4*ai*ai-ai*ai*ai-ai*ai*ai*ai)/6;
   wgh[4] = (-2*ai-ai*ai+2*ai*ai*ai+ai*ai*ai*ai)/24;
   wgh[5] = 0;

   // Derivatives of wgh wrt. ai:
   dwghda[0] = ( 1 - ai - 3*ai*ai+2*ai*ai*ai)/12;
   dwghda[1] = (-2+4*ai+1.5*ai*ai-2*ai*ai*ai)/3;
   dwghda[2] = -2.5*ai+ai*ai*ai;
   dwghda[3] = ( 2+4*ai-1.5*ai*ai-2*ai*ai*ai)/3;
   dwghda[4] = (-1 - ai + 3*ai*ai+2*ai*ai*ai)/12;
   dwghda[5] = 0;

   // Second derivatives of wgh wrt. ai:
   ddwghda[0] = -1.0/12 - 0.5*ai + 0.5*ai*ai;
   ddwghda[1] =  4.0/3 + ai - 2*ai*ai;
   ddwghda[2] = -2.5 + 3*ai*ai;
   ddwghda[3] =  4.0/3 - ai - 2*ai*ai;
   ddwghda[4] = -1.0/12 + 0.5*ai + 0.5*ai*ai;
   ddwghda[5] = 0;
}

//-----------------------------------------------------------------------
void Source::getsourcedwghlow( double ai, double wgh[6], double dwghda[6], double ddwghda[6] ) const
{
   // Lower component stencil, to use at lower boundaries, dirac derivative weights.
   // Moments k=0,1,2,3,4 exact, two cont. derivatives wrt. position

   // same as derivatives of dirac weights.
   wgh[0] = (-1 + ai + 3*ai*ai- 2*ai*ai*ai)/12;
   wgh[1] = ( 2 - 4*ai - 1.5*ai*ai + 2*ai*ai*ai)/3;
   wgh[2] = 2.5*ai-ai*ai*ai;
   wgh[3] = (-2 - 4*ai + 1.5*ai*ai + 2*ai*ai*ai)/3;
   wgh[4] = ( 1 + ai - 3*ai*ai- 2*ai*ai*ai)/12;
   wgh[5] = 0;

   // Derivatives of wgh wrt. ai:
   dwghda[0] =  1.0/12 + 0.5*ai - 0.5*ai*ai;
   dwghda[1] = -4.0/3 - ai + 2*ai*ai;
   dwghda[2] =   2.5 - 3*ai*ai;
   dwghda[3] = -4.0/3 + ai + 2*ai*ai;
   dwghda[4] =  1.0/12 - 0.5*ai - 0.5*ai*ai;
   dwghda[5] = 0;

   // Second derivatives of wgh wrt. ai:
   ddwghda[0] = 0.5 - ai;
   ddwghda[1] = -1 + 4*ai;
   ddwghda[2] =    -6*ai;
   ddwghda[3] =  1 + 4*ai;
   ddwghda[4] =  -0.5 - ai;
   ddwghda[5] = 0;
}

//-----------------------------------------------------------------------
void Source::getmetwgh( double ai, double wgh[8] ) const
{
   double pol = ai*ai*ai*ai*ai*ai*ai*(-251+135*ai+25*ai*ai-
                                      33*ai*ai*ai+6*ai*ai*ai*ai)/720;
   wgh[0] = -1.0/60*ai + 1.0/180*ai*ai + 1.0/48*ai*ai*ai + 23.0/144*ai*ai*ai*ai 
      - (17.0*ai + 223.0)*ai*ai*ai*ai*ai/720 - pol;
   wgh[1] = 3.0/20*ai -3.0/40*ai*ai -1.0/6*ai*ai*ai - 13.0/12*ai*ai*ai*ai + 
      97.0/45*ai*ai*ai*ai*ai + 1.0/6*ai*ai*ai*ai*ai*ai + 7*pol;
   wgh[2] = -0.75*ai +0.75*ai*ai+(13.0+155*ai)*ai*ai*ai/48 -103.0/16*ai*ai*ai*ai*ai
      - 121.0/240*ai*ai*ai*ai*ai*ai - 21*pol;
   wgh[3] = 1 - 49.0/36*ai*ai - 49.0/9*ai*ai*ai*ai+385.0/36*ai*ai*ai*ai*ai +
      61.0/72*ai*ai*ai*ai*ai*ai + 35*pol;
   wgh[4] = 0.75*ai + 0.75*ai*ai - 13.0/48*ai*ai*ai + 89.0/16*ai*ai*ai*ai - 
         1537.0/144*ai*ai*ai*ai*ai - 41.0/48*ai*ai*ai*ai*ai*ai - 35*pol;
   wgh[5] = -3.0/20*ai - 3.0/40*ai*ai + 1.0/6*ai*ai*ai - 41.0/12*ai*ai*ai*ai
      + 6.4*ai*ai*ai*ai*ai + 31.0/60*ai*ai*ai*ai*ai*ai + 21*pol;
   wgh[6] = 1.0/60*ai + 1.0/180*ai*ai - 1.0/48*ai*ai*ai + 167.0/144*ai*ai*ai*ai -
      1537.0/720*ai*ai*ai*ai*ai- 25.0/144*ai*ai*ai*ai*ai*ai - 7*pol;
   wgh[7] = -1.0/6*ai*ai*ai*ai + 11.0/36*ai*ai*ai*ai*ai + 1.0/40*ai*ai*ai*ai*ai*ai + pol;
}

//-----------------------------------------------------------------------
void Source::getmetdwgh( double ai, double wgh[8] ) const
{
   double pol = ai*ai*ai*ai*ai*ai*( -827 + 420*ai + 165*ai*ai - 180*ai*ai*ai
				    + 36*ai*ai*ai*ai)/720;

   wgh[0] = -1.0/60 + 1.0/90*ai + 1.0/16*ai*ai + 5.0/36*ai*ai*ai -
      55.0/144*ai*ai*ai*ai - 7.0/20*ai*ai*ai*ai*ai - pol;
   wgh[1] = 3.0/20*(1-ai) - 0.5*ai*ai - 5.0/6*ai*ai*ai + 47.0/18*ai*ai*ai*ai
      + 59.0/24*ai*ai*ai*ai*ai + 7*pol;
   wgh[2] = -0.75 + 1.5*ai + 13.0/16*ai*ai + 29.0/12*ai*ai*ai - 123.0/16*ai*ai*ai*ai -
      7.4*ai*ai*ai*ai*ai - 21*pol;
   wgh[3] = (-49.0*ai - 77.0*ai*ai*ai)/18 + 455.0/36*ai*ai*ai*ai + 99.0/8*ai*ai*ai*ai*ai 
           + 35*pol;
   wgh[4] = 0.75 + 1.5*ai - 13.0/16*ai*ai + 4.75*ai*ai*ai - 1805.0/144*ai*ai*ai*ai -
      149.0/12*ai*ai*ai*ai*ai - 35*pol;
   wgh[5] = -3.0/20*(1+ai) +0.5*ai*ai - 19.0/6*ai*ai*ai + 7.5*ai*ai*ai*ai +
      299.0/40*ai*ai*ai*ai*ai + 21*pol;
   wgh[6] = 1.0/60 + 1.0/90*ai - 1.0/16*ai*ai + 41.0/36*ai*ai*ai - 361.0/144*ai*ai*ai*ai -
      2.5*ai*ai*ai*ai*ai - 7*pol;
   wgh[7] = -1.0/6*ai*ai*ai + 13.0/36*ai*ai*ai*ai + 43.0/120*ai*ai*ai*ai*ai + pol;
}
//-----------------------------------------------------------------------
void Source::getmetwgh7( double ai, double wgh[7] ) const
{
   wgh[0] =ai*ai/180.0-ai*ai*ai*ai/144.0+ai*ai*ai*ai*ai*ai/720.0-ai/60.0+
	    ai*ai*ai/48.0-ai*ai*ai*ai*ai/240.0;
   wgh[1] =-3.0/40.0*ai*ai+ai*ai*ai*ai/12.0-ai*ai*ai*ai*ai*ai/120.0+
	    3.0/20.0*ai-ai*ai*ai/6.0+ai*ai*ai*ai*ai/60.0;
   wgh[2] =  3.0/4.0*ai*ai-13.0/48.0*ai*ai*ai*ai+ai*ai*ai*ai*ai*ai/48.0-
	    3.0/4.0*ai+13.0/48.0*ai*ai*ai-ai*ai*ai*ai*ai/48.0;
   wgh[3] =1.0-49.0/36.0*ai*ai+7.0/18.0*ai*ai*ai*ai-ai*ai*ai*ai*ai*ai/36.0;
   wgh[4] =3.0/4.0*ai*ai-13.0/48.0*ai*ai*ai*ai+3.0/4.0*ai-13.0/48.0*ai*ai*ai+
	    ai*ai*ai*ai*ai/48.0+ai*ai*ai*ai*ai*ai/48.0;
   wgh[5] =-3.0/40.0*ai*ai+ai*ai*ai*ai/12.0-3.0/20.0*ai+ai*ai*ai/6.0-
	    ai*ai*ai*ai*ai/60.0-ai*ai*ai*ai*ai*ai/120.0;
   wgh[6] = ai*ai/180.0-ai*ai*ai*ai/144.0+ai*ai*ai*ai*ai*ai/720.0+ai/60.0+
	    ai*ai*ai*ai*ai/240.0-ai*ai*ai/48.0;

}

//-----------------------------------------------------------------------
void Source::getmetdwgh7( double ai, double wgh[7] ) const
{
   wgh[0] =  -1.0/60.0+ai*ai/16.0-ai*ai*ai*ai/48.0+ai/90.0-
		 ai*ai*ai/36.0+ai*ai*ai*ai*ai/120.0;
   wgh[1] =3.0/20.0-ai*ai/2.0+ai*ai*ai*ai/12.0-3.0/20.0*ai+
		 ai*ai*ai/3.0-ai*ai*ai*ai*ai/20.0;
   wgh[2] =-3.0/4.0+13.0/16.0*ai*ai-5.0/48.0*ai*ai*ai*ai+
		 3.0/2.0*ai-13.0/12.0*ai*ai*ai+ai*ai*ai*ai*ai/8.0;
   wgh[3] =-49.0/18.0*ai+14.0/9.0*ai*ai*ai-ai*ai*ai*ai*ai/6.0;
   wgh[4] =3.0/4.0-13.0/16.0*ai*ai+3.0/2.0*ai-13.0/12.0*ai*ai*ai+
 		 5.0/48.0*ai*ai*ai*ai+ai*ai*ai*ai*ai/8.0;
   wgh[5] =-3.0/20.0+ai*ai/2.0-ai*ai*ai*ai/12.0-3.0/20.0*ai+
		 ai*ai*ai/3.0-ai*ai*ai*ai*ai/20.0;
   wgh[6] =1.0/60.0-ai*ai/16.0+ai*ai*ai*ai/48.0+ai/90.0-
		 ai*ai*ai/36.0+ai*ai*ai*ai*ai/120.0;
}

//-----------------------------------------------------------------------
void Source::set_grid_point_sources4( EW *a_EW, vector<GridPointSource*>& point_sources ) const
{
// note that this routine is called from all processors, for each input source 
   int i,j,k,g;
   a_EW->computeNearestGridPoint( i, j, k, g, mX0, mY0, mZ0 );
   double q, r, s;
   double h = a_EW->mGridSize[g];
   bool canBeInverted, curvilinear;
   double normwgh[4]={17.0/48.0, 59.0/48.0, 43.0/48.0, 49.0/48.0 };

   if( g == a_EW->mNumberOfGrids-1 && a_EW->topographyExists() )
   {
// Curvilinear
// Problem when the curvilinear mapping is NOT analytic:
// This routine can only compute the 's' coordinate if (mX0, mY0) is owned by this processor
      canBeInverted = a_EW->invert_curvilinear_grid_mapping( mX0, mY0, mZ0, q, r, s );
// if s < 0, the source is located above the grid and the call to
// find_curvilinear_derivatives_at_point will fail
      if (canBeInverted && s<0.)
      {
	 double xTop, yTop, zTop;
	 a_EW->curvilinear_grid_mapping(q, r, 0., xTop, yTop, zTop);
	 double lat, lon;
	 a_EW->computeGeographicCoord(mX0, mY0, lon, lat);
	 printf("Found a source above the curvilinear grid! Lat=%e, Lon=%e, source Z-level = %e, grid boundary Z = %e\n",
		lat, lon, mZ0, zTop);
	 MPI_Abort(MPI_COMM_WORLD,1);
      }
      curvilinear   = true;
   }
   else
   {
// Cartesian case
      q = mX0/h+1;
      r = mY0/h+1;
      s = (mZ0-a_EW->m_zmin[g])/h+1;
      canBeInverted = true;
      curvilinear   = false;
   }

   int Ni = a_EW->m_global_nx[g];
   int Nj = a_EW->m_global_ny[g];
   int Nz = a_EW->m_global_nz[g];

   int ic, jc, kc;
   bool upperbndry, lowerbndry, ccbndry;
   double ai, bi, ci;

// Delta distribution
   double wghi[6], wghj[6], wghk[6], wghix[6], wghjy[6], wghkz[6];
   double wghixx[6], wghjyy[6], wghkzz[6];
// Delta' distribution
   double dwghi[6], dwghj[6], dwghk[6], dwghix[6], dwghjy[6], dwghkz[6];
   double dwghixx[6], dwghjyy[6], dwghkzz[6];


   if( canBeInverted )
   {
 // Compute source location and weights in source discretization
      ic = static_cast<int>(floor(q));
      jc = static_cast<int>(floor(r));
      kc = static_cast<int>(floor(s));

// Bias stencil away from boundary, no source at ghost/padding points
      if( ic <= 2 )    ic = 3;
      if( ic >= Ni-2 ) ic = Ni-3;
      if( jc <= 2 )    jc = 3;
      if( jc >= Nj-2 ) jc = Nj-3;

// Six point stencil, with points, kc-2,..kc+3, Interior in domain
// if kc-2>=1, kc+3 <= Nz --> kc >= 3 and kc <= Nz-3
// Can evaluate with two ghost points if kc-2>=-1 and kc+3 <= Nz+2
//  --->  kc >=1 and kc <= Nz-1
//
      if( kc >= Nz )
	 kc = Nz-1;
      if( kc < 1 )
	 kc = 1;

// upper(surface) and lower boundaries , when the six point stencil kc-2,..kc+3
// make use of the first (k=1) or the last (k=Nz) interior point.
      upperbndry = (kc == 1    || kc == 2    || kc == 3  );
      lowerbndry = (kc == Nz-1 || kc == Nz-2 || kc == Nz-3 );

// ccbndry=true if at the interface between the curvilinear grid and the cartesian grid. 
// Defined as the six point stencil uses values from both grids.
      ccbndry = a_EW->topographyExists() &&  ( (upperbndry && g == a_EW->mNumberOfGrids-2) ||
					       (lowerbndry && g == a_EW->mNumberOfGrids-1)    );
// If not at the interface between curvilinear and Cartesian grids, bias stencil away 
// from the boundary.
      if( !ccbndry )
      {
	 if( kc <= 2 )    kc = 3;
	 if( kc >= Nz-2 ) kc = Nz-3;
      }
      ai=q-ic, bi=r-jc, ci=s-kc;

   // Delta distribution
      getsourcewgh( ai, wghi, wghix, wghixx );
      getsourcewgh( bi, wghj, wghjy, wghjyy );
      getsourcewgh( ci, wghk, wghkz, wghkzz );

// Delta' distribution
      getsourcedwgh( ai, dwghi, dwghix, dwghixx );
      getsourcedwgh( bi, dwghj, dwghjy, dwghjyy );
      getsourcedwgh( ci, dwghk, dwghkz, dwghkzz );

   // Special boundary stencil at free surface
      if( !ccbndry && (kc == 3 && ci <= 0) )
      {
	 getsourcewghlow( ci, wghk, wghkz, wghkzz );
	 getsourcedwghlow( ci, dwghk, dwghkz, dwghkzz );
      }

// Boundary correction, at upper boundary, but only if SBP operators are used there
      if( (g == a_EW->mNumberOfGrids-1) && a_EW->is_onesided(g,4)  )
      {
	 for( int k=0 ; k <= 5 ; k++ )
	 {
	    if( ( 1 <= k+kc-2) && ( k+kc-2 <= 4 ) )
	    {
	       wghk[k]    /= normwgh[k+kc-3];
	       dwghk[k]   /= normwgh[k+kc-3];
	       wghkz[k]   /= normwgh[k+kc-3];
	       dwghkz[k]  /= normwgh[k+kc-3];
	       wghkzz[k]  /= normwgh[k+kc-3];
	       dwghkzz[k] /= normwgh[k+kc-3];
	    }
	 }
      }
   }
   //   int myid;
   //   MPI_Comm_rank(MPI_COMM_WORLD, &myid );
   //   cout << myid << " SOURCE at " << ic << " " << jc << " "  << kc ;
   //   if( canBeInverted )
   //      cout << " can be inverted";
   //   else
   //      cout << " can not be inverted";

   // Point source. NOTE: Derivatives needed for source inversion not implemented for this case.
   if( !mIsMomentSource && canBeInverted )
   {
      for( int k=kc-2 ; k <= kc+3 ; k++ )
	 for( int j=jc-2 ; j <= jc+3 ; j++ )
	    for( int i=ic-2 ; i <= ic+3 ; i++ )
	    {
	       double wF = wghi[i-ic+2]*wghj[j-jc+2]*wghk[k-kc+2];
	       if( (wF != 0) && (mForces[0] != 0 || mForces[1] != 0 || mForces[2] != 0) 
		   && a_EW->interior_point_in_proc(i,j,g) ) // checks if (i,j) belongs to this processor
	       {
		  if( curvilinear )
		     wF /= a_EW->mJ(i,j,k);
		  else
		     wF /= h*h*h;

		  if( 1 <= k && k <= Nz )
		  {
		     GridPointSource* sourcePtr = new GridPointSource( mFreq, mT0,
								      i,j,k,g,
								      wF*mForces[0], wF*mForces[1], wF*mForces[2],
								      mTimeDependence, mNcyc, 
								      mPar, mNpar, mIpar, mNipar );
		     point_sources.push_back(sourcePtr);
		  }
		  if( k <= 1 && ccbndry && upperbndry )
		  {
		     int Nzp =a_EW->m_global_nz[g+1];
		     int kk = Nzp - 1 + k;

		     GridPointSource* sourcePtr = new GridPointSource( mFreq, mT0,
								      i,j,kk,g+1,
								      wF*mForces[0], wF*mForces[1], wF*mForces[2],
								      mTimeDependence, mNcyc, 
								      mPar, mNpar, mIpar, mNipar );
		     point_sources.push_back(sourcePtr);
		  }
		  if( k >= Nz && ccbndry && lowerbndry )
		  {
		     int kk = k-Nz + 1;
		     GridPointSource* sourcePtr = new GridPointSource( mFreq, mT0,
									  i,j,kk,g-1,
									  wF*mForces[0], wF*mForces[1], wF*mForces[2],
									  mTimeDependence, mNcyc, 
									  mPar, mNpar, mIpar, mNipar );
		     point_sources.push_back(sourcePtr);
		  }
	       }
	    }
   } 
   // Moment source.
   else if( mIsMomentSource )
   {
      double qX0[3], rX0[3], sX0[3];
      if( !curvilinear )
      {
	 // Cartesian case, constant metric
	 qX0[0] = 1/h;qX0[1]=0;  qX0[2]=0;
	 rX0[0] = 0;  rX0[1]=1/h;rX0[2]=0;
	 sX0[0] = 0;  sX0[1]=0;  sX0[2]=1/h;
      }	 
      else
      {
	 // Compute the curvilinear metric in the processor that owns the source.
	 //   (ic, jc are undefined if canBeInverted is false.)
	 double zdertmp[3]={0,0,0}, zq, zr, zs;
	 if( a_EW->interior_point_in_proc(ic,jc,g) && canBeInverted )
	 {
	    compute_metric_at_source( a_EW, q, r, s, ic, jc, kc, g, zq, zr, zs );
	    zdertmp[0] = zq;
	    zdertmp[1] = zr;
	    zdertmp[2] = zs;
	 }

	 //	 // Broadcast the computed metric to all processors. 
	 //	 // First find out the ID of the processor that computed the metric...
	 //	 int owner = -1;
	 //	 if( a_EW->interior_point_in_proc( ic, jc, g ) && canBeInverted )
	 //            MPI_Comm_rank(MPI_COMM_WORLD, &owner );
	 //	 int owntmp = owner;
	 //	 MPI_Allreduce( &owntmp, &owner, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
	 //	 // ...then broadcast the derivatives
	 //	 MPI_Bcast( zder, 3, MPI_DOUBLE, owner, MPI_COMM_WORLD );

         // Simpler solution
         double zder[3];
	 MPI_Allreduce( zdertmp, zder, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	 zq = zder[0];
	 zr = zder[1];
	 zs = zder[2];
	 
	 //	 cout << "zq = " << zq << " zr = " << zr << " zs = " << zs << endl;
	 //	 if( owntmp != -1 )
	 //	    cout << " *" ;
	 //	 cout << endl;
	 qX0[0] = 1/h;
	 qX0[1] = 0;
	 qX0[2] = 0;
	 rX0[0] = 0;
	 rX0[1] = 1/h;
	 rX0[2] = 0;
	 sX0[0] = -zq/(h*zs);
	 sX0[1] = -zr/(h*zs);
	 sX0[2] =   1/zs;
      }

   //         cout << "Recomputed metric at pt " << qX0[0] << " " << sX0[0] << " " << sX0[1] << " " << sX0[2] << endl;

	    // Gradients of sX0[0]=sX, sX0[1]=sY, and sX0[2]=sZ wrt. (q,r,s)
	    // NYI
      //            double dsX0[3], dsY0[3], dsZ0[3], d2sX0[6], d2sY0[6], d2sZ0[6];
      //	    dsX0[0] = 0;
      //	    dsX0[1] = 0;
      //	    dsX0[2] = 0;

      //	    dsY0[0] = 0; 
      //	    dsY0[1] = 0;
      //	    dsY0[2] = 0;

      //	    dsZ0[0] = 0;
      //	    dsZ0[1] = 0;
      //	    dsZ0[2] = 0;
	    // Hessians of sX0[0]=sX, sX0[1]=sY, and sX0[2]=sZ wrt. (q,r,s), in order qq,qr,qs,rr,rs,ss
      //	    d2sX0[0] = 0;
      //	    d2sX0[1] =0;
      //	    d2sX0[2] =0;
      //	    d2sX0[3] =0;
      //	    d2sX0[4] =0;
      //	    d2sX0[5] =0;

      //	    d2sY0[0] =0;
      //	    d2sY0[1] =0;
      //	    d2sY0[2] =0;
      //	    d2sY0[3] =0;
      //	    d2sY0[4] =0;
      //	    d2sY0[5] =0;

      //	    d2sZ0[0] =0;
      //	    d2sZ0[1] =0;
      //	    d2sZ0[2] =0;
      //	    d2sZ0[3] =0;
      //	    d2sZ0[4] =0;
      //	    d2sZ0[5] =0;

      if( canBeInverted )
      {
	 for( int k=kc-2 ; k <= kc+3 ; k++ )
	    for( int j=jc-2 ; j <= jc+3 ; j++ )
	       for( int i=ic-2 ; i <= ic+3 ; i++ )
	       {
		  double wFx=0, wFy=0, wFz=0, dsdp[27];
		  if( a_EW->interior_point_in_proc(i,j,g) ) 
		  {
		     wFx += qX0[0]*dwghi[i-ic+2]* wghj[j-jc+2]* wghk[k-kc+2];
		  //		  wFy += qX0[1]*dwghi[i-ic+2]* wghj[j-jc+2]* wghk[k-kc+2]; 
		  //		  wFz += qX0[2]*dwghi[i-ic+2]* wghj[j-jc+2]* wghk[k-kc+2];

		  //		  wFx +=  wghi[i-ic+2]*rX0[0]*dwghj[j-jc+2]* wghk[k-kc+2];
		     wFy +=  wghi[i-ic+2]*rX0[1]*dwghj[j-jc+2]* wghk[k-kc+2];
		  //		  wFz +=  wghi[i-ic+2]*rX0[2]*dwghj[j-jc+2]* wghk[k-kc+2];

		     wFx +=  wghi[i-ic+2]* wghj[j-jc+2]*sX0[0]*dwghk[k-kc+2];
		     wFy +=  wghi[i-ic+2]* wghj[j-jc+2]*sX0[1]*dwghk[k-kc+2];
		     wFz +=  wghi[i-ic+2]* wghj[j-jc+2]*sX0[2]*dwghk[k-kc+2];
		  // NOTE:  Source derivatives wrt. (x0,y0,z0) currently not implemented 
		  // for curvilinear grids.
		     double hi=1.0/h;
		     double hi2=hi*hi;
		     double wFxdx0 =dwghix[i-ic+2]*  wghj[j-jc+2]*  wghk[k-kc+2]*hi*qX0[0];
		     double wFxdy0 = dwghi[i-ic+2]* wghjy[j-jc+2]*  wghk[k-kc+2]*hi*rX0[1];
		     double wFxdz0 = dwghi[i-ic+2]*  wghj[j-jc+2]* wghkz[k-kc+2]*hi*sX0[2];
		     double wFydx0 = wghix[i-ic+2]* dwghj[j-jc+2]*  wghk[k-kc+2]*hi*qX0[0];
		     double wFydy0 =  wghi[i-ic+2]*dwghjy[j-jc+2]*  wghk[k-kc+2]*hi*rX0[1];
		     double wFydz0 =  wghi[i-ic+2]* dwghj[j-jc+2]* wghkz[k-kc+2]*hi*sX0[2];
		     double wFzdx0 = wghix[i-ic+2]*  wghj[j-jc+2]* dwghk[k-kc+2]*hi*qX0[0];
		     double wFzdy0 =  wghi[i-ic+2]* wghjy[j-jc+2]* dwghk[k-kc+2]*hi*rX0[1];
		     double wFzdz0 =  wghi[i-ic+2]*  wghj[j-jc+2]*dwghkz[k-kc+2]*hi*sX0[2];
		  //                  if( curvilinear && kc <= Nz-3 )
		  //		  {
		     // dsx0dx0[i+3*j] = d(sX0[j])/dxi
		  //                     wFxdx0 += wghi[i-ic+2]*wghj[j-jc+2]*dwghk[k-kc+2]*dsX0dX0[0];
		  //                     wFxdy0 += wghi[i-ic+2]*wghj[j-jc+2]*dwghk[k-kc+2]*dsX0dX0[1];
		  //		     wFxdz0 += wghi[i-ic+2]*wghj[j-jc+2]*dwghk[k-kc+2]*dxX0dX0[2];
		  //                     wFydx0 += wghi[i-ic+2]*wghj[j-jc+2]*dwghk[k-kc+2]*dsX0dX0[3];
		  //                     wFydy0 += wghi[i-ic+2]*wghj[j-jc+2]*dwghk[k-kc+2]*dsX0dX0[4];
		  //		     wFydz0 += wghi[i-ic+2]*wghj[j-jc+2]*dwghk[k-kc+2]*dxX0dX0[5];
		  //                     wFzdx0 += wghi[i-ic+2]*wghj[j-jc+2]*dwghk[k-kc+2]*dsX0dX0[6];
		  //                     wFzdy0 += wghi[i-ic+2]*wghj[j-jc+2]*dwghk[k-kc+2]*dsX0dX0[7];
		  //		     wFzdz0 += wghi[i-ic+2]*wghj[j-jc+2]*dwghk[k-kc+2]*dxX0dX0[8];
		  //		  }

		  // Second derivatives

		     double wFxdx0dx0 = dwghixx[i-ic+2]*  wghj[j-jc+2]*  wghk[k-kc+2]*hi2*hi;
		     double wFxdx0dy0 = dwghix[i-ic+2]*  wghjy[j-jc+2]*  wghk[k-kc+2]*hi2*hi;
		     double wFxdx0dz0 = dwghix[i-ic+2]*  wghj[j-jc+2]*  wghkz[k-kc+2]*hi2*hi;
		     double wFxdy0dy0 = dwghi[i-ic+2]* wghjyy[j-jc+2]*  wghk[k-kc+2]*hi2*hi;
		     double wFxdy0dz0 = dwghi[i-ic+2]* wghjy[j-jc+2]*  wghkz[k-kc+2]*hi2*hi;
		     double wFxdz0dz0 = dwghi[i-ic+2]*  wghj[j-jc+2]* wghkzz[k-kc+2]*hi2*hi;

		     double wFydx0dx0 = wghixx[i-ic+2]* dwghj[j-jc+2]*  wghk[k-kc+2]*hi2*hi;
		     double wFydx0dy0 = wghix[i-ic+2]* dwghjy[j-jc+2]*  wghk[k-kc+2]*hi2*hi;
		     double wFydx0dz0 = wghix[i-ic+2]* dwghj[j-jc+2]*  wghkz[k-kc+2]*hi2*hi;
		     double wFydy0dy0 = wghi[i-ic+2]*dwghjyy[j-jc+2]*  wghk[k-kc+2]*hi2*hi;
		     double wFydy0dz0 = wghi[i-ic+2]*dwghjy[j-jc+2]*  wghkz[k-kc+2]*hi2*hi;
		     double wFydz0dz0 = wghi[i-ic+2]* dwghj[j-jc+2]* wghkzz[k-kc+2]*hi2*hi;

		     double wFzdx0dx0 = wghixx[i-ic+2]*  wghj[j-jc+2]* dwghk[k-kc+2]*hi2*hi;
		     double wFzdx0dy0 = wghix[i-ic+2]*  wghjy[j-jc+2]* dwghk[k-kc+2]*hi2*hi;
		     double wFzdx0dz0 = wghix[i-ic+2]*  wghj[j-jc+2]* dwghkz[k-kc+2]*hi2*hi;
		     double wFzdy0dy0 = wghi[i-ic+2]* wghjyy[j-jc+2]* dwghk[k-kc+2]*hi2*hi;
		     double wFzdy0dz0 = wghi[i-ic+2]* wghjy[j-jc+2]* dwghkz[k-kc+2]*hi2*hi;
		     double wFzdz0dz0 = wghi[i-ic+2]*  wghj[j-jc+2]*dwghkzz[k-kc+2]*hi2*hi;
		     //                  if( curvilinear && kc <= Nz-3 )
		     //		  {

		     //		  }
               
		     double jaci;
		     if( curvilinear )
			jaci = 1/a_EW->mJ(i,j,k);
		     else
			jaci = 1.0/(h*h*h);

		     double fx = -(mForces[0]*wFx+mForces[1]*wFy+mForces[2]*wFz)*jaci;
		     double fy = -(mForces[1]*wFx+mForces[3]*wFy+mForces[4]*wFz)*jaci;
		     double fz = -(mForces[2]*wFx+mForces[4]*wFy+mForces[5]*wFz)*jaci;

		     // Derivatives with respect to (x0,y0,z0,mxx,mxy,mxz,myy,myz,mzz)
		     dsdp[0] = -(mForces[0]*wFxdx0+mForces[1]*wFydx0 + mForces[2]*wFzdx0)*jaci;
		     dsdp[1] = -(mForces[1]*wFxdx0+mForces[3]*wFydx0 + mForces[4]*wFzdx0)*jaci;
		     dsdp[2] = -(mForces[2]*wFxdx0+mForces[4]*wFydx0 + mForces[5]*wFzdx0)*jaci;
		     dsdp[3] = -(mForces[0]*wFxdy0+mForces[1]*wFydy0 + mForces[2]*wFzdy0)*jaci;
		     dsdp[4] = -(mForces[1]*wFxdy0+mForces[3]*wFydy0 + mForces[4]*wFzdy0)*jaci;
		     dsdp[5] = -(mForces[2]*wFxdy0+mForces[4]*wFydy0 + mForces[5]*wFzdy0)*jaci;
		     dsdp[6] = -(mForces[0]*wFxdz0+mForces[1]*wFydz0 + mForces[2]*wFzdz0)*jaci;
		     dsdp[7] = -(mForces[1]*wFxdz0+mForces[3]*wFydz0 + mForces[4]*wFzdz0)*jaci;
		     dsdp[8] = -(mForces[2]*wFxdz0+mForces[4]*wFydz0 + mForces[5]*wFzdz0)*jaci;
		     dsdp[9]  = -wFx*jaci;
		     dsdp[10] =  0;
		     dsdp[11] =  0;
		     dsdp[12]  =-wFy*jaci;
		     dsdp[13] = -wFx*jaci;
		     dsdp[14] =  0;
		     dsdp[15] = -wFz*jaci;
		     dsdp[16] =  0;
		     dsdp[17] = -wFx*jaci;
		     dsdp[18] =  0;
		     dsdp[19] = -wFy*jaci;
		     dsdp[20] =  0;
		     dsdp[21] =  0;
		     dsdp[22] = -wFz*jaci;
		     dsdp[23] = -wFy*jaci;
		     dsdp[24] =  0;
		     dsdp[25] =  0;
		     dsdp[26] = -wFz*jaci;

		     // Matrices needed for computing the Hessian wrt (x0,y0,z0,mxx,mxy,mxz,myy,myz,mzz)
		     double dddp[9], dh1[9], dh2[9], dh3[9];
		     dddp[0]  =-wFxdx0*jaci;
		     dddp[1]  =-wFxdy0*jaci;
		     dddp[2]  =-wFxdz0*jaci;
		     dddp[3]  =-wFydx0*jaci;
		     dddp[4]  =-wFydy0*jaci;
		     dddp[5]  =-wFydz0*jaci;
		     dddp[6]  =-wFzdx0*jaci;
		     dddp[7]  =-wFzdy0*jaci;
		     dddp[8]  =-wFzdz0*jaci;


		     // derivative of (dsdp[0],dsdp[3],dsdp[6]) (first component)
		     dh1[0] = -(mForces[0]*wFxdx0dx0 + mForces[1]*wFydx0dx0+mForces[2]*wFzdx0dx0)*jaci;
		     dh1[1] = -(mForces[0]*wFxdx0dy0 + mForces[1]*wFydx0dy0+mForces[2]*wFzdx0dy0)*jaci;
		     dh1[2] = -(mForces[0]*wFxdx0dz0 + mForces[1]*wFydx0dz0+mForces[2]*wFzdx0dz0)*jaci;

		     dh1[3] = dh1[1];
		     dh1[4] = -(mForces[0]*wFxdy0dy0 + mForces[1]*wFydy0dy0+mForces[2]*wFzdy0dy0)*jaci;
		     dh1[5] = -(mForces[0]*wFxdy0dz0 + mForces[1]*wFydy0dz0+mForces[2]*wFzdy0dz0)*jaci;

		     dh1[6] = dh1[2];
		     dh1[7] = dh1[5];
		     dh1[8] = -(mForces[0]*wFxdz0dz0 + mForces[1]*wFydz0dz0+mForces[2]*wFzdz0dz0)*jaci;

		     // derivative of (dsdp[1],dsdp[4],dsdp[7]) (second component)
		     dh2[0] = -(mForces[1]*wFxdx0dx0 + mForces[3]*wFydx0dx0+mForces[4]*wFzdx0dx0)*jaci;
		     dh2[1] = -(mForces[1]*wFxdx0dy0 + mForces[3]*wFydx0dy0+mForces[4]*wFzdx0dy0)*jaci;
		     dh2[2] = -(mForces[1]*wFxdx0dz0 + mForces[3]*wFydx0dz0+mForces[4]*wFzdx0dz0)*jaci;

		     dh2[3] = dh2[1];
		     dh2[4] = -(mForces[1]*wFxdy0dy0 + mForces[3]*wFydy0dy0+mForces[4]*wFzdy0dy0)*jaci;
		     dh2[5] = -(mForces[1]*wFxdy0dz0 + mForces[3]*wFydy0dz0+mForces[4]*wFzdy0dz0)*jaci;

		     dh2[6] = dh2[2];
		     dh2[7] = dh2[5];
		     dh2[8] = -(mForces[1]*wFxdz0dz0 + mForces[3]*wFydz0dz0+mForces[4]*wFzdz0dz0)*jaci;

		     // derivative of (dsdp[2],dsdp[5],dsdp[8]) (third component)
		     dh3[0] = -(mForces[2]*wFxdx0dx0 + mForces[4]*wFydx0dx0+mForces[5]*wFzdx0dx0)*jaci;
		     dh3[1] = -(mForces[2]*wFxdx0dy0 + mForces[4]*wFydx0dy0+mForces[5]*wFzdx0dy0)*jaci;
		     dh3[2] = -(mForces[2]*wFxdx0dz0 + mForces[4]*wFydx0dz0+mForces[5]*wFzdx0dz0)*jaci;

		     dh3[3] = dh3[1];
		     dh3[4] = -(mForces[2]*wFxdy0dy0 + mForces[4]*wFydy0dy0+mForces[5]*wFzdy0dy0)*jaci;
		     dh3[5] = -(mForces[2]*wFxdy0dz0 + mForces[4]*wFydy0dz0+mForces[5]*wFzdy0dz0)*jaci;

		     dh3[6] = dh3[2];
		     dh3[7] = dh3[5];
		     dh3[8] = -(mForces[2]*wFxdz0dz0 + mForces[4]*wFydz0dz0+mForces[5]*wFzdz0dz0)*jaci;

		     //                  if( i==42 && j==55 && k==39 )
		     //		  {
		     //		     cout.precision(16);
		     //                  cout << "-----------------------------------------------------------------------\n";
		     //		  cout << "     " << i <<  " " << j << " " << k << endl;
		     //                  cout << "dsp = " << dsdp[2] << " " << dsdp[5] << "  " << dsdp[8] << endl;
		     //                  cout << "dh  = " << dh3[0] << " " << dh3[1] << "  " << dh3[2] << endl;
		     //                  cout << "      " << dh3[3] << " " << dh3[4] << "  " << dh3[5] << endl;
		     //                  cout << "      " << dh3[6] << " " << dh3[7] << "  " << dh3[8] << endl;
		     //                  cout << "-----------------------------------------------------------------------" << endl;
		     //		  }

		     //		  if( mAmp != 0 && (fx != 0 || fy != 0 || fz != 0) )
		     if( 1 <= k && k <= Nz )
		     {
			GridPointSource* sourcePtr = new GridPointSource( mFreq, mT0, i, j, k, g, 
									  fx, fy, fz, mTimeDependence, mNcyc,
									  mPar, mNpar, mIpar, mNipar,
									  dsdp, dddp, dh1, dh2, dh3 );
			if( m_derivative >= 0 )
			   sourcePtr->set_derivative(m_derivative,m_dir);
			point_sources.push_back(sourcePtr);
		     }
		     if( k <= 1 && ccbndry && upperbndry )
		     {
			int Nzp =a_EW->m_global_nz[g+1];
			int kk = Nzp - 1 + k;
			GridPointSource* sourcePtr = new GridPointSource( mFreq, mT0, i, j, kk, g+1, 
									  fx, fy, fz, mTimeDependence, mNcyc,
									  mPar, mNpar, mIpar, mNipar,
									  dsdp, dddp, dh1, dh2, dh3 );
			if( m_derivative >= 0 )
			   sourcePtr->set_derivative(m_derivative,m_dir);
			point_sources.push_back(sourcePtr);
		     }
		     if( k >= Nz && ccbndry && lowerbndry )
		     {
			int kk = k-Nz + 1;
			GridPointSource* sourcePtr = new GridPointSource( mFreq, mT0, i, j, kk, g-1, 
									  fx, fy, fz, mTimeDependence, mNcyc,
									  mPar, mNpar, mIpar, mNipar,
									  dsdp, dddp, dh1, dh2, dh3 );
			if( m_derivative >= 0 )
			   sourcePtr->set_derivative(m_derivative,m_dir);
			point_sources.push_back(sourcePtr);
		     }
		  }
	    }
      }
   }
}

//-----------------------------------------------------------------------
void Source::exact_testmoments( int kx[3], int ky[3], int kz[3], double momex[3] )
{
   // Integrals over the domain of a polynomial of degree (kx,ky,kz) times the source
   if( !mIsMomentSource )
   {
      double x1, y1, z1;
      for( int c = 0; c < 3 ; c++ )
      {
	 if( kx[c] == 0 )
	    x1 = 1;
	 else
	    x1 = pow(mX0,kx[c]);
         if( ky[c] == 0 )
	    y1 = 1;
	 else
	    y1 = pow(mY0,ky[c]);
         if( kz[c] == 0 )
	    z1 = 1;
	 else
	    z1 = pow(mZ0,kz[c]);
	 momex[c] = mForces[c]*x1*y1*z1;
      }
   }
   else
   {
      double x1, y1, z1, xp1, yp1, zp1;
      for( int c = 0; c < 3 ; c++ )
      {
	 if( kx[c] == 0 )
	    x1 = 1;
	 else
	    x1 = pow(mX0,kx[c]);
	 if( kx[c] == 0 )
	    xp1 = 0;
	 else if( kx[c] == 1 )
            xp1 = -1;
	 else
	    xp1 =-kx[c]*pow(mX0,(kx[c]-1));

         if( ky[c] == 0 )
	    y1 = 1;
	 else
	    y1 = pow(mY0,ky[c]);
	 if( ky[c] == 0 )
	    yp1 = 0;
	 else if( ky[c] == 1 )
            yp1 = -1;
	 else
	    yp1 =-ky[c]*pow(mY0,(ky[c]-1));

         if( kz[c] == 0 )
	    z1 = 1;
	 else
	    z1 = pow(mZ0,kz[c]);
	 if( kz[c] == 0 )
	    zp1 = 0;
	 else if( kz[c] == 1 )
            zp1 = -1;
	 else
	    zp1 =-kz[c]*pow(mZ0,(kz[c]-1));
         if( c == 0 )
	    momex[c] = -(mForces[0]*xp1*y1*z1+mForces[1]*x1*yp1*z1+mForces[2]*x1*y1*zp1);
	 else if( c== 1 )
	    momex[c] = -(mForces[1]*xp1*y1*z1+mForces[3]*x1*yp1*z1+mForces[4]*x1*y1*zp1);
         else
	    momex[c] = -(mForces[2]*xp1*y1*z1+mForces[4]*x1*yp1*z1+mForces[5]*x1*y1*zp1);
      }
   }
}




//-----------------------------------------------------------------------
void Source::perturb( double h, int comp )
{
   if( comp == 0 )
      mX0 += h;
   else if( comp == 1 )
      mY0 += h;
   else if( comp == 2 )
      mZ0 += h;
   else if( comp >= 3 && comp <= 8 )
      mForces[comp-3] += h;
   else if( comp == 9 )
      mT0 += h;
   else
      mFreq += h;
}

//-----------------------------------------------------------------------
void Source::filter_timefunc( Filter* filter_ptr, double tstart, double dt, int nsteps )
{
   if( !m_is_filtered )
   {
      double (*timeFunc)(double f, double t,double* par, int npar, int* ipar, int nipar );
      switch( mTimeDependence )
      {
      case iRicker:
	 timeFunc = RickerWavelet;
	 break;
      case iGaussian :
	 timeFunc   = Gaussian;
	 break;
      case iRamp :
	 timeFunc = Ramp;
	 break;
      case iTriangle :
	 timeFunc = Triangle;
	 break;
      case iSawtooth :
	 timeFunc = Sawtooth;
	 break;
      case iSmoothWave :
	 timeFunc = SmoothWave;
	 break;
      case iErf :
	 timeFunc = Erf;
	 break;
      case iVerySmoothBump :
	 timeFunc = VerySmoothBump;
	 break;
      case iC6SmoothBump :
	 timeFunc = C6SmoothBump;
	 break;
      case iRickerInt :
	 timeFunc = RickerInt;
	 break;
      case iBrune :
	 timeFunc = Brune;
	 break;
      case iBruneSmoothed :
	 timeFunc = BruneSmoothed;
	 break;
      case iDBrune :
	 timeFunc = DBrune;
	 break;
      case iGaussianWindow :
	 timeFunc = GaussianWindow;
	 break;
      case iLiu :
	 timeFunc = Liu;
	 break;
      case iDirac :
	 timeFunc = Dirac;
	 break;
      case iDiscrete :
	 timeFunc = Discrete;
	 break;
      default:
	 cout << "ERROR in Source::filter_timefunc, source type not recoginzed" << endl;
      }

      // Convert to discrete representation
      double *discfunc = new double[nsteps];
      for (int k=0; k < nsteps; k++ )
	 discfunc[k] = timeFunc( mFreq, tstart+k*dt-mT0, mPar, mNpar, mIpar, mNipar );
      mTimeDependence = iDiscrete;

// Filter the discretized function 
      filter_ptr->evaluate( nsteps, &discfunc[0], &discfunc[0] );

// Give the source time function a smooth start if this is a 2-pass (forward + backward) bandpass filter
      if( filter_ptr->get_passes() == 2 && filter_ptr->get_type() == bandPass )
      {    
	 double wghv, xi;
	 int p0=3, p=20 ; // First non-zero time level, and number of points in ramp;

	 for( int i=1 ; i<=p0-1 ; i++ )
	 {
	    discfunc[i-1] = 0;
	 }
	 for( int i=p0 ; i<=p0+p ; i++ )
	 {
	    wghv = 0;
	    xi = (i-p0)/((double) p);
	 // polynomial P(xi), P(0) = 0, P(1)=1
	    wghv = xi*xi*xi*xi*(35-84*xi+70*xi*xi-20*xi*xi*xi);
	    discfunc[i-1] *=wghv;
	 }
      }

   // Save discrete function
      mNipar = 1;
      mIpar = new int[mNipar];
      mIpar[0] = nsteps;

      mFreq = 1./dt;   
      delete[] mPar;
      mNpar = nsteps+1;
      mPar = new double[mNpar];
      //      mPar[0] = tstart;
      mPar[0] = tstart-mT0;
      for( int i=0 ; i < nsteps; i++ )
         mPar[i+1] = discfunc[i];
      delete[] discfunc;

   // Build the spline representation
      spline_interpolation();
      m_is_filtered = true;
   }
}

//-----------------------------------------------------------------------
int Source::spline_interpolation( )
{
   // Assume mPar[1], to mPar[npts] contain the function
   // Assume mIpar[0] contains npts
   // Assume mFreq contains 1/dt, and mPar[0] is tstart.
   // Compute the six spline coefficients for each interval and return in mPar[1],to mPar[6*(npts-1)]
   if( mTimeDependence == iDiscrete )
   {
      int npts = mIpar[0];
      //            cout << "before spline interp" << endl;
      //            cout << "npts = " << npts << " t0 = " << mPar[0] << " dt= " << 1/mFreq << endl;
      //            for( int i=0 ; i < npts ; i++ )
      //      	 cout << "fun[" << i << "] = "<< mPar[i+1] << endl;

      Qspline quinticspline( npts, &mPar[1], mPar[0], 1/mFreq );
      double tstart = mPar[0];
      delete[] mPar;
      mPar = new double[6*(npts-1)+1];
      mNpar = 6*(npts-1)+1;
      mPar[0] = tstart;
      double* qsppt = quinticspline.get_polycof_ptr();
      for( int i=0 ; i < 6*(npts-1) ; i++ )
         mPar[i+1] = qsppt[i];
      //      cout << "after spline interp" << endl;
      //      for( int i=0 ; i < npts ; i++ )
      //	 cout << "fun[" << i << "] = "<< mPar[6*i+1] << endl;
      
      return 1;
   }
   else
      return 0;
}

//-----------------------------------------------------------------------
Source* Source::copy( std::string a_name )
{
   if( a_name == " " )
      a_name = mName;

   Source* retval = new Source();
   retval->m_i0 = m_i0;
   retval->m_j0 = m_j0;
   retval->m_k0 = m_k0;
   retval->m_grid = m_grid;
   retval->mName = a_name;
   retval->mIsMomentSource = mIsMomentSource;
   retval->mForces.push_back(mForces[0]);
   retval->mForces.push_back(mForces[1]);
   retval->mForces.push_back(mForces[2]);
   if( mIsMomentSource )
   {
      retval->mForces.push_back(mForces[3]);
      retval->mForces.push_back(mForces[4]);
      retval->mForces.push_back(mForces[5]);
   }
   retval->mFreq = mFreq;
   retval->mT0 = mT0;
   retval->mGridPointSet = mGridPointSet;
   retval->m_zRelativeToTopography = m_zRelativeToTopography;
   retval->mX0 = mX0;
   retval->mY0 = mY0;
   retval->mZ0 = mZ0;

   retval->mNpar = mNpar;
   retval->mPar = new double[mNpar];
   for( int i=0 ; i < mNpar ; i++ )
      retval->mPar[i] = mPar[i];

   retval->mNipar = mNipar;
   retval->mIpar = new int[mNipar];
   for( int i=0 ; i < mNipar ; i++ )
      retval->mIpar[i] = mIpar[i];

   retval->mNcyc = mNcyc;
   retval->m_derivative = m_derivative;
   retval->mTimeDependence = mTimeDependence;   
   for( int i=0 ; i < 11 ; i++ )
      retval->m_dir[i] = m_dir[i];
   retval->m_is_filtered = m_is_filtered;
   return retval;
}

//-----------------------------------------------------------------------
double Source::find_min_exponent() const
{
   // smallest number x, such that exp(x) does not cause underflow
  return -700.0;
}


//-----------------------------------------------------------------------
void Source::compute_metric_at_source( EW* a_EW, double q, double r, double s, int ic,
				       int jc, int kc, int g, double& zq, double& zr,
				       double& zs ) const
{
   int Nz = a_EW->m_global_nz[g];
   double h = a_EW->mGridSize[g];
   if( kc > Nz-3 )
   {
      // Treat downmost grid lines as cartesian
      zq = zr = 0;
      zs = h;
   }
   else
   {
      // Derivative of metric wrt. source position. Not yet fully implemented.
      double zqdx0, zqdy0, zqsz0, zrdx0, zrdy0, zrdz0;
	       
      // 3. Recompute metric to sixth order accuracy. Increased accuracy needed because
      //    of the multiplication with a singular (Dirac) function.
      //    compute only in the processor where the point is interior

      bool eightptstencil = true;
      double ai=q-ic;
      double bi=r-jc;
      double ci=s-kc;

      double d6cofi[8], d6cofj[8], d6cofk[8];
      if( eightptstencil )
      {
	 // Eight point stencil, smooth wrt. source position, 
	 // needed for source optimization
	 getmetdwgh( ai, d6cofi );
	 getmetdwgh( bi, d6cofj );
	 if( kc <= 3 && ci < 0 )
	 {
	    getmetdwgh7( ci, d6cofk );
	    d6cofk[7] = 0;
	 }
	 else
	    getmetdwgh( ci, d6cofk );
      }
      else
      {
	 // Seven point stencil, ok for forward solver.
	 getmetdwgh7( ai, d6cofi );
	 d6cofi[7] = 0;
	 getmetdwgh7( bi, d6cofj );
	 d6cofj[7] = 0;
	 getmetdwgh7( ci, d6cofk );
	 d6cofk[7] = 0;
      }

      double a6cofi[8], a6cofj[8], a6cofk[8];
      if( eightptstencil )
      {
	 getmetwgh( ai, a6cofi );
	 getmetwgh( bi, a6cofj );
	 if( kc <= 3 && ci < 0 )
	 {
	    getmetwgh7( ci, a6cofk );
	    a6cofk[7] = 0;
	 }
	 else
	    getmetwgh( ci, a6cofk );
      }
      else
      {
	 getmetwgh7( ai, a6cofi );
	 a6cofi[7] = 0;
	 getmetwgh7( bi, a6cofj );
	 a6cofj[7] = 0;
	 getmetwgh7( ci, a6cofk );
	 a6cofk[7] = 0;
      }

      // Assume grid uniform in x and y, compute metric with z=z(q,r,s)
      zq = zr = zs = 0;
      int order;
      double zetaBreak;
      a_EW->get_gridgen_info( order, zetaBreak );
      double zpar = (s-1)/(zetaBreak*(Nz-1));
      double kBreak = 1 + zetaBreak*(Nz-1);
      if( zpar >= 1 )
      {
	 zq = 0;
	 zr = 0;
	 //	       zqdx0 = 0;
	 //	       zqdy0 = 0;
	 //	       zrdx0 = 0;
	 //	       zrdy0 = 0;
	 //               zqdz0 = 1/h;
	 //	       zrdz0 = 1/h;
      }
      else 
      {
	 double tauq=0, taur=0;
	 double tauqdx0=0, tauqdy0=0, taurdx0=0, taurdy0=0;
	 for( int j=jc-3; j <= jc+4 ; j++ )
	    for( int i=ic-3; i <= ic+4 ; i++ )
	    {
	       tauq += d6cofi[i-(ic-3)]*a6cofj[j-(jc-3)]*a_EW->mTopoGridExt(i,j,1);
	       taur += a6cofi[i-(ic-3)]*d6cofj[j-(jc-3)]*a_EW->mTopoGridExt(i,j,1);
	       //                  tauqdx0 += dd6cofi[i-(ic-3)]* a6cofj[j-(jc-3)]*a_EW->mTopoGridExt(i,j,1);
	       //                  tauqdy0 +=  d6cofi[i-(ic-3)]*da6cofj[j-(jc-3)]*a_EW->mTopoGridExt(i,j,1);
	       //                  taurdx0 += da6cofi[i-(ic-3)]* d6cofj[j-(jc-3)]*a_EW->mTopoGridExt(i,j,1);
	       //                  taurdy0 +=  a6cofi[i-(ic-3)]*dd6cofj[j-(jc-3)]*a_EW->mTopoGridExt(i,j,1);
	    }
	 double powo = pow(1-zpar,order);
	 zq = (-tauq)*powo;
	 zr = (-taur)*powo;
	 //               zqdx0 = (-tauqdx0)*powo;
	 //               zqdy0 = (-tauqdy0)*powo;
	 //               zrdx0 = (-taurdx0)*powo;
	 //               zrdy0 = (-taurdy0)*powo;
      }
      double tauavg = 0;
      for( int j=jc-3; j <= jc+4 ; j++ )
	 for( int i=ic-3; i <= ic+4 ; i++ )
	    tauavg += a6cofi[i-(ic-3)]*a6cofj[j-(jc-3)]*a_EW->mTopoGridExt(i,j,1);

      // Compute dz/ds directly from the grid mapping, the explicit expression here
      // should be the same as in EW::curvilinear_grid_mapping
      double zMax = a_EW->m_zmin[a_EW->mNumberOfCartesianGrids-1] - (Nz-kBreak)*h;
      double c1 = zMax + tauavg - h*(kBreak-1);
      double z1d;
      for( int k=kc-3 ; k <= kc+ 4; k++ ) 
      {
	 zpar = (k-1)/(zetaBreak*(Nz-1));
	 if( zpar >= 1 )
	    z1d = zMax + (k-kBreak)*h;
	 else
	 {
	    z1d = (1-zpar)*(-tauavg) + zpar*(zMax + c1*(1-zpar));
	    for( int o=2 ; o < order ; o++ )
	       z1d += zpar*c1*pow(1-zpar,o);
	 }
	 zs += d6cofk[k-(kc-3)]*z1d;
      }
      //	    for( int k=kc-3; k <= kc+3 ; k++ )
      //	       for( int j=jc-3; j <= jc+3 ; j++ )
      //		  for( int i=ic-3; i <= ic+3 ; i++ )
      //		  {
      //		  zq += d6cofi[i-(ic-3)]*a6cofj[j-(jc-3)]*a6cofk[k-(kc-3)]*a_EW->mZ(i,j,k);
      //		  zr += a6cofi[i-(ic-3)]*d6cofj[j-(jc-3)]*a6cofk[k-(kc-3)]*a_EW->mZ(i,j,k);
      //		zs += a6cofi[i-(ic-3)]*a6cofj[j-(jc-3)]*d6cofk[k-(kc-3)]*a_EW->mZ(i,j,k);
      //		  }

      //      cout << " zpar at source = " << (s-1)/(zetaBreak*(Nz-1)) << endl;
   }
}
