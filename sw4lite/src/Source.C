#include "mpi.h"
#include "GridPointSource.h"
#include "Source.h"
#include "Require.h"

#include <fenv.h>
#include <cmath>

#include  "EW.h"
//#include "Filter.h"
//#include "Qspline.h"

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
	       float_sw4 frequency, 
	       float_sw4 t0,
	       float_sw4 x0, 
	       float_sw4 y0, 
	       float_sw4 z0,
	       float_sw4 Mxx,
	       float_sw4 Mxy,
	       float_sw4 Mxz,
	       float_sw4 Myy,
	       float_sw4 Myz,
	       float_sw4 Mzz,
	       timeDep tDep,
	       const char *name,
	       bool topodepth, 
	       int ncyc, 
	       float_sw4* pars, int npar, int* ipars, int nipar, bool correctForMu):
  mIsMomentSource(true),
  mFreq(frequency),
  mT0(t0),
  mX0(x0), mY0(y0), mZ0(z0), m_zTopo(-1e38),
  mIgnore(false),
  mTimeDependence(tDep),
  mNcyc(ncyc),
  m_zRelativeToTopography(topodepth),
  mShearModulusFactor(correctForMu),
  m_derivative(-1),
  m_is_filtered(false),
  m_myPoint(false)
{
   mForces.resize(6);
   mForces[0] = Mxx;
   mForces[1] = Mxy;
   mForces[2] = Mxz;
   mForces[3] = Myy;
   mForces[4] = Myz;
   mForces[5] = Mzz;
   mName = name;

   mNpar = npar;
   if( mNpar > 0 )
   {
      mPar   = new float_sw4[mNpar];
      for( int i= 0 ; i < mNpar ; i++ )
	 mPar[i] = pars[i];
   }
   else
   {
      mNpar = 2;
      mPar = new float_sw4[2];
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

   //   if( mTimeDependence == iDiscrete || mTimeDependence == iDiscrete6moments )
   //      spline_interpolation();
   //   else
   {
      mPar[0] = find_min_exponent();
      mPar[1] = mNcyc;
   }
   adjust_zcoord(a_ew);
   compute_mapped_coordinates(a_ew);
   m_i0 = static_cast<int>(round(mQ0));
   m_j0 = static_cast<int>(round(mR0));
   m_k0 = static_cast<int>(round(mS0));
   m_myPoint = a_ew->interior_point_in_proc( m_i0, m_j0, m_grid );
}

//-----------------------------------------------------------------------
Source::Source(EW *a_ew, float_sw4 frequency, float_sw4 t0,
	       float_sw4 x0, float_sw4 y0, float_sw4 z0,
	       float_sw4 Fx,
	       float_sw4 Fy,
	       float_sw4 Fz,
	       timeDep tDep,
	       const char *name, 
	       bool topodepth,
	       int ncyc, 
	       float_sw4* pars, int npar, int* ipars, int nipar, bool correctForMu ):
  mIsMomentSource(false),
  mFreq(frequency),
  mT0(t0),
  mX0(x0), mY0(y0), mZ0(z0), m_zTopo(-1e38),
  mIgnore(false),
  mTimeDependence(tDep),
  mNcyc(ncyc),
  m_zRelativeToTopography(topodepth),
  m_derivative(-1),
  m_is_filtered(false),
  mShearModulusFactor(correctForMu),
  m_myPoint(false)
{
  mForces.resize(3);
  mForces[0] = Fx;
  mForces[1] = Fy;
  mForces[2] = Fz;
  mName = name;

  mNpar = npar;
  if( mNpar > 0 )
  {
     mPar   = new float_sw4[mNpar];
     for( int i= 0 ; i < mNpar ; i++ )
	mPar[i] = pars[i];
  }
  else
  {
     mNpar = 2;
     mPar = new float_sw4[2];
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
  //  if( mTimeDependence == iDiscrete || mTimeDependence == iDiscrete6moments )
  //     spline_interpolation();
  //  else
  {
     mPar[0] = find_min_exponent();
     mPar[1] = mNcyc;
  }
   adjust_zcoord(a_ew);
   compute_mapped_coordinates(a_ew);
   m_i0 = static_cast<int>(round(mQ0));
   m_j0 = static_cast<int>(round(mR0));
   m_k0 = static_cast<int>(round(mS0));
   m_myPoint = a_ew->interior_point_in_proc( m_i0, m_j0, m_grid );
}

//-----------------------------------------------------------------------
Source::Source()
{
   mNpar = 0;
   mNipar = 0;
}

//-----------------------------------------------------------------------
Source::~Source()
{
   if( mNpar > 0 )
      delete[] mPar;
   if( mNipar > 0 )
      delete[] mIpar;
}

//-----------------------------------------------------------------------
float_sw4 Source::getX0() const
{
  return mX0;
}

//-----------------------------------------------------------------------
float_sw4 Source::getY0() const
{
  return mY0;
}

//-----------------------------------------------------------------------
float_sw4 Source::getZ0() const
{
  return mZ0;
}


//-----------------------------------------------------------------------
float_sw4 Source::getDepth() const
{
  return mZ0-m_zTopo;
}

//-----------------------------------------------------------------------
float_sw4 Source::getOffset() const
{
  return mT0;
}

//-----------------------------------------------------------------------
float_sw4 Source::getFrequency() const
{
  return mFreq;
}

//-----------------------------------------------------------------------
void Source::setMaxFrequency(float_sw4 max_freq)
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
void Source::getForces( float_sw4& fx, float_sw4& fy, float_sw4& fz ) const
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
void Source::getMoments( float_sw4& mxx, float_sw4& mxy, float_sw4& mxz, float_sw4& myy, float_sw4& myz, float_sw4& mzz ) const
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
void Source::setMoments( float_sw4 mxx, float_sw4 mxy, float_sw4 mxz, float_sw4 myy, float_sw4 myz, float_sw4 mzz )
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
float_sw4 Source::getAmplitude() const
{
  float_sw4 amplitude=0;
  if (mIsMomentSource)
  {
    float_sw4 msqr=0;
    for (int q=0; q<6; q++)
      msqr += SQR(mForces[q]);
    //    amplitude = mAmp*sqrt(msqr/2.);
    msqr += SQR(mForces[1])+SQR(mForces[2])+SQR(mForces[4]);
    amplitude = sqrt(0.5*msqr);
  }
  else
  {
    float_sw4 fsqr=0;
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
void Source::limit_frequency( int ppw, float_sw4 minvsoh )
{
   float_sw4 freqlim = minvsoh/(ppw);

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
float_sw4 Source::compute_t0_increase(float_sw4 t0_min) const
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
void Source::adjust_t0( float_sw4 dt0 )
{
   if( dt0 > 0 && !m_is_filtered )
      mT0 += dt0;
}

//-----------------------------------------------------------------------
float_sw4 Source::dt_to_resolve( int ppw ) const
{
  float_sw4 dt_resolved = 0;
  if( mTimeDependence == iBrune || mTimeDependence == iBruneSmoothed ||  mTimeDependence == iDBrune)
    {
      const float_sw4 t95 = 4.744/mFreq;
      dt_resolved = t95/ppw;
    }
  else
    {

    }
  return dt_resolved;
}

//-----------------------------------------------------------------------
int Source::ppw_to_resolve( float_sw4 dt ) const
{
  int ppw = 1;
  if( mTimeDependence == iBrune || mTimeDependence == iBruneSmoothed ||  mTimeDependence == iDBrune)
    {
      const float_sw4 t95 = 4.744/mFreq;
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
void Source::set_dirderivative( float_sw4 dir[11] )
{
   for( int i=0 ; i < 11 ; i++ )
      m_dir[i] = dir[i];
   m_derivative = 11;
}

//-----------------------------------------------------------------------
void Source::set_parameters( float_sw4 x[11] )
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
void Source::get_parameters( float_sw4 x[11] ) const
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
void Source::setFrequency( float_sw4 freq )
{
   mFreq=freq;
}

//-----------------------------------------------------------------------
void Source::adjust_zcoord( EW* a_ew )
{
// z-adjustment for depth below the topography, if necessary.
// Also defines m_zTopo and checks that the source is below the surface
   a_ew->find_topo_zcoord_all( mX0, mY0, m_zTopo );
   if( m_zRelativeToTopography )
      mZ0 += m_zTopo;
   if ( mZ0 < m_zTopo-1e-9)
   {
      mIgnore = true;
      cout << "Ignoring Source at X= " << mX0 << " Y= " << mY0 << " Z= " << mZ0 <<
	 "because it is above the surface, at z= " << m_zTopo<< "\n";
      return;
   }
}

//-----------------------------------------------------------------------
void Source::compute_mapped_coordinates( EW* a_ew )
{
// Defines the grid (m_grid) and the mapped coordinates (mQ0,mR0,mS0) within the grid.
   int i, j, k;
   // computeNearestGridPoint called to find m_grid:
   a_ew->computeNearestGridPoint( i, j, k, m_grid, mX0, mY0, mZ0 );
   // Use m_grid to invert mapping on the correct grid
   float_sw4 s0;
   a_ew->invert_grid_mapping( m_grid, mX0, mY0, mZ0, mQ0, mR0, s0 );
   // The curvilinear mapping can only be inverted by the owning processor
   // Broadcast the s-parameter to all processors, if this is a curvilinear grid.
   // Note: invert_grid_mapping sets s0=-1e38 if inversion fails.
   if( a_ew->m_topography_exists && m_grid == a_ew->mNumberOfGrids-1 )
      MPI_Allreduce( &s0, &mS0, 1, a_ew->m_mpifloat, MPI_MAX, a_ew->m_cartesian_communicator );
   else
      mS0 = s0; // Cartesian grid, all processes know s0.
}

//-----------------------------------------------------------------------
void Source::getsourcewgh( float_sw4 ai, float_sw4 wgh[6], float_sw4 dwghda[6], float_sw4 ddwghda[6] ) const
{
   // Moments k=0,1,2,3,4 exact, two cont. derivatives wrt. position
   float_sw4 p5 = ai*ai*ai*ai*ai*(5.0/3-7.0/24*ai -17/12.0*ai*ai+1.125*ai*ai*ai-0.25*ai*ai*ai*ai);
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
void Source::getsourcedwgh( float_sw4 ai, float_sw4 wgh[6], float_sw4 dwghda[6], float_sw4 ddwghda[6] ) const
{
   // Moments k=0,1,2,3,4 exact, two cont. derivatives wrt. position
   float_sw4 p5 = ai*ai*ai*ai*(-25.0/12-0.75*ai + 59.0/12*ai*ai - 4*ai*ai*ai + ai*ai*ai*ai);
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
void Source::getsourcewghlow( float_sw4 ai, float_sw4 wgh[6], float_sw4 dwghda[6], float_sw4 ddwghda[6] ) const
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
void Source::getsourcedwghlow( float_sw4 ai, float_sw4 wgh[6], float_sw4 dwghda[6], float_sw4 ddwghda[6] ) const
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
void Source::getmetwgh( float_sw4 ai, float_sw4 wgh[8], float_sw4 dwgh[8],
			float_sw4 ddwgh[8], float_sw4 dddwgh[8] ) const
{
   float_sw4 pol = ai*ai*ai*ai*ai*ai*ai*(-251+135*ai+25*ai*ai-
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

   // Derivative wrt. ai
   pol = ai*ai*ai*ai*ai*ai*(-1757.0/720 + 1.5*ai + 0.31250*ai*ai - (1.375*ai*ai*ai-0.275*ai*ai*ai*ai)/3);
   dwgh[0] = -1.0/60 + 1.0/90*ai+ ai*ai/16 + 23.0/36*ai*ai*ai - 223.0/144*ai*ai*ai*ai -
      17.0/120*ai*ai*ai*ai*ai - pol;
   dwgh[1] = 3.0/20 - 3.0/20*ai - 0.5*ai*ai-13.0/3*ai*ai*ai + 97.0/9*ai*ai*ai*ai +
      ai*ai*ai*ai*ai + 7*pol;
   dwgh[2] = -0.75 + 1.5*ai + 13.0/16*ai*ai + 155.0*ai*ai*ai/12-103.0*5.0/16*ai*ai*ai*ai
      - 121.0/40*ai*ai*ai*ai*ai - 21*pol;
   dwgh[3] = -49.0/18*ai - 4*49.0/9.0*ai*ai*ai + 385.0*5.0/36*ai*ai*ai*ai +
      61.0/12*ai*ai*ai*ai*ai + 35*pol;
   dwgh[4] = 0.75 + 1.5*ai - 13.0/16*ai*ai + 89.0/4*ai*ai*ai - 1537.0*5/144.0*ai*ai*ai*ai -
      41.0/8*ai*ai*ai*ai*ai - 35*pol;
   dwgh[5] = -3.0/20 - 3.0/20*ai + 0.5*ai*ai-41.0/3*ai*ai*ai + 32*ai*ai*ai*ai +
      3.1*ai*ai*ai*ai*ai + 21*pol;
   dwgh[6] = 1.0/60 + 1.0/90*ai - 1.0/16*ai*ai + 167.0/36*ai*ai*ai - 1537.0/144*ai*ai*ai*ai -
      25.0/24*ai*ai*ai*ai*ai - 7*pol;
   dwgh[7] = -2.0/3*ai*ai*ai + 55.0/36*ai*ai*ai*ai + 3.0/20*ai*ai*ai*ai*ai + pol;

   // Second derivative wrt. ai
   pol = ai*ai*ai*ai*ai*(-1757.0/120 + 10.5*ai + 2.5*ai*ai - 4.125*ai*ai*ai + 11.0/12*ai*ai*ai*ai);
   ddwgh[0] = 1.0/90 + 0.125*ai + 23.0/12*ai*ai - 223.0/36*ai*ai*ai - 17.0/24*ai*ai*ai*ai - pol;
   ddwgh[1] = -3.0/20 - ai - 13.0*ai*ai + 4*97.0/9.0*ai*ai*ai + 5*ai*ai*ai*ai + 7*pol;
   ddwgh[2] = 1.5 + 13.0/8*ai + 155.0/4*ai*ai - 103.0*5.0/4*ai*ai*ai - 121.0/8*ai*ai*ai*ai - 21*pol;
   ddwgh[3] = -49.0/18 - 4*49.0/3.0*ai*ai + 385.0*5.0/9.0*ai*ai*ai + 5*61.0/12*ai*ai*ai*ai + 35*pol;
   ddwgh[4] = 1.5 -13.0/8*ai+89.0*3.0/4*ai*ai - 1537.0*5.0/36*ai*ai*ai - 205.0/8*ai*ai*ai*ai - 35*pol;
   ddwgh[5] = -3.0/20 + ai - 41.0*ai*ai + 128*ai*ai*ai + 15.5*ai*ai*ai*ai + 21*pol;
   ddwgh[6] = 1.0/90 - 0.125*ai + 167.0/12*ai*ai - 1537.0/36*ai*ai*ai - 125.0/24*ai*ai*ai*ai - 7*pol;
   ddwgh[7] = -2*ai*ai + 220.0/36*ai*ai*ai + 0.75*ai*ai*ai*ai + pol;

   // Third derivative wrt. ai
   pol = ai*ai*ai*ai*(-1757.0/24 + 63*ai + 17.5*ai*ai - 33*ai*ai*ai + 8.25*ai*ai*ai*ai);
   dddwgh[0] = 0.125 + 23.0/6*ai-223.0/12*ai*ai-17.0/6*ai*ai*ai - pol;
   dddwgh[1] = -1 - 26.0*ai + 4*97.0/3*ai*ai + 20*ai*ai*ai + 7*pol;
   dddwgh[2] =  1.625 + 77.5*ai - 386.25*ai*ai -60.5*ai*ai*ai - 21*pol;
   dddwgh[3] = -392.0/3*ai + 1925.0/3*ai*ai + 305.0/3*ai*ai*ai + 35*pol;
   dddwgh[4] = -1.625 + 133.5*ai-7685.0/12*ai*ai - 102.5*ai*ai*ai - 35*pol;
   dddwgh[5] = 1 - 82.0*ai + 384.0*ai*ai + 62.0*ai*ai*ai + 21*pol;
   dddwgh[6] = -0.125 + 167.0/6*ai - 1537.0/12*ai*ai - 125.0/6*ai*ai*ai - 7*pol;
   dddwgh[7] = -4*ai + 220.0/12*ai*ai + 3*ai*ai*ai + pol;
}

//-----------------------------------------------------------------------
void Source::getmetdwgh( float_sw4 ai, float_sw4 wgh[8] ) const
{
   float_sw4 pol = ai*ai*ai*ai*ai*ai*( -827 + 420*ai + 165*ai*ai - 180*ai*ai*ai
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
void Source::getmetwgh7( float_sw4 ai, float_sw4 wgh[7] ) const
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
void Source::getmetdwgh7( float_sw4 ai, float_sw4 wgh[7] ) const
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
// Sources at grid refinement boundary
void Source::getsourcewghNM2sm6(  float_sw4 alph,  float_sw4 wghk[6] ) const
{
      wghk[0] = (1.0-10.0*alph*alph*alph+15.0*alph*alph*alph*alph-6.0*alph*
alph*alph*alph*alph)*(alph*alph*alph*alph/24.0+alph/12.0-alph*alph/24.0-alph*
alph*alph/12.0);
      wghk[1] = (1.0-10.0*alph*alph*alph+15.0*alph*alph*alph*alph-6.0*alph*
alph*alph*alph*alph)*(-alph*alph*alph*alph/6.0-2.0/3.0*alph+2.0/3.0*alph*alph+
alph*alph*alph/6.0)+(10.0*alph*alph*alph-15.0*alph*alph*alph*alph+6.0*alph*alph
*alph*alph*alph)*(-pow(alph-1.0,2.0)/30.0+alph/10.0-1.0/10.0+pow(alph-1.0,4.0)/
30.0-pow(alph-1.0,3.0)/10.0);
      wghk[2] = (1.0-10.0*alph*alph*alph+15.0*alph*alph*alph*alph-6.0*alph*
alph*alph*alph*alph)*(1.0+alph*alph*alph*alph/4.0-5.0/4.0*alph*alph)+(10.0*alph
*alph*alph-15.0*alph*alph*alph*alph+6.0*alph*alph*alph*alph*alph)*(5.0/8.0*pow(
alph-1.0,2.0)-3.0/4.0*alph+3.0/4.0-pow(alph-1.0,4.0)/8.0+pow(alph-1.0,3.0)/4.0)
;
      wghk[3] = (1.0-10.0*alph*alph*alph+15.0*alph*alph*alph*alph-6.0*alph*
alph*alph*alph*alph)*(-alph*alph*alph/6.0-alph*alph*alph*alph/6.0+2.0/3.0*alph+
2.0/3.0*alph*alph)+(10.0*alph*alph*alph-15.0*alph*alph*alph*alph+6.0*alph*alph*
alph*alph*alph)*(5.0/6.0-7.0/6.0*pow(alph-1.0,2.0)+alph/6.0+pow(alph-1.0,4.0)/
6.0-pow(alph-1.0,3.0)/6.0);
      wghk[4] = (1.0-10.0*alph*alph*alph+15.0*alph*alph*alph*alph-6.0*alph*
alph*alph*alph*alph)*(alph*alph*alph*alph/24.0-alph/12.0-alph*alph/24.0+alph*
alph*alph/12.0)+(10.0*alph*alph*alph-15.0*alph*alph*alph*alph+6.0*alph*alph*
alph*alph*alph)*(7.0/12.0*pow(alph-1.0,2.0)+alph/2.0-1.0/2.0-pow(alph-1.0,4.0)/
12.0);
      wghk[5] = (10.0*alph*alph*alph-15.0*alph*alph*alph*alph+6.0*alph*alph*
alph*alph*alph)*(-pow(alph-1.0,2.0)/120.0-alph/60.0+1.0/60.0+pow(alph-1.0,4.0)/
120.0+pow(alph-1.0,3.0)/60.0);
}

//-----------------------------------------------------------------------
void Source::getsourcedwghNM2sm6(  float_sw4 alph,  float_sw4 dwghk[6] ) const
{
      dwghk[0] = (1.0-10.0*alph*alph*alph+15.0*alph*alph*alph*alph-6.0*alph*
alph*alph*alph*alph)*(-1.0/12.0-alph*alph*alph/6.0+alph*alph/4.0+alph/12.0);
      dwghk[1] = (1.0-10.0*alph*alph*alph+15.0*alph*alph*alph*alph-6.0*alph*
alph*alph*alph*alph)*(2.0/3.0+2.0/3.0*alph*alph*alph-alph*alph/2.0-4.0/3.0*alph
)+(10.0*alph*alph*alph-15.0*alph*alph*alph*alph+6.0*alph*alph*alph*alph*alph)*(
-1.0/6.0-2.0/15.0*pow(alph-1.0,3.0)+3.0/10.0*pow(alph-1.0,2.0)+alph/15.0);
      dwghk[2] = (1.0-10.0*alph*alph*alph+15.0*alph*alph*alph*alph-6.0*alph*
alph*alph*alph*alph)*(-alph*alph*alph+5.0/2.0*alph)+(10.0*alph*alph*alph-15.0*
alph*alph*alph*alph+6.0*alph*alph*alph*alph*alph)*(2.0+pow(alph-1.0,3.0)/2.0
-3.0/4.0*pow(alph-1.0,2.0)-5.0/4.0*alph);
      dwghk[3] = (1.0-10.0*alph*alph*alph+15.0*alph*alph*alph*alph-6.0*alph*
alph*alph*alph*alph)*(-2.0/3.0+2.0/3.0*alph*alph*alph+alph*alph/2.0-4.0/3.0*
alph)+(10.0*alph*alph*alph-15.0*alph*alph*alph*alph+6.0*alph*alph*alph*alph*
alph)*(-5.0/2.0-2.0/3.0*pow(alph-1.0,3.0)+7.0/3.0*alph+pow(alph-1.0,2.0)/2.0);
      dwghk[4] = (1.0-10.0*alph*alph*alph+15.0*alph*alph*alph*alph-6.0*alph*
alph*alph*alph*alph)*(1.0/12.0-alph*alph*alph/6.0-alph*alph/4.0+alph/12.0)+(
10.0*alph*alph*alph-15.0*alph*alph*alph*alph+6.0*alph*alph*alph*alph*alph)*(2.0
/3.0+pow(alph-1.0,3.0)/3.0-7.0/6.0*alph);
      dwghk[5] = (10.0*alph*alph*alph-15.0*alph*alph*alph*alph+6.0*alph*alph*
alph*alph*alph)*(-pow(alph-1.0,3.0)/30.0-pow(alph-1.0,2.0)/20.0+alph/60.0);
}

//-----------------------------------------------------------------------
void Source::getsourcewghNM1sm6(  float_sw4 alph,  float_sw4 wghk[6] ) const
{
      wghk[0] = (1.0-10.0*alph*alph*alph+15.0*alph*alph*alph*alph-6.0*alph*
alph*alph*alph*alph)*(-alph*alph/30.0+alph/10.0+alph*alph*alph*alph/30.0-alph*
alph*alph/10.0);
      wghk[1] = (1.0-10.0*alph*alph*alph+15.0*alph*alph*alph*alph-6.0*alph*
alph*alph*alph*alph)*(5.0/8.0*alph*alph-3.0/4.0*alph-alph*alph*alph*alph/8.0+
alph*alph*alph/4.0)+(10.0*alph*alph*alph-15.0*alph*alph*alph*alph+6.0*alph*alph
*alph*alph*alph)*(-5.0/48.0*pow(alph-1.0,3.0)+pow(alph-1.0,4.0)/48.0+alph/6.0
-1.0/6.0+pow(alph-1.0,2.0)/24.0);
      wghk[2] = (1.0-10.0*alph*alph*alph+15.0*alph*alph*alph*alph-6.0*alph*
alph*alph*alph*alph)*(1.0-7.0/6.0*alph*alph+alph/6.0+alph*alph*alph*alph/6.0-
alph*alph*alph/6.0)+(10.0*alph*alph*alph-15.0*alph*alph*alph*alph+6.0*alph*alph
*alph*alph*alph)*(4.0/15.0*pow(alph-1.0,3.0)-pow(alph-1.0,4.0)/15.0-16.0/15.0*
alph+16.0/15.0+4.0/15.0*pow(alph-1.0,2.0));
      wghk[3] = (1.0-10.0*alph*alph*alph+15.0*alph*alph*alph*alph-6.0*alph*
alph*alph*alph*alph)*(7.0/12.0*alph*alph+alph/2.0-alph*alph*alph*alph/12.0)+(
10.0*alph*alph*alph-15.0*alph*alph*alph*alph+6.0*alph*alph*alph*alph*alph)*(1.0
/4.0+3.0/4.0*alph-3.0/16.0*pow(alph-1.0,3.0)-pow(alph-1.0,2.0)/2.0+pow(alph-1.0
,4.0)/16.0);
      wghk[4] = (1.0-10.0*alph*alph*alph+15.0*alph*alph*alph*alph-6.0*alph*
alph*alph*alph*alph)*(-alph*alph/120.0-alph/60.0+alph*alph*alph*alph/120.0+alph
*alph*alph/60.0)+(10.0*alph*alph*alph-15.0*alph*alph*alph*alph+6.0*alph*alph*
alph*alph*alph)*(pow(alph-1.0,3.0)/48.0-pow(alph-1.0,4.0)/48.0+alph/6.0-1.0/6.0
+5.0/24.0*pow(alph-1.0,2.0));
      wghk[5] = (10.0*alph*alph*alph-15.0*alph*alph*alph*alph+6.0*alph*alph*
alph*alph*alph)*(pow(alph-1.0,3.0)/240.0+pow(alph-1.0,4.0)/240.0-alph/60.0+1.0/
60.0-pow(alph-1.0,2.0)/60.0);

}

//-----------------------------------------------------------------------
void Source::getsourcedwghNM1sm6(  float_sw4 alph,  float_sw4 dwghk[6] ) const
{
      dwghk[0] = (1.0-10.0*alph*alph*alph+15.0*alph*alph*alph*alph-6.0*alph*
alph*alph*alph*alph)*(-1.0/10.0-2.0/15.0*alph*alph*alph+3.0/10.0*alph*alph+alph
/15.0);
      dwghk[1] = (1.0-10.0*alph*alph*alph+15.0*alph*alph*alph*alph-6.0*alph*
alph*alph*alph*alph)*(3.0/4.0+alph*alph*alph/2.0-3.0/4.0*alph*alph-5.0/4.0*alph
)+(10.0*alph*alph*alph-15.0*alph*alph*alph*alph+6.0*alph*alph*alph*alph*alph)*(
-alph/12.0-1.0/12.0+5.0/16.0*pow(alph-1.0,2.0)-pow(alph-1.0,3.0)/12.0);
      dwghk[2] = (1.0-10.0*alph*alph*alph+15.0*alph*alph*alph*alph-6.0*alph*
alph*alph*alph*alph)*(-1.0/6.0-2.0/3.0*alph*alph*alph+7.0/3.0*alph+alph*alph/
2.0)+(10.0*alph*alph*alph-15.0*alph*alph*alph*alph+6.0*alph*alph*alph*alph*alph
)*(-8.0/15.0*alph+8.0/5.0-4.0/5.0*pow(alph-1.0,2.0)+4.0/15.0*pow(alph-1.0,3.0))
;
      dwghk[3] = (1.0-10.0*alph*alph*alph+15.0*alph*alph*alph*alph-6.0*alph*
alph*alph*alph*alph)*(-1.0/2.0+alph*alph*alph/3.0-7.0/6.0*alph)+(10.0*alph*alph
*alph-15.0*alph*alph*alph*alph+6.0*alph*alph*alph*alph*alph)*(-7.0/4.0+9.0/16.0
*pow(alph-1.0,2.0)+alph-pow(alph-1.0,3.0)/4.0);
      dwghk[4] = (1.0-10.0*alph*alph*alph+15.0*alph*alph*alph*alph-6.0*alph*
alph*alph*alph*alph)*(1.0/60.0-alph*alph*alph/30.0-alph*alph/20.0+alph/60.0)+(
10.0*alph*alph*alph-15.0*alph*alph*alph*alph+6.0*alph*alph*alph*alph*alph)*(1.0
/4.0-pow(alph-1.0,2.0)/16.0-5.0/12.0*alph+pow(alph-1.0,3.0)/12.0);
      dwghk[5] = (10.0*alph*alph*alph-15.0*alph*alph*alph*alph+6.0*alph*alph*
alph*alph*alph)*(-1.0/60.0-pow(alph-1.0,2.0)/80.0+alph/30.0-pow(alph-1.0,3.0)/
60.0);
}

//-----------------------------------------------------------------------
void Source::getsourcewghNsm6(  float_sw4 alph,  float_sw4 wghk[6] ) const
{
      wghk[0] = (1.0-5.0/4.0*alph*alph*alph+15.0/16.0*alph*alph*alph*alph-3.0/
16.0*alph*alph*alph*alph*alph)*(-5.0/48.0*alph*alph*alph+alph*alph*alph*alph/
48.0+alph/6.0+alph*alph/24.0);
      wghk[1] = (1.0-5.0/4.0*alph*alph*alph+15.0/16.0*alph*alph*alph*alph-3.0/
16.0*alph*alph*alph*alph*alph)*(4.0/15.0*alph*alph*alph-alph*alph*alph*alph/
15.0-16.0/15.0*alph+4.0/15.0*alph*alph)+(5.0/4.0*alph*alph*alph-15.0/16.0*alph*
alph*alph*alph+3.0/16.0*alph*alph*alph*alph*alph)*(-4.0/105.0*pow(alph-2.0,2.0)
+pow(alph-2.0,4.0)/105.0+16.0/105.0*alph-32.0/105.0-4.0/105.0*pow(alph-2.0,3.0)
);
      wghk[2] = (1.0-5.0/4.0*alph*alph*alph+15.0/16.0*alph*alph*alph*alph-3.0/
16.0*alph*alph*alph*alph*alph)*(1.0+3.0/4.0*alph-3.0/16.0*alph*alph*alph-alph*
alph/2.0+alph*alph*alph*alph/16.0)+(5.0/4.0*alph*alph*alph-15.0/16.0*alph*alph*
alph*alph+3.0/16.0*alph*alph*alph*alph*alph)*(5.0/24.0*pow(alph-2.0,2.0)-pow(
alph-2.0,4.0)/48.0-alph/2.0+1.0+pow(alph-2.0,3.0)/16.0);
      wghk[3] = (1.0-5.0/4.0*alph*alph*alph+15.0/16.0*alph*alph*alph*alph-3.0/
16.0*alph*alph*alph*alph*alph)*(alph*alph*alph/48.0-alph*alph*alph*alph/48.0+
alph/6.0+5.0/24.0*alph*alph)+(5.0/4.0*alph*alph*alph-15.0/16.0*alph*alph*alph*
alph+3.0/16.0*alph*alph*alph*alph*alph)*(5.0/6.0+pow(alph-2.0,4.0)/48.0+alph/
12.0-pow(alph-2.0,3.0)/48.0-pow(alph-2.0,2.0)/3.0);
      wghk[4] = (1.0-5.0/4.0*alph*alph*alph+15.0/16.0*alph*alph*alph*alph-3.0/
16.0*alph*alph*alph*alph*alph)*(alph*alph*alph/240.0+alph*alph*alph*alph/240.0-
alph/60.0-alph*alph/60.0)+(5.0/4.0*alph*alph*alph-15.0/16.0*alph*alph*alph*alph
+3.0/16.0*alph*alph*alph*alph*alph)*(-pow(alph-2.0,3.0)/80.0+7.0/40.0*pow(alph
-2.0,2.0)-pow(alph-2.0,4.0)/80.0+3.0/10.0*alph-3.0/5.0);
      wghk[5] = (5.0/4.0*alph*alph*alph-15.0/16.0*alph*alph*alph*alph+3.0/16.0*
alph*alph*alph*alph*alph)*(-pow(alph-2.0,2.0)/84.0+pow(alph-2.0,4.0)/336.0-alph
/28.0+1.0/14.0+pow(alph-2.0,3.0)/112.0);
}

//-----------------------------------------------------------------------
void Source::getsourcedwghNsm6(  float_sw4 alph,  float_sw4 dwghk[6] ) const
{
   //      dwghk[0] = (1.0-5.0/4.0*alph*alph*alph+15.0/16.0*alph*alph*alph*alph-3.0/
   //16.0*alph*alph*alph*alph*alph)*(-alph/12.0-1.0/6.0+5.0/16.0*alph*alph-alph*alph
   //*alph/12.0);
   //      dwghk[1] = (1.0-5.0/4.0*alph*alph*alph+15.0/16.0*alph*alph*alph*alph-3.0/
   //16.0*alph*alph*alph*alph*alph)*(-8.0/15.0*alph+16.0/15.0-4.0/5.0*alph*alph+4.0/
   //15.0*alph*alph*alph)+(5.0/4.0*alph*alph*alph-15.0/16.0*alph*alph*alph*alph+3.0/
   //16.0*alph*alph*alph*alph*alph)*(-alph/12.0+1.0/12.0+5.0/16.0*pow(alph-3.0,2.0)-
   //pow(alph-3.0,3.0)/12.0);
   //      dwghk[2] = (1.0-5.0/4.0*alph*alph*alph+15.0/16.0*alph*alph*alph*alph-3.0/
   //16.0*alph*alph*alph*alph*alph)*(-3.0/4.0+9.0/16.0*alph*alph+alph-alph*alph*alph
   ///4.0)+(5.0/4.0*alph*alph*alph-15.0/16.0*alph*alph*alph*alph+3.0/16.0*alph*alph*
   //alph*alph*alph)*(-8.0/15.0*alph+8.0/3.0-4.0/5.0*pow(alph-3.0,2.0)+4.0/15.0*pow(
   //alph-3.0,3.0));
   //      dwghk[3] = (1.0-5.0/4.0*alph*alph*alph+15.0/16.0*alph*alph*alph*alph-3.0/
   //16.0*alph*alph*alph*alph*alph)*(-1.0/6.0-alph*alph/16.0-5.0/12.0*alph+alph*alph
   //*alph/12.0)+(5.0/4.0*alph*alph*alph-15.0/16.0*alph*alph*alph*alph+3.0/16.0*alph
   //*alph*alph*alph*alph)*(-15.0/4.0+9.0/16.0*pow(alph-3.0,2.0)+alph-pow(alph-3.0,
   //3.0)/4.0);
   //      dwghk[4] = (1.0-5.0/4.0*alph*alph*alph+15.0/16.0*alph*alph*alph*alph-3.0/
   //16.0*alph*alph*alph*alph*alph)*(1.0/60.0-alph*alph/80.0+alph/30.0-alph*alph*
   //alph/60.0)+(5.0/4.0*alph*alph*alph-15.0/16.0*alph*alph*alph*alph+3.0/16.0*alph*
   //alph*alph*alph*alph)*(13.0/12.0-pow(alph-3.0,2.0)/16.0-5.0/12.0*alph+pow(alph
   //-3.0,3.0)/12.0);
   //      dwghk[5] = (5.0/4.0*alph*alph*alph-15.0/16.0*alph*alph*alph*alph+3.0/16.0*
   //alph*alph*alph*alph*alph)*(-1.0/12.0-pow(alph-3.0,2.0)/80.0+alph/30.0-pow(alph
   //-3.0,3.0)/60.0);

      dwghk[0] = (1.0-5.0/4.0*alph*alph*alph+15.0/16.0*alph*alph*alph*alph-3.0/
16.0*alph*alph*alph*alph*alph)*(-alph/12.0-1.0/6.0+5.0/16.0*alph*alph-alph*alph
*alph/12.0);
      dwghk[1] = (1.0-5.0/4.0*alph*alph*alph+15.0/16.0*alph*alph*alph*alph-3.0/
16.0*alph*alph*alph*alph*alph)*(-8.0/15.0*alph+16.0/15.0-4.0/5.0*alph*alph+4.0/
15.0*alph*alph*alph)+(5.0/4.0*alph*alph*alph-15.0/16.0*alph*alph*alph*alph+3.0/
16.0*alph*alph*alph*alph*alph)*(-32.0/105.0-4.0/105.0*pow(alph-2.0,3.0)+4.0/
35.0*pow(alph-2.0,2.0)+8.0/105.0*alph);
      dwghk[2] = (1.0-5.0/4.0*alph*alph*alph+15.0/16.0*alph*alph*alph*alph-3.0/
16.0*alph*alph*alph*alph*alph)*(-3.0/4.0+9.0/16.0*alph*alph+alph-alph*alph*alph
/4.0)+(5.0/4.0*alph*alph*alph-15.0/16.0*alph*alph*alph*alph+3.0/16.0*alph*alph*
alph*alph*alph)*(4.0/3.0+pow(alph-2.0,3.0)/12.0-3.0/16.0*pow(alph-2.0,2.0)-5.0/
12.0*alph);
      dwghk[3] = (1.0-5.0/4.0*alph*alph*alph+15.0/16.0*alph*alph*alph*alph-3.0/
16.0*alph*alph*alph*alph*alph)*(-1.0/6.0-alph*alph/16.0-5.0/12.0*alph+alph*alph
*alph/12.0)+(5.0/4.0*alph*alph*alph-15.0/16.0*alph*alph*alph*alph+3.0/16.0*alph
*alph*alph*alph*alph)*(-17.0/12.0-pow(alph-2.0,3.0)/12.0+pow(alph-2.0,2.0)/16.0
+2.0/3.0*alph);
      dwghk[4] = (1.0-5.0/4.0*alph*alph*alph+15.0/16.0*alph*alph*alph*alph-3.0/
16.0*alph*alph*alph*alph*alph)*(1.0/60.0-alph*alph/80.0+alph/30.0-alph*alph*
alph/60.0)+(5.0/4.0*alph*alph*alph-15.0/16.0*alph*alph*alph*alph+3.0/16.0*alph*
alph*alph*alph*alph)*(2.0/5.0+pow(alph-2.0,3.0)/20.0+3.0/80.0*pow(alph-2.0,2.0)
-7.0/20.0*alph);
      dwghk[5] = (5.0/4.0*alph*alph*alph-15.0/16.0*alph*alph*alph*alph+3.0/16.0*
alph*alph*alph*alph*alph)*(-1.0/84.0-pow(alph-2.0,3.0)/84.0-3.0/112.0*pow(alph
-2.0,2.0)+alph/42.0);

}

//-----------------------------------------------------------------------
void Source::getsourcewghP1sm6(  float_sw4 alph,  float_sw4 wghk[6] ) const
{
      wghk[0] = (1.0-5.0/4.0*alph*alph*alph+15.0/16.0*alph*alph*alph*alph-3.0/
16.0*alph*alph*alph*alph*alph)*(-4.0/105.0*alph*alph+alph*alph*alph*alph/105.0+
16.0/105.0*alph-4.0/105.0*alph*alph*alph);
      wghk[1] = (1.0-5.0/4.0*alph*alph*alph+15.0/16.0*alph*alph*alph*alph-3.0/
16.0*alph*alph*alph*alph*alph)*(5.0/24.0*alph*alph-alph*alph*alph*alph/48.0-
alph/2.0+alph*alph*alph/16.0)+(5.0/4.0*alph*alph*alph-15.0/16.0*alph*alph*alph*
alph+3.0/16.0*alph*alph*alph*alph*alph)*(-pow(alph-2.0,2.0)/96.0+pow(alph-2.0,
4.0)/384.0+alph/24.0-1.0/12.0-pow(alph-2.0,3.0)/96.0);
      wghk[2] = (1.0-5.0/4.0*alph*alph*alph+15.0/16.0*alph*alph*alph*alph-3.0/
16.0*alph*alph*alph*alph*alph)*(1.0+alph*alph*alph*alph/48.0+alph/12.0-alph*
alph*alph/48.0-alph*alph/3.0)+(5.0/4.0*alph*alph*alph-15.0/16.0*alph*alph*alph*
alph+3.0/16.0*alph*alph*alph*alph*alph)*(pow(alph-2.0,2.0)/6.0-pow(alph-2.0,4.0
)/96.0-alph/3.0+2.0/3.0+pow(alph-2.0,3.0)/48.0);
      wghk[3] = (1.0-5.0/4.0*alph*alph*alph+15.0/16.0*alph*alph*alph*alph-3.0/
16.0*alph*alph*alph*alph*alph)*(-alph*alph*alph/80.0+7.0/40.0*alph*alph-alph*
alph*alph*alph/80.0+3.0/10.0*alph)+(5.0/4.0*alph*alph*alph-15.0/16.0*alph*alph*
alph*alph+3.0/16.0*alph*alph*alph*alph*alph)*(1.0-5.0/16.0*pow(alph-2.0,2.0)+
pow(alph-2.0,4.0)/64.0);
      wghk[4] = (1.0-5.0/4.0*alph*alph*alph+15.0/16.0*alph*alph*alph*alph-3.0/
16.0*alph*alph*alph*alph*alph)*(-alph*alph/84.0+alph*alph*alph*alph/336.0-alph/
28.0+alph*alph*alph/112.0)+(5.0/4.0*alph*alph*alph-15.0/16.0*alph*alph*alph*
alph+3.0/16.0*alph*alph*alph*alph*alph)*(-pow(alph-2.0,3.0)/48.0+pow(alph-2.0,
2.0)/6.0-pow(alph-2.0,4.0)/96.0+alph/3.0-2.0/3.0);
      wghk[5] = (5.0/4.0*alph*alph*alph-15.0/16.0*alph*alph*alph*alph+3.0/16.0
*alph*alph*alph*alph*alph)*(-pow(alph-2.0,2.0)/96.0+pow(alph-2.0,4.0)/384.0-
alph/24.0+1.0/12.0+pow(alph-2.0,3.0)/96.0);
}

//-----------------------------------------------------------------------
void Source::getsourcedwghP1sm6(  float_sw4 alph,  float_sw4 dwghk[6] ) const
{
      dwghk[0] = (1.0-5.0/4.0*alph*alph*alph+15.0/16.0*alph*alph*alph*alph-3.0/
16.0*alph*alph*alph*alph*alph)*(-16.0/105.0-4.0/105.0*alph*alph*alph+4.0/35.0*
alph*alph+8.0/105.0*alph);
      dwghk[1] = (1.0-5.0/4.0*alph*alph*alph+15.0/16.0*alph*alph*alph*alph-3.0/
16.0*alph*alph*alph*alph*alph)*(1.0/2.0+alph*alph*alph/12.0-3.0/16.0*alph*alph
-5.0/12.0*alph)+(5.0/4.0*alph*alph*alph-15.0/16.0*alph*alph*alph*alph+3.0/16.0*
alph*alph*alph*alph*alph)*(-1.0/12.0-pow(alph-2.0,3.0)/96.0+pow(alph-2.0,2.0)/
32.0+alph/48.0);
      dwghk[2] = (1.0-5.0/4.0*alph*alph*alph+15.0/16.0*alph*alph*alph*alph-3.0/
16.0*alph*alph*alph*alph*alph)*(-1.0/12.0-alph*alph*alph/12.0+alph*alph/16.0+
2.0/3.0*alph)+(5.0/4.0*alph*alph*alph-15.0/16.0*alph*alph*alph*alph+3.0/16.0*
alph*alph*alph*alph*alph)*(1.0+pow(alph-2.0,3.0)/24.0-pow(alph-2.0,2.0)/16.0-
alph/3.0);
      dwghk[3] = (1.0-5.0/4.0*alph*alph*alph+15.0/16.0*alph*alph*alph*alph-3.0/
16.0*alph*alph*alph*alph*alph)*(-3.0/10.0+alph*alph*alph/20.0+3.0/80.0*alph*
alph-7.0/20.0*alph)+(5.0/4.0*alph*alph*alph-15.0/16.0*alph*alph*alph*alph+3.0/
16.0*alph*alph*alph*alph*alph)*(5.0/8.0*alph-5.0/4.0-pow(alph-2.0,3.0)/16.0);
      dwghk[4] = (1.0-5.0/4.0*alph*alph*alph+15.0/16.0*alph*alph*alph*alph-3.0/
16.0*alph*alph*alph*alph*alph)*(1.0/28.0-alph*alph*alph/84.0-3.0/112.0*alph*
alph+alph/42.0)+(5.0/4.0*alph*alph*alph-15.0/16.0*alph*alph*alph*alph+3.0/16.0*
alph*alph*alph*alph*alph)*(1.0/3.0+pow(alph-2.0,3.0)/24.0+pow(alph-2.0,2.0)/
16.0-alph/3.0);
      dwghk[5] = (5.0/4.0*alph*alph*alph-15.0/16.0*alph*alph*alph*alph+3.0/16.0
*alph*alph*alph*alph*alph)*(-pow(alph-2.0,3.0)/96.0-pow(alph-2.0,2.0)/32.0+alph
/48.0);
}

//-----------------------------------------------------------------------
void Source::set_grid_point_sources4( EW *a_EW, vector<GridPointSource*>& point_sources ) const
{
// note that this routine is called from all processors, for each input source 
//   int i,j,k,g;
//   a_EW->computeNearestGridPoint( i, j, k, g, mX0, mY0, mZ0 );
//   float_sw4 q, r, s;
//   float_sw4 h = a_EW->mGridSize[g];
//   bool canBeInverted, curvilinear;
   float_sw4 normwgh[4]={17.0/48.0, 59.0/48.0, 43.0/48.0, 49.0/48.0 };

   //   if( g == a_EW->mNumberOfGrids-1 && a_EW->topographyExists() )
   //   {
   //// Curvilinear
   //// Problem when the curvilinear mapping is NOT analytic:
   //// This routine can only compute the 's' coordinate if (mX0, mY0) is owned by this processor
   //      canBeInverted = a_EW->invert_curvilinear_grid_mapping( g, mX0, mY0, mZ0, q, r, s );

   //      // Broadcast the computed s to all processors. 
   //      // First find out the ID of a processor that defines s ...
   //      int s_owner = -1;
   //      if( canBeInverted )
   //         MPI_Comm_rank(MPI_COMM_WORLD, &s_owner );
   //      int s_owner_tmp = s_owner;
   //      MPI_Allreduce( &s_owner_tmp, &s_owner, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
   //      // ...then broadcast s
   //      if( s_owner > -1 )
   //	 MPI_Bcast( &s, 1, MPI_FLOAT_SW4, s_owner, MPI_COMM_WORLD );
   //      else
   //      {
   //	 printf("ERROR in Source::set_grid_point_sources4, no processor could invert the grid mapping \n");
   //	 MPI_Abort(MPI_COMM_WORLD,1);
   //      }

   //// if s < 0, the source is located above the grid and the call to
   //// find_curvilinear_derivatives_at_point will fail
   //      if (s<0.)
   //      {
   //	 float_sw4 xTop, yTop, zTop;
   //	 a_EW->curvilinear_grid_mapping(q, r, 0., xTop, yTop, zTop);
   //	 float_sw4 lat, lon;
   //	 a_EW->computeGeographicCoord(mX0, mY0, lon, lat);
   //	 printf("Found a source above the curvilinear grid! Lat=%e, Lon=%e, source Z-level = %e, grid boundary Z = %e\n",
   //		lat, lon, mZ0, zTop);
   //	 MPI_Abort(MPI_COMM_WORLD,1);
   //      }
   //      curvilinear   = true;
   //      canBeInverted = true;
   //   }
   //   else
   //   {
// Cartesian case
//      q = mX0/h+1;
//      r = mY0/h+1;
//      s = (mZ0-a_EW->m_zmin[g])/h+1;
//      canBeInverted = true;
//      curvilinear   = false;
//   }
   bool canBeInverted = true;
   bool curvilinear = a_EW->topographyExists() && m_grid == a_EW->mNumberOfGrids-1;

   int Ni = a_EW->m_global_nx[m_grid];
   int Nj = a_EW->m_global_ny[m_grid];
   int Nz = a_EW->m_global_nz[m_grid];

   int ic, jc, kc;
   bool upperbndry, lowerbndry, ccbndry, gridrefbndry;
   float_sw4 ai, bi, ci;

// Delta distribution
   float_sw4 wghi[6], wghj[6], wghk[6], wghix[6], wghjy[6], wghkz[6];
   float_sw4 wghixx[6], wghjyy[6], wghkzz[6];
// Delta' distribution
   float_sw4 dwghi[6], dwghj[6], dwghk[6], dwghix[6], dwghjy[6], dwghkz[6];
   float_sw4 dwghixx[6], dwghjyy[6], dwghkzz[6];

   // k-weights across mesh refinement boundary
   float_sw4 wghkref[6], dwghkref[6], wghrefkz[6], wghrefkzz[6];

   if( canBeInverted )
   {
 // Compute source location and weights in source discretization
      ic = static_cast<int>(floor(mQ0));
      jc = static_cast<int>(floor(mR0));
      kc = static_cast<int>(floor(mS0));

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
      ccbndry = a_EW->topographyExists() &&  ( (upperbndry && m_grid == a_EW->mNumberOfGrids-2) ||
                                               (lowerbndry && m_grid == a_EW->mNumberOfGrids-1)    );
   //      ccbndry = false;

// gridrefbndry=true if at the interface between two cartesian grids of different refinements.
      gridrefbndry = (upperbndry && m_grid < a_EW->mNumberOfGrids-1 && !ccbndry) ||
	 (lowerbndry && m_grid>0 && !ccbndry );

// If not at the interface between different grids, bias stencil away 
// from the boundary.
      if( !ccbndry && !gridrefbndry )
      {
	 if( kc <= 2 )    kc = 3;
	 if( kc >= Nz-2 ) kc = Nz-3;
      }
      ai=mQ0-ic, bi=mR0-jc, ci=mS0-kc;

   // Delta distribution
      getsourcewgh( ai, wghi, wghix, wghixx );
      getsourcewgh( bi, wghj, wghjy, wghjyy );
      getsourcewgh( ci, wghk, wghkz, wghkzz );

// Delta' distribution
      getsourcedwgh( ai, dwghi, dwghix, dwghixx );
      getsourcedwgh( bi, dwghj, dwghjy, dwghjyy );
      getsourcedwgh( ci, dwghk, dwghkz, dwghkzz );

   // Special boundary stencil at free surface
      if( !ccbndry && !gridrefbndry && (kc == 3 && ci <= 0) )
      {
	 getsourcewghlow( ci, wghk, wghkz, wghkzz );
	 getsourcedwghlow( ci, dwghk, dwghkz, dwghkzz );
      }

  // Special source discretization across grid refinement boundary
      //      cout << "grid ref bndr " << gridrefbndry << " kc = " << kc << endl;
      if( gridrefbndry )
      {
	 if( lowerbndry )
	 {
	    if( kc == Nz-1 )
	    {
	       getsourcedwghNM1sm6( ci, dwghk );
	       getsourcewghNM1sm6(  ci,  wghk );
	       wghkref[3]  = wghk[3]*0.5;
	       wghk[3]     = wghk[3]*0.5;
	       wghkref[4]  = wghk[4];
	       wghkref[5]  = wghk[5];

	       dwghkref[3] = dwghk[3]*0.5;
	       dwghk[3]    = dwghk[3]*0.5;
	       dwghkref[4] = dwghk[4];
	       dwghkref[5] = dwghk[5];	       

	       wghkref[3] /= normwgh[0];
	       wghkref[4] /= normwgh[1];
	       wghkref[5] /= normwgh[2];
	       wghk[3]    /= normwgh[0];
	       wghk[2]    /= normwgh[1];
	       wghk[1]    /= normwgh[2];
	       wghk[0]    /= normwgh[3];
	       dwghkref[3] /= normwgh[0];
	       dwghkref[4] /= normwgh[1];
	       dwghkref[5] /= normwgh[2];
	       dwghk[3]    /= normwgh[0];
	       dwghk[2]    /= normwgh[1];
	       dwghk[1]    /= normwgh[2];
	       dwghk[0]    /= normwgh[3];
	    }
	    else if( kc == Nz-2 )
	    {
	       getsourcedwghNM2sm6( ci, dwghk );
	       getsourcewghNM2sm6(  ci,  wghk );
	       wghkref[4]  = wghk[4]*0.5;
	       wghk[4]     = wghk[4]*0.5;
	       wghkref[5]  = wghk[5];

	       dwghkref[4] = dwghk[4]*0.5;
	       dwghk[4]    = dwghk[4]*0.5;
	       dwghkref[5] = dwghk[5];

	       wghkref[4] /= normwgh[0];
	       wghkref[5] /= normwgh[1];
	       wghk[4]    /= normwgh[0];
	       wghk[3]    /= normwgh[1];
	       wghk[2]    /= normwgh[2];
	       wghk[1]    /= normwgh[3];
	       dwghkref[4] /= normwgh[0];
	       dwghkref[5] /= normwgh[1];
	       dwghk[4]    /= normwgh[0];
	       dwghk[3]    /= normwgh[1];
	       dwghk[2]    /= normwgh[2];
	       dwghk[1]    /= normwgh[3];
	    }
	    else if( kc == Nz-3 )
	    {
	       getsourcedwgh( ci, dwghk, wghrefkz, wghrefkzz );
	       getsourcewgh(  ci,  wghk, wghrefkz, wghrefkzz );
	       wghkref[5]  = wghk[5]*0.5;
	       wghk[5]     = wghk[5]*0.5;
	       dwghkref[5] = dwghk[5]*0.5;
	       dwghk[5]    = dwghk[5]*0.5;
 //	       cout << " sumwgh = " << dwghk[0]+dwghk[1]+dwghk[2]+dwghk[3]+dwghk[4]+dwghk[5]+dwghkref[5] << endl;

	       wghkref[5] /= normwgh[0];
	       wghk[5]    /= normwgh[0];
	       wghk[4]    /= normwgh[1];
	       wghk[3]    /= normwgh[2];
	       wghk[2]    /= normwgh[3];
	       dwghkref[5] /= normwgh[0];
	       dwghk[5]    /= normwgh[0];
	       dwghk[4]    /= normwgh[1];
	       dwghk[3]    /= normwgh[2];
	       dwghk[2]    /= normwgh[3];

	    }
	 }
	 else
	 {
	    if( kc == 1 )
	    {
	       getsourcedwghNsm6( 2*ci, dwghk );
	       getsourcewghNsm6(  2*ci,  wghk );
	       wghkref[0]  = wghk[0];
	       wghkref[1]  = wghk[1];
	       wghkref[2]  = wghk[2]*0.5;
	       wghk[2]     = wghk[2]*0.5;

	       dwghkref[0] = dwghk[0];
	       dwghkref[1] = dwghk[1];
	       dwghkref[2] = dwghk[2]*0.5;
	       dwghk[2]    = dwghk[2]*0.5;

//	       cout << "kc = 1  ref:  " << wghkref[0] << " " << wghkref[1] << " " << wghkref[2] << endl;
//	       cout << "        this: " << wghk[2] << " " << wghk[3] << " " << wghk[4] << " " << wghk[5] << endl;
//	       cout << "  middle sum: " << wghk[2]+wghkref[2] << endl;
//	       cout << " 2*ci = " << 2*ci << endl;

	       wghkref[0] /= normwgh[2];
	       wghkref[1] /= normwgh[1];
	       wghkref[2] /= normwgh[0];
	       wghk[2]    /= normwgh[0];
	       wghk[3]    /= normwgh[1];
	       wghk[4]    /= normwgh[2];
	       wghk[5]    /= normwgh[3];

	       dwghkref[0] /= normwgh[2];
	       dwghkref[1] /= normwgh[1];
	       dwghkref[2] /= normwgh[0];
	       dwghk[2]    /= normwgh[0];
	       dwghk[3]    /= normwgh[1];
	       dwghk[4]    /= normwgh[2];
	       dwghk[5]    /= normwgh[3];

	       dwghk[2] *= 2;
	       dwghk[3] *= 2;
	       dwghk[4] *= 2;
	       dwghk[5] *= 2;
	    }
	    else if( kc == 2 )
	    {
	       getsourcedwghP1sm6( 2*ci, dwghk );
	       getsourcewghP1sm6(  2*ci,  wghk );

	       wghkref[0]  = wghk[0];
	       wghkref[1]  = wghk[1]*0.5;
	       wghk[1]     = wghk[1]*0.5;

	       dwghkref[0] = dwghk[0];
	       dwghkref[1] = dwghk[1]*0.5;
	       dwghk[1]    = dwghk[1]*0.5;

	       wghkref[0] /= normwgh[1];
	       wghkref[1] /= normwgh[0];
	       wghk[1]    /= normwgh[0];
	       wghk[2]    /= normwgh[1];
	       wghk[3]    /= normwgh[2];
	       wghk[4]    /= normwgh[3];

	       dwghkref[0] /= normwgh[1];
	       dwghkref[1] /= normwgh[0];
	       dwghk[1]    /= normwgh[0];
	       dwghk[2]    /= normwgh[1];
	       dwghk[3]    /= normwgh[2];
	       dwghk[4]    /= normwgh[3];

	       dwghk[1] *= 2;
	       dwghk[2] *= 2;
	       dwghk[3] *= 2;
	       dwghk[4] *= 2;
	       dwghk[5] *= 2;
	       
	    }
	    else if( kc == 3 )
	    {
	       getsourcedwgh( ci, dwghk, wghrefkz, wghrefkzz );
	       for( int k=0 ; k <= 5 ;k++ )
		  dwghk[k] *= 0.5;
	       getsourcewgh(  ci,  wghk, wghrefkz, wghrefkzz );
	       wghkref[0]  = wghk[0]*0.5;
	       wghk[0]     = wghk[0]*0.5;
	       dwghkref[0] = dwghk[0]*0.5;
	       dwghk[0]    = dwghk[0]*0.5;

	       wghkref[0] /= normwgh[0];
	       wghk[0]    /= normwgh[0];
	       wghk[1]    /= normwgh[1];
	       wghk[2]    /= normwgh[2];
	       wghk[3]    /= normwgh[3];
	       dwghkref[0] /= normwgh[0];
	       dwghk[0]    /= normwgh[0];
	       dwghk[1]    /= normwgh[1];
	       dwghk[2]    /= normwgh[2];
	       dwghk[3]    /= normwgh[3];
	       dwghk[0] *= 2;
	       dwghk[1] *= 2;
	       dwghk[2] *= 2;
	       dwghk[3] *= 2;
	       dwghk[4] *= 2;
	       dwghk[5] *= 2;

	    }

	 }
      }

// Boundary correction, at upper boundary, but only if SBP operators are used there
//      if( !gridrefbndry && (g == a_EW->mNumberOfGrids-1) && a_EW->is_onesided(g,4)  )
      if( !gridrefbndry && a_EW->is_onesided(m_grid,4)  )
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
      if( !gridrefbndry && a_EW->is_onesided(m_grid,5)  )
      {
	 for( int k=0 ; k <= 5 ; k++ )
	 {
	    if( ( Nz-3 <= k+kc-2) && ( k+kc-2 <= Nz ) )
	    {
	       wghk[k]    /= normwgh[Nz-k-kc+2];
	       dwghk[k]   /= normwgh[Nz-k-kc+2];
	       wghkz[k]   /= normwgh[Nz-k-kc+2];
	       dwghkz[k]  /= normwgh[Nz-k-kc+2];
	       wghkzz[k]  /= normwgh[Nz-k-kc+2];
	       dwghkzz[k] /= normwgh[Nz-k-kc+2];
	    }
	 }
      }
   }
   int myid;
   MPI_Comm_rank(MPI_COMM_WORLD, &myid );
   //   cout << myid << " SOURCE at " << ic << " " << jc << " "  << kc ;
   //   if( canBeInverted )
   //      cout << " can be inverted";
   //   else
   //      cout << " can not be inverted";


  // If source at grid refinement interface, set up variables for 
  // discretization on grid on the other side of the interface
   int icref, jcref;
   float_sw4 h=a_EW->mGridSize[m_grid];
   float_sw4 airef, biref, wghiref[6], wghirefx[6], wghirefxx[6];
   float_sw4 wghjref[6], wghjrefy[6], wghjrefyy[6];
   float_sw4 dwghiref[6], dwghjref[6];
   if( gridrefbndry )
   {
      int Niref, Njref;
      float_sw4 qref, rref;
      if( kc-1 < Nz-kc )
      {
     // kc closer to upper boundary. Source spread to finer grid above.
	 qref = mX0/(0.5*h)+1;
	 rref = mY0/(0.5*h)+1;
	 Niref = a_EW->m_global_nx[m_grid+1];
	 Njref = a_EW->m_global_ny[m_grid+1];
      }
      else
      {
     // kc closer to lower boundary. Source spread to coarser grid below.
	 qref = mX0/(2*h)+1;
	 rref = mY0/(2*h)+1;
	 Niref = a_EW->m_global_nx[m_grid-1];
	 Njref = a_EW->m_global_ny[m_grid-1];
      }
      icref = static_cast<int>(qref);
      jcref = static_cast<int>(rref);
      if( icref <= 2 ) icref = 3;
      if( icref >= Niref-2 ) icref = Niref-3;
      if( jcref <= 2 ) jcref = 3;
      if( jcref >= Njref-2 ) jcref = Njref-3;
      airef = qref-icref;
      biref = rref-jcref;
      getsourcewgh( airef, wghiref, wghirefx, wghirefxx );
      getsourcewgh( biref, wghjref, wghjrefy, wghjrefyy );
      // reuse wghirefx,wghirefxx, these are assumed not to be used with grid refinement.
      getsourcedwgh( airef, dwghiref, wghirefx, wghirefxx );
      getsourcedwgh( biref, dwghjref, wghjrefy, wghjrefyy );
   }

   // Point source. NOTE: Derivatives needed for source inversion not implemented for this case.
   if( !mIsMomentSource )
   {
      for( int k=kc-2 ; k <= kc+3 ; k++ )
	 for( int j=jc-2 ; j <= jc+3 ; j++ )
	    for( int i=ic-2 ; i <= ic+3 ; i++ )
	    {
	       float_sw4 wF = wghi[i-ic+2]*wghj[j-jc+2]*wghk[k-kc+2];
	       if( (wF != 0) && (mForces[0] != 0 || mForces[1] != 0 || mForces[2] != 0) 
		   && a_EW->interior_point_in_proc(i,j,m_grid) ) // checks if (i,j) belongs to this processor
	       {
		  if( curvilinear )
		     wF /= a_EW->mJ(i,j,k);
		  else
		     wF /= h*h*h;

		  if( 1 <= k && k <= Nz )
		  {
		     GridPointSource* sourcePtr = new GridPointSource( mFreq, mT0,
								      i,j,k,m_grid,
								      wF*mForces[0], wF*mForces[1], wF*mForces[2],
								      mTimeDependence, mNcyc, 
								      mPar, mNpar, mIpar, mNipar );
		     point_sources.push_back(sourcePtr);
		  }
		  if( k <= 1 && ccbndry && upperbndry )
		  {
		     int Nzp =a_EW->m_global_nz[m_grid+1];
		     int kk = Nzp - 1 + k;

		     GridPointSource* sourcePtr = new GridPointSource( mFreq, mT0,
								      i,j,kk,m_grid+1,
								      wF*mForces[0], wF*mForces[1], wF*mForces[2],
								      mTimeDependence, mNcyc, 
								      mPar, mNpar, mIpar, mNipar );
		     point_sources.push_back(sourcePtr);
		  }
		  if( k >= Nz && ccbndry && lowerbndry )
		  {
		     int kk = k-Nz + 1;
		     GridPointSource* sourcePtr = new GridPointSource( mFreq, mT0,
									  i,j,kk,m_grid-1,
									  wF*mForces[0], wF*mForces[1], wF*mForces[2],
									  mTimeDependence, mNcyc, 
									  mPar, mNpar, mIpar, mNipar );
		     point_sources.push_back(sourcePtr);
		  }
	       }
	    }
      if( gridrefbndry )
      {
	 for( int k=kc-2 ; k <= kc+3 ; k++ )
	 {
	    if( k <= 1 )
	    {
	       // Finer grid above
	       for( int j=jcref-2 ; j <= jcref+3 ; j++ )
		  for( int i=icref-2 ; i <= icref+3 ; i++ )
		  {
		     float_sw4 wF = wghiref[i-icref+2]*wghjref[j-jcref+2]*wghkref[k-kc+2];
		     if( (wF != 0) && (mForces[0] != 0 || mForces[1] != 0 || mForces[2] != 0) 
			 && a_EW->interior_point_in_proc(i,j,m_grid+1) ) // checks if (i,j) belongs to this processor
		     {
			wF /= 0.125*h*h*h;
			int Nzp =a_EW->m_global_nz[m_grid+1];
			int kk = Nzp - 1 + k;
			GridPointSource* sourcePtr = new GridPointSource( mFreq, mT0,
								      i,j,kk,m_grid+1,
								      wF*mForces[0], wF*mForces[1], wF*mForces[2],
								      mTimeDependence, mNcyc, 
								      mPar, mNpar, mIpar, mNipar );
			point_sources.push_back(sourcePtr);
		     }
		  }
	    }
	    if( k >= Nz )
	    {
	       // Coarser grid below
	       for( int j=jcref-2 ; j <= jcref+3 ; j++ )
		  for( int i=icref-2 ; i <= icref+3 ; i++ )
		  {
		     float_sw4 wF = wghiref[i-icref+2]*wghjref[j-jcref+2]*wghkref[k-kc+2];
		     if( (wF != 0) && (mForces[0] != 0 || mForces[1] != 0 || mForces[2] != 0) 
			 && a_EW->interior_point_in_proc(i,j,m_grid-1) ) // checks if (i,j) belongs to this processor
		     {
			wF /= 8*h*h*h;
			int kk = k-Nz + 1;
			GridPointSource* sourcePtr = new GridPointSource( mFreq, mT0,
									  i,j,kk,m_grid-1,
									  wF*mForces[0], wF*mForces[1], wF*mForces[2],
									  mTimeDependence, mNcyc, 
									  mPar, mNpar, mIpar, mNipar );
			point_sources.push_back(sourcePtr);
		     }
		  }
	    }
	 }
      }// Grid refinement boundary
   } 
   // Moment source.
   else if( mIsMomentSource )
   {
      float_sw4 qX0[3], rX0[3], sX0[3];

      // Gradients of sX0[0]=sX, sX0[1]=sY, and sX0[2]=sZ wrt. (q,r,s)
      float_sw4 dsX0[3], dsY0[3], dsZ0[3];
      // Hessians of sX0[0]=sX, sX0[1]=sY, and sX0[2]=sZ wrt. (q,r,s), in order qq,qr,qs,rr,rs,ss
      //      float_sw4 d2sX0[6], d2sY0[6], d2sZ0[6];
      if( !curvilinear )
      {
	 // Cartesian case, constant metric
	 qX0[0] = 1/h;qX0[1]=0;  qX0[2]=0;
	 rX0[0] = 0;  rX0[1]=1/h;rX0[2]=0;
	 sX0[0] = 0;  sX0[1]=0;  sX0[2]=1/h;
	 dsX0[0] = dsX0[1] = dsX0[2] = 0;
	 dsY0[0] = dsY0[1] = dsY0[2] = 0;
	 dsZ0[0] = dsZ0[1] = dsZ0[2] = 0;
	 //         for( int i=0 ; i < 6  ; i++ )
	 //	    d2sX0[i] = d2sY0[i] = d2sZ0[i] = 0;
      }	 
      else
      {
	 // Compute the curvilinear metric in the processor that owns the source.
	 //   (ic, jc are undefined if canBeInverted is false.)
	 float_sw4 zdertmp[9]={0,0,0,0,0,0,0,0,0}, zq, zr, zs;
	 float_sw4 zqq, zqr, zqs, zrr, zrs, zss;
	 if( a_EW->interior_point_in_proc(ic,jc,m_grid) && canBeInverted )
	 {
	    compute_metric_at_source( a_EW, mQ0, mR0, mS0, ic, jc, kc, m_grid, zq, zr, zs,
				      zqq, zqr, zqs, zrr, zrs, zss );
	    zdertmp[0] = zq;
	    zdertmp[1] = zr;
	    zdertmp[2] = zs;
	    zdertmp[3] = zqq;
	    zdertmp[4] = zqr;
	    zdertmp[5] = zqs;
	    zdertmp[6] = zrr;
	    zdertmp[7] = zrs;
	    zdertmp[8] = zss;
	 }
         float_sw4 zder[9];
	 MPI_Allreduce( zdertmp, zder, 9, a_EW->m_mpifloat, MPI_SUM, MPI_COMM_WORLD );
	 zq  = zder[0];
	 zr  = zder[1];
	 zs  = zder[2];
	 zqq = zder[3];
	 zqr = zder[4];
	 zqs = zder[5];
	 zrr = zder[6];
	 zrs = zder[7];
	 zss = zder[8];

	 qX0[0] = 1/h;
	 qX0[1] = 0;
	 qX0[2] = 0;
	 rX0[0] = 0;
	 rX0[1] = 1/h;
	 rX0[2] = 0;
	 sX0[0] = -zq/(h*zs);
	 sX0[1] = -zr/(h*zs);
	 sX0[2] =   1/zs;

         float_sw4 deni = 1/(h*zs*zs);
	 dsX0[0] = -(zqq*zs-zqs*zq)*deni;
	 dsX0[1] = -(zqr*zs-zrs*zq)*deni;
	 dsX0[2] = -(zqs*zs-zss*zq)*deni;

	 dsY0[0] = -(zqr*zs-zqs*zr)*deni;
	 dsY0[1] = -(zrr*zs-zrs*zr)*deni;
	 dsY0[2] = -(zrs*zs-zss*zr)*deni;

         deni *= h;
	 dsZ0[0] = -zqs*deni;
	 dsZ0[1] = -zrs*deni;
	 dsZ0[2] = -zss*deni;
      }

	    // Gradients of sX0[0]=sX, sX0[1]=sY, and sX0[2]=sZ wrt. (q,r,s)
	    // NYI
      //            float_sw4 dsX0[3], dsY0[3], dsZ0[3], d2sX0[6], d2sY0[6], d2sZ0[6];
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

      //      if( canBeInverted )
      {
	 for( int k=kc-2 ; k <= kc+3 ; k++ )
	    for( int j=jc-2 ; j <= jc+3 ; j++ )
	       for( int i=ic-2 ; i <= ic+3 ; i++ )
	       {
		  float_sw4 wFx=0, wFy=0, wFz=0, dsdp[27];
		  if( a_EW->interior_point_in_proc(i,j,m_grid) ) 
		  {
		     //                     cout << " src at " << i << " " << j << " " << k << endl;
		     wFx += qX0[0]*dwghi[i-ic+2]* wghj[j-jc+2]* wghk[k-kc+2];
		  //		  wFy += qX0[1]*dwghi[i-ic+2]* wghj[j-jc+2]* wghk[k-kc+2]; 
		  //		  wFz += qX0[2]*dwghi[i-ic+2]* wghj[j-jc+2]* wghk[k-kc+2];

		  //		  wFx +=  wghi[i-ic+2]*rX0[0]*dwghj[j-jc+2]* wghk[k-kc+2];
		     wFy +=  wghi[i-ic+2]*rX0[1]*dwghj[j-jc+2]* wghk[k-kc+2];
		  //		  wFz +=  wghi[i-ic+2]*rX0[2]*dwghj[j-jc+2]* wghk[k-kc+2];

		     wFx +=  wghi[i-ic+2]* wghj[j-jc+2]*sX0[0]*dwghk[k-kc+2];
		     //		     wFx +=  sX0[0];
		     wFy +=  wghi[i-ic+2]* wghj[j-jc+2]*sX0[1]*dwghk[k-kc+2];
		     wFz +=  wghi[i-ic+2]* wghj[j-jc+2]*sX0[2]*dwghk[k-kc+2];

		     float_sw4 hi=1.0/h;
		     float_sw4 hi2=hi*hi;
		     float_sw4 wFxdx0 =dwghix[i-ic+2]*  wghj[j-jc+2]*  wghk[k-kc+2]*hi*qX0[0];
		     float_sw4 wFxdy0 = dwghi[i-ic+2]* wghjy[j-jc+2]*  wghk[k-kc+2]*hi*rX0[1];
		     //                     float_sw4 wFxdy0 = 0;
		     float_sw4 wFxdz0 = dwghi[i-ic+2]*  wghj[j-jc+2]* wghkz[k-kc+2]*hi*sX0[2];

		     float_sw4 wFydx0 = wghix[i-ic+2]* dwghj[j-jc+2]*  wghk[k-kc+2]*hi*qX0[0];
		     float_sw4 wFydy0 =  wghi[i-ic+2]*dwghjy[j-jc+2]*  wghk[k-kc+2]*hi*rX0[1];
		     float_sw4 wFydz0 =  wghi[i-ic+2]* dwghj[j-jc+2]* wghkz[k-kc+2]*hi*sX0[2];

		     float_sw4 wFzdx0 = wghix[i-ic+2]*  wghj[j-jc+2]* dwghk[k-kc+2]*sX0[2]*qX0[0];
		     float_sw4 wFzdy0 =  wghi[i-ic+2]* wghjy[j-jc+2]* dwghk[k-kc+2]*sX0[2]*rX0[1];
		     float_sw4 wFzdz0 =  wghi[i-ic+2]*  wghj[j-jc+2]*dwghkz[k-kc+2]*sX0[2]*sX0[2];

		     if( curvilinear && kc <= Nz-3 )
		     {
			wFxdx0 += dwghi[i-ic+2]*wghj[j-jc+2]* wghkz[k-kc+2]*hi*sX0[0] +
                                  wghix[i-ic+2]*wghj[j-jc+2]* dwghk[k-kc+2]*hi*sX0[0] +
                                   wghi[i-ic+2]*wghj[j-jc+2]*dwghkz[k-kc+2]*sX0[0]*sX0[0] +
			    wghi[i-ic+2]*  wghj[j-jc+2]* dwghk[k-kc+2]*(hi*dsX0[0]+sX0[0]*dsX0[2]);

                        wFxdy0 += dwghi[i-ic+2]*  wghj[j-jc+2]* wghkz[k-kc+2]*hi*sX0[1] +
                                  wghi[i-ic+2]*  wghjy[j-jc+2]* dwghk[k-kc+2]*hi*sX0[0] +
                                  wghi[i-ic+2]*  wghj[j-jc+2]* dwghkz[k-kc+2]*sX0[0]*sX0[1] +
 			    wghi[i-ic+2]*  wghj[j-jc+2]* dwghk[k-kc+2]*(hi*dsX0[1]+sX0[1]*dsX0[2]);
			//			wFxdy0 += (hi*dsX0[1]+sX0[1]*dsX0[2]);

                        wFxdz0 += wghi[i-ic+2]*  wghj[j-jc+2]* dwghkz[k-kc+2]*sX0[0]*sX0[2] +
			          wghi[i-ic+2]*  wghj[j-jc+2]*  dwghk[k-kc+2]*sX0[2]*dsX0[2];

                        wFydx0 += wghi[i-ic+2]* dwghj[j-jc+2]* wghkz[k-kc+2]*hi*sX0[0] +
			          wghix[i-ic+2]* wghj[j-jc+2]* dwghk[k-kc+2]*hi*sX0[1] +
                                  wghi[i-ic+2]*  wghj[j-jc+2]* dwghkz[k-kc+2]*sX0[0]*sX0[1] +
			          wghi[i-ic+2]*  wghj[j-jc+2]* dwghk[k-kc+2]*(hi*dsY0[0]+sX0[0]*dsY0[2]);

                        wFydy0 += wghi[i-ic+2]* dwghj[j-jc+2]* wghkz[k-kc+2]*hi*sX0[1] +
			          wghi[i-ic+2]* wghjy[j-jc+2]* dwghk[k-kc+2]*hi*sX0[1] +
                                  wghi[i-ic+2]*  wghj[j-jc+2]* dwghkz[k-kc+2]*sX0[1]*sX0[1] +
			          wghi[i-ic+2]*  wghj[j-jc+2]* dwghk[k-kc+2]*(hi*dsY0[1]+sX0[1]*dsY0[2]);

                        wFydz0 += wghi[i-ic+2]*  wghj[j-jc+2]* dwghkz[k-kc+2]*sX0[1]*sX0[2] +
 			          wghi[i-ic+2]*  wghj[j-jc+2]* dwghk[k-kc+2]*sX0[2]*dsY0[2];

			wFzdx0 += wghi[i-ic+2]*  wghj[j-jc+2]* dwghkz[k-kc+2]*sX0[0]*sX0[2] +
 			          wghi[i-ic+2]*  wghj[j-jc+2]* dwghk[k-kc+2]*(hi*dsZ0[0] + sX0[0]*dsZ0[2]);

			wFzdy0 += wghi[i-ic+2]*  wghj[j-jc+2]* dwghkz[k-kc+2]*sX0[1]*sX0[2] +
 			          wghi[i-ic+2]*  wghj[j-jc+2]* dwghk[k-kc+2]*(hi*dsZ0[1] + sX0[1]*dsZ0[2]);

                        wFzdz0 += wghi[i-ic+2]*  wghj[j-jc+2]* dwghk[k-kc+2]*(sX0[2]*dsZ0[2]);

		     }
		  // NOTE:  Source second derivatives wrt. (x0,y0,z0) currently not yet implemented 
		  // for curvilinear grids.

		  // Second derivatives

		     float_sw4 wFxdx0dx0 = dwghixx[i-ic+2]*  wghj[j-jc+2]*  wghk[k-kc+2]*hi2*hi;
		     float_sw4 wFxdx0dy0 = dwghix[i-ic+2]*  wghjy[j-jc+2]*  wghk[k-kc+2]*hi2*hi;
		     float_sw4 wFxdx0dz0 = dwghix[i-ic+2]*  wghj[j-jc+2]*  wghkz[k-kc+2]*hi2*hi;
		     float_sw4 wFxdy0dy0 = dwghi[i-ic+2]* wghjyy[j-jc+2]*  wghk[k-kc+2]*hi2*hi;
		     float_sw4 wFxdy0dz0 = dwghi[i-ic+2]* wghjy[j-jc+2]*  wghkz[k-kc+2]*hi2*hi;
		     float_sw4 wFxdz0dz0 = dwghi[i-ic+2]*  wghj[j-jc+2]* wghkzz[k-kc+2]*hi2*hi;

		     float_sw4 wFydx0dx0 = wghixx[i-ic+2]* dwghj[j-jc+2]*  wghk[k-kc+2]*hi2*hi;
		     float_sw4 wFydx0dy0 = wghix[i-ic+2]* dwghjy[j-jc+2]*  wghk[k-kc+2]*hi2*hi;
		     float_sw4 wFydx0dz0 = wghix[i-ic+2]* dwghj[j-jc+2]*  wghkz[k-kc+2]*hi2*hi;
		     float_sw4 wFydy0dy0 = wghi[i-ic+2]*dwghjyy[j-jc+2]*  wghk[k-kc+2]*hi2*hi;
		     float_sw4 wFydy0dz0 = wghi[i-ic+2]*dwghjy[j-jc+2]*  wghkz[k-kc+2]*hi2*hi;
		     float_sw4 wFydz0dz0 = wghi[i-ic+2]* dwghj[j-jc+2]* wghkzz[k-kc+2]*hi2*hi;

		     float_sw4 wFzdx0dx0 = wghixx[i-ic+2]*  wghj[j-jc+2]* dwghk[k-kc+2]*hi2*hi;
		     float_sw4 wFzdx0dy0 = wghix[i-ic+2]*  wghjy[j-jc+2]* dwghk[k-kc+2]*hi2*hi;
		     float_sw4 wFzdx0dz0 = wghix[i-ic+2]*  wghj[j-jc+2]* dwghkz[k-kc+2]*hi2*hi;
		     float_sw4 wFzdy0dy0 = wghi[i-ic+2]* wghjyy[j-jc+2]* dwghk[k-kc+2]*hi2*hi;
		     float_sw4 wFzdy0dz0 = wghi[i-ic+2]* wghjy[j-jc+2]* dwghkz[k-kc+2]*hi2*hi;
		     float_sw4 wFzdz0dz0 = wghi[i-ic+2]*  wghj[j-jc+2]*dwghkzz[k-kc+2]*hi2*hi;
               
		     float_sw4 jaci;
		     if( curvilinear )
			jaci = 1/a_EW->mJ(i,j,k);
		     else
			jaci = 1.0/(h*h*h);

		     float_sw4 fx = -(mForces[0]*wFx+mForces[1]*wFy+mForces[2]*wFz)*jaci;
		     float_sw4 fy = -(mForces[1]*wFx+mForces[3]*wFy+mForces[4]*wFz)*jaci;
		     float_sw4 fz = -(mForces[2]*wFx+mForces[4]*wFy+mForces[5]*wFz)*jaci;

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
		     float_sw4 dddp[9], dh1[9], dh2[9], dh3[9];
		     dddp[0]  =-wFxdx0*jaci;
		     dddp[1]  =-wFxdy0*jaci;
		     dddp[2]  =-wFxdz0*jaci;
		     dddp[3]  =-wFydx0*jaci;
		     dddp[4]  =-wFydy0*jaci;
		     dddp[5]  =-wFydz0*jaci;
		     dddp[6]  =-wFzdx0*jaci;
		     dddp[7]  =-wFzdy0*jaci;
		     dddp[8]  =-wFzdz0*jaci;

		     //                     if( i == ic && j == jc && k == kc )
		     //		     {
		     //                        cout.precision(16);
		     //                        cout << "forcing " << endl;
		     //			cout << fx << " " << fy << " " << fz << endl;
		     //			cout << "gradient dsdp- = " << endl;
		     //			for( int dd=0 ; dd < 9 ; dd++ )
		     //			   cout << dsdp[dd] << endl;
			//                        cout << "wFz and dwFzdx0 at " << ic << " " << jc << " " << kc << endl;
			//			cout << wFz << endl;
			//			cout << wFzdx0 << endl;
			//			cout << wFzdy0 << endl;
			//			cout << wFzdz0 << endl;
		     //		     }

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
			GridPointSource* sourcePtr = new GridPointSource( mFreq, mT0, i, j, k, m_grid, 
									  fx, fy, fz, mTimeDependence, mNcyc,
									  mPar, mNpar, mIpar, mNipar,
									  dsdp, dddp, dh1, dh2, dh3 );
			if( m_derivative >= 0 )
			   sourcePtr->set_derivative(m_derivative,m_dir);
			point_sources.push_back(sourcePtr);
		     }
		     if( k <= 1 && ccbndry && upperbndry )
		     {
			int Nzp =a_EW->m_global_nz[m_grid+1];
			int kk = Nzp - 1 + k;
			GridPointSource* sourcePtr = new GridPointSource( mFreq, mT0, i, j, kk, m_grid+1, 
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
			GridPointSource* sourcePtr = new GridPointSource( mFreq, mT0, i, j, kk, m_grid-1, 
									  fx, fy, fz, mTimeDependence, mNcyc,
									  mPar, mNpar, mIpar, mNipar,
									  dsdp, dddp, dh1, dh2, dh3 );
			if( m_derivative >= 0 )
			   sourcePtr->set_derivative(m_derivative,m_dir);
			point_sources.push_back(sourcePtr);
		     }
		  }
	       }
	 if( gridrefbndry )
	 {
	    // Tese arrays are currently undefined across the mesh refinement boundary.
	    // --> Source inversion can not be done if the source is located at the interface.
	    float_sw4 dddp[9], dh1[9], dh2[9], dh3[9], dsdp[27];
	    for( int k=kc-2 ; k <= kc+3 ; k++ )
	    {
	       if( k <= 1 )
	       {
		  float_sw4 hi = 1.0/(0.5*h);
		  float_sw4 jaci = 1.0/(0.125*h*h*h);
		  for( int j=jcref-2 ; j <= jcref+3 ; j++ )
		     for( int i=icref-2 ; i <= icref+3 ; i++ )
		     {
			float_sw4 wFx=0, wFy=0, wFz=0;
			if( a_EW->interior_point_in_proc(i,j,m_grid+1) ) 
			{
			   wFx = dwghiref[i-icref+2]* wghjref[j-jcref+2]* wghkref[k-kc+2]*hi;
			   wFy =  wghiref[i-icref+2]*dwghjref[j-jcref+2]* wghkref[k-kc+2]*hi;
			   wFz =  wghiref[i-icref+2]* wghjref[j-jcref+2]*dwghkref[k-kc+2]*hi;
			   float_sw4 fx = -(mForces[0]*wFx+mForces[1]*wFy+mForces[2]*wFz)*jaci;
			   float_sw4 fy = -(mForces[1]*wFx+mForces[3]*wFy+mForces[4]*wFz)*jaci;
			   float_sw4 fz = -(mForces[2]*wFx+mForces[4]*wFy+mForces[5]*wFz)*jaci;
			   int Nzp =a_EW->m_global_nz[m_grid+1];
			   int kk = Nzp - 1 + k;
			   GridPointSource* sourcePtr = new GridPointSource( mFreq, mT0, i, j, kk, m_grid+1, 
									  fx, fy, fz, mTimeDependence, mNcyc,
									  mPar, mNpar, mIpar, mNipar,
									  dsdp, dddp, dh1, dh2, dh3 );
			   point_sources.push_back(sourcePtr);
			}
		     }
	       }
	       if( k >= Nz )
	       {
		  float_sw4 jaci = 1.0/(8*h*h*h);
		  float_sw4 hi = 1.0/(2*h);
		  for( int j=jcref-2 ; j <= jcref+3 ; j++ )
		     for( int i=icref-2 ; i <= icref+3 ; i++ )
		     {
			float_sw4 wFx=0, wFy=0, wFz=0;
			if( a_EW->interior_point_in_proc(i,j,m_grid-1) ) 
			{
			   wFx = dwghiref[i-icref+2]* wghjref[j-jcref+2]* wghkref[k-kc+2]*hi;
			   wFy =  wghiref[i-icref+2]*dwghjref[j-jcref+2]* wghkref[k-kc+2]*hi;
			   wFz =  wghiref[i-icref+2]* wghjref[j-jcref+2]*dwghkref[k-kc+2]*2*hi;
			   float_sw4 fx = -(mForces[0]*wFx+mForces[1]*wFy+mForces[2]*wFz)*jaci;
			   float_sw4 fy = -(mForces[1]*wFx+mForces[3]*wFy+mForces[4]*wFz)*jaci;
			   float_sw4 fz = -(mForces[2]*wFx+mForces[4]*wFy+mForces[5]*wFz)*jaci;
			   int kk = k-Nz + 1;
			   GridPointSource* sourcePtr = new GridPointSource( mFreq, mT0, i, j, kk, m_grid-1, 
									  fx, fy, fz, mTimeDependence, mNcyc,
									  mPar, mNpar, mIpar, mNipar,
									  dsdp, dddp, dh1, dh2, dh3 );
			   point_sources.push_back(sourcePtr);
			}
		     }
		  
	       }
	    }
	 }
      }
   }
}


//-----------------------------------------------------------------------
void Source::exact_testmoments( int kx[3], int ky[3], int kz[3], float_sw4 momex[3] )
{
   // Integrals over the domain of a polynomial of degree (kx,ky,kz) times the source
   if( !mIsMomentSource )
   {
      float_sw4 x1, y1, z1;
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
      float_sw4 x1, y1, z1, xp1, yp1, zp1;
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
void Source::perturb( float_sw4 h, int comp )
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
Source* Source::copy( std::string a_name )
{
   if( a_name == " " )
      a_name = mName;

   Source* retval = new Source();
   retval->m_i0 = m_i0;
   retval->m_j0 = m_j0;
   retval->m_k0 = m_k0;
   retval->mQ0 = mQ0;
   retval->mR0 = mR0;
   retval->mS0 = mS0;
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
//   retval->mGridPointSet = mGridPointSet;
   retval->m_myPoint = m_myPoint;
   retval->m_zRelativeToTopography = m_zRelativeToTopography;
   retval->mX0 = mX0;
   retval->mY0 = mY0;
   retval->mZ0 = mZ0;

   retval->mNpar = mNpar;
   retval->mPar = new float_sw4[mNpar];
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

   retval->m_zTopo = m_zTopo;
   retval->mIgnore = mIgnore;
   retval->mShearModulusFactor = mShearModulusFactor;

   return retval;
}

//-----------------------------------------------------------------------
float_sw4 Source::find_min_exponent() const
{
   // smallest number x, such that exp(x) does not cause underflow
  return -700.0;
}


//-----------------------------------------------------------------------
void Source::compute_metric_at_source( EW* a_EW, float_sw4 q, float_sw4 r, float_sw4 s, int ic,
				       int jc, int kc, int g, float_sw4& zq, float_sw4& zr,
				       float_sw4& zs, float_sw4& zqq, float_sw4& zqr, float_sw4& zqs,
				       float_sw4& zrr, float_sw4& zrs, float_sw4& zss ) const
{
   int Nz = a_EW->m_global_nz[g];
   float_sw4 h = a_EW->mGridSize[g];
   if( kc > Nz-3 )
   {
      // Treat downmost grid lines as cartesian
      zq = zr = 0;
      zs = h;
      zqq = zqr = zqs = zrr = zrs = zss = 0;
   }
   else
   {
      bool analytic_derivative = true;
      // Derivative of metric wrt. source position. Not yet fully implemented.
      //      double zqdx0, zqdy0, zqsz0, zrdx0, zrdy0, zrdz0;
	       
      // 3. Recompute metric to sixth order accuracy. Increased accuracy needed because
      //    of the multiplication with a singular (Dirac) function.
      //    compute only in the processor where the point is interior

      //      bool eightptstencil = true;
      float_sw4 ai=q-ic;
      float_sw4 bi=r-jc;
      float_sw4 ci=s-kc;

      //      double d6cofi[8], d6cofj[8], d6cofk[8];
      //      double dd6cofi[8], dd6cofj[8], dd6cofk[8];
      //      if( eightptstencil )
      //      {
      //	 // Eight point stencil, smooth wrt. source position, 
      //	 // needed for source optimization
      //	 getmetdwgh( ai, d6cofi );
      //	 getmetdwgh( bi, d6cofj );
      //	 if( kc <= 3 && ci < 0 )
      //	 {
      //	    getmetdwgh7( ci, d6cofk );
      ///	    d6cofk[7] = 0;
      //	 }
      //	 else
      //	    getmetdwgh( ci, d6cofk );
      //      }
      //      else
      //      {
      //	 // Seven point stencil, ok for forward solver.
      //	 getmetdwgh7( ai, d6cofi );
      //	 d6cofi[7] = 0;
      //	 getmetdwgh7( bi, d6cofj );
      //	 d6cofj[7] = 0;
      //	 getmetdwgh7( ci, d6cofk );
      //	 d6cofk[7] = 0;
      //      }

      float_sw4 a6cofi[8], a6cofj[8], a6cofk[8];
      float_sw4 d6cofi[8], d6cofj[8], d6cofk[8];
      float_sw4 dd6cofi[8], dd6cofj[8], dd6cofk[8];
      float_sw4 ddd6cofi[8], ddd6cofj[8], ddd6cofk[8];
      //      if( eightptstencil )
      //      {
	 getmetwgh( ai, a6cofi, d6cofi, dd6cofi, ddd6cofi );
	 getmetwgh( bi, a6cofj, d6cofj, dd6cofj, ddd6cofj );
	 getmetwgh( ci, a6cofk, d6cofk, dd6cofk, ddd6cofk );
      //	 if( kc <= 3 && ci < 0 )
      //	 {
      //	    getmetwgh7( ci, a6cofk );
      //	    a6cofk[7] = 0;
      //	 }
      //	 else
      //	    getmetwgh( ci, a6cofk );
   //      }
   //      else
   //      {
   //	 getmetwgh7( ai, a6cofi );
   //	 a6cofi[7] = 0;
   //	 getmetwgh7( bi, a6cofj );
   //	 a6cofj[7] = 0;
   //	 getmetwgh7( ci, a6cofk );
   //	 a6cofk[7] = 0;
   //      }

      // Assume grid uniform in x and y, compute metric with z=z(q,r,s)
      zq = zr = zs = 0;
      int order=a_EW->m_grid_interpolation_order;
      float_sw4 zetaBreak=a_EW->m_zetaBreak;
      
      float_sw4 zpar = (s-1)/(zetaBreak*(Nz-1));
      float_sw4 kBreak = 1 + zetaBreak*(Nz-1);

      if( zpar >= 1 )
      {
	 zq = 0;
	 zr = 0;
         zs = h;
	 zqq = zqr = zqs = zrr = zrs = zss = 0;
      }
      else 
      {
	 float_sw4 pp = pow(1-zpar,order-1);
	 float_sw4 powo = (1-zpar)*pp;
	 float_sw4 dpowo = -order*pp/zetaBreak;
	 float_sw4 tauavg = 0;
	 float_sw4 tauq=0, taur=0;
	 float_sw4 tauqq=0, tauqr=0, taurr=0;
	 for( int j=jc-3; j <= jc+4 ; j++ )
	    for( int i=ic-3; i <= ic+4 ; i++ )
	    {
	       tauavg += a6cofi[i-(ic-3)]* a6cofj[j-(jc-3)]*a_EW->mTopoGridExt(i,j,1);
	       tauq  +=  d6cofi[i-(ic-3)]* a6cofj[j-(jc-3)]*a_EW->mTopoGridExt(i,j,1);
	       taur  +=  a6cofi[i-(ic-3)]* d6cofj[j-(jc-3)]*a_EW->mTopoGridExt(i,j,1);
               tauqq += dd6cofi[i-(ic-3)]* a6cofj[j-(jc-3)]*a_EW->mTopoGridExt(i,j,1);
	       tauqr +=  d6cofi[i-(ic-3)]* d6cofj[j-(jc-3)]*a_EW->mTopoGridExt(i,j,1);
               taurr +=  a6cofi[i-(ic-3)]*dd6cofj[j-(jc-3)]*a_EW->mTopoGridExt(i,j,1);
	    }
	 //	 double powo = pow(1-zpar,order);
	 //      double dpowo = -order*pow(1-zpar,order-1);
	 zq = (-tauq)*powo;
	 zr = (-taur)*powo;
	 zqq = (-tauqq)*powo;
	 zqr = (-tauqr)*powo;
	 zrr = (-taurr)*powo;
         zqs = (-tauq)*dpowo/(Nz-1);
         zrs = (-taur)*dpowo/(Nz-1);

      // Compute dz/ds directly from the grid mapping, the explicit expression here
      // should be the same as in EW::curvilinear_grid_mapping
	 float_sw4 zMax = a_EW->m_zmin[a_EW->mNumberOfCartesianGrids-1] - (Nz-kBreak)*h;
	 float_sw4 c1 = zMax + tauavg - h*(kBreak-1);

      // Divide by Nz-1 to make consistent with undivided differences
         if( analytic_derivative )
	 {
	    zs  = h + c1*(-dpowo)/(Nz-1);
	    zss = -c1*order*(order-1)*pow(1-zpar,order-2)/(zetaBreak*zetaBreak*(Nz-1)*(Nz-1));
	 //         cout << "AN: zs = " << zs << " zss= " << zss << endl;
	 }
	 else
	 {
	    zs = 0;
	    zss= 0;
	    float_sw4 z1d=0;
	    for( int k=kc-3 ; k <= kc+4; k++ ) 
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
	       zs  += d6cofk[k-(kc-3)]*z1d;
	       zss += dd6cofk[k-(kc-3)]*z1d;
	    } 
	 //         cout << "NU: zs = " << zs << " zss= " << zss << endl;
	 }
      }
   }
}
