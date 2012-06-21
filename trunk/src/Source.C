#include "mpi.h"
#include "GridPointSource.h"
#include "Source.h"
#include "Require.h"

#include <fenv.h>
#include <cmath>

#include  "EW.h"

using namespace std;

#define SQR(x) ((x)*(x))

Source::Source(EW *a_wpp, 
	       double amplitude, 
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
	       int ncyc):
//  mEW(a_wpp),
  mAmp(amplitude),
  mIsMomentSource(true),
  mFreq(frequency),
  mT0(t0),
  m_zRelativeToTopography(false),
  mX0(x0), mY0(y0), mZ0(z0),
  mGridPointSet(false),
  mTimeDependence(tDep),
  mIgnore(false),
  mNcyc(ncyc),
  m_derivative(-1)
{
  mForces.resize(6);
  mForces[0] = Mxx;
  mForces[1] = Mxy;
  mForces[2] = Mxz;
  mForces[3] = Myy;
  mForces[4] = Myz;
  mForces[5] = Mzz;
  mPar  = new double[2];
  mName = name;

  a_wpp->computeNearestGridPoint(m_i0,m_j0,m_k0,m_grid,mX0,mY0,mZ0);
}


Source::Source(EW *a_wpp,double amplitude, double frequency, double t0,
	       double x0, double y0, double z0,
	       double Fx,
	       double Fy,
	       double Fz,
	       timeDep tDep,
	       const char *name, int ncyc):
//  mEW(a_wpp),
  mAmp(amplitude),
  mIsMomentSource(false),
  mFreq(frequency),
  mT0(t0),
  m_zRelativeToTopography(false),
  mX0(x0), mY0(y0), mZ0(z0),
  mGridPointSet(false),
  mTimeDependence(tDep),
  mIgnore(false),
  mNcyc(ncyc),
  m_derivative(-1)
{
  mPar = new double[2];
  mForces.resize(3);
  mForces[0] = Fx;
  mForces[1] = Fy;
  mForces[2] = Fz;
  mName = name;

  a_wpp->computeNearestGridPoint(m_i0,m_j0,m_k0,m_grid,mX0,mY0,mZ0);
}


Source::~Source()
{
  //   if( mPar != NULL )
  delete[] mPar;
}

double Source::getX0() const
{
  return mX0;
}

double Source::getY0() const
{
  return mY0;
}

double Source::getZ0() const
{
  return mZ0;
}

bool Source::ignore() const
{
  return mIgnore;
}

double Source::getOffset() const
{
  return mT0;
}

double Source::getFrequency() const
{
  return mFreq;
}

void Source::setMaxFrequency(double max_freq)
{
  if (mFreq > max_freq)
    mFreq=max_freq;
}

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
void Source::getMoments( double& mxx, double& myy, double& mzz, double& mxy, double& mxz, double& myz ) const
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
double Source::getAmplitude() const
{
  double amplitude=0;
  if (mIsMomentSource)
  {
    double msqr=0;
    for (int q=0; q<6; q++)
      msqr += SQR(mForces[q]);
    amplitude = mAmp*sqrt(msqr/2.);
  }
  else
  {
    double fsqr=0;
    for (int q=0; q<3; q++)
      fsqr += SQR(mForces[q]);
    amplitude = mAmp*sqrt(fsqr);
  }
  
  return amplitude;
}

//-----------------------------------------------------------------------
ostream& operator<<( ostream& output, const Source& s )
{
  output << s.mName << (s.isMomentSource()? " moment":" force") << " source term" << endl;
   output << "   Location (X,Y,Z) = " << s.mX0 << "," << s.mY0 << "," << s.mZ0 << " in grid no " << s.m_grid << endl;
   output << "   Strength " << s.mAmp;
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
       mTimeDependence == iLiu )
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
double Source::compute_t0_increase(double t0_min)
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
   if( dt0 > 0 )
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
void Source::get_parameters( double x[11] )
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
void Source::getsourcewgh( double ai, double wgh[6], double dwghda[6], double ddwghda[6] )
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
void Source::getsourcedwgh( double ai, double wgh[6], double dwghda[6], double ddwghda[6] )
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
void Source::set_grid_point_sources4( EW *a_EW, vector<GridPointSource*>& point_sources )
{
   int i,j,k,g;
   a_EW->computeNearestGridPoint( i, j, k, g, mX0, mY0, mZ0 );
   double q, r, s;
   double h = a_EW->mGridSize[g];
   double normwgh[4]={17.0/48.0, 59.0/48.0, 43.0/48.0, 49.0/48.0 };

// Cartesian case
   q = mX0/h+1;
   r = mY0/h+1;
   s = (mZ0-a_EW->m_zmin[g])/h+1;

   int Ni = a_EW->m_global_nx[g];
   int Nj = a_EW->m_global_ny[g];
   int Nz = a_EW->m_global_nz[g];

   int ic = static_cast<int>(floor(q));
   int jc = static_cast<int>(floor(r));
   int kc = static_cast<int>(floor(s));

// Bias stencil away from boundary, no source at ghost/padding points
   if( ic <= 2 )    ic = 3;
   if( ic >= Ni-2 ) ic = Ni-3;
   if( jc <= 2 )    jc = 3;
   if( jc >= Nj-2 ) jc = Nj-3;
   if( kc <= 2 )    kc = 3;
   if( kc >= Nz-2 ) kc = Nz-3;

   double ai=q-ic, bi=r-jc, ci=s-kc;
// Delta distribution
   double wghi[6], wghj[6], wghk[6], wghix[6], wghjy[6], wghkz[6];
   double wghixx[6], wghjyy[6], wghkzz[6];
   getsourcewgh( ai, wghi, wghix, wghixx );
   getsourcewgh( bi, wghj, wghjy, wghjyy );
   getsourcewgh( ci, wghk, wghkz, wghkzz );

// Delta' distribution
   double dwghi[6], dwghj[6], dwghk[6], dwghix[6], dwghjy[6], dwghkz[6];
   double dwghixx[6], dwghjyy[6], dwghkzz[6];
   getsourcedwgh( ai, dwghi, dwghix, dwghixx );
   getsourcedwgh( bi, dwghj, dwghjy, dwghjyy );
   getsourcedwgh( ci, dwghk, dwghkz, dwghkzz );

// Boundary correction
   for( int k=0 ; k <= 5 ; k++ )
   {
      //      if( ( 1 <= k+ic-2) && ( k+ic-2 <= 4 ) )
      //      {
      //         wghi[k]  /= normwgh[k+ic-3];
      //         dwghi[k] /= normwgh[k+ic-3];
      //      }
      //      if( ( Ni-3 <= k+ic-2) && ( k+ic-2 <= Ni ) )
      //      {
      //         wghi[k]  /= normwgh[Ni-k-ic+2];
      //         dwghi[k] /= normwgh[Ni-k-ic+2];
      //      }
      //      if( ( 1 <= k+jc-2) && ( k+jc-2 <= 4 ) )
      //      {
      //         wghj[k]  /= normwgh[k+jc-3];
      //         dwghj[k] /= normwgh[k+jc-3];
      //      }
      //      if( ( Nj-3 <= k+jc-2) && ( k+jc-2 <= Nj ) )
      //      {
      //         wghj[k]  /= normwgh[Nj-k-jc+2];
      //         dwghj[k] /= normwgh[Nj-k-jc+2];
      //      }
      if( ( 1 <= k+kc-2) && ( k+kc-2 <= 4 ) )
      {
         wghk[k]  /= normwgh[k+kc-3];
         dwghk[k] /= normwgh[k+kc-3];
	 wghkz[k] /= normwgh[k+kc-3];
	 dwghkz[k] /= normwgh[k+kc-3];
      }
      if( ( Nz-3 <= k+kc-2) && ( k+kc-2 <= Nz ) )
      {
         wghk[k]  /= normwgh[Nz-k-kc+2];
         dwghk[k] /= normwgh[Nz-k-kc+2];
         wghkz[k] /= normwgh[Nz-k-kc+2];
         dwghkz[k] /= normwgh[Nz-k-kc+2];
      }
   }
   if( !mIsMomentSource )
   {
      for( int k=kc-2 ; k <= kc+3 ; k++ )
	 for( int j=jc-2 ; j <= jc+3 ; j++ )
	    for( int i=ic-2 ; i <= ic+3 ; i++ )
	    {
	       double wF = wghi[i-ic+2]*wghj[j-jc+2]*wghk[k-kc+2];
	       if( (wF*mAmp != 0) && (mForces[0] != 0 || mForces[1] != 0 || mForces[2] != 0) 
		   && a_EW->interior_point_in_proc(i,j,g) )
	       {
		  wF /= h*h*h;
		  {
		     GridPointSource* sourcePtr = new GridPointSource(
								      mAmp*wF, mFreq, mT0,
								      i,j,k,g,
								      mForces[0], mForces[1], mForces[2],
								      mTimeDependence, mNcyc );
		     point_sources.push_back(sourcePtr);
		  }
	       }
	    }
   }
   else
   {
      for( int k=kc-2 ; k <= kc+3 ; k++ )
	 for( int j=jc-2 ; j <= jc+3 ; j++ )
	    for( int i=ic-2 ; i <= ic+3 ; i++ )
	    {
	       double wFx=0, wFy=0, wFz=0, dsdp[27];
	       if( a_EW->interior_point_in_proc(i,j,g) ) 
	       {
		  wFx = dwghi[i-ic+2]* wghj[j-jc+2]* wghk[k-kc+2];
		  wFy =  wghi[i-ic+2]*dwghj[j-jc+2]* wghk[k-kc+2];
		  wFz =  wghi[i-ic+2]* wghj[j-jc+2]*dwghk[k-kc+2];
                  double hi=1.0/h;
                  double hi2=hi*hi;
                  double wFxdx0 =dwghix[i-ic+2]*  wghj[j-jc+2]*  wghk[k-kc+2]*hi;
                  double wFxdy0 = dwghi[i-ic+2]* wghjy[j-jc+2]*  wghk[k-kc+2]*hi;
                  double wFxdz0 = dwghi[i-ic+2]*  wghj[j-jc+2]* wghkz[k-kc+2]*hi;
                  double wFydx0 = wghix[i-ic+2]* dwghj[j-jc+2]*  wghk[k-kc+2]*hi;
                  double wFydy0 =  wghi[i-ic+2]*dwghjy[j-jc+2]*  wghk[k-kc+2]*hi;
                  double wFydz0 =  wghi[i-ic+2]* dwghj[j-jc+2]* wghkz[k-kc+2]*hi;
                  double wFzdx0 = wghix[i-ic+2]*  wghj[j-jc+2]* dwghk[k-kc+2]*hi;
                  double wFzdy0 =  wghi[i-ic+2]* wghjy[j-jc+2]* dwghk[k-kc+2]*hi;
                  double wFzdz0 =  wghi[i-ic+2]*  wghj[j-jc+2]*dwghkz[k-kc+2]*hi;

		  // Second derivatives
                  double wFxdx0dx0 = dwghixx[i-ic+2]*  wghj[j-jc+2]*  wghk[k-kc+2]*hi2;
                  double wFxdx0dy0 = dwghix[i-ic+2]*  wghjy[j-jc+2]*  wghk[k-kc+2]*hi2;
                  double wFxdx0dz0 = dwghix[i-ic+2]*  wghj[j-jc+2]*  wghkz[k-kc+2]*hi2;
		  double wFxdy0dy0 = dwghi[i-ic+2]* wghjyy[j-jc+2]*  wghk[k-kc+2]*hi2;
                  double wFxdy0dz0 = dwghi[i-ic+2]* wghjy[j-jc+2]*  wghkz[k-kc+2]*hi2;
		  double wFxdz0dz0 = dwghi[i-ic+2]*  wghj[j-jc+2]* wghkzz[k-kc+2]*hi2;

                  double wFydx0dx0 = wghixx[i-ic+2]* dwghj[j-jc+2]*  wghk[k-kc+2]*hi2;
                  double wFydx0dy0 = wghix[i-ic+2]* dwghjy[j-jc+2]*  wghk[k-kc+2]*hi2;
                  double wFydx0dz0 = wghix[i-ic+2]* dwghj[j-jc+2]*  wghkz[k-kc+2]*hi2;
                  double wFydy0dy0 = wghi[i-ic+2]*dwghjyy[j-jc+2]*  wghk[k-kc+2]*hi2;
                  double wFydy0dz0 = wghi[i-ic+2]*dwghjy[j-jc+2]*  wghkz[k-kc+2]*hi2;
                  double wFydz0dz0 = wghi[i-ic+2]* dwghj[j-jc+2]* wghkzz[k-kc+2]*hi2;

                  double wFzdx0dx0 = wghixx[i-ic+2]*  wghj[j-jc+2]* dwghk[k-kc+2]*hi2;
                  double wFzdx0dy0 = wghix[i-ic+2]*  wghjy[j-jc+2]* dwghk[k-kc+2]*hi2;
                  double wFzdx0dz0 = wghix[i-ic+2]*  wghj[j-jc+2]* dwghkz[k-kc+2]*hi2;
                  double wFzdy0dy0 = wghi[i-ic+2]* wghjyy[j-jc+2]* dwghk[k-kc+2]*hi2;
                  double wFzdy0dz0 = wghi[i-ic+2]* wghjy[j-jc+2]* dwghkz[k-kc+2]*hi2;
                  double wFzdz0dz0 = wghi[i-ic+2]*  wghj[j-jc+2]*dwghkzz[k-kc+2]*hi2;

		  double jaci = 1.0/(h*h*h*h);

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
		  {
		     GridPointSource* sourcePtr = new GridPointSource( mAmp, mFreq, mT0, i, j, k, g, 
								       fx, fy, fz, mTimeDependence, mNcyc,
								       dsdp, dddp, dh1, dh2, dh3 );
                     if( m_derivative >= 0 )
			sourcePtr->set_derivative(m_derivative,m_dir);
		     point_sources.push_back(sourcePtr);
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
      for( int c = 0; c <= 3 ; c++ )
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
      for( int c = 0; c <= 3 ; c++ )
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
void Source::set_grid_point_sources( EW *a_EW, vector<GridPointSource*>& point_sources )
{
   int i,j,k,g;
   a_EW->computeNearestGridPoint( i, j, k, g, mX0, mY0, mZ0 );
   double q, r, s;
   double h = a_EW->mGridSize[g];
   bool canBeInverted, curvilinear;

   if( g == a_EW->mNumberOfGrids-1 && a_EW->topographyExists() )
   {
// Curvilinear
      canBeInverted = a_EW->invert_curvilinear_grid_mapping( mX0, mY0, mZ0, q, r, s );
// if s < 0, the source is located above the grid and the call to
// find_curvilinear_derivatives_at_point will fail
      if (canBeInverted && s<0.)
      {
	double xTop, yTop, zTop;
	a_EW->curvilinear_grid_mapping(q, r, 0., xTop, yTop, zTop);
	double lat, lon;
	a_EW->computeGeographicCoord(mX0, mY0, lon, lat);
	printf("Found a source above the curvilinear grid! Lat=%e, Lon=%e, source Z-level = %e, grid boundary Z = %e\n", lat, lon, mZ0, zTop);
		
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
   if( canBeInverted )
   {

      int Ni = a_EW->m_global_nx[g];
      int Nj = a_EW->m_global_ny[g];
      int Nz = a_EW->m_global_nz[g];
//      int Nz = a_EW->m_kEnd[g]- a_EW->m_ghost_points;

      int ic3 = static_cast<int>(round(q));
      int jc3 = static_cast<int>(round(r));
   // Bias stencil away from boundary
      if( ic3 <= 1 ) ic3 = 2;
      if( ic3 >= Ni ) ic3 = Ni-1;
      if( jc3 <= 1 ) jc3 = 2;
      if( jc3 >= Nj ) jc3 = Nj-1;

      int ic4 = static_cast<int>(floor(q));
      int jc4 = static_cast<int>(floor(r));
      // Bias stencil away from boundary
      if( ic4 <= 1 ) ic4 = 2;
      if( ic4 >= Ni-1 ) ic4 = Ni-2;
      if( jc4 <= 1 ) jc4 = 2;
      if( jc4 >= Nj-1 ) jc4 = Nj-2;

      int kc3 = static_cast<int>(round(s));
      int kc4 = static_cast<int>(floor(s));

      // if kc4=Nz, point is on boundary, move stencil to interior
      if( kc4 >= Nz ) kc4 = Nz-1;

      bool upperbndry = (kc3 == 1  || kc3 == 2);
      bool lowerbndry = (kc3 == Nz || kc3 == Nz-1);
      if( mIsMomentSource )
      {
         upperbndry  = upperbndry || kc4 == 1    || kc4 == 2;
         lowerbndry  = lowerbndry || kc4 == Nz-1 || kc4 == Nz-2;
      }

   // ccbndry : source on Cartesian/Curvliniear boundary
      bool ccbndry  =  (upperbndry && g == a_EW->mNumberOfGrids-2 && a_EW->topographyExists()) ||
	 (lowerbndry && g == a_EW->mNumberOfGrids-1 && a_EW->topographyExists());

   // refbndry : source on Cartesian/Cartesian grid refinement bndry
      bool refbndry = (upperbndry && g < a_EW->mNumberOfGrids-1 && !ccbndry) ||
	 (lowerbndry && g>0 && !ccbndry );
   
      if( !refbndry && !ccbndry )
      {
	 // If not interpolation boundary, bias stencil away from boundary, 
	 if( kc3 <= 2 ) kc3 = 2;
	 if( kc3 >= Nz ) kc3 = Nz-1;
	 if( kc4 <= 2 ) kc4 = 2;
	 if( kc4 >= Nz-1 ) kc4 = Nz-2;
      }
      //      cout << "refbndry " << refbndry << endl;
      //      cout << "grid     " << g << endl;
      //      cout << "kc3      " << kc3 << endl;
      //      cout << "Nz       " << Nz << endl;
      //      cout << "beta     " << s-kc3 << endl;

      if( refbndry )
      {
	// Source on grid refinement boundary 
	 if( !mIsMomentSource )
	 {
       // Point source
	    //            cout << "Point boundary g= " << g << " kc3 = " << kc3 << endl;
	    double ci=s-kc3;
	    double wghk[3];
	    if( kc3 == 2 || kc3 == Nz-1 )
	    {
	       wghk[0] = 0.5*(ci*ci-ci);
	       wghk[1] = (1-ci*ci);
	       wghk[2] = 0.5*(ci*ci+ci);
	       for( int k= kc3-1 ; k<= kc3+1 ; k++ )
		 distribute_source_xyplane( a_EW, point_sources, g, k, wghk[k-kc3+1] );
	       if( kc3 == 2 )
	       {
//		  int Nzp = a_EW->m_kEnd[g+1]- a_EW->m_ghost_points;
		 int Nzp = a_EW->m_global_nz[g+1];
		 distribute_source_xyplane( a_EW, point_sources, g+1, Nzp, wghk[0] );
	       }
	       if( kc3 == Nz-1 )
		 distribute_source_xyplane( a_EW, point_sources, g-1, 1, wghk[2] );
	    }
	    else if( kc3 == 1 )
	    {
               wghk[0] = (4*ci*ci-4*ci)/3;
	       wghk[1] = (1+ci-2*ci*ci);
	       wghk[2] = (ci*ci*4+ci*2)/6;
	       for( int k= 1 ; k<= 2 ; k++ )
		 distribute_source_xyplane( a_EW, point_sources, g, k, wghk[k] );
//	       int Nzp = a_EW->m_kEnd[g+1]- a_EW->m_ghost_points;
	       int Nzp = a_EW->m_global_nz[g+1];
	       for( int k= Nzp-1 ; k<= Nzp ; k++ )
		 distribute_source_xyplane( a_EW, point_sources, g+1, k, wghk[k-Nzp+1] );
	    }
	    else if( kc3 == Nz )
	    {
	       wghk[0] = (ci*ci-2*ci)/3;
	       wghk[1] = (1-0.5*ci*ci+0.5*ci);
	       wghk[2] = (ci*ci+ci)/6;
	       for( int k= Nz-1 ; k<= Nz ; k++ )
		 distribute_source_xyplane( a_EW, point_sources, g, k, wghk[k-Nz+1] );
	       //	       wghk[1] = wghk[1]/2;
	       for( int k= 1 ; k<= 2 ; k++ )
		 distribute_source_xyplane( a_EW, point_sources, g-1, k, wghk[k] );
	    }
	 }
	 else
	 {
 // Moment source on grid refinement boundary
	    double ci=s-kc3;
	    double wghk[3];
            double cid = s-kc4;
	    double dwghk[4];
            double wgz, dwgz;

	    // 8 different cases:
	    if( kc3 == Nz-2 && kc4 == Nz-2 )
	    {
	       wghk[0] = 0.5*(ci*ci-ci);
	       wghk[1] = (1-ci*ci);
	       wghk[2] = 0.5*(ci*ci+ci);

               dwghk[0] = 1.0/3-cid+0.5*cid*cid;
               dwghk[1] = 0.5+2*cid-1.5*cid*cid;
               dwghk[2] = -1-cid+1.5*cid*cid;
	       dwghk[3] = 1.0/6 -0.5*cid*cid;
               for( int k = kc4-1 ; k<= Nz ; k++ )
	       {
                  if( kc3-1 <= k )
		     wgz = wghk[k-kc3+1];
		  else
		     wgz = 0;
		  dwgz = dwghk[k-kc4+1];
                  distribute_source_xyplane_mom( a_EW, point_sources, g, k, wgz, dwgz );
	       }
               wgz = 0;
	       dwgz = dwghk[3]*2;
	       distribute_source_xyplane_mom( a_EW, point_sources, g-1, 1, wgz, dwgz );
	    }
            else if( kc3 == Nz-1 && kc4 == Nz-2 )
	    {
	       wghk[0] = 0.5*(ci*ci-ci);
	       wghk[1] = (1-ci*ci);
	       wghk[2] = 0.5*(ci*ci+ci);
  
               dwghk[0] = 1.0/3-cid+0.5*cid*cid;
               dwghk[1] = 0.5+2*cid-1.5*cid*cid;
               dwghk[2] = -1-cid+1.5*cid*cid;
	       dwghk[3] = 1.0/6 -0.5*cid*cid;
               for( int k = kc4-1 ; k<= Nz ; k++ )
	       {
                  if( kc3-1 <= k )
		     wgz = wghk[k-kc3+1];
		  else
		     wgz = 0;
		  dwgz = dwghk[k-kc4+1];
                  distribute_source_xyplane_mom( a_EW, point_sources, g, k, wgz, dwgz );
	       }
               wgz = wghk[2];
	       dwgz = dwghk[3]*2;
	       distribute_source_xyplane_mom( a_EW, point_sources, g-1, 1, wgz, dwgz );
	    }
            else if( kc3 == Nz-1 && kc4 == Nz-1 )
	    {
	       wghk[0] = 0.5*(ci*ci-ci);
	       wghk[1] = (1-ci*ci);
	       wghk[2] = 0.5*(ci*ci+ci);

               dwghk[0] = 0.375-cid+0.375*cid*cid;
               dwghk[1] = 1.0/3 + 2*cid - cid*cid;
               dwghk[2] = -0.75-cid+0.75*cid*cid;
	       dwghk[3] = 1.0/12-cid*cid*0.25;
               for( int k = kc4-1 ; k<= Nz ; k++ )
	       {
                  if( kc3-1 <= k )
		     wgz = wghk[k-kc3+1];
		  else
		     wgz = 0;
		  dwgz = dwghk[k-kc4+1];
                  distribute_source_xyplane_mom( a_EW, point_sources, g, k, wgz, dwgz );
	       }
               wgz = wghk[2];               
	       dwgz = dwghk[2]*2;
	       distribute_source_xyplane_mom( a_EW, point_sources, g-1, 1, wgz, dwgz );
               wgz = 0;
	       dwgz = dwghk[3];
	       distribute_source_xyplane_mom( a_EW, point_sources, g-1, 2, wgz, dwgz );
	    }
            else if( kc3 == Nz && kc4 == Nz-1 )
	    {
	       wghk[0] = (ci*ci-2*ci)/3;
	       wghk[1] = (1-0.5*ci*ci+0.5*ci);
	       wghk[2] = (ci*ci+ci)/6;

               dwghk[0] = 0.375-cid+0.375*cid*cid;
               dwghk[1] = 1.0/3 + 2*cid - cid*cid;
               dwghk[2] = -0.75-cid+0.75*cid*cid;
	       dwghk[3] = 1.0/12-cid*cid*0.25;
               for( int k = kc4-1 ; k<= Nz ; k++ )
	       {
                  if( kc3-1 <= k )
		     wgz = wghk[k-kc3+1];
		  else
		     wgz = 0;
		  dwgz = dwghk[k-kc4+1];
                  distribute_source_xyplane_mom( a_EW, point_sources, g, k, wgz, dwgz );
	       }
               wgz = wghk[1];               
	       dwgz = dwghk[2]*2;
	       distribute_source_xyplane_mom( a_EW, point_sources, g-1, 1, wgz, dwgz );
               wgz = wghk[2];
	       dwgz = dwghk[3];
	       distribute_source_xyplane_mom( a_EW, point_sources, g-1, 2, wgz, dwgz );

	    }
            else if( kc3 == 1 && kc4 == 1 )
	    {
               wghk[0] = (4*ci*ci-4*ci)/3;
	       wghk[1] = (1+ci-2*ci*ci);
	       wghk[2] = (ci*ci*4+ci*2)/6;

               dwghk[0] = 8.0/15-1.6*cid+0.8*cid*cid;
               dwghk[1] = -0.5+5*cid-3*cid*cid;
               dwghk[2] = -2.0/3-2*cid+2*cid*cid;
	       dwghk[3] = 0.1+0.2*cid-0.6*cid*cid;

               for( int k=1 ; k <= kc4+2 ; k++ )
	       {
                  if( k <= kc3+1 )
		     wgz = wghk[k-kc3+1];
		  else
		     wgz = 0;
		  dwgz = dwghk[k-kc4+1];
                  distribute_source_xyplane_mom( a_EW, point_sources, g, k, wgz, dwgz );
	       }
//	       int Nzp = a_EW->m_kEnd[g+1]- a_EW->m_ghost_points;
	       int Nzp = a_EW->m_global_nz[g+1];
               wgz  = wghk[0];
	       dwgz = dwghk[0];
               distribute_source_xyplane_mom( a_EW, point_sources, g+1, Nzp-1, wgz, dwgz );
               wgz  = wghk[1];
	       dwgz = dwghk[1]/2;
               distribute_source_xyplane_mom( a_EW, point_sources, g+1, Nzp, wgz, dwgz );
	    }
            else if( kc3 == 2 && kc4 == 1 )
	    {
	       wghk[0] = 0.5*(ci*ci-ci);
	       wghk[1] = (1-ci*ci);
	       wghk[2] = 0.5*(ci*ci+ci);

               dwghk[0] = 8.0/15-1.6*cid+0.8*cid*cid;
               dwghk[1] = -0.5+5*cid-3*cid*cid;
               dwghk[2] = -2.0/3-2*cid+2*cid*cid;
	       dwghk[3] = 0.1+0.2*cid-0.6*cid*cid;

               for( int k=1 ; k <= kc4+2 ; k++ )
	       {
                  if( k <= kc3+1 )
		     wgz = wghk[k-kc3+1];
		  else
		     wgz = 0;
		  dwgz = dwghk[k-kc4+1];
                  distribute_source_xyplane_mom( a_EW, point_sources, g, k, wgz, dwgz );
	       }
//	       int Nzp = a_EW->m_kEnd[g+1]- a_EW->m_ghost_points;
	       int Nzp = a_EW->m_global_nz[g+1];
               wgz  = 0;
	       dwgz = dwghk[0];
               distribute_source_xyplane_mom( a_EW, point_sources, g+1, Nzp-1, wgz, dwgz );
               wgz  = wghk[0];
	       dwgz = dwghk[1]/2;
               distribute_source_xyplane_mom( a_EW, point_sources, g+1, Nzp, wgz, dwgz );
	    }
            else if( kc3 == 2 && kc4 == 2 )
	    {
	       wghk[0] = 0.5*(ci*ci-ci);
	       wghk[1] = (1-ci*ci);
	       wghk[2] = 0.5*(ci*ci+ci);

               dwghk[0] = 1.0/3-cid+0.5*cid*cid;
               dwghk[1] = 0.5+2*cid-1.5*cid*cid;
               dwghk[2] = -1-cid+1.5*cid*cid;
	       dwghk[3] = 1.0/6 -0.5*cid*cid;
               for( int k=1 ; k <= kc4+2 ; k++ )
	       {
                  if( k <= kc3+1 )
		     wgz = wghk[k-kc3+1];
		  else
		     wgz = 0;
		  dwgz = dwghk[k-kc4+1];
                  distribute_source_xyplane_mom( a_EW, point_sources, g, k, wgz, dwgz );
	       }
//	       int Nzp = a_EW->m_kEnd[g+1]- a_EW->m_ghost_points;
	       int Nzp = a_EW->m_global_nz[g+1];
               wgz  = wghk[0];
	       dwgz = dwghk[0]/2;
               distribute_source_xyplane_mom( a_EW, point_sources, g+1, Nzp, wgz, dwgz );

	    }
            else if( kc3 == 3 && kc4 == 2 )
	    {
	       wghk[0] = 0.5*(ci*ci-ci);
	       wghk[1] = (1-ci*ci);
	       wghk[2] = 0.5*(ci*ci+ci);

               dwghk[0] = 1.0/3-cid+0.5*cid*cid;
               dwghk[1] = 0.5+2*cid-1.5*cid*cid;
               dwghk[2] = -1-cid+1.5*cid*cid;
	       dwghk[3] = 1.0/6 -0.5*cid*cid;
               for( int k=1 ; k <= kc4+2 ; k++ )
	       {
                  if( k <= kc3+1 )
		     wgz = wghk[k-kc3+1];
		  else
		     wgz = 0;
		  dwgz = dwghk[k-kc4+1];
                  distribute_source_xyplane_mom( a_EW, point_sources, g, k, wgz, dwgz );
	       }
//	       int Nzp = a_EW->m_kEnd[g+1]- a_EW->m_ghost_points;
	       int Nzp = a_EW->m_global_nz[g+1];
               wgz  = 0;
	       dwgz = dwghk[0]/2;
               distribute_source_xyplane_mom( a_EW, point_sources, g+1, Nzp, wgz, dwgz );
	    }
	    else
	    {
	       cout << "Undefined configuration for moment source at a grid refinement boundary"
		    << endl;
	       MPI_Abort(MPI_COMM_WORLD,1);
	    }
	 }
      }
      else
      {
      // Interior to a grid, or on Cartesian/Curvilinear boundary

// Point source:
	 double ai=q-ic3, bi=r-jc3, ci=s-kc3;

// Delta distribution
	 double wghi[3] = {0.5*(ai*ai-ai),1-ai*ai,0.5*(ai*ai+ai)};
	 double wghj[3] = {0.5*(bi*bi-bi),1-bi*bi,0.5*(bi*bi+bi)};
	 double wghk[3] = {0.5*(ci*ci-ci),1-ci*ci,0.5*(ci*ci+ci)};

// boundary correction. Note: on the curvilinear grid, k=1 is always the real
// boundary, and k=Nk is always an interpolation boundary.

	 if( ic3 == 2 ) wghi[0] *= 2;
	 if( ic3 == Ni-1 ) wghi[2] *= 2;
	 if( jc3 == 2 ) wghj[0] *= 2;
	 if( jc3 == Nj-1 ) wghj[2] *= 2;
	 if( kc3 == 2    && !(ccbndry && upperbndry) ) wghk[0] *= 2;
	 if( kc3 == Nz-1 && !(ccbndry && lowerbndry) ) wghk[2] *= 2;

	 if( !mIsMomentSource )
	 {
	    for( int k=kc3-1 ; k <= kc3+1 ; k++ )
	       for( int j=jc3-1 ; j <= jc3+1 ; j++ )
		  for( int i=ic3-1 ; i <= ic3+1 ; i++ )
		  {
		     //		     if( 0 <= k-kc3+1 <= 2 )
		     //		     {
			double wF = wghi[i-ic3+1]*wghj[j-jc3+1]*wghk[k-kc3+1];

			if( (wF*mAmp != 0) && (mForces[0] != 0 || mForces[1] != 0 || mForces[2] != 0) 
			    && a_EW->point_in_proc(i,j,g) )
			{
			   if( curvilinear )
			      wF /= a_EW->mJ(i,j,k);
			   else
			      wF /= h*h*h;

			   if( 1 <= k && k <= Nz )
			   {
			      GridPointSource* sourcePtr = new GridPointSource(
							     mAmp*wF, mFreq, mT0,
							     i,j,k,g,
							     mForces[0], mForces[1], mForces[2],
							     mTimeDependence, mNcyc );
			      point_sources.push_back(sourcePtr);
			   }
			   if( k <= 1 && ccbndry && upperbndry )
			   {
			     int Nzp = a_EW->m_global_nz[g+1];
//			      int Nzp = a_EW->m_kEnd[g+1]- a_EW->m_ghost_points;
			      int kk = Nzp - 1 + k;
			      GridPointSource* sourcePtr = new GridPointSource(
							     mAmp*wF, mFreq, mT0,
							     i,j,kk,g+1,
							     mForces[0], mForces[1], mForces[2],
							     mTimeDependence, mNcyc );
			      point_sources.push_back(sourcePtr);
			   }
			   if( k >= Nz && ccbndry && lowerbndry )
			   {
			      int kk = k - Nz + 1;
			      GridPointSource* sourcePtr = new GridPointSource(
							     mAmp*wF, mFreq, mT0,
							     i,j,kk,g-1,
							     mForces[0], mForces[1], mForces[2],
							     mTimeDependence, mNcyc );
			      point_sources.push_back(sourcePtr);
			   }
			   //			}
		     }
		  }
	 } // end if !momentSource (i.e., point force)
	 else
         {
// Moment source:
	    ai = q-ic4; bi = r-jc4; ci=s-kc4;

// Delta distribution derivatives (this is a 3rd order accurate approx. of the derivative)
	    double dwghi[4] = {1.0/3-ai+0.5*ai*ai, 0.5+2*ai-1.5*ai*ai, -1-ai+1.5*ai*ai, 1.0/6-0.5*ai*ai};
	    double dwghj[4] = {1.0/3-bi+0.5*bi*bi, 0.5+2*bi-1.5*bi*bi, -1-bi+1.5*bi*bi, 1.0/6-0.5*bi*bi};
	    double dwghk[4] = {1.0/3-ci+0.5*ci*ci, 0.5+2*ci-1.5*ci*ci, -1-ci+1.5*ci*ci, 1.0/6-0.5*ci*ci};
         
// boundary correction. 
	    if( ic4 == 2 )    dwghi[0] *= 2;
	    if( ic4 == Ni-2 ) dwghi[3] *= 2;
	    if( jc4 == 2 )    dwghj[0] *= 2;
	    if( jc4 == Nj-2 ) dwghj[3] *= 2;
	    if( kc4 == 2    && !(ccbndry && upperbndry) ) dwghk[0] *= 2;
	    if( kc4 == Nz-2 && !(ccbndry && lowerbndry) ) dwghk[3] *= 2;

	    double qX0[3], rX0[3], sX0[3];// qX0[0] = dq/dx, qX0[1] = dq/dy, qX0[2] = dq/dz, etc
            if( curvilinear )
	    {
	       //	      if (!a_EW->find_curvilinear_derivatives_at_point( q, r, s, qX0, rX0, sX0 ))
	       //	      {
	       //		 cout << "Unable to obtain derivatives of the curvilinear grid for source:" << endl << *this << endl;
	       //		 MPI_Abort(MPI_COMM_WORLD,1);
	       //
	       //	      }

	      double d4cofi[5], d4cofj[5], d4cofk[5];
	      // redefine ai,bi,ci, old def not needed below

              int kc3p = kc3;
	      if( kc3p >= Nz ) kc3p = Nz-1;

	      ai=q-ic3;
	      bi=r-jc3;
	      ci=s-kc3p;
	      d4cofi[0] = (1.0-ai)/12-ai*ai/4+ai*ai*ai/6;
	      d4cofi[1] =(-1.0+2*ai)*2.0/3+ai*ai/2-2*ai*ai*ai/3;
	      d4cofi[2] =-5*ai/2+ai*ai*ai;
	      d4cofi[3] = (1.0+2*ai)*2.0/3-ai*ai/2-2*ai*ai*ai/3;
	      d4cofi[4] = (-1.0-ai)/12+ai*ai/4+ai*ai*ai/6;

	      d4cofj[0] = (1.0-bi)/12-bi*bi/4+bi*bi*bi/6;
	      d4cofj[1] =(-1.0+2*bi)*2.0/3+bi*bi/2-2*bi*bi*bi/3;
	      d4cofj[2] =-5*bi/2+bi*bi*bi;
	      d4cofj[3] = (1.0+2*bi)*2.0/3-bi*bi/2-2*bi*bi*bi/3;
	      d4cofj[4] = (-1.0-bi)/12+bi*bi/4+bi*bi*bi/6;

	      d4cofk[0] = (1.0-ci)/12-ci*ci/4+ci*ci*ci/6;
	      d4cofk[1] =(-1.0+2*ci)*2.0/3+ci*ci/2-2*ci*ci*ci/3;
	      d4cofk[2] =-5*ci/2+ci*ci*ci;
	      d4cofk[3] = (1.0+2*ci)*2.0/3-ci*ci/2-2*ci*ci*ci/3;
	      d4cofk[4] = (-1.0-ci)/12+ci*ci/4+ci*ci*ci/6;

              double a4cofi[5], a4cofj[5], a4cofk[5];
	      a4cofi[0] = ai*(1-ai*ai)/12;
	      a4cofi[1] = -2*ai/3+ai*ai/2+ai*ai*ai/6;
	      a4cofi[2] = 1-ai*ai;
	      a4cofi[3] =  2*ai/3+ai*ai/2-ai*ai*ai/6;
	      a4cofi[4] = ai*(-1+ai*ai)/12;

	      a4cofj[0] = bi*(1-bi*bi)/12;
	      a4cofj[1] = -2*bi/3+bi*bi/2+bi*bi*bi/6;
	      a4cofj[2] = 1-bi*bi;
	      a4cofj[3] =  2*bi/3+bi*bi/2-bi*bi*bi/6;
	      a4cofj[4] = bi*(-1+bi*bi)/12;

	      a4cofk[0] = ci*(1-ci*ci)/12;
	      a4cofk[1] = -2*ci/3+ci*ci/2+ci*ci*ci/6;
	      a4cofk[2] = 1-ci*ci;
	      a4cofk[3] =  2*ci/3+ci*ci/2-ci*ci*ci/6;
	      a4cofk[4] = ci*(-1+ci*ci)/12;

	      //	      cout << "metric at pt " << qX0[0] << " " << rX0[0] << " " << sX0[0] << endl;
	      //	      cout << "metric at pt " << qX0[1] << " " << rX0[1] << " " << sX0[1] << endl;
	      //	      cout << "metric at pt " << qX0[2] << " " << rX0[2] << " " << sX0[2] << endl;

	      // Assume grid uniform in x and y, compute metric with z=z(q,r,s)

              double zq = 0, zr = 0, zs=0;
              for( int k=kc3p-2; k <= kc3p+2 ; k++ )
		 for( int j=jc3-2; j <= jc3+2 ; j++ )
		    for( int i=ic3-2; i <= ic3+2 ; i++ )
		    {
		       zq += d4cofi[i-(ic3-2)]*a4cofj[j-(jc3-2)]*a4cofk[k-(kc3p-2)]*a_EW->mZ(i,j,k);
		       zr += a4cofi[i-(ic3-2)]*d4cofj[j-(jc3-2)]*a4cofk[k-(kc3p-2)]*a_EW->mZ(i,j,k);
		       zs += a4cofi[i-(ic3-2)]*a4cofj[j-(jc3-2)]*d4cofk[k-(kc3p-2)]*a_EW->mZ(i,j,k);
		    }
	      qX0[0] = 1/h;
	      qX0[1] = 0;
	      qX0[2] = 0;
	      rX0[0] = 0;
	      rX0[1] = 1/h;
	      rX0[2] = 0;
	      sX0[0] = -zq/(h*zs);
	      sX0[1] = -zr/(h*zs);
	      sX0[2] = 1/zs;
	      //              cout << "Recomputed metric:" << endl;
	      //	      cout << "metric at pt " << qX0[0] << " " << rX0[0] << " " << sX0[0] << endl;
	      //	      cout << "metric at pt " << qX0[1] << " " << rX0[1] << " " << sX0[1] << endl;
	      //	      cout << "metric at pt " << qX0[2] << " " << rX0[2] << " " << sX0[2] << endl;

	    } // end if curvilinear
	    else
	    {
               qX0[0] = 1/h;qX0[1]=0;  qX0[2]=0;
	       rX0[0] = 0;  rX0[1]=1/h;rX0[2]=0;
	       sX0[0] = 0;  sX0[1]=0;  sX0[2]=1/h;
	    }

	    for( int k=kc4-1 ; k <= kc4+2 ; k++ )
	       for( int j=jc4-1 ; j <= jc4+2 ; j++ )
		  for( int i=ic4-1 ; i <= ic4+2 ; i++ )
		  {
		     double wFx=0, wFy=0, wFz=0;
// Note that point_in_proc can not happen sooner because a Moment source gets distributed over a 4x4x4 stencil, so some 
// points might be on a neighboring processor
		     if( a_EW->point_in_proc(i,j,g) ) 
		     {
			if( 0 <= j-jc3+1 && j-jc3+1 <= 2 && 0 <= k-kc3+1 && k-kc3+1 <= 2 )
			{
			   wFx += qX0[0]*dwghi[i-ic4+1]*wghj[j-jc3+1]*wghk[k-kc3+1];
			   wFy += qX0[1]*dwghi[i-ic4+1]*wghj[j-jc3+1]*wghk[k-kc3+1];
			   wFz += qX0[2]*dwghi[i-ic4+1]*wghj[j-jc3+1]*wghk[k-kc3+1];
			}
			if(  0 <= i-ic3+1 && i-ic3+1 <= 2 && 0 <= k-kc3+1 && k-kc3+1 <= 2 )
			{
			   wFx += wghi[i-ic3+1]*rX0[0]*dwghj[j-jc4+1]*wghk[k-kc3+1];
			   wFy += wghi[i-ic3+1]*rX0[1]*dwghj[j-jc4+1]*wghk[k-kc3+1];
			   wFz += wghi[i-ic3+1]*rX0[2]*dwghj[j-jc4+1]*wghk[k-kc3+1];
			}
			if(  0 <= i-ic3+1 && i-ic3+1 <= 2 && 0 <= j-jc3+1 && j-jc3+1 <= 2 )
			{
			   wFx += wghi[i-ic3+1]*wghj[j-jc3+1]*sX0[0]*dwghk[k-kc4+1];
			   wFy += wghi[i-ic3+1]*wghj[j-jc3+1]*sX0[1]*dwghk[k-kc4+1];
			   wFz += wghi[i-ic3+1]*wghj[j-jc3+1]*sX0[2]*dwghk[k-kc4+1];
			}
			double jaci;
			if( curvilinear )
			   jaci = 1.0/a_EW->mJ(i,j,k);
			else
			   jaci = 1.0/(h*h*h);

			double fx = -(mForces[0]*wFx+mForces[1]*wFy+mForces[2]*wFz)*jaci;
			double fy = -(mForces[1]*wFx+mForces[3]*wFy+mForces[4]*wFz)*jaci;
			double fz = -(mForces[2]*wFx+mForces[4]*wFy+mForces[5]*wFz)*jaci;

			if( mAmp != 0 && (fx != 0 || fy != 0 || fz != 0) )
			{
			   if( 1 <= k && k <= Nz )
			   {
			      GridPointSource* sourcePtr = new GridPointSource( mAmp, mFreq, mT0, i, j, k, g, 
							      fx, fy, fz, mTimeDependence, mNcyc );
			      point_sources.push_back(sourcePtr);
			   }
			   if( k <= 1 && ccbndry && upperbndry )
			   {
			     int Nzp = a_EW->m_global_nz[g+1];
//			      int Nzp = a_EW->m_kEnd[g+1]- a_EW->m_ghost_points;
			      int kk = Nzp - 1 + k;
			      GridPointSource* sourcePtr = new GridPointSource( mAmp, mFreq, mT0, i, j, kk, g+1, 
							      fx, fy, fz, mTimeDependence, mNcyc );
			      point_sources.push_back(sourcePtr);
			   }
			   if( k >= Nz && ccbndry && lowerbndry )
			   {
			      int kk = k - Nz + 1;
			      GridPointSource* sourcePtr = new GridPointSource( mAmp, mFreq, mT0, i, j, kk, g-1, 
							     fx, fy, fz, mTimeDependence, mNcyc );
			      point_sources.push_back(sourcePtr);

			   }
			}
		     }
		  } // end for i...
	    
	 } // end momentSource
	 //	 double m0[3]={0,0,0};
	 //	 for( unsigned int s=0 ; s < point_sources.size() ; s++ )
	 //	 {
	 //	    double jac, x, y, z;
	 //	    int i0= point_sources[s]->m_i0;
	 //	    int j0= point_sources[s]->m_j0;
	 //            int k0= point_sources[s]->m_k0;
	 //	    if( curvilinear )
	 //	    {
	 //	       jac = a_EW->mJ(i0,j0,k0);
	 //               x   = a_EW->mX(i0,j0,k0);
	 //               y   = a_EW->mY(i0,j0,k0);
	 //               z   = a_EW->mZ(i0,j0,k0);
	 //	    }
	 //	    else
	 //	    {
	 //	       jac = (h*h*h);
	 //               x   = h*(i0-1);
	 //               y   = h*(j0-1);
	 //               z   = h*(k0-1);
	 //	    }
	 //	    //            cout << *point_sources[s] << endl;
	 //           double wgh = point_sources[s]->mAmp;
	 //	    //            double wgh = 1;
	 //	    m0[0] += point_sources[s]->mForces[0]*wgh*jac*x*x*x;
	 //	    m0[1] += point_sources[s]->mForces[1]*wgh*jac*y*y*y;
	 //	    m0[2] += point_sources[s]->mForces[2]*wgh*jac*z*z*z;
	 //	 }
	 //	 //	 cout << "Zero moments are " << m0[0] << " " << m0[1] << " " << m0[2] << endl;
	 //	 //	 cout << "Error " << m0[0]-1 << " " << m0[1]-1 << " " << m0[2]-1 << endl;
	 //	 //	 cout << "Error " << m0[0]-2*mX0 << " " << m0[1]-2*mY0 << " " << m0[2]-2*mZ0 << endl;
	 //	 //	 cout << "Error " << m0[0]-3*mX0*mX0 << " " << m0[1]-3*mY0*mY0 << " " << m0[2]-3*mZ0*mZ0 << endl;
	 //	 //	 cout << "Error " << m0[0]-1 << " " << m0[1]-1 << " " << m0[2]-1 << endl;

      } // end, interior to a grid, or on Cartesian/Curvilinear boundary      
   } // end if canBeInverted
} // end set_grid_point_sources



double Source::distributedhat(double x) {
  if (x <= 1 && x > -1){
    double xi=x*x;
    // return 3.0/4.0*(1-xi);
    // return 15.0/16.0*(1-2.0*xi+4.0*pow(xi,2));
    return 315.0/512.0*(3.0-20.0*xi+42.0*pow(xi,2)-36.0*pow(xi,3)+11.0*pow(xi,4)); 
    // This polynomial satisfies 4 moment conditions and has 2 continous derivatives 
  }  
  else
    return 0.0; 
}

double Source::dist_d_dx_dirac(double x) {
  if (x < 1 && x > -1){
    double xi=x*x;
    //  return 0.5*(-659.8388671878608*x*x*x*x*x*x*x*x*x
    // 		+1935.527343751010*x*x*x*x*x*x*x
    // 		-2009.970703125978*x*x*x*x*x
    // 		+852.7148437503693*x*x*x
    // // 		-118.4326171875402*x);
    return (5.608630371051406e+02*x*x*x*x*x*x*x*x*x*x*x
 	    -2.144476318343753e+03*x*x*x*x*x*x*x*x*x
 	    +3.145231933571939e+03*x*x*x*x*x*x*x
 	    -2.177468261704705e+03*x*x*x*x*x
 	    +6.928308105429047e+02*x*x*x
 	    -7.698120117152610e+01*x);
    // m=6 k=6
    //     return x*(-4.811325073226033e+01*xi*xi*xi*xi*xi*xi
    // 	      +2.598115539554080e+02*xi*xi*xi*xi*xi
    // 	      -6.804588317906174e+02*xi*xi*xi*xi
    // 	      +9.828849792561246e+02*xi*xi*xi
    // 	      -8.041786193933241e+02*xi*xi
    // 	      +3.505393981977732e+02*xi
    // 	      -6.343093872160293e+01);
    
    //return (35.0/32.0)*(-6*x+12*x*x*x-6*x*x*x*x*x);
    //    return x*(-16699.0/282.0*2.0
    //	      +4.0*5969.0/28.0*xi    
    //	      -6.0*68674.0/205.0*xi*xi    
    //	      +8.0*53227.0/220.0*xi*xi*xi       
    //	      -10.0*4091.0/62.0*xi*xi*xi*xi);            
    // This is the derivative of a polynomial which
    // satisfies 6 moment conditions and has 2 continous derivatives 
  }  
  else
    return 0.0;
}



// first we define the hat function needed to approximate
// the dirac distribution...
double Source::hat(double x)
{
  if (x <= 0 && x > -1)
    return 1 + x;
  else if (x > 0 && x <= 1)
    return 1 - x;
  else
    return 0.0;
}
// ...and then the derivative of thepiecewise cubic
// function used to approximate the dreivative of said function
double Source::d_qubic_dx(double x)
{
  if (x >= -2 && x < -1)
    return 1./6.*(11. + 12.*x + 3.*pow(x,2));
  else if (x >= -1 && x < 0)
    return 1./2.*(1. - 4.*x - 3.*pow(x,2));
  else if (x == 0)
    return 0.0;
  else if (x > 0 && x <= 1)
    return 1./2.*(-1. - 4.*x + 3.*pow(x,2));
  else if (x > 1 && x <= 2)
    return 1./6.*(-11. + 12.*x - 3.*pow(x,2));
  else 
    return 0.0;
}
// ...and also the lower order non-symmetric approximation used close to
// the surface
double Source::d_quadric_dx(double x)
{
  if (x >= -1 && x < 0)
    return 1./2.*(3. + 2.*x);
  else if (x >= 0 && x <= 1)
    return -2.*x;
  else if (x > 1 && x <= 2)
    return 1./2.*(-3. + 2.*x);
  else
    return 0.0;
}
// // we use the derivative of the hat
// // function close to free surfaces
// double Source::d_hat_dx(double x)
// {
//   if (x <= 0 && x > -1)
//     return 1;
//   else if (x > 0 && x <= 1)
//     return -1;
//   else
//     return 0.0;
// }

int Source::sgn(double arg)
{
  if (arg < 0)
    return -1;
  else 
    return 1;
}

//-----------------------------------------------------------------------
double Source::find_min_exponent()
{
  return -700.0;
}

//-----------------------------------------------------------------------
void Source::distribute_source_xyplane( EW *a_EW, vector<GridPointSource*>& point_sources, 
				        int g, int k, double wghz )
{
   int Ni = a_EW->m_global_nx[g];
   int Nj = a_EW->m_global_ny[g];
   int Nz = a_EW->m_global_nz[g];
//   int Nz = a_EW->m_kEnd[g]- a_EW->m_ghost_points;

   double h = a_EW->mGridSize[g];
   double q = mX0/h+1;
   double r = mY0/h+1;

   int ic3 = static_cast<int>(round(q));
   int jc3 = static_cast<int>(round(r));

// Bias stencil away from boundary
   if( ic3 <= 1 )  ic3 = 2;
   if( ic3 >= Ni ) ic3 = Ni-1;
   if( jc3 <= 1 )  jc3 = 2;
   if( jc3 >= Nj ) jc3 = Nj-1;

   double ai=q-ic3, bi=r-jc3;
// Delta distribution
   double wghi[3] = {0.5*(ai*ai-ai),1-ai*ai,0.5*(ai*ai+ai)};
   double wghj[3] = {0.5*(bi*bi-bi),1-bi*bi,0.5*(bi*bi+bi)};

// Boundary correction
   if( ic3 == 2 )    wghi[0] *= 2;
   if( ic3 == Ni-1 ) wghi[2] *= 2;
   if( jc3 == 2 )    wghj[0] *= 2;
   if( jc3 == Nj-1 ) wghj[2] *= 2;

   for( int j=jc3-1 ; j <= jc3+1 ; j++ )
      for( int i=ic3-1 ; i <= ic3+1 ; i++ )
      {
	 double wF = wghi[i-ic3+1]*wghj[j-jc3+1]*wghz;
	 if( (wF*mAmp != 0) && (mForces[0] != 0 || mForces[1] != 0 || mForces[2] != 0) 
	     && a_EW->point_in_proc(i,j,g) )
	 {
            wF /= h*h*h;
	    if( 1 <= k && k <= Nz )
	    {
	       GridPointSource* sourcePtr = new GridPointSource(
					      mAmp*wF, mFreq, mT0,
					      i,j,k,g,
					      mForces[0], mForces[1], mForces[2],
					      mTimeDependence, mNcyc );
	       point_sources.push_back(sourcePtr);
	    }
	 }
      }
}       

//-----------------------------------------------------------------------
void Source::distribute_source_xyplane_mom( EW *a_EW, vector<GridPointSource*>& point_sources,
					    int g, int k, double wghz, double dwghz )
{
   int Ni = a_EW->m_global_nx[g];
   int Nj = a_EW->m_global_ny[g];
   int Nz = a_EW->m_global_nz[g];
//   int Nz = a_EW->m_kEnd[g]- a_EW->m_ghost_points;

   double h = a_EW->mGridSize[g];
   double q = mX0/h+1;
   double r = mY0/h+1;

// Stencil for point source
   int ic3 = static_cast<int>(round(q));
   int jc3 = static_cast<int>(round(r));

// Bias stencil away from boundary
   if( ic3 <= 1 )  ic3 = 2;
   if( ic3 >= Ni ) ic3 = Ni-1;
   if( jc3 <= 1 )  jc3 = 2;
   if( jc3 >= Nj ) jc3 = Nj-1;

   double ai=q-ic3, bi=r-jc3;
// Delta distribution
   double wghi[3] = {0.5*(ai*ai-ai),1-ai*ai,0.5*(ai*ai+ai)};
   double wghj[3] = {0.5*(bi*bi-bi),1-bi*bi,0.5*(bi*bi+bi)};

// Boundary correction
   if( ic3 == 2 )    wghi[0] *= 2;
   if( ic3 == Ni-1 ) wghi[2] *= 2;
   if( jc3 == 2 )    wghj[0] *= 2;
   if( jc3 == Nj-1 ) wghj[2] *= 2;

// Stencil for moment source

   int ic4 = static_cast<int>(floor(q));
   int jc4 = static_cast<int>(floor(r));

// Bias stencil away from boundary
   if( ic4 <= 1 )    ic4 = 2;
   if( ic4 >= Ni-1 ) ic4 = Ni-2;
   if( jc4 <= 1 )    jc4 = 2;
   if( jc4 >= Nj-1 ) jc4 = Nj-2;

// Bias stencil away from boundary
   ai = q-ic4; bi = r-jc4;
// Delta distribution derivatives
   double dwghi[4] = {1.0/3-ai+0.5*ai*ai, 0.5+2*ai-1.5*ai*ai, -1-ai+1.5*ai*ai, 1.0/6-0.5*ai*ai};
   double dwghj[4] = {1.0/3-bi+0.5*bi*bi, 0.5+2*bi-1.5*bi*bi, -1-bi+1.5*bi*bi, 1.0/6-0.5*bi*bi};

// Boundary correction
   if( ic4 == 2 )    dwghi[0] *= 2;
   if( ic4 == Ni-2 ) dwghi[3] *= 2;
   if( jc4 == 2 )    dwghj[0] *= 2;
   if( jc4 == Nj-2 ) dwghj[3] *= 2;

   for( int j=jc4-1 ; j <= jc4+2 ; j++ )
      for( int i=ic4-1 ; i <= ic4+2 ; i++ )
      {
	 double wFx=0, wFy=0, wFz=0;
	 if( a_EW->point_in_proc(i,j,g) ) 
	 {
	    if( 0 <= j-jc3+1 && j-jc3+1 <= 2  )
	    {
	       wFx += dwghi[i-ic4+1]*wghj[j-jc3+1]*wghz;
	    }
	    if(  0 <= i-ic3+1 && i-ic3+1 <= 2 )
	    {
	       wFy += wghi[i-ic3+1]*dwghj[j-jc4+1]*wghz;
	    }
	    if(  0 <= i-ic3+1 && i-ic3+1 <= 2 && 0 <= j-jc3+1 && j-jc3+1 <= 2 )
	    {
	       wFz += wghi[i-ic3+1]*wghj[j-jc3+1]*dwghz;
	    }
	    double jaci=1/(h*h*h*h);
	    if( 1 <= k && k <= Nz )
	    {
	       double fx = -(mForces[0]*wFx+mForces[1]*wFy+mForces[2]*wFz)*jaci;
	       double fy = -(mForces[1]*wFx+mForces[3]*wFy+mForces[4]*wFz)*jaci;
	       double fz = -(mForces[2]*wFx+mForces[4]*wFy+mForces[5]*wFz)*jaci;
	       if( mAmp != 0 && (fx != 0 || fy != 0 || fz != 0) )
	       {
		  GridPointSource* sourcePtr = new GridPointSource( mAmp, mFreq, mT0, i, j, k, g, 
								    fx, fy, fz, mTimeDependence, 
								    mNcyc );
		  point_sources.push_back(sourcePtr);
	       }
	    }
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
Source* Source::copy( EW* a_ew )
{

   Source* retval;
   if( !mIsMomentSource )
      retval = new Source( a_ew, mAmp, mFreq, mT0, mX0, mY0, mZ0, mForces[0],
			   mForces[1], mForces[2], mTimeDependence, mName.c_str(), mNcyc );
   else
      retval = new Source( a_ew, mAmp, mFreq, mT0, mX0, mY0, mZ0, mForces[0],
			   mForces[1], mForces[2], mForces[3], mForces[4], mForces[5],
			   mTimeDependence, mName.c_str(), mNcyc );
   retval->m_derivative = m_derivative;
   return retval;
}
