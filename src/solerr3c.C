#if defined(SW4_USE_RAJA) // SW4 backend
#include "EW.h"
#else // SW4 backend
#endif // SW4 backend
#include "sw4.h"
#if defined(SW4_USE_RAJA) // SW4 backend
#else // SW4 backend
#include <iostream>
#include "EW.h"
//-----------------------------------------------------------------------
#endif // SW4 backend
void EW::solerr3_ci( int ib, int ie, int jb, int je, int kb, int ke,
                     float_sw4 h, float_sw4* __restrict__ uex,
                     float_sw4* __restrict__ u, float_sw4& li,
                     float_sw4& l2, float_sw4& xli, float_sw4 zmin, float_sw4 x0,
                     float_sw4 y0, float_sw4 z0, float_sw4 radius,
                     int imin, int imax, int jmin, int jmax,
#if defined(SW4_USE_RAJA) // SW4 backend
int kmin, int kmax) {

#else // SW4 backend
int kmin, int kmax,
                     int geocube, int i0, int i1, int j0, int j1, int k0, int k1 )
{

#endif // SW4 backend
li = 0;
   l2 = 0;
   xli= 0;
   float_sw4 sradius2 = radius*radius;
   if( radius < 0 )
      sradius2 = -sradius2;

#if defined(SW4_USE_RAJA) // SW4 backend
#else // SW4 backend
int imxerr, jmxerr, kmxerr;

#endif // SW4 backend
const size_t ni = ie-ib+1;
   const size_t nij = ni*(je-jb+1);
   const float_sw4 h3 = h*h*h;
   const size_t nijk = nij*(ke-kb+1);

for( int c=0 ; c<3 ;c++)
   {
      float_sw4 liloc=0, l2loc=0, xliloc=0;
#pragma omp parallel for reduction(max:liloc,xliloc) reduction(+:l2loc)
         for( size_t k=kmin; k <= kmax ; k++ )
            for( size_t j=jmin; j <= jmax ; j++ )
               for( size_t i=imin; i <= imax ; i++ )
               {

#if defined(SW4_USE_RAJA) // SW4 backend
#else // SW4 backend
bool inside =  i<i0 || i> i1 || j<j0 || j>j1 || k<k0 || k>k1 ;

#endif // SW4 backend
if( ((i-1)*h-x0)*((i-1)*h-x0)+((j-1)*h-y0)*((j-1)*h-y0)+
                      ((k-1)*h+zmin-z0)*((k-1)*h+zmin-z0) >
#if defined(SW4_USE_RAJA) // SW4 backend
sradius2) {

#else // SW4 backend
sradius2 &&
                      ( geocube != 1 ||( geocube==1 && inside ) ) )
                  {

#endif // SW4 backend
size_t ind = i-ib+ni*(j-jb)+nij*(k-kb)+nijk*c;
                     float_sw4 err=fabs(u[ind]-uex[ind]);
             //		     if( fabs(u[ind]-uex[ind])>liloc )
                     //			liloc = fabs(u[ind]-uex[ind]);

#if defined(SW4_USE_RAJA) // SW4 backend
if (err > liloc) liloc = err;

#else // SW4 backend
if( err > liloc )
                     {
                        liloc = err;
                        imxerr = i;
                        jmxerr = j;
                        kmxerr = k;
                     }

#endif // SW4 backend
if( fabs(uex[ind])>xliloc )
                        xliloc = fabs(uex[ind]);
                     //		     l2loc += h3*(u[ind]-uex[ind])*(u[ind]-uex[ind]);

l2loc += h3*err*err;
                  }
               }
         l2 += l2loc;
         li  =  liloc >  li ?  liloc: li;
         xli = xliloc > xli ? xliloc:xli;

}
}

//-----------------------------------------------------------
void EW::solerrgp_ci( int ifirst, int ilast, int jfirst, int jlast,
                      int kfirst, int klast, float_sw4 h,
                      float_sw4* __restrict__ uex, float_sw4* __restrict__ u,
                      float_sw4& li, float_sw4& l2 )
{
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   const float_sw4 h3 = h*h*h;
   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);
   li = 0;
   l2 = 0;
   //   int k=kfirst+1;
   for( int c=0 ; c<3 ;c++)
   {
      float_sw4 liloc=0, l2loc=0;
#pragma omp parallel for reduction(max:liloc) reduction(+:l2loc)
      for( int j=jfirst+2; j<= jlast-2 ; j++ )
         for( int i=ifirst+2; i<= ilast-2 ;i++ )
         {
  // exact solution in array 'uex'
            size_t ind = base + i+ ni*j + nij*(kfirst+1)+nijk*c;
            float_sw4 err = abs( u[ind]-uex[ind] );
            if( liloc < err )
            {
               liloc = err;
               //	       cout << " num err= " << err << " at " << i << " " << j << " " << kfirst+1 << ", c= " << c << endl;
               //	       cout << "u, uex = " << u[ind] << " " << uex[ind] << endl;
            }
                //	    liloc = liloc > err ? liloc : err;
            l2loc += h3*err*err;
         }
      li  = li>liloc?li:liloc;
      l2 += l2loc;
   }
   //      k=klast-1;
   for( int c=0 ; c<3 ;c++)
   {
      float_sw4 liloc=0, l2loc=0;
#pragma omp parallel for reduction(max:liloc) reduction(+:l2loc)
      for( int j=jfirst+2; j<= jlast-2 ; j++ )
         for( int i=ifirst+2; i<= ilast-2 ;i++ )
         {
            size_t ind = base + i+ ni*j + nij*(klast-1)+nijk*c;
            float_sw4 err = abs( u[ind]-uex[ind] );
            if( liloc < err )
            {
               liloc = err;
               //	       cout << " num err= " << err << " at " << i << " " << j << " " << klast-1 << ", c= " << c <<  endl;
               //	       cout << "u, uex = " << u[ind] << " " << uex[ind] << endl;
            }
               //	    liloc = liloc > err ? liloc : err;
            l2loc += h3*err*err;
         }
      li  = li>liloc?li:liloc;
      l2 += l2loc;
   }
}

//-----------------------------------------------------------------------
void EW::solerr3c_ci( int ib, int ie, int jb, int je, int kb, int ke,
                      float_sw4* __restrict__ uex, float_sw4* __restrict__ u,
                      float_sw4* __restrict__ x, float_sw4* __restrict__ y,
                      float_sw4* __restrict__ z, float_sw4* __restrict__ jac,
                      float_sw4& li, float_sw4& l2, float_sw4& xli, float_sw4 x0,
                      float_sw4 y0, float_sw4 z0, float_sw4 radius,
                      int imin, int imax, int jmin, int jmax, int kmin, int kmax,
                      int usesg, float_sw4* __restrict__ strx, float_sw4* __restrict__ stry )
{
   li = 0;
   l2 = 0;
   xli= 0;
   float_sw4 sradius2 = radius*radius;
   if( radius < 0 )
      sradius2 = -sradius2;
   const size_t ni = ie-ib+1;
   const size_t nij = ni*(je-jb+1);
   const size_t nijk = nij*(ke-kb+1);
   const size_t base  = -(ib+ni*jb+nij*kb);

#if defined(SW4_USE_RAJA) // SW4 backend
#else // SW4 backend
int imxerr=0, jmxerr=0, kmxerr=0;

#endif // SW4 backend
for( int c=0 ; c<3 ;c++)
   {
      float_sw4 liloc=0, xliloc=0, l2loc=0;
#pragma omp parallel for reduction(max:liloc,xliloc) reduction(+:l2loc)
         for( size_t k=kmin; k <= kmax ; k++ )
            for( size_t j=jmin; j <= jmax ; j++ )
               for( size_t i=imin; i <= imax ; i++ )
               {
                  size_t ind = base+i+ni*j+nij*k;
                  float_sw4 dist = (x[ind]-x0)*(x[ind]-x0)+(y[ind]-y0)*(y[ind]-y0)+
                     (z[ind]-z0)*(z[ind]-z0);
                  if( dist > sradius2 )
                  {
                     size_t ind3 = ind+nijk*c;
                     float_sw4 err = fabs(u[ind3]-uex[ind3]);
                     //		     liloc = liloc > err ? liloc:err;

#if defined(SW4_USE_RAJA) // SW4 backend
liloc = liloc > err ? liloc : err;
            xliloc = xliloc > fabs(uex[ind3]) ? xliloc : fabs(uex[ind3]);
            // xliloc = xliloc > uex[ind3] ? xliloc : uex[ind3]; // ORG PRE
            // CURVIMR std::cout<<"SOL"<<ind3<<" "<<xliloc<<"
            // "<<uex[ind3]<<"\n";

#else // SW4 backend
#endif // SW4 backend
#if defined(SW4_USE_RAJA) // SW4 backend
#else // SW4 backend
if( liloc < err )
                     {
                        liloc = err;
                        imxerr = i;
                        jmxerr = j;
                        kmxerr = k;
                     }
                     xliloc = xliloc > fabs(uex[ind3]) ? xliloc:fabs(uex[ind3]);

#endif // SW4 backend
if( usesg != 1 )
                        l2loc += jac[ind]*err*err;
                     else
                        l2loc += jac[ind]*err*err/(strx[i-ib]*stry[j-jb]);
                  }
               }
          li =  liloc> li? liloc: li;

#if defined(SW4_USE_RAJA) // SW4 backend
l2 += l2loc;

#else // SW4 backend
#endif // SW4 backend
xli = xliloc>xli?xliloc:xli;

#if defined(SW4_USE_RAJA) // SW4 backend
}
#else // SW4 backend
l2 += l2loc;

#endif // SW4 backend
}
   //   std::cout << "Max error on grid: " << li << " at (i,j,k)= " << imxerr << " " << jmxerr << " " << kmxerr << std::endl;
#if defined(SW4_USE_RAJA) // SW4 backend
#ifdef XL_BUG_152435FIXED
//-----------------------------------------------------------------------
// This version does not compile due to a bug in XL :
//
// LLNL: SW4 ICE in clangtana with -qsmp=omp (152435).
// Use workaround below until this is fixed.
void EW::meterr4c_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                     int klast, float_sw4* __restrict__ met,
                     float_sw4* __restrict__ metex, float_sw4* __restrict__ jac,
                     float_sw4* __restrict__ jacex, float_sw4 li[5],
                     float_sw4 l2[5], int imin, int imax, int jmin, int jmax,
                     int kmin, int kmax, float_sw4 h) {
  const size_t ni = ilast - ifirst + 1;
  const size_t nij = ni * (jlast - jfirst + 1);
  const size_t nijk = nij * (klast - kfirst + 1);
  const size_t base = -(ifirst + ni * jfirst + nij * kfirst);

  for (int c = 0; c < 5; c++) li[c] = l2[c] = 0;

  const float_sw4 isqh = 1 / sqrt(h);
  const float_sw4 ih3 = 1 / (h * h * h);
#pragma omp parallel
  for (int c = 0; c < 5; c++)
#pragma omp for reduction(max : li[:5]) reduction(+ : l2[:5])
    for (size_t k = kmin; k <= kmax; k++)
      for (size_t j = jmin; j <= jmax; j++)
        for (size_t i = imin; i <= imax; i++) {
          size_t ind = base + i + ni * j + nij * k;
          float_sw4 err;
          if (c < 4)
            err = fabs(met[ind + c * nijk] - metex[ind + c * nijk]) * isqh;
          else
            err = fabs(jac[ind] - jacex[ind]) * ih3;
          li[c] = li[c] > err ? li[c] : err;
          l2[c] += jacex[ind] * err * err;
        }
}
#else
//-----------------------------------------------------------------------
void EW::meterr4c_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                     int klast, float_sw4* __restrict__ met,
                     float_sw4* __restrict__ metex, float_sw4* __restrict__ jac,
                     float_sw4* __restrict__ jacex, float_sw4 li[5],
                     float_sw4 l2[5], int imin, int imax, int jmin, int jmax,
                     int kmin, int kmax, float_sw4 h) {
  const size_t ni = ilast - ifirst + 1;
  const size_t nij = ni * (jlast - jfirst + 1);
  const size_t nijk = nij * (klast - kfirst + 1);
  const size_t base = -(ifirst + ni * jfirst + nij * kfirst);

  for (int c = 0; c < 5; c++) li[c] = l2[c] = 0;

  float_sw4 tmp_li;
  float_sw4 tmp_l2;
  const float_sw4 isqh = 1 / sqrt(h);
  const float_sw4 ih3 = 1 / (h * h * h);
  for (int c = 0; c < 5; c++) {
    tmp_li = li[c];
    tmp_l2 = l2[c];
#pragma omp parallel for reduction(max : tmp_li) reduction(+ : tmp_l2)
    for (size_t k = kmin; k <= kmax; k++)
      for (size_t j = jmin; j <= jmax; j++)
        for (size_t i = imin; i <= imax; i++) {
          size_t ind = base + i + ni * j + nij * k;
          float_sw4 err;
          if (c < 4)
            err = fabs(met[ind + c * nijk] - metex[ind + c * nijk]) * isqh;
          else
            err = fabs(jac[ind] - jacex[ind]) * ih3;
          li[c] = li[c] > err ? li[c] : err;
          l2[c] += jacex[ind] * err * err;
        }
    li[c] = tmp_li;
    l2[c] = tmp_l2;
  }
}
#endif
#else // SW4 backend
}

//-----------------------------------------------------------------------
void EW::meterr4c_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                     int klast, float_sw4* __restrict__ met, float_sw4* __restrict__ metex,
                     float_sw4* __restrict__ jac, float_sw4* __restrict__ jacex,
                     float_sw4 li[5], float_sw4 l2[5], int imin, int imax, int jmin,
                     int jmax, int kmin, int kmax, float_sw4 h )
{
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);

   for( int c=0; c< 5;c ++ )
      li[c] = l2[c] = 0;

   const float_sw4 isqh = 1/sqrt(h);
   const float_sw4 ih3  = 1/(h*h*h);
#pragma omp parallel
   for( int c=0 ; c < 5 ; c++ )
#pragma omp for reduction(max:li[c]) reduction(+:l2[c])
      for( size_t k=kmin; k <= kmax ; k++ )
         for( size_t j=jmin; j <= jmax ; j++ )
            for( size_t i=imin; i <= imax ; i++ )
            {
               size_t ind = base+i+ni*j+nij*k;
               float_sw4 err;
               if( c < 4 )
                  err = fabs( met[ind+c*nijk]- metex[ind+c*nijk] )*isqh;
               else
                  err = fabs(jac[ind]-jacex[ind])*ih3;
               li[c] = li[c]>err?li[c]:err;
               l2[c] += jacex[ind]*err*err;
            }
}
#endif // SW4 backend
