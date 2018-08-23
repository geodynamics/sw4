#include "sw4.h"
#include "cf_interface.h"
#include <cstddef>
#include <cstdio>
#include "Mspace.h"
#include "policies.h"
#include "caliper.h"
#include "foralls.h"

//-----------------------------------------------------------------------
void update_unext( int ib, int ie, int jb, int je, int kb, int ke,
		   float_sw4* __restrict__ a_unext, float_sw4* __restrict__ a_up,
		   float_sw4* __restrict__ a_lutt, float_sw4* __restrict__ a_force,
		   float_sw4* __restrict__ a_rho, float_sw4 cof, int kic )
{
  SW4_MARK_FUNCTION;
// Reversed indexation
#define Lutt(c,i,j,k) a_lutt[-base3_lutt+i+ni*(j)+nij*(k)+nijk_lutt*(c)]   
#define Unext(c,i,j,k) a_unext[-base3_unext+i+ni*(j)+nij*(k)+nijk_unext*(c)]   
#define force(c,i,j,k) a_force[-base3_force+i+ni*(j)+nij*(k)+nijk_force*(c)]   

// up and rho have different dimensions in the k-direction   
#define up(c,i,j,k)   a_up[-base3_u+i+ni*(j)+nij*(k)+nijk_u*(c)]   
#define rho(i,j,k)   a_rho[-base_rho+i+ni*(j)+nij*(k)]   

  const long int ni    = ie-ib+1;
  const long int nij   = ni*(je-jb+1);
  const long int base_rho = (ib+ni*jb+nij*kb);

  
  const long int nijk_u  = nij*(ke-kb+1);
  const long int base3_u = (ib+ni*jb+nij*kb+nijk_u);

  const long int nijk_unext = nij*(1);
  const long int base3_unext = (ib+ni*jb+nij*kic+nijk_unext);

  const long int nijk_lutt = nij*(1);
  const long int base3_lutt = (ib+ni*jb+nij*kic+nijk_lutt);

  const long int nijk_force = nij*(1);
  const long int base3_force = (ib+ni*jb+nij*kic+nijk_force);

  long int c, j, i;
#pragma omp parallel private(i,j,c)
{
  for (c=1; c<=3; c++)
  {
#pragma omp for
    for(j=jb+2; j <= je-2 ; j++ )
    {
#pragma ivdep
      for(i=ib+2; i <= ie-2 ; i++ )
      {
	Unext(c,i,j,kic) = up(c,i,j,kic) + cof*(Lutt(c,i,j,kic)+force(c,i,j,kic))/rho(i,j,kic);
      } // end for i
    } // end for j
  } // end for c
 
} // end omp parallel

#undef Unext
#undef up
#undef force
#undef Lutt
} // end dpdmt_wind

//-----------------------------------------------------------------------
void dpdmt_wind( int ib, int ie, int jb, int je, int kb_tt, int ke_tt, int kb_u, int ke_u,
		 float_sw4* __restrict__ a_up, float_sw4* __restrict__ a_u,
		 float_sw4* __restrict__ a_um, float_sw4* __restrict__ a_utt,
		 float_sw4 dt2i )
{
  SW4_MARK_FUNCTION;
   // Reversed indexation
#define up(c,i,j,k)   a_up[-base3_u+i+ni*(j)+nij*(k)+nijk_u*(c)]   
#define u(c,i,j,k)    a_u[-base3_u+i+ni*(j)+nij*(k)+nijk_u*(c)]   
#define um(c,i,j,k)   a_um[-base3_u+i+ni*(j)+nij*(k)+nijk_u*(c)]
// u_tt has different dimensions in the k-direction   
#define u_tt(c,i,j,k) a_utt[-base3_tt+i+ni*(j)+nij*(k)+nijk_tt*(c)]   
  const long int ni    = ie-ib+1;
  const long int nij   = ni*(je-jb+1);
  
  const long int nijk_u  = nij*(ke_u-kb_u+1);
  const long int nijk_tt  = nij*(ke_tt-kb_tt+1);

  const long int base3_u = (ib+ni*jb+nij*kb_u+nijk_u);
  const long int base3_tt = (ib+ni*jb+nij*kb_tt+nijk_tt);

//   long int c, k, j, i;
// #pragma omp parallel private(k,i,j,c)
// {
//   for (c=1; c<=3; c++)
//   {
//     for(k= kb_tt; k <= ke_tt ; k++ )
//     {
// #pragma omp for
//       for(j=jb; j <= je ; j++ )
//       {
// #pragma simd
// #pragma ivdep
// 	for(i=ib; i <= ie ; i++ )
// 	{
// 	  u_tt(c,i,j,k) = dt2i*( up(c,i,j,k)-2*u(c,i,j,k)+um(c,i,j,k) );
// 	}
//       }
//     }
//   }
#define NO_COLLAPSE
#if defined(NO_COLLAPSE)
  Range<16> I(ib,ie+1);
  Range<4>J(jb,je+1);
  Range<4>K(kb_tt,ke_tt+1);
  forall3async(I,J,K, [=]RAJA_DEVICE(int i,int j,int k){
#pragma unroll
      for(int c=1;c<4;c++)
	u_tt(c,i,j,k) = dt2i*( up(c,i,j,k)-2*u(c,i,j,k)+um(c,i,j,k) );
    });
  
#else
  RAJA::RangeSegment i_range(ib,ie+1);
  RAJA::RangeSegment j_range(jb,je+1);
  RAJA::RangeSegment k_range(kb_tt,ke_tt+1);
  RAJA::RangeSegment c_range(1,4);
  RAJA::kernel<DPDMT_WIND_LOOP_POL_ASYNC>(
		       RAJA::make_tuple(i_range,j_range,k_range,c_range),
		       [=]RAJA_DEVICE (long int i,long int j, long int k,long int c) {
			 u_tt(c,i,j,k) = dt2i*( up(c,i,j,k)-2*u(c,i,j,k)+um(c,i,j,k) );
		       });
#endif
  //SYNC_STREAM;
#undef up
#undef u
#undef um
#undef u_tt
} // end dpdmt_wind


