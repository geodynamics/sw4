#include "sw4.h"
#include "EW.h"
#include "Require.h"
#include <iostream>

//#include <cuda_runtime.h>
#include "RAJA/RAJA.hxx"
using namespace RAJA;
#include "mynvtx.h"
#ifdef CUDA_CODE
typedef NestedPolicy<ExecList<cuda_threadblock_x_exec<4>,cuda_threadblock_y_exec<4>,
			      cuda_threadblock_z_exec<32>>>
  EXEC;
typedef NestedPolicy<ExecList<cuda_threadblock_x_exec<16>,cuda_threadblock_y_exec<16>,
			      cuda_threadblock_z_exec<16>>>
  EXEC_BC;

typedef NestedPolicy<ExecList<cuda_threadblock_x_exec<32>,cuda_threadblock_y_exec<32>>>
  EXEC_BC2;


typedef NestedPolicy<ExecList<seq_exec, seq_exec, seq_exec>>
  CEXEC_BC;

typedef RAJA::cuda_exec<1024> FEXEC;

#define SYNC_DEVICE cudaDeviceSynchronize();
#else

typedef NestedPolicy<ExecList<omp_parallel_for_exec,omp_parallel_for_exec,
			      omp_parallel_for_exec>>
  EXEC_BC;
typedef NestedPolicy<ExecList<omp_parallel_for_exec,omp_parallel_for_exec,
			      omp_parallel_for_exec>>
  EXEC;

typedef RAJA::omp_parallel_for_exec FEXEC;
#define SYNC_DEVICE 
#endif

void PrintPointerAttributes(void *ptr);

using namespace std;
//-----------------------------------------------------------------------
void EW::corrfort( int ib, int ie, int jb, int je, int kb, int ke, float_sw4* up,
		   float_sw4* lu, float_sw4* fo, float_sw4* rho, float_sw4 dt4 )
{
   const float_sw4 dt4i12= dt4/12;
   const size_t npts = static_cast<size_t>((ie-ib+1))*(je-jb+1)*(ke-kb+1);
   // cudaMemPrefetchAsync(up, 
   // 			3*npts*sizeof(float_sw4),
   // 			0,
   // 			0);
   // cudaMemPrefetchAsync(lu, 
   // 			3*npts*sizeof(float_sw4),
   // 			0,
   // 			0);
   // cudaMemPrefetchAsync(fo, 
   // 			3*npts*sizeof(float_sw4),
   // 			0,
   // 			0);
   // cudaMemPrefetchAsync(rho, 
   // 			npts*sizeof(float_sw4),
   // 			0,
   // 			0);
   // SYNC_DEVICE
   if( m_corder )
   {
     //for( size_t i=0 ; i < npts ; i++ )
     forall< FEXEC> (0,npts,[=] RAJA_DEVICE(size_t i)
      {
	 float_sw4 dt4i12orh = dt4i12/rho[i];
	 up[i  ]      += dt4i12orh*(lu[i  ]     +fo[i  ]);
	 up[i+npts]   += dt4i12orh*(lu[i+npts]  +fo[i+npts]);
	 up[i+2*npts] += dt4i12orh*(lu[i+2*npts]+fo[i+2*npts]);
      });
     SYNC_DEVICE
   }
   else
   {
     //for( size_t i=0 ; i < npts ; i++ )
     forall< FEXEC> (0,npts,[=] RAJA_DEVICE(size_t i)
      {
	 float_sw4 dt4i12orh = dt4i12/rho[i];
	 up[3*i  ] += dt4i12orh*(lu[3*i  ]+fo[3*i  ]);
	 up[3*i+1] += dt4i12orh*(lu[3*i+1]+fo[3*i+1]);
	 up[3*i+2] += dt4i12orh*(lu[3*i+2]+fo[3*i+2]);
      });
     SYNC_DEVICE
   }
}

//-----------------------------------------------------------------------
void EW::predfort( int ib, int ie, int jb, int je, int kb, int ke, float_sw4* up,
		   float_sw4* u, float_sw4* um, float_sw4* lu, float_sw4* fo,
		   float_sw4* rho, float_sw4 dt2 )
{
   const size_t npts = static_cast<size_t>((ie-ib+1))*(je-jb+1)*(ke-kb+1);
   // cudaMemPrefetchAsync(up, 
   // 			3*npts*sizeof(float_sw4),
   // 			0,
   // 			0);
   // cudaMemPrefetchAsync(u, 
   // 			3*npts*sizeof(float_sw4),
   // 			0,
   // 			0);
   // cudaMemPrefetchAsync(um, 
   // 			3*npts*sizeof(float_sw4),
   // 			0,
   // 			0);
   // cudaMemPrefetchAsync(lu, 
   // 			3*npts*sizeof(float_sw4),
   // 			0,
   // 			0);
   // cudaMemPrefetchAsync(fo, 
   // 			3*npts*sizeof(float_sw4),
   // 			0,
   // 			0);
   // cudaMemPrefetchAsync(rho, 
   // 			npts*sizeof(float_sw4),
   // 			0,
   // 			0);
   // SYNC_DEVICE
   if( m_corder )
   {
      // Like this ?
     //for( size_t i=0 ; i < npts ; i++ )
     forall< FEXEC> (0,npts,[=] RAJA_DEVICE(size_t i)
      {
	 float_sw4 dt2orh = dt2/rho[i];
	 up[i  ]      = 2*u[i  ]     -um[i  ]      + dt2orh*(lu[i  ]     +fo[i  ]);
	 up[i+npts]   = 2*u[i+npts]  -um[i+npts]   + dt2orh*(lu[i+npts]  +fo[i+npts]);
	 up[i+2*npts] = 2*u[i+2*npts]-um[i+2*npts] + dt2orh*(lu[i+2*npts]+fo[i+2*npts]);
      });
     SYNC_DEVICE
// Alternatives:
//      // Or this ?
//      for( int c=0 ; c < 3 ; c++ )
//	 for( size_t i=npts*c ; i < (c+1)*npts ; i++ )
//	    up[i] = 2*u[i] - um[i] + dt2/rho[i-c*npts]*(lu[i]+fo[i]);
//      // Or, perhaps, unrolled ?
//      for( size_t i=0 ; i < npts ; i++ )
//	 up[i] = 2*u[i] - um[i] + dt2/rho[i]*(lu[i]+fo[i]);
//      for( size_t i=npts ; i < 2*npts ; i++ )
//	 up[i] = 2*u[i] - um[i] + dt2/rho[i-npts]*(lu[i]+fo[i]);
//      for( size_t i=2*npts ; i < 3*npts ; i++ )
//	 up[i] = 2*u[i] - um[i] + dt2/rho[i-2*npts]*(lu[i]+fo[i]);
//    // Or, as one long array ?
//      for( size_t i=0 ; i < 3*npts ; i++ )
//	 up[i] = 2*u[i] - um[i] + dt2/rho[i%npts]*(lu[i]+fo[i]);
   }
   else
   {
     //for( size_t i=0 ; i < npts ; i++ )
forall<FEXEC > (0,npts,[=] RAJA_DEVICE(size_t i)
      {
	 float_sw4 dt2orh = dt2/rho[i];
	 up[3*i  ] = 2*u[3*i  ]-um[3*i  ] + dt2orh*(lu[3*i  ]+fo[3*i  ]);
	 up[3*i+1] = 2*u[3*i+1]-um[3*i+1] + dt2orh*(lu[3*i+1]+fo[3*i+1]);
	 up[3*i+2] = 2*u[3*i+2]-um[3*i+2] + dt2orh*(lu[3*i+2]+fo[3*i+2]);
      });
SYNC_DEVICE
   }
}

//-----------------------------------------------------------------------
void EW::dpdmtfort( int ib, int ie, int jb, int je, int kb, int ke, const float_sw4* __restrict__ up,
		    const float_sw4* __restrict__ u, const float_sw4* __restrict__ um,
		    float_sw4* __restrict__ u2, float_sw4 dt2i )
{
   const size_t npts = static_cast<size_t>((ie-ib+1))*(je-jb+1)*(ke-kb+1);
   // cudaMemPrefetchAsync(up, 
   // 			3*npts*sizeof(float_sw4),
   // 			0,
   // 			0);
   // cudaMemPrefetchAsync(u, 
   // 			3*npts*sizeof(float_sw4),
   // 			0,
   // 			0);
   // cudaMemPrefetchAsync(um, 
   // 			3*npts*sizeof(float_sw4),
   // 			0,
   // 			0);
   // cudaMemPrefetchAsync(u2, 
   // 			3*npts*sizeof(float_sw4),
   // 			0,
   // 			0);
   //SYNC_DEVICE
   // Code is slower with this RAJA loop due to the cost of copying data. Prefecth version is the fastest:: 2.70s on 20 procs. Without this loop it is around 2.4s
   // To fix. either port more loops to device or wait for the UM copies to reach 30GB/s
   forall< FEXEC> (0,3*npts,[=] RAJA_DEVICE(size_t i){
       //   for( size_t i = 0 ; i < 3*npts ; i++ )
      u2[i] = dt2i*(up[i]-2*u[i]+um[i]);
     });
   SYNC_DEVICE
   //   if( m_corder )
   //   {
   //      for( size_t i=0 ; i < npts ; i++ )
   //      {
   //	 u2[i  ]      = dt2i*(up[i  ]      - 2*u[i  ]      + um[i  ]);
   //	 u2[i+npts]   = dt2i*(up[i+npts]   - 2*u[i+npts]   + um[i+npts]);
   //	 u2[i+2*npts] = dt2i*(up[i+2*npts] - 2*u[i+2*npts] + um[i+2*npts]);
   //      }
   //   }
   //   else
   //   {
   //      for( size_t i=0 ; i < npts ; i++ )
   //      {
   //	 u2[3*i  ] = dt2i*(up[3*i  ] - 2*u[3*i  ] + um[3*i  ]);
   //	 u2[3*i+1] = dt2i*(up[3*i+1] - 2*u[3*i+1] + um[3*i+1]);
   //	 u2[3*i+2] = dt2i*(up[3*i+2] - 2*u[3*i+2] + um[3*i+2]);
   //      }
   //   }
}

//-----------------------------------------------------------------------
void EW::solerr3fort( int ib, int ie, int jb, int je, int kb, int ke,
		      float_sw4 h, float_sw4* uex, float_sw4* u, float_sw4& li,
		      float_sw4& l2, float_sw4& xli, float_sw4 zmin, float_sw4 x0,
		      float_sw4 y0, float_sw4 z0, float_sw4 radius,
		      int imin, int imax, int jmin, int jmax, int kmin, int kmax )
{

   li = 0;
   l2 = 0;
   xli= 0;
   float_sw4 sradius2 = radius*radius;
   if( radius < 0 )
      sradius2 = -sradius2;
   const size_t ni = ie-ib+1;
   const size_t nij = ni*(je-jb+1);
   const double h3 = h*h*h;
   if( m_corder )
   {
      const size_t nijk = nij*(ke-kb+1);
      for( int c=0 ; c<3 ;c++)
	 for( size_t k=kmin; k <= kmax ; k++ )
	    for( size_t j=jmin; j <= jmax ; j++ )
	       for( size_t i=imin; i <= imax ; i++ )
	       {
		  if( ((i-1)*h-x0)*((i-1)*h-x0)+((j-1)*h-y0)*((j-1)*h-y0)+
		      ((k-1)*h+zmin-z0)*((k-1)*h+zmin-z0) > sradius2 )
		  {
		     size_t ind = i-ib+ni*(j-jb)+nij*(k-kb)+nijk*c;
		     {
			if( fabs(u[ind]-uex[ind])>li )
			   li = fabs(u[ind]-uex[ind]);
			if( uex[ind]>xli )
			   xli = uex[ind];
			l2 += h3*(u[ind]-uex[ind])*(u[ind]-uex[ind]);
		     }
		  }
	       }
   }
   else
   {
      for( size_t k=kmin; k <= kmax ; k++ )
	 for( size_t j=jmin; j <= jmax ; j++ )
	    for( size_t i=imin; i <= imax ; i++ )
	    {
	       if( ((i-1)*h-x0)*((i-1)*h-x0)+((j-1)*h-y0)*((j-1)*h-y0)+
		   ((k-1)*h+zmin-z0)*((k-1)*h+zmin-z0) > sradius2 )
	       {
		  size_t ind = i-ib+ni*(j-jb)+nij*(k-kb);
		  for( int c=0 ; c<3 ;c++)
		  {
		     if( fabs(u[3*ind+c]-uex[3*ind+c])>li )
			li = fabs(u[3*ind+c]-uex[3*ind+c]);
		     if( uex[3*ind+c]>xli )
			xli = uex[3*ind+c];
		     l2 += h3*(u[3*ind+c]-uex[3*ind+c])*(u[3*ind+c]-uex[3*ind+c]);
		  }
	       }
	    }
   }
}

//-----------------------------------------------------------------------
void EW::bcfortsg( int ib, int ie, int jb, int je, int kb, int ke, int wind[36], 
		   int nx, int ny, int nz, float_sw4* u, float_sw4 h, boundaryConditionType bccnd[6],
		   float_sw4 sbop[5], float_sw4* mu, float_sw4* la, float_sw4 t,
		   float_sw4* bforce1, float_sw4* bforce2, float_sw4* bforce3, 
		   float_sw4* bforce4, float_sw4* bforce5, float_sw4* bforce6,
		   float_sw4 om, float_sw4 ph, float_sw4 cv,
		   float_sw4* strx, float_sw4* stry )
{
   const float_sw4 d4a = 2.0/3.0;
   const float_sw4 d4b = -1.0/12.0;
   const size_t ni  = ie-ib+1;
   const size_t nij = ni*(je-jb+1);
   for( int s=0 ; s < 6 ; s++ )
   {
//      printf("in bcfortsg, s=%d bccnd[s]=%d\n",s,bccnd[s]); fflush(stdout);
      if( bccnd[s]==1 || bccnd[s]==2 )
      {
         // Would it make any sense to line up pointers to the bforceN vectors in
         // a vector and collapse these 6 blocks? I think they're all the same.
         // This would make this code simpler and speed parameter passing.
         // OTOH there is some benefit in having the compiler know what s is.
         size_t idel = 1+wind[1+6*s]-wind[6*s];
         size_t ijdel = idel * (1+wind[3+6*s]-wind[2+6*s]);
	 if( s== 0 )
	 {
            // Note - don't use collapse(2) since loop indices are used in loop body;
            // it's hard to get them back after collapse
	   //for( k=wind[4+6*s]; k <= wind[5+6*s] ; k++ ) {
	   //   int qq = (k-wind[4+6*s])*ijdel;
	   //   for( int j=wind[2+6*s]; j <= wind[3+6*s] ; j++ ) {
	   //	  for( int i=wind[6*s]; i <= wind[1+6*s] ; i++ ) {
	   int is = wind[0+6*s];
	   int ie = wind[1+6*s];
	   int js = wind[2+6*s];
	   int je = wind[3+6*s];
	   int ks = wind[4+6*s];
	   int ke = wind[5+6*s];
	   //PrintPointerAttributes((void*)bforce1);
	   forallN<EXEC_BC, int, int,int>(
					  RangeSegment(ks,ke+1),
					  RangeSegment(js,je+1),
					  RangeSegment(is,ie+1),
					  [=] RAJA_DEVICE(int k, int j, int i) {
					 size_t ind = i-ib+ni*(j-jb)+nij*(k-kb);
					 int qq=(k-ks)*ijdel+(i-is)+(j-js)*(ie-is+1);
					 u[3*ind  ] = bforce1[  3*qq];
					 u[3*ind+1] = bforce1[1+3*qq];
					 u[3*ind+2] = bforce1[2+3*qq];
				       });
      
	 }
	 else if( s== 1 )
	 {
	   int is = wind[0+6*s];
	   int ie = wind[1+6*s];
	   int js = wind[2+6*s];
	   int je = wind[3+6*s];
	   int ks = wind[4+6*s];
	   int ke = wind[5+6*s];
	   //PrintPointerAttributes((void*)bforce1);
	   forallN<EXEC_BC, int, int,int>(
					  RangeSegment(ks,ke+1),
					  RangeSegment(js,je+1),
					  RangeSegment(is,ie+1),
					  [=] RAJA_DEVICE(int k, int j, int i) {
					    size_t ind = i-ib+ni*(j-jb)+nij*(k-kb);
					    int qq=(k-ks)*ijdel+(i-is)+(j-js)*(ie-is+1);
					    u[3*ind  ] = bforce2[  3*qq];
					    u[3*ind+1] = bforce2[1+3*qq];
					    u[3*ind+2] = bforce2[2+3*qq];
					  });
	  
	 }
	 else if( s==2 )
	 {
	   int is = wind[0+6*s];
	   int ie = wind[1+6*s];
	   int js = wind[2+6*s];
	   int je = wind[3+6*s];
	   int ks = wind[4+6*s];
	   int ke = wind[5+6*s];
	   //PrintPointerAttributes((void*)bforce1);
	   forallN<EXEC_BC, int, int,int>(
					  RangeSegment(ks,ke+1),
					  RangeSegment(js,je+1),
					  RangeSegment(is,ie+1),
					  [=] RAJA_DEVICE(int k, int j, int i) {
					    size_t ind = i-ib+ni*(j-jb)+nij*(k-kb);
					    int qq=(k-ks)*ijdel+(i-is)+(j-js)*(ie-is+1);
					    u[3*ind  ] = bforce3[  3*qq];
					    u[3*ind+1] = bforce3[1+3*qq];
					    u[3*ind+2] = bforce3[2+3*qq];
					  });
	   
	 }
	 else if( s==3 )
	 {
	   int is = wind[0+6*s];
	   int ie = wind[1+6*s];
	   int js = wind[2+6*s];
	   int je = wind[3+6*s];
	   int ks = wind[4+6*s];
	   int ke = wind[5+6*s];
	   //PrintPointerAttributes((void*)bforce1);
	   forallN<EXEC_BC, int, int,int>(
					  RangeSegment(ks,ke+1),
					  RangeSegment(js,je+1),
					  RangeSegment(is,ie+1),
					  [=] RAJA_DEVICE(int k, int j, int i) {
					    size_t ind = i-ib+ni*(j-jb)+nij*(k-kb);
					    int qq=(k-ks)*ijdel+(i-is)+(j-js)*(ie-is+1);
					    u[3*ind  ] = bforce4[  3*qq];
					    u[3*ind+1] = bforce4[1+3*qq];
					    u[3*ind+2] = bforce4[2+3*qq];
					  });
	  
	 }
	 else if( s==4 )
	 {
	   int is = wind[0+6*s];
	   int ie = wind[1+6*s];
	   int js = wind[2+6*s];
	   int je = wind[3+6*s];
	   int ks = wind[4+6*s];
	   int ke = wind[5+6*s];
	   //PrintPointerAttributes((void*)bforce1);
	   forallN<EXEC_BC, int, int,int>(
					  RangeSegment(ks,ke+1),
					  RangeSegment(js,je+1),
					  RangeSegment(is,ie+1),
					  [=] RAJA_DEVICE(int k, int j, int i) {
					    size_t ind = i-ib+ni*(j-jb)+nij*(k-kb);
					    int qq=(k-ks)*ijdel+(i-is)+(j-js)*(ie-is+1);
					    u[3*ind  ] = bforce5[  3*qq];
					    u[3*ind+1] = bforce5[1+3*qq];
					    u[3*ind+2] = bforce5[2+3*qq];
					  });
	   
	 }
	 else if( s==5 )
	 {
	   
	   int is = wind[0+6*s];
	   int ie = wind[1+6*s];
	   int js = wind[2+6*s];
	   int je = wind[3+6*s];
	   int ks = wind[4+6*s];
	   int ke = wind[5+6*s];
	   //PrintPointerAttributes((void*)bforce1);
	   forallN<EXEC_BC, int, int,int>(
					  RangeSegment(ks,ke+1),
					  RangeSegment(js,je+1),
					  RangeSegment(is,ie+1),
					  [=] RAJA_DEVICE(int k, int j, int i) {
					    size_t ind = i-ib+ni*(j-jb)+nij*(k-kb);
					    int qq=(k-ks)*ijdel+(i-is)+(j-js)*(ie-is+1);
					    u[3*ind  ] = bforce6[  3*qq];
					    u[3*ind+1] = bforce6[1+3*qq];
					    u[3*ind+2] = bforce6[2+3*qq];
					  });
	   
	   
	 }
      }
      else if( bccnd[s]==3 )
      {
	 if( s==0 )
	 {
            #pragma omp parallel for 
	    for( int k=wind[4+6*s]; k <= wind[5+6*s] ; k++ )
	       for( int j=wind[2+6*s]; j <= wind[3+6*s] ; j++ )
		  for( int i=wind[6*s]; i <= wind[1+6*s] ; i++ )
		  {
		     size_t ind = i-ib+ni*(j-jb)+nij*(k-kb);
		     size_t indp= ind+nx;
		     u[3*ind  ] = u[3*indp];
		     u[3*ind+1] = u[3*indp+1];
		     u[3*ind+2] = u[3*indp+2];
		  }
	 }
	 else if( s==1 )
	 {
            #pragma omp parallel for
	    for( int k=wind[4+6*s]; k <= wind[5+6*s] ; k++ )
	       for( int j=wind[2+6*s]; j <= wind[3+6*s] ; j++ )
		  for( int i=wind[6*s]; i <= wind[1+6*s] ; i++ )
		  {
		     size_t ind = i-ib+ni*(j-jb)+nij*(k-kb);
		     size_t indp= ind-nx;
		     u[3*ind  ] = u[3*indp];
		     u[3*ind+1] = u[3*indp+1];
		     u[3*ind+2] = u[3*indp+2];
		  }
	 }
	 else if( s==2 )
	 {
            #pragma omp parallel for
	    for( int k=wind[4+6*s]; k <= wind[5+6*s] ; k++ )
	       for( int j=wind[2+6*s]; j <= wind[3+6*s] ; j++ )
		  for( int i=wind[6*s]; i <= wind[1+6*s] ; i++ )
		  {
		     size_t ind = i-ib+ni*(j-jb)+nij*(k-kb);
		     size_t indp= ind+ni*ny;
		     u[3*ind  ] = u[3*indp];
		     u[3*ind+1] = u[3*indp+1];
		     u[3*ind+2] = u[3*indp+2];
		  }
	 }
	 else if( s==3 )
	 {
            #pragma omp parallel for
	    for( int k=wind[4+6*s]; k <= wind[5+6*s] ; k++ )
	       for( int j=wind[2+6*s]; j <= wind[3+6*s] ; j++ )
		  for( int i=wind[6*s]; i <= wind[1+6*s] ; i++ )
		  {
		     size_t ind = i-ib+ni*(j-jb)+nij*(k-kb);
		     size_t indp= ind-ni*ny;
		     u[3*ind  ] = u[3*indp];
		     u[3*ind+1] = u[3*indp+1];
		     u[3*ind+2] = u[3*indp+2];
		  }
	 }
	 else if( s==4 )
	 {
            #pragma omp parallel for
	    for( int k=wind[4+6*s]; k <= wind[5+6*s] ; k++ )
	       for( int j=wind[2+6*s]; j <= wind[3+6*s] ; j++ )
		  for( int i=wind[6*s]; i <= wind[1+6*s] ; i++ )
		  {
		     size_t ind = i-ib+ni*(j-jb)+nij*(k-kb);
		     size_t indp= ind+nij*nz;
		     u[3*ind  ] = u[3*indp];
		     u[3*ind+1] = u[3*indp+1];
		     u[3*ind+2] = u[3*indp+2];
		  }
	 }
	 else if( s==5 )
	 {
            #pragma omp parallel for 
	    for( int k=wind[4+6*s]; k <= wind[5+6*s] ; k++ )
	       for( int j=wind[2+6*s]; j <= wind[3+6*s] ; j++ )
		  for( int i=wind[6*s]; i <= wind[1+6*s] ; i++ )
		  {
		     size_t ind = i-ib+ni*(j-jb)+nij*(k-kb);
		     size_t indp= ind-nij*nz;
		     u[3*ind  ] = u[3*indp];
		     u[3*ind+1] = u[3*indp+1];
		     u[3*ind+2] = u[3*indp+2];
		  }
	 }
      }
      else if( bccnd[s]==0 )
      {
	 REQUIRE2( s == 4 || s == 5, "EW::bcfortsg,  ERROR: Free surface condition"
		  << " not implemented for side " << s << endl);
	 if( s==4 )
	 {
	    int k=1, kl=1;
            #pragma omp parallel for   
	    for( int j=jb+2 ; j <= je-2 ; j++ )
	       for( int i=ib+2 ; i <= ie-2 ; i++ )
	       {
		  size_t qq = i-ib+ni*(j-jb);
		  size_t ind = i-ib+ni*(j-jb)+nij*(k-kb);
		  float_sw4 wx = strx[i-ib]*(d4a*(u[2+3*ind+3]-u[2+3*ind-3])+d4b*(u[2+3*ind+6]-u[2+3*ind-6]));
		  float_sw4 ux = strx[i-ib]*(d4a*(u[  3*ind+3]-u[  3*ind-3])+d4b*(u[  3*ind+6]-u[  3*ind-6]));
		  float_sw4 wy = stry[j-jb]*(d4a*(u[2+3*ind+3*ni]-u[2+3*ind-3*ni])+
					     d4b*(u[2+3*ind+6*ni]-u[2+3*ind-6*ni]));
		  float_sw4 vy = stry[j-jb]*(d4a*(u[1+3*ind+3*ni]-u[1+3*ind-3*ni])+
					     d4b*(u[1+3*ind+6*ni]-u[1+3*ind-6*ni]));
		  float_sw4 uz=0, vz=0, wz=0;
		  for( int w=1 ; w <= 4 ; w++ )
		  {
		     uz += sbop[w]*u[  3*ind+3*nij*kl*(w-1)];
		     vz += sbop[w]*u[1+3*ind+3*nij*kl*(w-1)];
		     wz += sbop[w]*u[2+3*ind+3*nij*kl*(w-1)];
		  }
		  u[  3*ind-3*nij*kl] = (-uz-kl*wx+kl*h*bforce5[  3*qq]/mu[ind])/sbop[0];
		  u[1+3*ind-3*nij*kl] = (-vz-kl*wy+kl*h*bforce5[1+3*qq]/mu[ind])/sbop[0];
		  u[2+3*ind-3*nij*kl] = (-wz + (-kl*la[ind]*(ux+vy)+kl*h*bforce5[2+3*qq])/
					 (2*mu[ind]+la[ind]))/sbop[0];
	       }
	 }
	 else
	 {
	    int k=nz, kl=-1;
            #pragma omp parallel for   
	    for( int j=jb+2 ; j <= je-2 ; j++ )
	       for( int i=ib+2 ; i <= ie-2 ; i++ )
	       {
		  size_t qq = i-ib+ni*(j-jb);
		  size_t ind = i-ib+ni*(j-jb)+nij*(k-kb);
		  float_sw4 wx = strx[i-ib]*(d4a*(u[2+3*ind+3]-u[2+3*ind-3])+d4b*(u[2+3*ind+6]-u[2+3*ind-6]));
		  float_sw4 ux = strx[i-ib]*(d4a*(u[  3*ind+3]-u[  3*ind-3])+d4b*(u[  3*ind+6]-u[  3*ind-6]));
		  float_sw4 wy = stry[j-jb]*(d4a*(u[2+3*ind+3*ni]-u[2+3*ind-3*ni])+
					     d4b*(u[2+3*ind+6*ni]-u[2+3*ind-6*ni]));
		  float_sw4 vy = stry[j-jb]*(d4a*(u[1+3*ind+3*ni]-u[1+3*ind-3*ni])+
					     d4b*(u[1+3*ind+6*ni]-u[1+3*ind-6*ni]));
		  float_sw4 uz=0, vz=0, wz=0;
		  for( int w=1 ; w <= 4 ; w++ )
		  {
		     uz += sbop[w]*u[  3*ind+3*nij*kl*(w-1)];
		     vz += sbop[w]*u[1+3*ind+3*nij*kl*(w-1)];
		     wz += sbop[w]*u[2+3*ind+3*nij*kl*(w-1)];
		  }
		  u[  3*ind-3*nij*kl] = (-uz-kl*wx+kl*h*bforce6[  3*qq]/mu[ind])/sbop[0];
		  u[1+3*ind-3*nij*kl] = (-vz-kl*wy+kl*h*bforce6[1+3*qq]/mu[ind])/sbop[0];
		  u[2+3*ind-3*nij*kl] = (-wz+(-kl*la[ind]*(ux+vy)+kl*h*bforce6[2+3*qq])/
					 (2*mu[ind]+la[ind]))/sbop[0];
	       }
	 }
      }
   }
   SYNC_DEVICE;
}

//-----------------------------------------------------------------------
void EW::bcfortsg_indrev( int ib, int ie, int jb, int je, int kb, int ke, int wind[36], 
		   int nx, int ny, int nz, float_sw4* u, float_sw4 h, boundaryConditionType bccnd[6],
		   float_sw4 sbop[5], float_sw4* mu, float_sw4* la, float_sw4 t,
		   float_sw4* bforce1, float_sw4* bforce2, float_sw4* bforce3, 
		   float_sw4* bforce4, float_sw4* bforce5, float_sw4* bforce6,
		   float_sw4 om, float_sw4 ph, float_sw4 cv,
		   float_sw4* strx, float_sw4* stry )
{
  PUSH_RANGE("EW::bcfortsg_indrev",5);
   const float_sw4 d4a = 2.0/3.0;
   const float_sw4 d4b = -1.0/12.0;
   const size_t ni  = ie-ib+1;
   const size_t nij = ni*(je-jb+1);
   const size_t npts = static_cast<size_t>((ie-ib+1))*(je-jb+1)*(ke-kb+1);
   for( int s=0 ; s < 6 ; s++ )
   {
     int is = wind[0+6*s];
     int ie = wind[1+6*s];
     int js = wind[2+6*s];
     int je = wind[3+6*s];
     int ks = wind[4+6*s];
     int ke = wind[5+6*s];
      if( bccnd[s]==1 || bccnd[s]==2 )
      {
	PUSH_RANGE_PAYLOAD("bcfortsg_indrev_RAJA",4,s);
         size_t idel = 1+wind[1+6*s]-wind[6*s];
         size_t ijdel = idel * (1+wind[3+6*s]-wind[2+6*s]);
	 if( s== 0 )
	 {
	   PUSH_RANGE("LOOP-0",0);
	   forallN<EXEC_BC, int, int,int>(
					  RangeSegment(ks,ke+1),
					  RangeSegment(js,je+1),
					  RangeSegment(is,ie+1),
					  [=] RAJA_DEVICE(int k, int j, int i) {
					    size_t ind = i-ib+ni*(j-jb)+nij*(k-kb);
					    size_t qq=(k-ks)*ijdel+(i-is)+(j-js)*(ie-is+1);
					    u[ind  ]      = bforce1[  3*qq];
					    u[ind+npts]   = bforce1[1+3*qq];
					    u[ind+2*npts] = bforce1[2+3*qq];
					  });
	   SYNC_DEVICE;POP_RANGE;
	 }
	 else if( s== 1 )
	 {
	   PUSH_RANGE("LOOP-1",1);
	   forallN<EXEC_BC, int, int,int>(
					  RangeSegment(ks,ke+1),
					  RangeSegment(js,je+1),
					  RangeSegment(is,ie+1),
					  [=] RAJA_DEVICE(int k, int j, int i) {
					    size_t ind = i-ib+ni*(j-jb)+nij*(k-kb);
					    size_t qq=(k-ks)*ijdel+(i-is)+(j-js)*(ie-is+1);
					    u[ind]        = bforce2[  3*qq];
					    u[ind+npts]   = bforce2[1+3*qq];
					    u[ind+2*npts] = bforce2[2+3*qq];
		    
					  });
	   SYNC_DEVICE;POP_RANGE;
	 }
	 else if( s==2 )
	 {
	   forallN<EXEC_BC, int, int,int>(
					  RangeSegment(ks,ke+1),
					  RangeSegment(js,je+1),
					  RangeSegment(is,ie+1),
					  [=] RAJA_DEVICE(int k, int j, int i){
					    size_t ind = i-ib+ni*(j-jb)+nij*(k-kb);
					    size_t qq=(k-ks)*ijdel+(i-is)+(j-js)*(ie-is+1);
					    u[ind  ] = bforce3[  3*qq];
					    u[ind+npts] = bforce3[1+3*qq];
					    u[ind+2*npts] = bforce3[2+3*qq];
					  });
	 }
	 else if( s==3 )
	 {
	   forallN<EXEC_BC, int, int,int>(
					  RangeSegment(ks,ke+1),
					  RangeSegment(js,je+1),
					  RangeSegment(is,ie+1),
					  [=] RAJA_DEVICE(int k, int j, int i){
					    size_t ind = i-ib+ni*(j-jb)+nij*(k-kb);
					    size_t qq=(k-ks)*ijdel+(i-is)+(j-js)*(ie-is+1);
					    u[ind  ] = bforce4[  3*qq];
					    u[ind+npts] = bforce4[1+3*qq];
					    u[ind+2*npts] = bforce4[2+3*qq];
					  });
	 }
	 else if( s==4 )
	 {
	   forallN<EXEC_BC, int, int,int>(
					  RangeSegment(ks,ke+1),
					  RangeSegment(js,je+1),
					  RangeSegment(is,ie+1),
					  [=] RAJA_DEVICE(int k, int j, int i) {
					    size_t ind = i-ib+ni*(j-jb)+nij*(k-kb);
					    size_t qq=(k-ks)*ijdel+(i-is)+(j-js)*(ie-is+1);
					    u[ind  ] = bforce5[  3*qq];
					    u[ind+npts] = bforce5[1+3*qq];
					    u[ind+2*npts] = bforce5[2+3*qq];
					  });

	 }
	 else if( s==5 )
	   forallN<EXEC_BC, int, int,int>(
					  RangeSegment(ks,ke+1),
					  RangeSegment(js,je+1),
					  RangeSegment(is,ie+1),
					  [=] RAJA_DEVICE(int k, int j, int i) {
					    size_t ind = i-ib+ni*(j-jb)+nij*(k-kb);
					    size_t qq=(k-ks)*ijdel+(i-is)+(j-js)*(ie-is+1);
					    u[ind  ] = bforce6[  3*qq];
					    u[ind+npts] = bforce6[1+3*qq];
					    u[ind+2*npts] = bforce6[2+3*qq];
					  });
	 POP_RANGE;
      }
      else if( bccnd[s]==3 )
      {
	std::cout<<"ENABLE RAJA IN THIS SECTION 1\n";
	 if( s==0 )
	 {
            #pragma omp parallel for
	    for( int k=wind[4+6*s]; k <= wind[5+6*s] ; k++ )
	       for( int j=wind[2+6*s]; j <= wind[3+6*s] ; j++ )
		  for( int i=wind[6*s]; i <= wind[1+6*s] ; i++ )
		  {
		     size_t ind = i-ib+ni*(j-jb)+nij*(k-kb);
		     size_t indp= ind+nx;
		     u[ind  ] = u[indp];
		     u[ind+npts] = u[indp+npts];
		     u[ind+2*npts] = u[indp+2*npts];
		  }
	 }
	 else if( s==1 )
	 {
            #pragma omp parallel for 
	    for( int k=wind[4+6*s]; k <= wind[5+6*s] ; k++ )
	       for( int j=wind[2+6*s]; j <= wind[3+6*s] ; j++ )
		  for( int i=wind[6*s]; i <= wind[1+6*s] ; i++ )
		  {
		     size_t ind = i-ib+ni*(j-jb)+nij*(k-kb);
		     size_t indp= ind-nx;
		     u[ind  ] = u[indp];
		     u[ind+npts] = u[indp+npts];
		     u[ind+2*npts] = u[indp+2*npts];
		  }
	 }
	 else if( s==2 )
	 {
            #pragma omp parallel for
	    for( int k=wind[4+6*s]; k <= wind[5+6*s] ; k++ )
	       for( int j=wind[2+6*s]; j <= wind[3+6*s] ; j++ )
		  for( int i=wind[6*s]; i <= wind[1+6*s] ; i++ )
		  {
		     size_t ind = i-ib+ni*(j-jb)+nij*(k-kb);
		     size_t indp= ind+ni*ny;
		     u[ind  ] = u[indp];
		     u[ind+npts] = u[indp+npts];
		     u[ind+2*npts] = u[indp+2*npts];
		  }
	 }
	 else if( s==3 )
	 {
            #pragma omp parallel for 
	    for( int k=wind[4+6*s]; k <= wind[5+6*s] ; k++ )
	       for( int j=wind[2+6*s]; j <= wind[3+6*s] ; j++ )
		  for( int i=wind[6*s]; i <= wind[1+6*s] ; i++ )
		  {
		     size_t ind = i-ib+ni*(j-jb)+nij*(k-kb);
		     size_t indp= ind-ni*ny;
		     u[ind  ] = u[indp];
		     u[ind+npts] = u[indp+npts];
		     u[ind+2*npts] = u[indp+2*npts];
		  }
	 }
	 else if( s==4 )
	 {
            #pragma omp parallel for
	    for( int k=wind[4+6*s]; k <= wind[5+6*s] ; k++ )
	       for( int j=wind[2+6*s]; j <= wind[3+6*s] ; j++ )
		  for( int i=wind[6*s]; i <= wind[1+6*s] ; i++ )
		  {
		     size_t ind = i-ib+ni*(j-jb)+nij*(k-kb);
		     size_t indp= ind+nij*nz;
		     u[ind  ] = u[indp];
		     u[ind+npts] = u[indp+npts];
		     u[ind+2*npts] = u[indp+2*npts];
		  }
	 }
	 else if( s==5 )
	 {
            #pragma omp parallel for 
	    for( int k=wind[4+6*s]; k <= wind[5+6*s] ; k++ )
	       for( int j=wind[2+6*s]; j <= wind[3+6*s] ; j++ )
		  for( int i=wind[6*s]; i <= wind[1+6*s] ; i++ )
		  {
		     size_t ind = i-ib+ni*(j-jb)+nij*(k-kb);
		     size_t indp= ind-nij*nz;
		     u[ind  ] = u[indp];
		     u[ind+npts] = u[indp+npts];
		     u[ind+2*npts] = u[indp+2*npts];
		  }
	 }
      }
      else if( bccnd[s]==0 )
      {
	REQUIRE2( s == 4 || s == 5, "EW::bcfortsg_indrev,  ERROR: Free surface condition"
		  << " not implemented for side " << s << endl);
	 if( s==4 )
	 {
	    int k=1, kl=1;
	    forallN<EXEC_BC2, int, int>(
					  RangeSegment(jb+2,je-1),
					  RangeSegment(ib+2,ie-1),
					  [=] RAJA_DEVICE(int j, int i) {
		  size_t qq = i-ib+ni*(j-jb);
		  size_t ind = i-ib+ni*(j-jb)+nij*(k-kb);
		  float_sw4 wx = strx[i-ib]*(d4a*(u[2*npts+ind+1]-u[2*npts+ind-1])+d4b*(u[2*npts+ind+2]-u[2*npts+ind-2]));
		  float_sw4 ux = strx[i-ib]*(d4a*(u[ind+1]-u[ind-1])+d4b*(u[ind+2]-u[ind-2]));
		  float_sw4 wy = stry[j-jb]*(d4a*(u[2*npts+ind+ni  ]-u[2*npts+ind-ni  ])+
					     d4b*(u[2*npts+ind+2*ni]-u[2*npts+ind-2*ni]));
		  float_sw4 vy = stry[j-jb]*(d4a*(u[npts+ind+ni]  -u[npts+ind-ni])+
					     d4b*(u[npts+ind+2*ni]-u[npts+ind-2*ni]));
		  float_sw4 uz=0, vz=0, wz=0;
		  for( int w=1 ; w <= 4 ; w++ )
		  {
		     uz += sbop[w]*u[       ind+nij*kl*(w-1)];
		     vz += sbop[w]*u[npts  +ind+nij*kl*(w-1)];
		     wz += sbop[w]*u[2*npts+ind+nij*kl*(w-1)];
		  }
		  u[       ind-nij*kl] = (-uz-kl*wx+kl*h*bforce5[  3*qq]/mu[ind])/sbop[0];
		  u[npts  +ind-nij*kl] = (-vz-kl*wy+kl*h*bforce5[1+3*qq]/mu[ind])/sbop[0];
		  u[2*npts+ind-nij*kl] = (-wz + (-kl*la[ind]*(ux+vy)+kl*h*bforce5[2+3*qq])/
					 (2*mu[ind]+la[ind]))/sbop[0];
					  });
	 }
	 else
	 {
	    int k=nz, kl=-1;

	       forallN<EXEC_BC2, int, int>(
					  RangeSegment(jb+2,je-1),
					  RangeSegment(ib+2,ie-1),
					  [=] RAJA_DEVICE(int j, int i) {
		  size_t qq = i-ib+ni*(j-jb);
		  size_t ind = i-ib+ni*(j-jb)+nij*(k-kb);
		  float_sw4 wx = strx[i-ib]*(d4a*(u[2*npts+ind+1]-u[2*npts+ind-1])+d4b*(u[2*npts+ind+2]-u[2*npts+ind-2]));
		  float_sw4 ux = strx[i-ib]*(d4a*(u[ind+1]-u[ind-1])+d4b*(u[ind+2]-u[ind-2]));
		  float_sw4 wy = stry[j-jb]*(d4a*(u[2*npts+ind+ni  ]-u[2*npts+ind-ni  ])+
					     d4b*(u[2*npts+ind+2*ni]-u[2*npts+ind-2*ni]));
		  float_sw4 vy = stry[j-jb]*(d4a*(u[npts+ind+ni]  -u[npts+ind-ni])+
					     d4b*(u[npts+ind+2*ni]-u[npts+ind-2*ni]));
		  float_sw4 uz=0, vz=0, wz=0;
		  for( int w=1 ; w <= 4 ; w++ )
		  {
		     uz += sbop[w]*u[       ind+nij*kl*(w-1)];
		     vz += sbop[w]*u[npts  +ind+nij*kl*(w-1)];
		     wz += sbop[w]*u[2*npts+ind+nij*kl*(w-1)];
		  }
		  u[       ind-nij*kl] = (-uz-kl*wx+kl*h*bforce6[  3*qq]/mu[ind])/sbop[0];
		  u[npts  +ind-nij*kl] = (-vz-kl*wy+kl*h*bforce6[1+3*qq]/mu[ind])/sbop[0];
		  u[2*npts+ind-nij*kl] = (-wz+(-kl*la[ind]*(ux+vy)+kl*h*bforce6[2+3*qq])/
					 (2*mu[ind]+la[ind]))/sbop[0];
					  });
	 }
      }
   }
   SYNC_DEVICE;
   POP_RANGE;
}

//----------------------------------------------------------------------- 15%
void EW::addsgd4fort( int ifirst, int ilast, int jfirst, int jlast,
		      int kfirst, int klast,
		      float_sw4* a_up, const float_sw4* __restrict__ a_u, const float_sw4* __restrict__ a_um, const float_sw4* __restrict__ a_rho,
		      const float_sw4* __restrict__ a_dcx,  const float_sw4* __restrict__ a_dcy,  const float_sw4* __restrict__ a_dcz,
		      const float_sw4* __restrict__ a_strx, const float_sw4* __restrict__ a_stry, const float_sw4* __restrict__ a_strz,
		      const float_sw4* __restrict__ a_cox,  const float_sw4* __restrict__ a_coy,  const float_sw4* __restrict__ a_coz,
		      float_sw4 beta )
{
   if( beta != 0 )
   {
#define rho(i,j,k) a_rho[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)]
#define up(c,i,j,k) a_up[c + 3*(i-ifirst)+3*ni*(j-jfirst)+3*nij*(k-kfirst)]
#define u(c,i,j,k)   a_u[c + 3*(i-ifirst)+3*ni*(j-jfirst)+3*nij*(k-kfirst)]
#define um(c,i,j,k) a_um[c + 3*(i-ifirst)+3*ni*(j-jfirst)+3*nij*(k-kfirst)]
#define strx(i) a_strx[(i-ifirst)]
#define dcx(i) a_dcx[(i-ifirst)]
#define cox(i) a_cox[(i-ifirst)]
#define stry(j) a_stry[(j-jfirst)]
#define dcy(j) a_dcy[(j-jfirst)]
#define coy(j) a_coy[(j-jfirst)]
#define strz(k) a_strz[(k-kfirst)]
#define dcz(k) a_dcz[(k-kfirst)]
#define coz(k) a_coz[(k-kfirst)]

     //prefetch(const_cast<void*>(static_cast<const void*>((a_rho))));
     // prefetch(const_cast<void*>(static_cast<const void*>((a_up))));
     //prefetch(const_cast<void*>(static_cast<const void*>((a_u))));
     //prefetch(const_cast<void*>(static_cast<const void*>((a_um))));
     //prefetch(const_cast<void*>(static_cast<const void*>((a_strx))));
     //prefetch(const_cast<void*>(static_cast<const void*>((a_dcx))));
     //prefetch(const_cast<void*>(static_cast<const void*>((a_cox))));
     //prefetch(const_cast<void*>(static_cast<const void*>((a_stry))));
     //prefetch(const_cast<void*>(static_cast<const void*>((a_dcy))));
     //prefetch(const_cast<void*>(static_cast<const void*>((a_coy))));
     //prefetch(const_cast<void*>(static_cast<const void*>((a_strz))));
     //prefetch(const_cast<void*>(static_cast<const void*>((a_dcz))));
     //prefetch(const_cast<void*>(static_cast<const void*>((a_coz))));
     
     const size_t ni = ilast-ifirst+1;
      const size_t nij = ni*(jlast-jfirst+1);

      forallN<EXEC, int, int,int>(
				  RangeSegment(kfirst+2,klast-1),
				  RangeSegment(jfirst+2,jlast-1),
				  RangeSegment(ifirst+2,ilast-1),
				  [=] RAJA_DEVICE(int k, int j, int i) {

	       float_sw4 birho=beta/rho(i,j,k);
#pragma unroll 
	       for( int c=0 ; c < 3 ; c++ )
	       {
		  up(c,i,j,k) -= birho*( 
		  // x-differences
		   strx(i)*coy(j)*coz(k)*(
       rho(i+1,j,k)*dcx(i+1)*
                   ( u(c,i+2,j,k) -2*u(c,i+1,j,k)+ u(c,i,  j,k))
      -2*rho(i,j,k)*dcx(i)  *
                   ( u(c,i+1,j,k) -2*u(c,i,  j,k)+ u(c,i-1,j,k))
      +rho(i-1,j,k)*dcx(i-1)*
                   ( u(c,i,  j,k) -2*u(c,i-1,j,k)+ u(c,i-2,j,k)) 
      -rho(i+1,j,k)*dcx(i+1)*
                   (um(c,i+2,j,k)-2*um(c,i+1,j,k)+um(c,i,  j,k)) 
      +2*rho(i,j,k)*dcx(i)  *
                   (um(c,i+1,j,k)-2*um(c,i,  j,k)+um(c,i-1,j,k)) 
      -rho(i-1,j,k)*dcx(i-1)*
                   (um(c,i,  j,k)-2*um(c,i-1,j,k)+um(c,i-2,j,k)) ) +
// y-differences
      stry(j)*cox(i)*coz(k)*(
      +rho(i,j+1,k)*dcy(j+1)*
                   ( u(c,i,j+2,k) -2*u(c,i,j+1,k)+ u(c,i,j,  k)) 
      -2*rho(i,j,k)*dcy(j)  *
                   ( u(c,i,j+1,k) -2*u(c,i,j,  k)+ u(c,i,j-1,k))
      +rho(i,j-1,k)*dcy(j-1)*
                   ( u(c,i,j,  k) -2*u(c,i,j-1,k)+ u(c,i,j-2,k)) 
      -rho(i,j+1,k)*dcy(j+1)*
                   (um(c,i,j+2,k)-2*um(c,i,j+1,k)+um(c,i,j,  k)) 
      +2*rho(i,j,k)*dcy(j)  *
                   (um(c,i,j+1,k)-2*um(c,i,j,  k)+um(c,i,j-1,k)) 
      -rho(i,j-1,k)*dcy(j-1)*
                   (um(c,i,j,  k)-2*um(c,i,j-1,k)+um(c,i,j-2,k)) ) +
       strz(k)*cox(i)*coy(j)*(
// z-differences
      +rho(i,j,k+1)*dcz(k+1)* 
                 ( u(c,i,j,k+2) -2*u(c,i,j,k+1)+ u(c,i,j,k  )) 
      -2*rho(i,j,k)*dcz(k)  *
                 ( u(c,i,j,k+1) -2*u(c,i,j,k  )+ u(c,i,j,k-1))
      +rho(i,j,k-1)*dcz(k-1)*
                 ( u(c,i,j,k  ) -2*u(c,i,j,k-1)+ u(c,i,j,k-2)) 
      -rho(i,j,k+1)*dcz(k+1)*
                 (um(c,i,j,k+2)-2*um(c,i,j,k+1)+um(c,i,j,k  )) 
      +2*rho(i,j,k)*dcz(k)  *
                 (um(c,i,j,k+1)-2*um(c,i,j,k  )+um(c,i,j,k-1)) 
      -rho(i,j,k-1)*dcz(k-1)*
                 (um(c,i,j,k  )-2*um(c,i,j,k-1)+um(c,i,j,k-2)) ) 
					 );

	       }
				  });
      SYNC_DEVICE
 
#undef rho
#undef up
#undef u
#undef um
#undef strx
#undef dcx
#undef cox
#undef stry
#undef dcy
#undef coy
#undef strz
#undef dcz
#undef coz
   }
}

//-----------------------------------------------------------------------
void EW::addsgd6fort( int ifirst, int ilast, int jfirst, int jlast,
		      int kfirst, int klast,
		      float_sw4* a_up, float_sw4* a_u, float_sw4* a_um, float_sw4* a_rho,
		      float_sw4* a_dcx,  float_sw4* a_dcy,  float_sw4* a_dcz,
		      float_sw4* a_strx, float_sw4* a_stry, float_sw4* a_strz,
		      float_sw4* a_cox,  float_sw4* a_coy,  float_sw4* a_coz,
		      float_sw4 beta )
{
   if( beta != 0 )
   {
#define rho(i,j,k) a_rho[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)]
#define up(c,i,j,k) a_up[c + 3*(i-ifirst)+3*ni*(j-jfirst)+3*nij*(k-kfirst)]
#define u(c,i,j,k) a_u[c + 3*(i-ifirst)+3*ni*(j-jfirst)+3*nij*(k-kfirst)]
#define um(c,i,j,k) a_um[c + 3*(i-ifirst)+3*ni*(j-jfirst)+3*nij*(k-kfirst)]
#define strx(i) a_strx[(i-ifirst)]
#define dcx(i) a_dcx[(i-ifirst)]
#define cox(i) a_cox[(i-ifirst)]
#define stry(j) a_stry[(j-jfirst)]
#define dcy(j) a_dcy[(j-jfirst)]
#define coy(j) a_coy[(j-jfirst)]
#define strz(k) a_strz[(k-kfirst)]
#define dcz(k) a_dcz[(k-kfirst)]
#define coz(k) a_coz[(k-kfirst)]
      const size_t ni = ilast-ifirst+1;
      const size_t nij = ni*(jlast-jfirst+1);
#pragma omp parallel for
      for( int k=kfirst+3; k <= klast-3 ; k++ )
	 for( int j=jfirst+3; j <= jlast-3 ; j++ )
	    for( int i=ifirst+3; i <= ilast-3 ; i++ )
	    {
	       float_sw4 birho=0.5*beta/rho(i,j,k);
#pragma simd
#pragma ivdep
	       for( int c=0 ; c < 3 ; c++ )
	       {
		 up(c,i,j,k) += birho*( 
       strx(i)*coy(j)*coz(k)*(
// x-differences
         (rho(i+2,j,k)*dcx(i+2)+rho(i+1,j,k)*dcx(i+1))*(
         u(c,i+3,j,k) -3*u(c,i+2,j,k)+ 3*u(c,i+1,j,k)- u(c,i, j,k) 
      -(um(c,i+3,j,k)-3*um(c,i+2,j,k)+3*um(c,i+1,j,k)-um(c,i, j,k)) )
      -3*(rho(i+1,j,k)*dcx(i+1)+rho(i,j,k)*dcx(i))*(
         u(c,i+2,j,k)- 3*u(c,i+1,j,k)+ 3*u(c,i, j,k)- u(c,i-1,j,k)
      -(um(c,i+2,j,k)-3*um(c,i+1,j,k)+3*um(c,i, j,k)-um(c,i-1,j,k)) )
      +3*(rho(i,j,k)*dcx(i)+rho(i-1,j,k)*dcx(i-1))*(
         u(c,i+1,j,k)- 3*u(c,i,  j,k)+3*u(c,i-1,j,k)- u(c,i-2,j,k) 
      -(um(c,i+1,j,k)-3*um(c,i, j,k)+3*um(c,i-1,j,k)-um(c,i-2,j,k)) )
       - (rho(i-1,j,k)*dcx(i-1)+rho(i-2,j,k)*dcx(i-2))*(
         u(c,i, j,k)- 3*u(c,i-1,j,k)+ 3*u(c,i-2,j,k)- u(c,i-3,j,k) 
      -(um(c,i, j,k)-3*um(c,i-1,j,k)+3*um(c,i-2,j,k)-um(c,i-3,j,k)) )
                 ) +  stry(j)*cox(i)*coz(k)*(
// y-differences
         (rho(i,j+2,k)*dcy(j+2)+rho(i,j+1,k)*dcy(j+1))*(
         u(c,i,j+3,k) -3*u(c,i,j+2,k)+ 3*u(c,i,j+1,k)- u(c,i,  j,k)
      -(um(c,i,j+3,k)-3*um(c,i,j+2,k)+3*um(c,i,j+1,k)-um(c,i,  j,k)) )
      -3*(rho(i,j+1,k)*dcy(j+1)+rho(i,j,k)*dcy(j))*(
         u(c,i,j+2,k) -3*u(c,i,j+1,k)+ 3*u(c,i,  j,k)- u(c,i,j-1,k) 
      -(um(c,i,j+2,k)-3*um(c,i,j+1,k)+3*um(c,i,  j,k)-um(c,i,j-1,k)) )
      +3*(rho(i,j,k)*dcy(j)+rho(i,j-1,k)*dcy(j-1))*(
         u(c,i,j+1,k)- 3*u(c,i, j,k)+ 3*u(c,i,j-1,k)- u(c,i,j-2,k) 
      -(um(c,i,j+1,k)-3*um(c,i, j,k)+3*um(c,i,j-1,k)-um(c,i,j-2,k)) )
       - (rho(i,j-1,k)*dcy(j-1)+rho(i,j-2,k)*dcy(j-2))*(
         u(c,i, j,k)- 3*u(c,i,j-1,k)+  3*u(c,i,j-2,k)- u(c,i,j-3,k) 
      -(um(c,i, j,k)-3*um(c,i,j-1,k)+ 3*um(c,i,j-2,k)-um(c,i,j-3,k)) )
                 ) +  strz(k)*cox(i)*coy(j)*(
// z-differences
         ( rho(i,j,k+2)*dcz(k+2) + rho(i,j,k+1)*dcz(k+1) )*(
         u(c,i,j,k+3)- 3*u(c,i,j,k+2)+ 3*u(c,i,j,k+1)- u(c,i,  j,k) 
      -(um(c,i,j,k+3)-3*um(c,i,j,k+2)+3*um(c,i,j,k+1)-um(c,i,  j,k)) )
      -3*(rho(i,j,k+1)*dcz(k+1)+rho(i,j,k)*dcz(k))*(
         u(c,i,j,k+2) -3*u(c,i,j,k+1)+ 3*u(c,i,  j,k)- u(c,i,j,k-1) 
      -(um(c,i,j,k+2)-3*um(c,i,j,k+1)+3*um(c,i,  j,k)-um(c,i,j,k-1)) )
      +3*(rho(i,j,k)*dcz(k)+rho(i,j,k-1)*dcz(k-1))*(
         u(c,i,j,k+1)- 3*u(c,i,  j,k)+ 3*u(c,i,j,k-1)-u(c,i,j,k-2) 
      -(um(c,i,j,k+1)-3*um(c,i,  j,k)+3*um(c,i,j,k-1)-um(c,i,j,k-2)) )
       - (rho(i,j,k-1)*dcz(k-1)+rho(i,j,k-2)*dcz(k-2))*(
         u(c,i,  j,k) -3*u(c,i,j,k-1)+ 3*u(c,i,j,k-2)- u(c,i,j,k-3)
      -(um(c,i,  j,k)-3*um(c,i,j,k-1)+3*um(c,i,j,k-2)-um(c,i,j,k-3)) )
					     )  );
	       }
	    }
#undef rho
#undef up
#undef u
#undef um
#undef strx
#undef dcx
#undef cox
#undef stry
#undef dcy
#undef coy
#undef strz
#undef dcz
#undef coz
   }
}

//-----------------------------------------------------------------------
void EW::addsgd4fort_indrev( int ifirst, int ilast, int jfirst, int jlast,
		      int kfirst, int klast,
		      float_sw4* a_up, float_sw4* a_u, float_sw4* a_um, float_sw4* a_rho,
		      float_sw4* a_dcx,  float_sw4* a_dcy,  float_sw4* a_dcz,
		      float_sw4* a_strx, float_sw4* a_stry, float_sw4* a_strz,
		      float_sw4* a_cox,  float_sw4* a_coy,  float_sw4* a_coz,
		      float_sw4 beta )
{

  //std::cout<<"addsgd\n";
   if( beta != 0 )
   {
#define rho(i,j,k) a_rho[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)]
#define up(c,i,j,k) a_up[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define u(c,i,j,k)   a_u[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define um(c,i,j,k) a_um[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define strx(i) a_strx[(i-ifirst)]
#define dcx(i) a_dcx[(i-ifirst)]
#define cox(i) a_cox[(i-ifirst)]
#define stry(j) a_stry[(j-jfirst)]
#define dcy(j) a_dcy[(j-jfirst)]
#define coy(j) a_coy[(j-jfirst)]
#define strz(k) a_strz[(k-kfirst)]
#define dcz(k) a_dcz[(k-kfirst)]
#define coz(k) a_coz[(k-kfirst)]

     //prefetch(const_cast<void*>(static_cast<const void*>((a_rho))));
     //prefetch(const_cast<void*>(static_cast<const void*>((a_up))));
     //prefetch(const_cast<void*>(static_cast<const void*>((a_u))));
     //prefetch(const_cast<void*>(static_cast<const void*>((a_um))));
     //prefetch(const_cast<void*>(static_cast<const void*>((a_strx))));
     //prefetch(const_cast<void*>(static_cast<const void*>((a_dcx))));
     //prefetch(const_cast<void*>(static_cast<const void*>((a_cox))));
     //prefetch(const_cast<void*>(static_cast<const void*>((a_stry))));
     //prefetch(const_cast<void*>(static_cast<const void*>((a_dcy))));
     //prefetch(const_cast<void*>(static_cast<const void*>((a_coy))));
     //prefetch(const_cast<void*>(static_cast<const void*>((a_strz))));
     //prefetch(const_cast<void*>(static_cast<const void*>((a_dcz))));
     //prefetch(const_cast<void*>(static_cast<const void*>((a_coz))));

      const size_t ni = ilast-ifirst+1;
      const size_t nij = ni*(jlast-jfirst+1);
      const size_t npts = nij*(klast-kfirst+1);

// AP: The for c loop could be inside the for i loop. The simd, ivdep pragmas should be outside the inner-most loop
//      for( int c=0 ; c < 3 ; c++ )
	//#pragma omp parallel for
	//      for( int k=kfirst+2; k <= klast-2 ; k++ )
	//	 for( int j=jfirst+2; j <= jlast-2 ; j++ )
	//#pragma simd
	//#pragma ivdep
	//    for( int i=ifirst+2; i <= ilast-2 ; i++ )
	forallN<EXEC, int, int,int>(
				  RangeSegment(kfirst+2,klast-1),
				  RangeSegment(jfirst+2,jlast-1),
				  RangeSegment(ifirst+2,ilast-1),
				  [=] RAJA_DEVICE(int k, int j, int i){
	       float_sw4 birho=beta/rho(i,j,k);
	       for( int c=0 ; c < 3 ; c++ ){
		  up(c,i,j,k) -= birho*( 
		  // x-differences
		   strx(i)*coy(j)*coz(k)*(
       rho(i+1,j,k)*dcx(i+1)*
                   ( u(c,i+2,j,k) -2*u(c,i+1,j,k)+ u(c,i,  j,k))
      -2*rho(i,j,k)*dcx(i)  *
                   ( u(c,i+1,j,k) -2*u(c,i,  j,k)+ u(c,i-1,j,k))
      +rho(i-1,j,k)*dcx(i-1)*
                   ( u(c,i,  j,k) -2*u(c,i-1,j,k)+ u(c,i-2,j,k)) 
      -rho(i+1,j,k)*dcx(i+1)*
                   (um(c,i+2,j,k)-2*um(c,i+1,j,k)+um(c,i,  j,k)) 
      +2*rho(i,j,k)*dcx(i)  *
                   (um(c,i+1,j,k)-2*um(c,i,  j,k)+um(c,i-1,j,k)) 
      -rho(i-1,j,k)*dcx(i-1)*
                   (um(c,i,  j,k)-2*um(c,i-1,j,k)+um(c,i-2,j,k)) ) +
// y-differences
      stry(j)*cox(i)*coz(k)*(
      +rho(i,j+1,k)*dcy(j+1)*
                   ( u(c,i,j+2,k) -2*u(c,i,j+1,k)+ u(c,i,j,  k)) 
      -2*rho(i,j,k)*dcy(j)  *
                   ( u(c,i,j+1,k) -2*u(c,i,j,  k)+ u(c,i,j-1,k))
      +rho(i,j-1,k)*dcy(j-1)*
                   ( u(c,i,j,  k) -2*u(c,i,j-1,k)+ u(c,i,j-2,k)) 
      -rho(i,j+1,k)*dcy(j+1)*
                   (um(c,i,j+2,k)-2*um(c,i,j+1,k)+um(c,i,j,  k)) 
      +2*rho(i,j,k)*dcy(j)  *
                   (um(c,i,j+1,k)-2*um(c,i,j,  k)+um(c,i,j-1,k)) 
      -rho(i,j-1,k)*dcy(j-1)*
                   (um(c,i,j,  k)-2*um(c,i,j-1,k)+um(c,i,j-2,k)) ) +
       strz(k)*cox(i)*coy(j)*(
// z-differences
      +rho(i,j,k+1)*dcz(k+1)* 
                 ( u(c,i,j,k+2) -2*u(c,i,j,k+1)+ u(c,i,j,k  )) 
      -2*rho(i,j,k)*dcz(k)  *
                 ( u(c,i,j,k+1) -2*u(c,i,j,k  )+ u(c,i,j,k-1))
      +rho(i,j,k-1)*dcz(k-1)*
                 ( u(c,i,j,k  ) -2*u(c,i,j,k-1)+ u(c,i,j,k-2)) 
      -rho(i,j,k+1)*dcz(k+1)*
                 (um(c,i,j,k+2)-2*um(c,i,j,k+1)+um(c,i,j,k  )) 
      +2*rho(i,j,k)*dcz(k)  *
                 (um(c,i,j,k+1)-2*um(c,i,j,k  )+um(c,i,j,k-1)) 
      -rho(i,j,k-1)*dcz(k-1)*
                 (um(c,i,j,k  )-2*um(c,i,j,k-1)+um(c,i,j,k-2)) ) 
					 );

	       }
				  }); SYNC_DEVICE
#undef rho
#undef up
#undef u
#undef um
#undef strx
#undef dcx
#undef cox
#undef stry
#undef dcy
#undef coy
#undef strz
#undef dcz
#undef coz
   }
}

//-----------------------------------------------------------------------
void EW::addsgd6fort_indrev( int ifirst, int ilast, int jfirst, int jlast,
		      int kfirst, int klast,
		      float_sw4* a_up, float_sw4* a_u, float_sw4* a_um, float_sw4* a_rho,
		      float_sw4* a_dcx,  float_sw4* a_dcy,  float_sw4* a_dcz,
		      float_sw4* a_strx, float_sw4* a_stry, float_sw4* a_strz,
		      float_sw4* a_cox,  float_sw4* a_coy,  float_sw4* a_coz,
		      float_sw4 beta )
{
   if( beta != 0 )
   {
#define rho(i,j,k) a_rho[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)]
#define up(c,i,j,k) a_up[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define u(c,i,j,k)   a_u[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define um(c,i,j,k) a_um[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define strx(i) a_strx[(i-ifirst)]
#define dcx(i) a_dcx[(i-ifirst)]
#define cox(i) a_cox[(i-ifirst)]
#define stry(j) a_stry[(j-jfirst)]
#define dcy(j) a_dcy[(j-jfirst)]
#define coy(j) a_coy[(j-jfirst)]
#define strz(k) a_strz[(k-kfirst)]
#define dcz(k) a_dcz[(k-kfirst)]
#define coz(k) a_coz[(k-kfirst)]
      const size_t ni = ilast-ifirst+1;
      const size_t nij = ni*(jlast-jfirst+1);
      const size_t npts = nij*(klast-kfirst+1);
      for( int c=0 ; c < 3 ; c++ )
#pragma omp parallel for
      for( int k=kfirst+3; k <= klast-3 ; k++ )
	 for( int j=jfirst+3; j <= jlast-3 ; j++ )
#pragma simd
#pragma ivdep
	    for( int i=ifirst+3; i <= ilast-3 ; i++ )
	    {
	       float_sw4 birho=0.5*beta/rho(i,j,k);
	       {
		 up(c,i,j,k) += birho*( 
       strx(i)*coy(j)*coz(k)*(
// x-differences
         (rho(i+2,j,k)*dcx(i+2)+rho(i+1,j,k)*dcx(i+1))*(
         u(c,i+3,j,k) -3*u(c,i+2,j,k)+ 3*u(c,i+1,j,k)- u(c,i, j,k) 
      -(um(c,i+3,j,k)-3*um(c,i+2,j,k)+3*um(c,i+1,j,k)-um(c,i, j,k)) )
      -3*(rho(i+1,j,k)*dcx(i+1)+rho(i,j,k)*dcx(i))*(
         u(c,i+2,j,k)- 3*u(c,i+1,j,k)+ 3*u(c,i, j,k)- u(c,i-1,j,k)
      -(um(c,i+2,j,k)-3*um(c,i+1,j,k)+3*um(c,i, j,k)-um(c,i-1,j,k)) )
      +3*(rho(i,j,k)*dcx(i)+rho(i-1,j,k)*dcx(i-1))*(
         u(c,i+1,j,k)- 3*u(c,i,  j,k)+3*u(c,i-1,j,k)- u(c,i-2,j,k) 
      -(um(c,i+1,j,k)-3*um(c,i, j,k)+3*um(c,i-1,j,k)-um(c,i-2,j,k)) )
       - (rho(i-1,j,k)*dcx(i-1)+rho(i-2,j,k)*dcx(i-2))*(
         u(c,i, j,k)- 3*u(c,i-1,j,k)+ 3*u(c,i-2,j,k)- u(c,i-3,j,k) 
      -(um(c,i, j,k)-3*um(c,i-1,j,k)+3*um(c,i-2,j,k)-um(c,i-3,j,k)) )
                 ) +  stry(j)*cox(i)*coz(k)*(
// y-differences
         (rho(i,j+2,k)*dcy(j+2)+rho(i,j+1,k)*dcy(j+1))*(
         u(c,i,j+3,k) -3*u(c,i,j+2,k)+ 3*u(c,i,j+1,k)- u(c,i,  j,k)
      -(um(c,i,j+3,k)-3*um(c,i,j+2,k)+3*um(c,i,j+1,k)-um(c,i,  j,k)) )
      -3*(rho(i,j+1,k)*dcy(j+1)+rho(i,j,k)*dcy(j))*(
         u(c,i,j+2,k) -3*u(c,i,j+1,k)+ 3*u(c,i,  j,k)- u(c,i,j-1,k) 
      -(um(c,i,j+2,k)-3*um(c,i,j+1,k)+3*um(c,i,  j,k)-um(c,i,j-1,k)) )
      +3*(rho(i,j,k)*dcy(j)+rho(i,j-1,k)*dcy(j-1))*(
         u(c,i,j+1,k)- 3*u(c,i, j,k)+ 3*u(c,i,j-1,k)- u(c,i,j-2,k) 
      -(um(c,i,j+1,k)-3*um(c,i, j,k)+3*um(c,i,j-1,k)-um(c,i,j-2,k)) )
       - (rho(i,j-1,k)*dcy(j-1)+rho(i,j-2,k)*dcy(j-2))*(
         u(c,i, j,k)- 3*u(c,i,j-1,k)+  3*u(c,i,j-2,k)- u(c,i,j-3,k) 
      -(um(c,i, j,k)-3*um(c,i,j-1,k)+ 3*um(c,i,j-2,k)-um(c,i,j-3,k)) )
                 ) +  strz(k)*cox(i)*coy(j)*(
// z-differences
         ( rho(i,j,k+2)*dcz(k+2) + rho(i,j,k+1)*dcz(k+1) )*(
         u(c,i,j,k+3)- 3*u(c,i,j,k+2)+ 3*u(c,i,j,k+1)- u(c,i,  j,k) 
      -(um(c,i,j,k+3)-3*um(c,i,j,k+2)+3*um(c,i,j,k+1)-um(c,i,  j,k)) )
      -3*(rho(i,j,k+1)*dcz(k+1)+rho(i,j,k)*dcz(k))*(
         u(c,i,j,k+2) -3*u(c,i,j,k+1)+ 3*u(c,i,  j,k)- u(c,i,j,k-1) 
      -(um(c,i,j,k+2)-3*um(c,i,j,k+1)+3*um(c,i,  j,k)-um(c,i,j,k-1)) )
      +3*(rho(i,j,k)*dcz(k)+rho(i,j,k-1)*dcz(k-1))*(
         u(c,i,j,k+1)- 3*u(c,i,  j,k)+ 3*u(c,i,j,k-1)-u(c,i,j,k-2) 
      -(um(c,i,j,k+1)-3*um(c,i,  j,k)+3*um(c,i,j,k-1)-um(c,i,j,k-2)) )
       - (rho(i,j,k-1)*dcz(k-1)+rho(i,j,k-2)*dcz(k-2))*(
         u(c,i,  j,k) -3*u(c,i,j,k-1)+ 3*u(c,i,j,k-2)- u(c,i,j,k-3)
      -(um(c,i,  j,k)-3*um(c,i,j,k-1)+3*um(c,i,j,k-2)-um(c,i,j,k-3)) )
					     )  );
	       }
	    }
#undef rho
#undef up
#undef u
#undef um
#undef strx
#undef dcx
#undef cox
#undef stry
#undef dcy
#undef coy
#undef strz
#undef dcz
#undef coz
   }
}

//-----------------------------------------------------------------------
void EW::addsgd4cfort( int ifirst, int ilast, int jfirst, int jlast,
		       int kfirst, int klast,
		       float_sw4* __restrict__ a_up, float_sw4* __restrict__ a_u, 
		       float_sw4* __restrict__ a_um, float_sw4* __restrict__ a_rho,
		       float_sw4* __restrict__ a_dcx, float_sw4* __restrict__ a_dcy, 
		       float_sw4* __restrict__ a_strx, float_sw4* __restrict__ a_stry, 
		       float_sw4* __restrict__ a_jac, float_sw4* __restrict__ a_cox,  
		       float_sw4* __restrict__ a_coy, float_sw4 beta )
{
   if( beta != 0 )
   {
#define rho(i,j,k) a_rho[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)]
#define up(c,i,j,k) a_up[c + 3*(i-ifirst)+3*ni*(j-jfirst)+3*nij*(k-kfirst)]
#define u(c,i,j,k) a_u[c + 3*(i-ifirst)+3*ni*(j-jfirst)+3*nij*(k-kfirst)]
#define um(c,i,j,k) a_um[c + 3*(i-ifirst)+3*ni*(j-jfirst)+3*nij*(k-kfirst)]
#define jac(i,j,k) a_jac[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)]
#define strx(i) a_strx[(i-ifirst)]
#define dcx(i) a_dcx[(i-ifirst)]
#define cox(i) a_cox[(i-ifirst)]
#define stry(j) a_stry[(j-jfirst)]
#define dcy(j) a_dcy[(j-jfirst)]
#define coy(j) a_coy[(j-jfirst)]

      const size_t ni = ilast-ifirst+1;
      const size_t nij = ni*(jlast-jfirst+1);
#pragma omp parallel for
      for( int k=kfirst+2; k <= klast-2 ; k++ )
	 for( int j=jfirst+2; j <= jlast-2 ; j++ )
	    for( int i=ifirst+2; i <= ilast-2 ; i++ )
	    {
	       float_sw4 irhoj=beta/(rho(i,j,k)*jac(i,j,k));
#pragma simd
#pragma ivdep
	       for( int c=0 ; c < 3 ; c++ )
	       {
		  up(c,i,j,k) -= irhoj*( 
		  // x-differences
		   strx(i)*coy(j)*(
		   rho(i+1,j,k)*dcx(i+1)*jac(i+1,j,k)*
                   ( u(c,i+2,j,k) -2*u(c,i+1,j,k)+ u(c,i,  j,k))
		   -2*rho(i,j,k)*dcx(i)*jac(i,j,k)*
                   ( u(c,i+1,j,k) -2*u(c,i,  j,k)+ u(c,i-1,j,k))
		   +rho(i-1,j,k)*dcx(i-1)*jac(i-1,j,k)*
                   ( u(c,i,  j,k) -2*u(c,i-1,j,k)+ u(c,i-2,j,k)) 
		   -rho(i+1,j,k)*dcx(i+1)*jac(i+1,j,k)*
                   (um(c,i+2,j,k)-2*um(c,i+1,j,k)+um(c,i,  j,k)) 
		   +2*rho(i,j,k)*dcx(i)*jac(i,j,k)*
                   (um(c,i+1,j,k)-2*um(c,i,  j,k)+um(c,i-1,j,k)) 
		   -rho(i-1,j,k)*dcx(i-1)*jac(i-1,j,k)*
                   (um(c,i,  j,k)-2*um(c,i-1,j,k)+um(c,i-2,j,k)) ) +
// y-differences
		   stry(j)*cox(i)*(
		    +rho(i,j+1,k)*dcy(j+1)*jac(i,j+1,k)*
                   ( u(c,i,j+2,k) -2*u(c,i,j+1,k)+ u(c,i,j,  k)) 
		    -2*rho(i,j,k)*dcy(j)*jac(i,j,k)*
                   ( u(c,i,j+1,k) -2*u(c,i,j,  k)+ u(c,i,j-1,k))
		    +rho(i,j-1,k)*dcy(j-1)*jac(i,j-1,k)*
                   ( u(c,i,j,  k) -2*u(c,i,j-1,k)+ u(c,i,j-2,k)) 
		    -rho(i,j+1,k)*dcy(j+1)*jac(i,j+1,k)*
                   (um(c,i,j+2,k)-2*um(c,i,j+1,k)+um(c,i,j,  k)) 
		    +2*rho(i,j,k)*dcy(j)*jac(i,j,k)*
                   (um(c,i,j+1,k)-2*um(c,i,j,  k)+um(c,i,j-1,k)) 
		    -rho(i,j-1,k)*dcy(j-1)*jac(i,j-1,k)*
		    (um(c,i,j,  k)-2*um(c,i,j-1,k)+um(c,i,j-2,k)) ) );
	       }
	    } 
#undef rho
#undef up
#undef u
#undef um
#undef strx
#undef dcx
#undef cox
#undef stry
#undef dcy
#undef coy
#undef jac
   }
}

//-----------------------------------------------------------------------
void EW::addsgd6cfort( int ifirst, int ilast, int jfirst, int jlast,
		       int kfirst, int klast,
		       float_sw4* __restrict__ a_up, float_sw4* __restrict__ a_u,
		       float_sw4* __restrict__ a_um, float_sw4* __restrict__ a_rho,
		       float_sw4* __restrict__ a_dcx,  float_sw4* __restrict__ a_dcy,
		       float_sw4* __restrict__ a_strx, float_sw4* __restrict__ a_stry,
		       float_sw4* __restrict__ a_jac, float_sw4* __restrict__ a_cox,
		       float_sw4* __restrict__ a_coy,
		       float_sw4 beta )
{
   if( beta != 0 )
   {
#define rho(i,j,k) a_rho[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)]
#define up(c,i,j,k) a_up[c + 3*(i-ifirst)+3*ni*(j-jfirst)+3*nij*(k-kfirst)]
#define u(c,i,j,k) a_u[c + 3*(i-ifirst)+3*ni*(j-jfirst)+3*nij*(k-kfirst)]
#define um(c,i,j,k) a_um[c + 3*(i-ifirst)+3*ni*(j-jfirst)+3*nij*(k-kfirst)]
#define jac(i,j,k) a_jac[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)]
#define strx(i) a_strx[(i-ifirst)]
#define dcx(i) a_dcx[(i-ifirst)]
#define cox(i) a_cox[(i-ifirst)]
#define stry(j) a_stry[(j-jfirst)]
#define dcy(j) a_dcy[(j-jfirst)]
#define coy(j) a_coy[(j-jfirst)]
      const size_t ni = ilast-ifirst+1;
      const size_t nij = ni*(jlast-jfirst+1);
#pragma omp parallel for
      for( int k=kfirst+3; k <= klast-3 ; k++ )
	 for( int j=jfirst+3; j <= jlast-3 ; j++ )
	    for( int i=ifirst+3; i <= ilast-3 ; i++ )
	    {
	       float_sw4 birho=0.5*beta/(rho(i,j,k)*jac(i,j,k));
#pragma simd
#pragma ivdep
	       for( int c=0 ; c < 3 ; c++ )
	       {
		 up(c,i,j,k) += birho*( 
       strx(i)*coy(j)*(
// x-differences
      (rho(i+2,j,k)*dcx(i+2)*jac(i+2,j,k)+rho(i+1,j,k)*dcx(i+1)*jac(i+1,j,k))*(
         u(c,i+3,j,k) -3*u(c,i+2,j,k)+ 3*u(c,i+1,j,k)- u(c,i, j,k) 
      -(um(c,i+3,j,k)-3*um(c,i+2,j,k)+3*um(c,i+1,j,k)-um(c,i, j,k)) )
      -3*(rho(i+1,j,k)*dcx(i+1)*jac(i+1,j,k)+rho(i,j,k)*dcx(i)*jac(i,j,k))*(
         u(c,i+2,j,k)- 3*u(c,i+1,j,k)+ 3*u(c,i, j,k)- u(c,i-1,j,k)
      -(um(c,i+2,j,k)-3*um(c,i+1,j,k)+3*um(c,i, j,k)-um(c,i-1,j,k)) )
      +3*(rho(i,j,k)*dcx(i)*jac(i,j,k)+rho(i-1,j,k)*dcx(i-1)*jac(i-1,j,k))*(
         u(c,i+1,j,k)- 3*u(c,i,  j,k)+3*u(c,i-1,j,k)- u(c,i-2,j,k) 
      -(um(c,i+1,j,k)-3*um(c,i, j,k)+3*um(c,i-1,j,k)-um(c,i-2,j,k)) )
      - (rho(i-1,j,k)*dcx(i-1)*jac(i-1,j,k)+rho(i-2,j,k)*dcx(i-2)*jac(i-2,j,k))*(
         u(c,i, j,k)- 3*u(c,i-1,j,k)+ 3*u(c,i-2,j,k)- u(c,i-3,j,k) 
      -(um(c,i, j,k)-3*um(c,i-1,j,k)+3*um(c,i-2,j,k)-um(c,i-3,j,k)) )
		       ) +  stry(j)*cox(i)*( 
// y-differences
     (rho(i,j+2,k)*dcy(j+2)*jac(i,j+2,k)+rho(i,j+1,k)*dcy(j+1)*jac(i,j+1,k))*(
         u(c,i,j+3,k) -3*u(c,i,j+2,k)+ 3*u(c,i,j+1,k)- u(c,i,  j,k)
      -(um(c,i,j+3,k)-3*um(c,i,j+2,k)+3*um(c,i,j+1,k)-um(c,i,  j,k)) )
     -3*(rho(i,j+1,k)*dcy(j+1)*jac(i,j+1,k)+rho(i,j,k)*dcy(j)*jac(i,j,k))*(
         u(c,i,j+2,k) -3*u(c,i,j+1,k)+ 3*u(c,i,  j,k)- u(c,i,j-1,k) 
      -(um(c,i,j+2,k)-3*um(c,i,j+1,k)+3*um(c,i,  j,k)-um(c,i,j-1,k)) )
     +3*(rho(i,j,k)*dcy(j)*jac(i,j,k)+rho(i,j-1,k)*dcy(j-1)*jac(i,j-1,k))*(
         u(c,i,j+1,k)- 3*u(c,i, j,k)+ 3*u(c,i,j-1,k)- u(c,i,j-2,k) 
      -(um(c,i,j+1,k)-3*um(c,i, j,k)+3*um(c,i,j-1,k)-um(c,i,j-2,k)) )
     - (rho(i,j-1,k)*dcy(j-1)*jac(i,j-1,k)+rho(i,j-2,k)*dcy(j-2)*jac(i,j-2,k))*(
         u(c,i, j,k)- 3*u(c,i,j-1,k)+  3*u(c,i,j-2,k)- u(c,i,j-3,k) 
      -(um(c,i, j,k)-3*um(c,i,j-1,k)+ 3*um(c,i,j-2,k)-um(c,i,j-3,k)) )
					     )  );
	       }
	    }
#undef rho
#undef up
#undef u
#undef um
#undef strx
#undef dcx
#undef cox
#undef stry
#undef dcy
#undef coy
#undef jac
   }
}

//-----------------------------------------------------------------------
void EW::addsgd4cfort_indrev( int ifirst, int ilast, int jfirst, int jlast,
			      int kfirst, int klast,
			      float_sw4* __restrict__ a_up, float_sw4* __restrict__ a_u,
			      float_sw4* __restrict__ a_um, float_sw4* __restrict__ a_rho,
			      float_sw4* __restrict__ a_dcx, float_sw4* __restrict__ a_dcy,
			      float_sw4* __restrict__ a_strx, float_sw4* __restrict__ a_stry, 
			      float_sw4* __restrict__ a_jac, float_sw4* __restrict__ a_cox,
			      float_sw4* __restrict__ a_coy, float_sw4 beta )
{
   if( beta != 0 )
   {
#define rho(i,j,k) a_rho[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)]
#define up(c,i,j,k) a_up[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define u(c,i,j,k)   a_u[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define um(c,i,j,k) a_um[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define jac(i,j,k) a_jac[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)]
#define strx(i) a_strx[(i-ifirst)]
#define dcx(i) a_dcx[(i-ifirst)]
#define cox(i) a_cox[(i-ifirst)]
#define stry(j) a_stry[(j-jfirst)]
#define dcy(j) a_dcy[(j-jfirst)]
#define coy(j) a_coy[(j-jfirst)]

      const size_t ni   =     (ilast-ifirst+1);
      const size_t nij  =  ni*(jlast-jfirst+1);
      const size_t npts = nij*(klast-kfirst+1);
      for( int c=0 ; c < 3 ; c++ )
#pragma omp parallel for
      for( int k=kfirst+2; k <= klast-2 ; k++ )
	 for( int j=jfirst+2; j <= jlast-2 ; j++ )
#pragma simd
#pragma ivdep
	    for( int i=ifirst+2; i <= ilast-2 ; i++ )
	    {
	       float_sw4 irhoj=beta/(rho(i,j,k)*jac(i,j,k));
	       {
		  up(c,i,j,k) -= irhoj*( 
		  // x-differences
		   strx(i)*coy(j)*(
		   rho(i+1,j,k)*dcx(i+1)*jac(i+1,j,k)*
                   ( u(c,i+2,j,k) -2*u(c,i+1,j,k)+ u(c,i,  j,k))
		   -2*rho(i,j,k)*dcx(i)*jac(i,j,k)*
                   ( u(c,i+1,j,k) -2*u(c,i,  j,k)+ u(c,i-1,j,k))
		   +rho(i-1,j,k)*dcx(i-1)*jac(i-1,j,k)*
                   ( u(c,i,  j,k) -2*u(c,i-1,j,k)+ u(c,i-2,j,k)) 
		   -rho(i+1,j,k)*dcx(i+1)*jac(i+1,j,k)*
                   (um(c,i+2,j,k)-2*um(c,i+1,j,k)+um(c,i,  j,k)) 
		   +2*rho(i,j,k)*dcx(i)*jac(i,j,k)*
                   (um(c,i+1,j,k)-2*um(c,i,  j,k)+um(c,i-1,j,k)) 
		   -rho(i-1,j,k)*dcx(i-1)*jac(i-1,j,k)*
                   (um(c,i,  j,k)-2*um(c,i-1,j,k)+um(c,i-2,j,k)) ) +
// y-differences
		   stry(j)*cox(i)*(
		    +rho(i,j+1,k)*dcy(j+1)*jac(i,j+1,k)*
                   ( u(c,i,j+2,k) -2*u(c,i,j+1,k)+ u(c,i,j,  k)) 
		    -2*rho(i,j,k)*dcy(j)*jac(i,j,k)*
                   ( u(c,i,j+1,k) -2*u(c,i,j,  k)+ u(c,i,j-1,k))
		    +rho(i,j-1,k)*dcy(j-1)*jac(i,j-1,k)*
                   ( u(c,i,j,  k) -2*u(c,i,j-1,k)+ u(c,i,j-2,k)) 
		    -rho(i,j+1,k)*dcy(j+1)*jac(i,j+1,k)*
                   (um(c,i,j+2,k)-2*um(c,i,j+1,k)+um(c,i,j,  k)) 
		    +2*rho(i,j,k)*dcy(j)*jac(i,j,k)*
                   (um(c,i,j+1,k)-2*um(c,i,j,  k)+um(c,i,j-1,k)) 
		    -rho(i,j-1,k)*dcy(j-1)*jac(i,j-1,k)*
		    (um(c,i,j,  k)-2*um(c,i,j-1,k)+um(c,i,j-2,k)) ) );
	       }
	    } 
#undef rho
#undef up
#undef u
#undef um
#undef strx
#undef dcx
#undef cox
#undef stry
#undef dcy
#undef coy
#undef jac
   }
}

//-----------------------------------------------------------------------
void EW::addsgd6cfort_indrev(  int ifirst, int ilast, int jfirst, int jlast,
			       int kfirst, int klast,
			       float_sw4* __restrict__ a_up, float_sw4* __restrict__ a_u,
			       float_sw4* __restrict__ a_um, float_sw4* __restrict__ a_rho,
			       float_sw4* __restrict__ a_dcx, float_sw4* __restrict__ a_dcy,
			       float_sw4* __restrict__ a_strx, float_sw4* __restrict__ a_stry, 
			       float_sw4* __restrict__ a_jac, float_sw4* __restrict__ a_cox,
			       float_sw4* __restrict__ a_coy, float_sw4 beta )
{
   if( beta != 0 )
   {
#define rho(i,j,k) a_rho[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)]
#define up(c,i,j,k) a_up[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define u(c,i,j,k)   a_u[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define um(c,i,j,k) a_um[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define jac(i,j,k) a_jac[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)]
#define strx(i) a_strx[(i-ifirst)]
#define dcx(i) a_dcx[(i-ifirst)]
#define cox(i) a_cox[(i-ifirst)]
#define stry(j) a_stry[(j-jfirst)]
#define dcy(j) a_dcy[(j-jfirst)]
#define coy(j) a_coy[(j-jfirst)]
      const size_t ni = ilast-ifirst+1;
      const size_t nij = ni*(jlast-jfirst+1);
      const size_t npts = nij*(klast-kfirst+1);
	       for( int c=0 ; c < 3 ; c++ )
#pragma omp parallel for
      for( int k=kfirst+3; k <= klast-3 ; k++ )
	 for( int j=jfirst+3; j <= jlast-3 ; j++ )
#pragma simd
#pragma ivdep
	    for( int i=ifirst+3; i <= ilast-3 ; i++ )
	    {
	       float_sw4 birho=0.5*beta/(rho(i,j,k)*jac(i,j,k));
	       {
		 up(c,i,j,k) += birho*( 
       strx(i)*coy(j)*(
// x-differences
      (rho(i+2,j,k)*dcx(i+2)*jac(i+2,j,k)+rho(i+1,j,k)*dcx(i+1)*jac(i+1,j,k))*(
         u(c,i+3,j,k) -3*u(c,i+2,j,k)+ 3*u(c,i+1,j,k)- u(c,i, j,k) 
      -(um(c,i+3,j,k)-3*um(c,i+2,j,k)+3*um(c,i+1,j,k)-um(c,i, j,k)) )
      -3*(rho(i+1,j,k)*dcx(i+1)*jac(i+1,j,k)+rho(i,j,k)*dcx(i)*jac(i,j,k))*(
         u(c,i+2,j,k)- 3*u(c,i+1,j,k)+ 3*u(c,i, j,k)- u(c,i-1,j,k)
      -(um(c,i+2,j,k)-3*um(c,i+1,j,k)+3*um(c,i, j,k)-um(c,i-1,j,k)) )
      +3*(rho(i,j,k)*dcx(i)*jac(i,j,k)+rho(i-1,j,k)*dcx(i-1)*jac(i-1,j,k))*(
         u(c,i+1,j,k)- 3*u(c,i,  j,k)+3*u(c,i-1,j,k)- u(c,i-2,j,k) 
      -(um(c,i+1,j,k)-3*um(c,i, j,k)+3*um(c,i-1,j,k)-um(c,i-2,j,k)) )
      - (rho(i-1,j,k)*dcx(i-1)*jac(i-1,j,k)+rho(i-2,j,k)*dcx(i-2)*jac(i-2,j,k))*(
         u(c,i, j,k)- 3*u(c,i-1,j,k)+ 3*u(c,i-2,j,k)- u(c,i-3,j,k) 
      -(um(c,i, j,k)-3*um(c,i-1,j,k)+3*um(c,i-2,j,k)-um(c,i-3,j,k)) )
		       ) +  stry(j)*cox(i)*( 
// y-differences
     (rho(i,j+2,k)*dcy(j+2)*jac(i,j+2,k)+rho(i,j+1,k)*dcy(j+1)*jac(i,j+1,k))*(
         u(c,i,j+3,k) -3*u(c,i,j+2,k)+ 3*u(c,i,j+1,k)- u(c,i,  j,k)
      -(um(c,i,j+3,k)-3*um(c,i,j+2,k)+3*um(c,i,j+1,k)-um(c,i,  j,k)) )
     -3*(rho(i,j+1,k)*dcy(j+1)*jac(i,j+1,k)+rho(i,j,k)*dcy(j)*jac(i,j,k))*(
         u(c,i,j+2,k) -3*u(c,i,j+1,k)+ 3*u(c,i,  j,k)- u(c,i,j-1,k) 
      -(um(c,i,j+2,k)-3*um(c,i,j+1,k)+3*um(c,i,  j,k)-um(c,i,j-1,k)) )
     +3*(rho(i,j,k)*dcy(j)*jac(i,j,k)+rho(i,j-1,k)*dcy(j-1)*jac(i,j-1,k))*(
         u(c,i,j+1,k)- 3*u(c,i, j,k)+ 3*u(c,i,j-1,k)- u(c,i,j-2,k) 
      -(um(c,i,j+1,k)-3*um(c,i, j,k)+3*um(c,i,j-1,k)-um(c,i,j-2,k)) )
     - (rho(i,j-1,k)*dcy(j-1)*jac(i,j-1,k)+rho(i,j-2,k)*dcy(j-2)*jac(i,j-2,k))*(
         u(c,i, j,k)- 3*u(c,i,j-1,k)+  3*u(c,i,j-2,k)- u(c,i,j-3,k) 
      -(um(c,i, j,k)-3*um(c,i,j-1,k)+ 3*um(c,i,j-2,k)-um(c,i,j-3,k)) )
					     )  );
	       }
	    }
#undef rho
#undef up
#undef u
#undef um
#undef strx
#undef dcx
#undef cox
#undef stry
#undef dcy
#undef coy
#undef jac
   }
}

//-----------------------------------------------------------------------
void EW::GetStencilCoefficients( float_sw4* _acof, float_sw4* _ghcof,
				 float_sw4* _bope, float_sw4* _sbop )
{
#define acof(q,k,m) (_acof[q-1+6*(k-1)+48*(m-1)])
#define ghcof(k) (_ghcof[k-1])
#define bope(q,k) (_bope[q-1+6*(k-1)])

   ghcof(1) = 12.0/17;
   ghcof(2) = 0;
   ghcof(3) = 0;
   ghcof(4) = 0;
   ghcof(5) = 0;
   ghcof(6) = 0;

   acof(1,1,1) = 104.0/289.0;
   acof(1,1,2) = -2476335.0/2435692.0;
   acof(1,1,3) = -16189.0/84966.0;
   acof(1,1,4) = -9.0/3332.0;
   acof(1,1,5) = 0;
   acof(1,1,6) = 0;
   acof(1,1,7) = 0;
   acof(1,1,8) = 0;
   acof(1,2,1) = -516.0/289.0;
   acof(1,2,2) = 544521.0/1217846.0;
   acof(1,2,3) = 2509879.0/3653538.0;
   acof(1,2,4) = 0;
   acof(1,2,5) = 0;
   acof(1,2,6) = 0;
   acof(1,2,7) = 0;
   acof(1,2,8) = 0;
   acof(1,3,1) = 312.0/289.0;
   acof(1,3,2) = 1024279.0/2435692.0;
   acof(1,3,3) = -687797.0/1217846.0;
   acof(1,3,4) = 177.0/3332.0;
   acof(1,3,5) = 0;
   acof(1,3,6) = 0;
   acof(1,3,7) = 0;
   acof(1,3,8) = 0;
   acof(1,4,1) = -104.0/289.0;
   acof(1,4,2) = 181507.0/1217846.0;
   acof(1,4,3) = 241309.0/3653538.0;
   acof(1,4,4) = 0;
   acof(1,4,5) = 0;
   acof(1,4,6) = 0;
   acof(1,4,7) = 0;
   acof(1,4,8) = 0;
   acof(1,5,1) = 0;
   acof(1,5,2) = 0;
   acof(1,5,3) = 5.0/2193.0;
   acof(1,5,4) = -48.0/833.0;
   acof(1,5,5) = 0;
   acof(1,5,6) = 0;
   acof(1,5,7) = 0;
   acof(1,5,8) = 0;
   acof(1,6,1) = 0;
   acof(1,6,2) = 0;
   acof(1,6,3) = 0;
   acof(1,6,4) = 6.0/833.0;
   acof(1,6,5) = 0;
   acof(1,6,6) = 0;
   acof(1,6,7) = 0;
   acof(1,6,8) = 0;
   acof(1,7,1) = 0;
   acof(1,7,2) = 0;
   acof(1,7,3) = 0;
   acof(1,7,4) = 0;
   acof(1,7,5) = 0;
   acof(1,7,6) = 0;
   acof(1,7,7) = 0;
   acof(1,7,8) = 0;
   acof(1,8,1) = 0;
   acof(1,8,2) = 0;
   acof(1,8,3) = 0;
   acof(1,8,4) = 0;
   acof(1,8,5) = 0;
   acof(1,8,6) = 0;
   acof(1,8,7) = 0;
   acof(1,8,8) = 0;
   acof(2,1,1) = 12.0/17.0;
   acof(2,1,2) = 544521.0/4226642.0;
   acof(2,1,3) = 2509879.0/12679926.0;
   acof(2,1,4) = 0;
   acof(2,1,5) = 0;
   acof(2,1,6) = 0;
   acof(2,1,7) = 0;
   acof(2,1,8) = 0;
   acof(2,2,1) = -59.0/68.0;
   acof(2,2,2) = -1633563.0/4226642.0;
   acof(2,2,3) = -21510077.0/25359852.0;
   acof(2,2,4) = -12655.0/372939.0;
   acof(2,2,5) = 0;
   acof(2,2,6) = 0;
   acof(2,2,7) = 0;
   acof(2,2,8) = 0;
   acof(2,3,1) = 2.0/17.0;
   acof(2,3,2) = 1633563.0/4226642.0;
   acof(2,3,3) = 2565299.0/4226642.0;
   acof(2,3,4) = 40072.0/372939.0;
   acof(2,3,5) = 0;
   acof(2,3,6) = 0;
   acof(2,3,7) = 0;
   acof(2,3,8) = 0;
   acof(2,4,1) = 3.0/68.0;
   acof(2,4,2) = -544521.0/4226642.0;
   acof(2,4,3) = 987685.0/25359852.0;
   acof(2,4,4) = -14762.0/124313.0;
   acof(2,4,5) = 0;
   acof(2,4,6) = 0;
   acof(2,4,7) = 0;
   acof(2,4,8) = 0;
   acof(2,5,1) = 0;
   acof(2,5,2) = 0;
   acof(2,5,3) = 1630.0/372939.0;
   acof(2,5,4) = 18976.0/372939.0;
   acof(2,5,5) = 0;
   acof(2,5,6) = 0;
   acof(2,5,7) = 0;
   acof(2,5,8) = 0;
   acof(2,6,1) = 0;
   acof(2,6,2) = 0;
   acof(2,6,3) = 0;
   acof(2,6,4) = -1.0/177.0;
   acof(2,6,5) = 0;
   acof(2,6,6) = 0;
   acof(2,6,7) = 0;
   acof(2,6,8) = 0;
   acof(2,7,1) = 0;
   acof(2,7,2) = 0;
   acof(2,7,3) = 0;
   acof(2,7,4) = 0;
   acof(2,7,5) = 0;
   acof(2,7,6) = 0;
   acof(2,7,7) = 0;
   acof(2,7,8) = 0;
   acof(2,8,1) = 0;
   acof(2,8,2) = 0;
   acof(2,8,3) = 0;
   acof(2,8,4) = 0;
   acof(2,8,5) = 0;
   acof(2,8,6) = 0;
   acof(2,8,7) = 0;
   acof(2,8,8) = 0;
   acof(3,1,1) = -96.0/731.0;
   acof(3,1,2) = 1024279.0/6160868.0;
   acof(3,1,3) = -687797.0/3080434.0;
   acof(3,1,4) = 177.0/8428.0;
   acof(3,1,5) = 0;
   acof(3,1,6) = 0;
   acof(3,1,7) = 0;
   acof(3,1,8) = 0;
   acof(3,2,1) = 118.0/731.0;
   acof(3,2,2) = 1633563.0/3080434.0;
   acof(3,2,3) = 2565299.0/3080434.0;
   acof(3,2,4) = 40072.0/271803.0;
   acof(3,2,5) = 0;
   acof(3,2,6) = 0;
   acof(3,2,7) = 0;
   acof(3,2,8) = 0;
   acof(3,3,1) = -16.0/731.0;
   acof(3,3,2) = -5380447.0/6160868.0;
   acof(3,3,3) = -3569115.0/3080434.0;
   acof(3,3,4) = -331815.0/362404.0;
   acof(3,3,5) = -283.0/6321.0;
   acof(3,3,6) = 0;
   acof(3,3,7) = 0;
   acof(3,3,8) = 0;
   acof(3,4,1) = -6.0/731.0;
   acof(3,4,2) = 544521.0/3080434.0;
   acof(3,4,3) = 2193521.0/3080434.0;
   acof(3,4,4) = 8065.0/12943.0;
   acof(3,4,5) = 381.0/2107.0;
   acof(3,4,6) = 0;
   acof(3,4,7) = 0;
   acof(3,4,8) = 0;
   acof(3,5,1) = 0;
   acof(3,5,2) = 0;
   acof(3,5,3) = -14762.0/90601.0;
   acof(3,5,4) = 32555.0/271803.0;
   acof(3,5,5) = -283.0/2107.0;
   acof(3,5,6) = 0;
   acof(3,5,7) = 0;
   acof(3,5,8) = 0;
   acof(3,6,1) = 0;
   acof(3,6,2) = 0;
   acof(3,6,3) = 0;
   acof(3,6,4) = 9.0/2107.0;
   acof(3,6,5) = -11.0/6321.0;
   acof(3,6,6) = 0;
   acof(3,6,7) = 0;
   acof(3,6,8) = 0;
   acof(3,7,1) = 0;
   acof(3,7,2) = 0;
   acof(3,7,3) = 0;
   acof(3,7,4) = 0;
   acof(3,7,5) = 0;
   acof(3,7,6) = 0;
   acof(3,7,7) = 0;
   acof(3,7,8) = 0;
   acof(3,8,1) = 0;
   acof(3,8,2) = 0;
   acof(3,8,3) = 0;
   acof(3,8,4) = 0;
   acof(3,8,5) = 0;
   acof(3,8,6) = 0;
   acof(3,8,7) = 0;
   acof(3,8,8) = 0;
   acof(4,1,1) = -36.0/833.0;
   acof(4,1,2) = 181507.0/3510262.0;
   acof(4,1,3) = 241309.0/10530786.0;
   acof(4,1,4) = 0;
   acof(4,1,5) = 0;
   acof(4,1,6) = 0;
   acof(4,1,7) = 0;
   acof(4,1,8) = 0;
   acof(4,2,1) = 177.0/3332.0;
   acof(4,2,2) = -544521.0/3510262.0;
   acof(4,2,3) = 987685.0/21061572.0;
   acof(4,2,4) = -14762.0/103243.0;
   acof(4,2,5) = 0;
   acof(4,2,6) = 0;
   acof(4,2,7) = 0;
   acof(4,2,8) = 0;
   acof(4,3,1) = -6.0/833.0;
   acof(4,3,2) = 544521.0/3510262.0;
   acof(4,3,3) = 2193521.0/3510262.0;
   acof(4,3,4) = 8065.0/14749.0;
   acof(4,3,5) = 381.0/2401.0;
   acof(4,3,6) = 0;
   acof(4,3,7) = 0;
   acof(4,3,8) = 0;
   acof(4,4,1) = -9.0/3332.0;
   acof(4,4,2) = -181507.0/3510262.0;
   acof(4,4,3) = -2647979.0/3008796.0;
   acof(4,4,4) = -80793.0/103243.0;
   acof(4,4,5) = -1927.0/2401.0;
   acof(4,4,6) = -2.0/49.0;
   acof(4,4,7) = 0;
   acof(4,4,8) = 0;
   acof(4,5,1) = 0;
   acof(4,5,2) = 0;
   acof(4,5,3) = 57418.0/309729.0;
   acof(4,5,4) = 51269.0/103243.0;
   acof(4,5,5) = 1143.0/2401.0;
   acof(4,5,6) = 8.0/49.0;
   acof(4,5,7) = 0;
   acof(4,5,8) = 0;
   acof(4,6,1) = 0;
   acof(4,6,2) = 0;
   acof(4,6,3) = 0;
   acof(4,6,4) = -283.0/2401.0;
   acof(4,6,5) = 403.0/2401.0;
   acof(4,6,6) = -6.0/49.0;
   acof(4,6,7) = 0;
   acof(4,6,8) = 0;
   acof(4,7,1) = 0;
   acof(4,7,2) = 0;
   acof(4,7,3) = 0;
   acof(4,7,4) = 0;
   acof(4,7,5) = 0;
   acof(4,7,6) = 0;
   acof(4,7,7) = 0;
   acof(4,7,8) = 0;
   acof(4,8,1) = 0;
   acof(4,8,2) = 0;
   acof(4,8,3) = 0;
   acof(4,8,4) = 0;
   acof(4,8,5) = 0;
   acof(4,8,6) = 0;
   acof(4,8,7) = 0;
   acof(4,8,8) = 0;
   acof(5,1,1) = 0;
   acof(5,1,2) = 0;
   acof(5,1,3) = 5.0/6192.0;
   acof(5,1,4) = -1.0/49.0;
   acof(5,1,5) = 0;
   acof(5,1,6) = 0;
   acof(5,1,7) = 0;
   acof(5,1,8) = 0;
   acof(5,2,1) = 0;
   acof(5,2,2) = 0;
   acof(5,2,3) = 815.0/151704.0;
   acof(5,2,4) = 1186.0/18963.0;
   acof(5,2,5) = 0;
   acof(5,2,6) = 0;
   acof(5,2,7) = 0;
   acof(5,2,8) = 0;
   acof(5,3,1) = 0;
   acof(5,3,2) = 0;
   acof(5,3,3) = -7381.0/50568.0;
   acof(5,3,4) = 32555.0/303408.0;
   acof(5,3,5) = -283.0/2352.0;
   acof(5,3,6) = 0;
   acof(5,3,7) = 0;
   acof(5,3,8) = 0;
   acof(5,4,1) = 0;
   acof(5,4,2) = 0;
   acof(5,4,3) = 28709.0/151704.0;
   acof(5,4,4) = 51269.0/101136.0;
   acof(5,4,5) = 381.0/784.0;
   acof(5,4,6) = 1.0/6.0;
   acof(5,4,7) = 0;
   acof(5,4,8) = 0;
   acof(5,5,1) = 0;
   acof(5,5,2) = 0;
   acof(5,5,3) = -349.0/7056.0;
   acof(5,5,4) = -247951.0/303408.0;
   acof(5,5,5) = -577.0/784.0;
   acof(5,5,6) = -5.0/6.0;
   acof(5,5,7) = -1.0/24.0;
   acof(5,5,8) = 0;
   acof(5,6,1) = 0;
   acof(5,6,2) = 0;
   acof(5,6,3) = 0;
   acof(5,6,4) = 1135.0/7056.0;
   acof(5,6,5) = 1165.0/2352.0;
   acof(5,6,6) = 1.0/2.0;
   acof(5,6,7) = 1.0/6.0;
   acof(5,6,8) = 0;
   acof(5,7,1) = 0;
   acof(5,7,2) = 0;
   acof(5,7,3) = 0;
   acof(5,7,4) = 0;
   acof(5,7,5) = -1.0/8.0;
   acof(5,7,6) = 1.0/6.0;
   acof(5,7,7) = -1.0/8.0;
   acof(5,7,8) = 0;
   acof(5,8,1) = 0;
   acof(5,8,2) = 0;
   acof(5,8,3) = 0;
   acof(5,8,4) = 0;
   acof(5,8,5) = 0;
   acof(5,8,6) = 0;
   acof(5,8,7) = 0;
   acof(5,8,8) = 0;
   acof(6,1,1) = 0;
   acof(6,1,2) = 0;
   acof(6,1,3) = 0;
   acof(6,1,4) = 1.0/392.0;
   acof(6,1,5) = 0;
   acof(6,1,6) = 0;
   acof(6,1,7) = 0;
   acof(6,1,8) = 0;
   acof(6,2,1) = 0;
   acof(6,2,2) = 0;
   acof(6,2,3) = 0;
   acof(6,2,4) = -1.0/144.0;
   acof(6,2,5) = 0;
   acof(6,2,6) = 0;
   acof(6,2,7) = 0;
   acof(6,2,8) = 0;
   acof(6,3,1) = 0;
   acof(6,3,2) = 0;
   acof(6,3,3) = 0;
   acof(6,3,4) = 3.0/784.0;
   acof(6,3,5) = -11.0/7056.0;
   acof(6,3,6) = 0;
   acof(6,3,7) = 0;
   acof(6,3,8) = 0;
   acof(6,4,1) = 0;
   acof(6,4,2) = 0;
   acof(6,4,3) = 0;
   acof(6,4,4) = -283.0/2352.0;
   acof(6,4,5) = 403.0/2352.0;
   acof(6,4,6) = -1.0/8.0;
   acof(6,4,7) = 0;
   acof(6,4,8) = 0;
   acof(6,5,1) = 0;
   acof(6,5,2) = 0;
   acof(6,5,3) = 0;
   acof(6,5,4) = 1135.0/7056.0;
   acof(6,5,5) = 1165.0/2352.0;
   acof(6,5,6) = 1.0/2.0;
   acof(6,5,7) = 1.0/6.0;
   acof(6,5,8) = 0;
   acof(6,6,1) = 0;
   acof(6,6,2) = 0;
   acof(6,6,3) = 0;
   acof(6,6,4) = -47.0/1176.0;
   acof(6,6,5) = -5869.0/7056.0;
   acof(6,6,6) = -3.0/4.0;
   acof(6,6,7) = -5.0/6.0;
   acof(6,6,8) = -1.0/24.0;
   acof(6,7,1) = 0;
   acof(6,7,2) = 0;
   acof(6,7,3) = 0;
   acof(6,7,4) = 0;
   acof(6,7,5) = 1.0/6.0;
   acof(6,7,6) = 1.0/2.0;
   acof(6,7,7) = 1.0/2.0;
   acof(6,7,8) = 1.0/6.0;
   acof(6,8,1) = 0;
   acof(6,8,2) = 0;
   acof(6,8,3) = 0;
   acof(6,8,4) = 0;
   acof(6,8,5) = 0;
   acof(6,8,6) = -1.0/8.0;
   acof(6,8,7) = 1.0/6.0;
   acof(6,8,8) = -1.0/8.0;

   bope(1,1) = -24.0/17.0;
   bope(1,2) = 59.0/34.0;
   bope(1,3) = -4.0/17.0;
   bope(1,4) = -3.0/34.0;
   bope(1,5) = 0;
   bope(1,6) = 0;
   bope(1,7) = 0;
   bope(1,8) = 0;
   bope(2,1) = -1.0/2.0;
   bope(2,2) = 0;
   bope(2,3) = 1.0/2.0;
   bope(2,4) = 0;
   bope(2,5) = 0;
   bope(2,6) = 0;
   bope(2,7) = 0;
   bope(2,8) = 0;
   bope(3,1) = 4.0/43.0;
   bope(3,2) = -59.0/86.0;
   bope(3,3) = 0;
   bope(3,4) = 59.0/86.0;
   bope(3,5) = -4.0/43.0;
   bope(3,6) = 0;
   bope(3,7) = 0;
   bope(3,8) = 0;
   bope(4,1) = 3.0/98.0;
   bope(4,2) = 0;
   bope(4,3) = -59.0/98.0;
   bope(4,4) = 0;
   bope(4,5) = 32.0/49.0;
   bope(4,6) = -4.0/49.0;
   bope(4,7) = 0;
   bope(4,8) = 0;
   bope(5,1) = 0;
   bope(5,2) = 0;
   float_sw4 d4a = 2.0/3;
   float_sw4 d4b = -1.0/12;
   bope(5,3) = -d4b;
   bope(5,4) = -d4a;
   bope(5,5) = 0;
   bope(5,6) =  d4a;
   bope(5,7) =  d4b;
   bope(5,8) = 0;
   bope(6,1) = 0;
   bope(6,2) = 0;
   bope(6,3) = 0;
   bope(6,4) = -d4b;
   bope(6,5) = -d4a;
   bope(6,6) = 0;
   bope(6,7) =  d4a;
   bope(6,8) =  d4b;
#undef acof
#undef ghcof
#undef bope
   _sbop[0] = -1.0/4;
   _sbop[1] = -5.0/6;
   _sbop[2] =  3.0/2;
   _sbop[3] = -1.0/2;
   _sbop[4] =  1.0/12;
}
