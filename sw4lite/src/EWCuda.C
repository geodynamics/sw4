#include "EWCuda.h"
#include <iostream>
using namespace std;

EWCuda::EWCuda( int ndevice, int nstream )
{
   m_ndevice = ndevice;
   m_nstream = nstream;
   m_active_gpu = 0;
#ifdef SW4_CUDA
   if( nstream > 0 )
      m_stream  = new cudaStream_t[nstream];
   else
      m_stream = static_cast<cudaStream_t*>(0);
   for( int s=0 ; s < nstream ; s++ )
   {
      cudaError_t retcode;
      retcode = cudaStreamCreate(&m_stream[s]);
      if( retcode != cudaSuccess )
	 cout << "Error EWCuda::EWCuda, cudaStreamCreate no " << s << " returned " <<
	    cudaGetErrorString(retcode) << endl;
   }
#endif
}

void EWCuda::reset_gpu()
{
#ifdef SW4_CUDA
   if( m_ndevice > 0 )
      cudaDeviceReset();
#endif
}

void EWCuda::sync_stream( int st )
{
   cudaError_t retcode;
   retcode = cudaStreamSynchronize(m_stream[st]);
   if( retcode != cudaSuccess )
      cout << "Error EWCuda::EWCuda, cudaStreamSynchronize no " << st << " returned " <<
	 cudaGetErrorString(retcode) << endl;
}

EWCuda::~EWCuda()
{
   for( int s=0 ; s < m_nstream ; s++ )
      cudaStreamDestroy(m_stream[s]);
   if( m_nstream > 0 )
      delete[] m_stream;
}
