#include "EWCuda.h"

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
#endif
}

void EWCuda::reset_gpu()
{
#ifdef SW4_CUDA
   if( m_ndevice > 0 )
      cudaDeviceReset();
#endif
}
