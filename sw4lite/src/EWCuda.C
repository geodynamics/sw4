#include "EWCuda.h"
#include <iostream>
using namespace std;

//-----------------------------------------------------------------------
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

//-----------------------------------------------------------------------
void EWCuda::reset_gpu()
{
#ifdef SW4_CUDA
   if( m_ndevice > 0 )
      cudaDeviceReset();
#endif
}

//-----------------------------------------------------------------------
void EWCuda::sync_stream( int st )
{
#ifdef SW4_CUDA
   cudaError_t retcode;
   retcode = cudaStreamSynchronize(m_stream[st]);
   if( retcode != cudaSuccess )
      cout << "Error EWCuda::EWCuda, cudaStreamSynchronize no " << st << " returned " <<
	 cudaGetErrorString(retcode) << endl;
#endif
}

//-----------------------------------------------------------------------
void EWCuda::initialize_gpu(int myrank)
{
#ifdef SW4_CUDA
   if( m_ndevice > 0){

     cudaError_t retcode  = cudaSetDevice(0);
     if (retcode != cudaSuccess)
	cout << "Error cudaSetDevice: "  << cudaGetErrorString(retcode) << endl;
     else
        cudaDeviceReset();
        int myDevice = sched_getcpu()/(40);
        int cpu = sched_getcpu();
        cout << "myrank = " << myrank <<  "  cpuid = " << cpu << "  mydevice = " << myDevice  << endl;

        retcode  = cudaSetDevice(myDevice);
        if (retcode != cudaSuccess)
      	   cout << "Error cudaSetDevice: "  << cudaGetErrorString(retcode) << endl;
        else
        {
	   for (int i = 0; i < m_nstream; i++)
	   {
	      retcode = cudaStreamCreate(&m_stream[i]);
	      if( retcode != cudaSuccess )
	         cout << "Error EWCuda::EWCuda, cudaStreamCreate no " << i << " returned " <<
		          cudaGetErrorString(retcode) << endl;
	   }
        }
   }
#endif
}

//-----------------------------------------------------------------------
EWCuda::~EWCuda()
{
#ifdef SW4_CUDA
   for( int s=0 ; s < m_nstream ; s++ )
      cudaStreamDestroy(m_stream[s]);
   if( m_nstream > 0 )
      delete[] m_stream;
#endif
}
