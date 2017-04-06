#ifndef SW4_EWCUDA
#define SW4_EWCUDA

#ifdef SW4_CUDA
#include <cuda_runtime.h>
#endif

class EWCuda 
{
 public:
   int m_nstream, m_ndevice, m_active_gpu;
#ifdef SW4_CUDA
   cudaStream_t* m_stream;
#endif
   EWCuda( int ndevice, int nstream ); 
   ~EWCuda();
   bool has_gpu(){return m_ndevice>0;}
   void reset_gpu();
   void initialize_gpu( int myrank );
   void sync_stream( int st );
};

#define RADIUS 2
#define DIAMETER (2*RADIUS+1)
#define RHS4_BLOCKX 16
#define RHS4_BLOCKY 16
#define ADDSGD4_BLOCKX 24
#define ADDSGD4_BLOCKY 16

#endif
