#include "EWCuda.h"
#include <iostream>
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
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

int mystrcmp(void const *a, void const *b) { 
  char const *aa = (char const *)a;
  char const *bb = (char const *)b;

  return strcmp(aa, bb);
}

extern "C"
void setupgpu(int verbose)
{
#ifdef SW4_CUDA
  int rank, nrank;
  MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  char (*hosts)[MPI_MAX_PROCESSOR_NAME] = (char (*)[MPI_MAX_PROCESSOR_NAME])malloc(nrank*(sizeof *hosts));
  char hostname[MPI_MAX_PROCESSOR_NAME];
  int namelen;
  
  // each rank collects hostname of all nodes
  MPI_Get_processor_name(hostname, &namelen);
  strcpy(hosts[rank], hostname);
  for (int i = 0; i < nrank; i++) {
    MPI_Bcast(hosts[i], MPI_MAX_PROCESSOR_NAME, MPI_CHAR, i, MPI_COMM_WORLD);
  }

  // sort list of names
  qsort(hosts, nrank, MPI_MAX_PROCESSOR_NAME, mystrcmp);

  // assign same color to the same node
  int color = 0;
  for (int i = 0; i < nrank; i++) {
    if (i > 0) {
      if (strcmp(hosts[i-1], hosts[i]) != 0) color++;
    }
    if (strcmp(hosts[i], hostname) == 0) break;
  }
  
  MPI_Comm new_comm;
  MPI_Comm_split(MPI_COMM_WORLD, color, 0, &new_comm);
  int new_rank;
  MPI_Comm_rank(new_comm, &new_rank);

  int ngpu;
  cudaError_t ierr;
  ierr = cudaGetDeviceCount(&ngpu);
  if (ierr != cudaSuccess) {
    printf("%s\n", cudaGetErrorString(ierr));
    exit(0);
  }

  int igpu = new_rank % ngpu;
  if (verbose) 
    printf("P(%d): %s using gpu id %d\n", rank, hosts[rank], igpu);
  ierr = cudaSetDevice(igpu);
  if (ierr != cudaSuccess) {
    printf("%s\n", cudaGetErrorString(ierr));
    exit(0);
  }

  free(hosts);
#endif
}

//-----------------------------------------------------------------------
void EWCuda::initialize_gpu(int myrank)
{
#ifdef SW4_CUDA
   if( m_ndevice > 0){
     setupgpu(0);
     cudaError_t retcode;
     for (int i = 0; i < m_nstream; i++)
     {
       retcode = cudaStreamCreate(&m_stream[i]);
       if( retcode != cudaSuccess )
         cout << "Error EWCuda::EWCuda, cudaStreamCreate no " << i << " returned " <<
           cudaGetErrorString(retcode) << endl;
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
