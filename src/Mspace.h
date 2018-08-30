#include "sw4.h"
#ifndef __MSPACE_H__
#define __MSPACE_H__

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include <tuple>
#include <mpi.h>
#ifdef SW4_USE_UMPIRE
#include "umpire/ResourceManager.hpp"
//#include "umpire/Umpire.hpp"
//#include "umpire/Allocator.hpp"
#include "umpire/strategy/DynamicPool.hpp"
#endif
#if defined(ENABLE_CUDA)
#include "cuda_runtime.h"
#include <nvml.h>
void CheckError(cudaError_t const err, const char* file, char const* const fun, const int line);
void prefetch_to_device(const float_sw4 *ptr);
#define SW4_CheckDeviceError(err) CheckError(err,__FILE__, __FUNCTION__, __LINE__)
#define PROFILER_START SW4_CheckDeviceError(cudaProfilerStart())
#define PROFILER_STOP SW4_CheckDeviceError(cudaProfilerStop())
#else
#define SW4_CheckDeviceError(err) 
#endif
void check_mem();
enum Space { Host, Managed,Device,Pinned,Managed_temps};
void * operator new(std::size_t size,Space loc) throw(std::bad_alloc) ;
void operator delete(void *ptr, Space loc) throw();
void * operator new[](std::size_t size,Space loc) throw(std::bad_alloc) ;
void * operator new[](std::size_t size,Space loc,const char *file,int line);
void operator delete[](void *ptr, Space loc) throw();
void presetGPUID();
void print_hwm();
struct global_variable_holder_struct {
  size_t gpu_memory_hwm ;
  size_t curr_mem;
  size_t max_mem;
};

extern struct global_variable_holder_struct global_variables;
#ifdef ENABLE_CUDA
#define SW4_TRACK_MEMORY_ALLOCATIONS 1
#endif
#if defined(SW4_TRACK_MEMORY_ALLOCATIONS)
#define SW4_NEW(type, arg )				\
  ( new(type,__FILE__,__LINE__) arg)

#if defined(ENABLE_MEMORY_ASSERTS)
#define ASSERT_MANAGED(ptr)\
  ( assert_check_managed((ptr),__FILE__,__LINE__))

#define ASSERT_HOST(ptr)\
  ( assert_check_host((ptr),__FILE__,__LINE__))

#else

#define ASSERT_MANAGED(ptr) 

#define ASSERT_HOST(ptr) 
#endif

#define PREFETCH(ptr)\
  ( prefetch_to_device((ptr)))

void assert_check_host(void *ptr, const char *file, int line);
void assert_check_managed(void *ptr, const char *file, int line);

#define PTR_PUSH(type,ptr,size)			\
  ( ptr_push(ptr,type,size,__FILE__,__LINE__))

void ptr_push(void *ptr, Space type, size_t size,const char *file, int line);

#else

#define SW4_NEW(type, arg)				\
  ( new(type) arg)

#define ASSERT_MANAGED(ptr) 

#define ASSERT_HOST(ptr) 

#define PTR_PUSH(type,ptr,size) 

#define PREFETCH(ptr)

ssize_t getsize(const void *ptr);
#endif

// THIS WILL HAVE TO BE MODIFIED FOR NON_GPU MACHINES
class Managed {
public:
  //static size_t mem_total;
  static int host;
  static int device;
  Managed(){
  }
  ~Managed(){
    //mem_total=0;
  }
  void *operator new(size_t len) {
    void *ptr;
    //mem_total+=len;
    //std::cout<<"Total mem is now "<<mem_total/1024/1024<<" MB \n";
    // std::cout<<"Call to Managed class "<<len<<"\n";
    SW4_CheckDeviceError(cudaMallocManaged(&ptr, len));
    //SW4_CheckDeviceError(cudaDeviceSynchronize());
    return ptr;
  }

  void *operator new[](size_t len) {
    void *ptr;
    //mem_total+=len;
    //std::cout<<"Total mem is now "<<mem_total/1204/1024<<" MB \n";
    //std::cout<<"Call to Managed class "<<len<<"\n";
    SW4_CheckDeviceError(cudaMallocManaged(&ptr, len));
    //SW4_CheckDeviceError(cudaDeviceSynchronize());
    return ptr;
  }
  

  void operator delete(void *ptr) {
    //SW4_CheckDeviceError(cudaDeviceSynchronize());
    SW4_CheckDeviceError(cudaFree(ptr));
  }

  void operator delete[](void *ptr) {
    //SW4_CheckDeviceError(cudaDeviceSynchronize());
    SW4_CheckDeviceError(cudaFree(ptr));
  }
};
//size_t Managed::mem_total=0;

// AUTOPEEL CODE
class Apc{
public:
  std::ofstream ofile;
  int counter;
  Apc(char *s);
  ~Apc();
};

template <typename T,typename N>
  std::string line(T,N) = delete;

std::string line(int n,int C);
std::string line(int *n,int C);
std::string line(double n, int C);
std::string line(double *n,int C);
std::string line(char n,int C);

template<typename T>
void autopeel(Apc &apc,T only){
  apc.ofile<<line(only,apc.counter);
  return;
}

template<typename T,typename... Args>
void autopeel(Apc &apc, T first, Args&&... args){
  apc.ofile<<line(first,apc.counter);
  apc.counter++;
  autopeel(apc,args...);
}

// END AUTOPEEL CODE


#endif
