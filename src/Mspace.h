#ifndef __MSPACE_H__
#define __MSPACE_H__

#include <cstdlib>
#include <iostream>

#if defined(ENABLE_CUDA)
#include "cuda_runtime.h"
void CheckError(cudaError_t const err, const char* file, char const* const fun, const int line);
#define SW4_CheckDeviceError(err) CheckError(err,__FILE__, __FUNCTION__, __LINE__)
#else
#define SW4_CheckDeviceError(err) ()
#endif

enum Space { Host, Managed};
void * operator new(std::size_t size,Space loc) throw(std::bad_alloc) ;
void operator delete(void *ptr, Space loc) throw();
void * operator new[](std::size_t size,Space loc) throw(std::bad_alloc) ;
void * operator new[](std::size_t size,Space loc,const char *file,int line);
void operator delete[](void *ptr, Space loc) throw();




#define SW4_TRACK_MEMORY_ALLOCATIONS 1
#if defined(SW4_TRACK_MEMORY_ALLOCATIONS)
#define SW4_NEW(type, arg )				\
  ( new(type,__FILE__,__LINE__) arg)

#define ASSERT_MANAGED(ptr)\
  ( assert_check_managed((ptr),__FILE__,__LINE__))

#define ASSERT_HOST(ptr)\
  ( assert_check_host((ptr),__FILE__,__LINE__))

void assert_check_host(void *ptr, const char *file, int line);
void assert_check_managed(void *ptr, const char *file, int line);

#define PTR_PUSH(type,ptr)				\
  ( ptr_push(ptr,type,__FILE__,__LINE__))

void ptr_push(void *ptr, Space type, const char *file, int line);

#else

#define SW4_NEW(type, arg)				\
  ( new(type) arg)

#define ASSERT_MANAGED(ptr) ()

#define ASSERT_HOST(ptr) ()

#define PTR_PUSH(ptr) ()


#endif


#endif
