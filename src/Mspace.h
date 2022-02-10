#include "sw4.h"
#ifndef __MSPACE_H__
#define __MSPACE_H__

#include <mpi.h>

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "policies.h"
#ifdef SW4_USE_UMPIRE
#include "umpire/ResourceManager.hpp"
//#include "umpire/Umpire.hpp"
//#include "umpire/Allocator.hpp"
//#include "umpire/strategy/DynamicPool.hpp"
#include "umpire/strategy/QuickPool.hpp"
//#include "umpire/strategy/MixedPool.hpp"
//#include "umpire/util/StatisticsDatabase.hpp"
#include "umpire/strategy/AlignedAllocator.hpp"
#include "umpire/strategy/AllocationAdvisor.hpp"
#include "umpire/strategy/MonotonicAllocationStrategy.hpp"
#include "umpire/util/Macros.hpp"
#endif
#ifdef USE_MAGMA
#include "magma_v2.h"
#endif
#include <cstddef>

#if defined(ENABLE_CUDA)
#include <cuda_profiler_api.h>
#include <nvml.h>

#include "cuda_runtime.h"
bool mpi_supports_device_buffers();
void CheckError(cudaError_t const err, const char *file, char const *const fun,
                const int line);
void prefetch_to_device(const float_sw4 *ptr);
#define SW4_CheckDeviceError(err) \
  CheckError(err, __FILE__, __FUNCTION__, __LINE__)
#define PROFILER_START SW4_CheckDeviceError(cudaProfilerStart())
#define PROFILER_STOP SW4_CheckDeviceError(cudaProfilerStop())

#elif ENABLE_HIP
#ifndef SW4_NO_ROCTRACER
#include <roctracer/roctracer.h>
#include <roctracer/roctracer_ext.h>
#endif
void CheckError(hipError_t const err, const char *file, char const *const fun,
                const int line);
#define SW4_CheckDeviceError(err) \
  CheckError(err, __FILE__, __FUNCTION__, __LINE__)
//#define PROFILER_START SW4_CheckDeviceError(hipProfilerStart())
#ifndef SW4_NO_ROCTRACER
#define PROFILER_STOP SW4_CheckDeviceError(hipProfilerStop())
#define PROFILER_START roctracer_start()
#else
#define PROFILER_START
#define PROFILER_STOP
#endif
#else
// void CheckError(hipError_t const err, const char *file, char const *const
// fun,
//                const int line);
#define SW4_CheckDeviceError(err)
#define PROFILER_START
#define PROFILER_STOP
#endif

void check_mem();
void global_prefetch();
enum class Space : unsigned int {
  Host,
  Managed,
  Device,
  Pinned,
  Managed_temps,
  Space_Error
};
template <typename Enumeration>
auto as_int(Enumeration const value) ->
    typename std::underlying_type<Enumeration>::type {
  return static_cast<typename std::underlying_type<Enumeration>::type>(value);
}
Space GML(const void *ptr);
void *operator new(std::size_t size, Space loc) throw();
void operator delete(void *ptr, Space loc) throw();
void *operator new[](std::size_t size, Space loc) throw();
void *operator new[](std::size_t size, Space loc, const char *file, int line);
void operator delete[](void *ptr, Space loc) throw();
void operator delete(void *ptr, Space loc, const char *file, int line) throw();
void operator delete[](void *ptr, Space loc, const char *file,
                       int line) throw();
int presetGPUID(int mpi_rank, int local_rank, int local_size);
void print_hwm(int rank);
struct global_variable_holder_struct {
  size_t gpu_memory_hwm;
  size_t curr_mem;
  size_t max_mem;
  size_t host_mem_hwm;
  size_t host_curr_mem;
  size_t host_max_mem;
  int device;
  int num_devices;
  bool firstCycle;
  int current_step;
  size_t buffer_size;
  float_sw4 *device_buffer;
  std::vector<std::tuple<char *, size_t>> massprefetch;
};

extern struct global_variable_holder_struct global_variables;
#ifdef ENABLE_CUDA
//#define SW4_TRACK_MEMORY_ALLOCATIONS 1
#endif
#if defined(SW4_TRACK_MEMORY_ALLOCATIONS)
#define SW4_NEW(type, arg) (new (type, __FILE__, __LINE__) arg)

#if defined(ENABLE_MEMORY_ASSERTS)
#define ASSERT_MANAGED(ptr) (assert_check_managed((ptr), __FILE__, __LINE__))

#define ASSERT_HOST(ptr) (assert_check_host((ptr), __FILE__, __LINE__))

#else

#define ASSERT_MANAGED(ptr)

#define ASSERT_HOST(ptr)
#endif

#define PREFETCH(ptr) (prefetch_to_device((ptr)))

void assert_check_host(void *ptr, const char *file, int line);
void assert_check_managed(void *ptr, const char *file, int line);

#define PTR_PUSH(type, ptr, size) \
  (ptr_push(ptr, type, size, __FILE__, __LINE__))

void ptr_push(void *ptr, Space type, size_t size, const char *file, int line);

#else

#define SW4_NEW(type, arg) (new (type) arg)

#define ASSERT_MANAGED(ptr)

#define ASSERT_HOST(ptr)

#define PTR_PUSH(type, ptr, size)

#define PREFETCH(ptr)

ssize_t getsize(const void *ptr);
#endif

// THIS WILL HAVE TO BE MODIFIED FOR NON_GPU MACHINES
// Used in GridPointSource
class Managed {
 public:
  static size_t ocount;
  static size_t hwm;
  // static size_t mem_total;
  // static int host;
  // static int device;

  Managed() {}
  ~Managed() {
    // mem_total=0;
  }
#if defined(ENABLE_GPU)
  void *operator new(size_t len);

  void *operator new[](size_t len);

  void operator delete(void *ptr);
  void operator delete[](void *ptr);
#endif
};

#define SW4_CHRONO_NOW std::chrono::high_resolution_clock::now()

#define SW4_CHRONO_DURATION_MS(arg1, arg2) \
  std::chrono::duration_cast<std::chrono::milliseconds>((arg2) - (arg1)).count()
#define SW4_CHRONO_DURATION_US(arg1, arg2) \
  std::chrono::duration_cast<std::chrono::microseconds>((arg2) - (arg1)).count()

// AUTOPEEL CODE
class Apc {
 public:
  std::ofstream ofile;
  int counter;
  Apc(char *s);
  ~Apc();
};

template <typename T, typename N>
std::string line(T, N) = delete;

std::string line(int n, int C);
std::string line(int *n, int C);
std::string line(double n, int C);
std::string line(double *n, int C);
std::string line(char n, int C);

template <typename T>
void autopeel(Apc &apc, T only) {
  apc.ofile << line(only, apc.counter);
  return;
}

template <typename T, typename... Args>
void autopeel(Apc &apc, T first, Args &&... args) {
  apc.ofile << line(first, apc.counter);
  apc.counter++;
  autopeel(apc, args...);
}

// END AUTOPEEL CODE
template <class T>
void Write(T &t, std::string filename) {
  for (auto &i : t) {
    i.swrite(filename);
  }
}
void invert(float_sw4 *A, int N);

template <class T>
void spacecopy(T *&dst, T *&src, Space dst_space, Space src_space,
               size_t size) {
  if (src_space == Space::Host) {
#ifdef ENABLE_CUDA
    SW4_CheckDeviceError(
        cudaMemcpy(dst, src, size * sizeof(T), cudaMemcpyHostToDevice));
#elif ENABLE_HIP
    SW4_CheckDeviceError(
        hipMemcpy(dst, src, size * sizeof(T), hipMemcpyHostToDevice));
#else
    std::cout << "ERROR is Mspace:;copy:: undefined copy operation\n";
#endif
    return;
  }

  if (dst_space == Space::Host) {
#ifdef ENABLE_CUDA
    SW4_CheckDeviceError(
        cudaMemcpy(dst, src, size * sizeof(T), cudaMemcpyDeviceToHost));
#elif ENABLE_HIP
    SW4_CheckDeviceError(
        hipMemcpy(dst, src, size * sizeof(T), hipMemcpyDeviceToHost));
#else
    std::cout << "ERROR is Mspace:;copy:: undefined copy operation\n";
#endif
    return;
  }

  std::cout << "ERROR :: Undefined Mspace::copy operation from "
            << as_int(src_space) << " to " << as_int(dst_space) << "\n";
}
#endif
