#ifndef SW4_CUDA_COMPAT_H
#define SW4_CUDA_COMPAT_H

#include <cuda_runtime.h>

#include <cstddef>

namespace sw4 {
namespace cuda {

inline cudaError_t mem_prefetch_async(const void* ptr, std::size_t count,
                                      int device,
                                      cudaStream_t stream = nullptr) {
#if CUDART_VERSION >= 13000
  cudaMemLocation location{};
  if (device == cudaCpuDeviceId) {
    location.type = cudaMemLocationTypeHost;
  } else {
    location.type = cudaMemLocationTypeDevice;
    location.id = device;
  }
  return cudaMemPrefetchAsync(ptr, count, location, 0, stream);
#else
  return cudaMemPrefetchAsync(ptr, count, device, stream);
#endif
}

inline cudaError_t mem_advise(const void* ptr, std::size_t count,
                              cudaMemoryAdvise advice, int device) {
#if CUDART_VERSION >= 13000
  cudaMemLocation location{};
  if (device == cudaCpuDeviceId) {
    location.type = cudaMemLocationTypeHost;
  } else {
    location.type = cudaMemLocationTypeDevice;
    location.id = device;
  }
  return cudaMemAdvise(ptr, count, advice, location);
#else
  return cudaMemAdvise(ptr, count, advice, device);
#endif
}

}  // namespace cuda
}  // namespace sw4

#endif  // SW4_CUDA_COMPAT_H
