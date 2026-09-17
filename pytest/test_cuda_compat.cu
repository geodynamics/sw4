#include <iostream>

#include "CudaCompat.h"

namespace {

bool check(cudaError_t error, const char* operation) {
  if (error == cudaSuccess) return true;
  std::cerr << operation << " failed: " << cudaGetErrorString(error) << '\n';
  return false;
}

__global__ void set_value(int* value) { *value = 42; }

}  // namespace

int main() {
  int device = 0;
  if (!check(cudaGetDevice(&device), "cudaGetDevice")) return 1;

  int* value = nullptr;
  if (!check(cudaMallocManaged(&value, sizeof(*value)), "cudaMallocManaged"))
    return 1;

  bool passed =
      check(sw4::cuda::mem_advise(value, sizeof(*value),
                                  cudaMemAdviseSetPreferredLocation, device),
            "cudaMemAdvise") &&
      check(sw4::cuda::mem_prefetch_async(value, sizeof(*value), device),
            "cudaMemPrefetchAsync(device)");

  if (passed) {
    set_value<<<1, 1>>>(value);
    passed = check(cudaGetLastError(), "set_value launch") &&
             check(cudaDeviceSynchronize(), "set_value synchronization") &&
             check(sw4::cuda::mem_prefetch_async(value, sizeof(*value),
                                                 cudaCpuDeviceId),
                   "cudaMemPrefetchAsync(host)") &&
             check(cudaDeviceSynchronize(), "host prefetch synchronization");
  }

  if (passed && *value != 42) {
    std::cerr << "managed-memory value mismatch: expected 42, got " << *value
              << '\n';
    passed = false;
  }

  if (!check(cudaFree(value), "cudaFree")) passed = false;
  if (passed) std::cout << "CUDA compatibility test passed\n";
  return passed ? 0 : 1;
}
