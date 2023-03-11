// Compile using :
// hipcc -ftrapv -g -O3 --offload-arch=gfx90a:xnack-:sramecc+ overflow.C
#include <iostream>
#include "hip/hip_runtime.h"
#include "hip/hip_runtime_api.h"

template <typename Func>
__global__ void runOnGPUKernel(Func f) {
  f();
}

template <typename LoopBody>
void runOnGPU(LoopBody &&body) {
  int tpb = 1;
  int blocks=1;
  // printf("Launching the kernel blocks= %d tpb= %d \n",blocks,tpb);
  // forallkernel<<<blocks, tpb>>>(start, end, body);
  hipLaunchKernelGGL(runOnGPUKernel<LoopBody>, blocks, tpb, 0, 0,
                     body);
  hipStreamSynchronize(0);
}

int main(int argc, char **argv){
  int i = 123456;
  int j = 345678;

  int k=0;
  //int k = i*j;
  std::cout<<"K = "<<k<<"\n";

  runOnGPU([=]__device__(){
	     int k = i*j;
	   });
      
}
  
