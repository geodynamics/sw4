//hipcc -x hip --std=c++17 -O3 -ftrapv -g --offload-arch=gfx90a trap.C
// Code generates int overflows and host and device based on 1st argument
#include <iostream>

#include "hip/hip_runtime.h"
#include "hip/hip_runtime_api.h"


template <typename Func>
__global__ void forallkernel(Func f) {
  f();
}

template <typename LoopBody>
void RunOnGPU(LoopBody &&body) {
  int tpb = 1;
  int blocks = 1;
  hipLaunchKernelGGL(forallkernel<LoopBody>, blocks, tpb, 0, 0,
                     body);
}

int main(int argc, char **argv){

  int a = 1234578;
  int b = 5678901;

  int c=0;
  if ((argc==2 ) &&(strcmp(argv[1],"H") ==0)){
    std::cout<<"Host overflow\n";
    c = a*b;
    std::cout<<"c is "<<c<<"\n";
  } else {
    std::cout<<"Device overflow \n";
    RunOnGPU([=]__device__(){
	       int64_t c=a*b;
	     });
    hipStreamSynchronize(0);
  }
}
