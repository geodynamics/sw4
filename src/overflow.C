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

void set(int &in){
  in = -2;
}

int main(int argc, char **argv){
  int i = 123456;
  int j = 345678;

  //Dangers of mixing size_t and int
  int ii;
  set(ii);
  ii=std::stoi(argv[1]);
  //std::cout<<" II is "<<ii<<" "<<(size_t)ii<<"\n";
  size_t off = 2000;
  ssize_t off2 = off;

  int one = off/ii;
  int two = off2/ii;
  size_t iii=ii;
  int three = off/iii;
  std::cout<<"MULTIPLY = "<<std::is_signed_v<decltype(off*ii)><<" "<<
    std::is_signed_v<decltype(off2*ii)><<" "<<
    std::is_signed_v<decltype(off*iii)><<"\n";
  std::cout<<"VALUES = "<<off*ii<<" "<<off2*ii<<" "<<off*iii<<"\n\n\n";
  std::cout<<" DIVISION =  "<<std::is_signed_v<decltype(off/ii)><<" "<<
    std::is_signed_v<decltype(off2/ii)><<" "<<
    std::is_signed_v<decltype(off/iii)><<"\n";
  std::cout<<" VALUES = "<<one<<" "<<two<<" "<<three<<"\n";
  int k=0;
  //int k = i*j;
  std::cout<<"K = "<<k<<"\n";

  runOnGPU([=]__device__(){
	     ssize_t k = (long long)i*j;
	     int sk=k;
	   });
      
}
  
