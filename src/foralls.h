#ifdef ENABLE_CUDA
template<typename Func>
__global__ void forallkernel(int start,int N,Func f){
  int tid=start+threadIdx.x+blockIdx.x*blockDim.x;
  if  (tid<N) f(tid);
}
template<typename LoopBody>
void forall(int start, int end, LoopBody &&body){
  int tpb=32;
  int blocks=(end-start)/tpb;
  blocks=((end-start)%tpb==0)?blocks:blocks+1;
  //printf("Launching the kernel blocks= %d tpb= %d \n",blocks,tpb);
  forallkernel<<<blocks,tpb>>>(start,end,body);
  cudaDeviceSynchronize();
}
template<typename LoopBody>
void forallasync(int start, int end, LoopBody &&body){
  int tpb=1024;
  int blocks=(end-start)/tpb;
  blocks=((end-start)%tpb==0)?blocks:blocks+1;
  //printf("Launching the kernel blocks= %d tpb= %d \n",blocks,tpb);
  forallkernel<<<blocks,tpb>>>(start,end,body);
  //cudaDeviceSynchronize();
}

template<int N, typename LoopBody>
void forall(int start, int end, LoopBody &&body){
  int tpb=1024;
  int blocks=(end-start)/tpb;
  blocks=((end-start)%tpb==0)?blocks:blocks+1;
  printf("Launching the kernel blocks= %d tpb= %d on line %d\n",blocks,tpb,N);
  forallkernel<<<blocks,tpb>>>(start,end,body);
  cudaDeviceSynchronize();
}


template<typename Func>
__global__ void forallkernelB(int start,int N,Func f){
  int tid=start+threadIdx.x+blockIdx.x*blockDim.x;
  const int B=8;
  for(int i=tid;i<N;i+=B*blockDim.x * gridDim.x){
    f(i);
    //int ii=i+blockDim.x * gridDim.x;
#pragma unroll(B-1)
    for (int ii=1;ii<B;ii++){
      int iii=i+ii*blockDim.x*gridDim.x;
      if (iii<N) f(iii);
    }
  }
}
template<typename LoopBody>
void forallB(int start, int end, LoopBody &&body){
  int tpb=1024;
  int blocks=52;
  //blocks=((end-start)%tpb==0)?blocks:blocks+1;
  printf("Launching the kernel\n");
  forallkernelB<<<blocks,tpb>>>(start,end,body);
  cudaDeviceSynchronize();
}




template<int N>
class Range{
public:
  Range(int istart,int iend) : start(istart),end(iend),tpb(N)
{ blocks = (end-start)/N; blocks=((end-start)%N==0)?blocks:blocks+1;};
  int start;
  int end;
  int blocks;
  int tpb;
};

template<int N,int M>
class RangeGS{
public:
 RangeGS(int istart,int iend) : start(istart),end(iend),tpb(N),blocks(M)
{};
  int start;
  int end;
  int blocks;
  int tpb;
};


template<typename Func>
__global__ void forall3kernel(const int start0,const int N0,const int start1,const int N1, const int start2, const int N2, Func f){
  int tid0=start0+threadIdx.x+blockIdx.x*blockDim.x;
  int tid1=start1+threadIdx.y+blockIdx.y*blockDim.y;
  int tid2=start2+threadIdx.z+blockIdx.z*blockDim.z;
  if  ((tid0<N0) &&(tid1<N1) &&(tid2<N2)) f(tid0,tid1,tid2);
}

template<typename Func>
__global__ void forall3gskernel(int start0,int N0,int start1,int N1, int start2, int N2, Func f){
  for(int i=start0+threadIdx.x+blockIdx.x*blockDim.x;i<N0;i+=blockDim.x * gridDim.x)
    for(int j=start1+threadIdx.y+blockIdx.y*blockDim.y;j<N1;j+=blockDim.y * gridDim.y)
      for(int k=start2+threadIdx.z+blockIdx.z*blockDim.z;k<N2;k+=blockDim.z * gridDim.z)
	f(i,j,k);
}

template<typename LoopBody>
void forall3(int start0, int end0, int start1, int end1, int start2, int end2, LoopBody &&body){
  int tpb0 = 16;
  int tpb1 = 16 ;
  int tpb2 = 1024/(tpb0*tpb1);
  
  int block0=(end0-start0)/tpb0;
  block0=((end0-start0)%tpb0==0)?block0:block0+1;
  int block1=(end1-start1)/tpb1;
  block1=((end1-start1)%tpb1==0)?block1:block1+1;
  int block2=(end2-start2)/tpb2;
  block2=((end2-start2)%tpb2==0)?block2:block2+1;

  std::cout<<" BLOCKS "<<block0<<" "<<block1<<" "<<block2<<"\n";
  dim3 tpb(tpb0,tpb1,tpb2);
  dim3 blocks(block0,block1,block2);
  
  printf("Launching the kernel 3d \n");
  forall3kernel<<<blocks,tpb>>>(start0,end0,start1,end1,start2,end2,body);
  cudaDeviceSynchronize();
}

template<typename T1, typename T2, typename T3, typename LoopBody>
void forall3async(T1 &irange, T2 &jrange, T3 &krange, LoopBody &&body){
  
  dim3 tpb(irange.tpb,jrange.tpb,krange.tpb);
  dim3 blocks(irange.blocks,jrange.blocks,krange.blocks);
  
  forall3kernel<<<blocks,tpb>>>(irange.start,irange.end,jrange.start,jrange.end,krange.start,krange.end,body);
}
template<typename T1, typename T2, typename T3, typename LoopBody>
void forall3(T1 &irange, T2 &jrange, T3 &krange, LoopBody &&body){
  
  // dim3 tpb(irange.tpb,jrange.tpb,krange.tpb);
  // dim3 blocks(irange.blocks,jrange.blocks,krange.blocks);
  
  // forall3kernel<<<blocks,tpb>>>(irange.start,irange.end,jrange.start,jrange.end,krange.start,krange.end,body);
  forall3async(irange,jrange,krange,body);
  cudaStreamSynchronize(0);
}

template<typename T1, typename T2, typename T3, typename LoopBody>
void forall3GSasync(T1 &irange, T2 &jrange, T3 &krange, LoopBody &&body){
  
  dim3 tpb(irange.tpb,jrange.tpb,krange.tpb);
  dim3 blocks(irange.blocks,jrange.blocks,krange.blocks);
  
  forall3gskernel<<<blocks,tpb>>>(irange.start,irange.end,jrange.start,jrange.end,krange.start,krange.end,body);
}

template<typename T1, typename T2, typename T3, typename LoopBody>
void forall3GS(T1 &irange, T2 &jrange, T3 &krange, LoopBody &&body){
  
  // dim3 tpb(irange.tpb,jrange.tpb,krange.tpb);
  // dim3 blocks(irange.blocks,jrange.blocks,krange.blocks);
  
  // forall3gskernel<<<blocks,tpb>>>(irange.start,irange.end,jrange.start,jrange.end,krange.start,krange.end,body);
  forall3async(irange,jrange,krange,body);
  cudaStreamSynchronize(0);
}

// Forall2 

template<typename Func>
__global__ void forall2kernel(const int start0,const int N0,const int start1,const int N1,Func f){
  int tid0=start0+threadIdx.x+blockIdx.x*blockDim.x;
  int tid1=start1+threadIdx.y+blockIdx.y*blockDim.y;
  if  ((tid0<N0) &&(tid1<N1)) f(tid0,tid1);
}

template<typename T1, typename T2, typename LoopBody>
void forall2async(T1 &irange, T2 &jrange, LoopBody &&body){
  
  dim3 tpb(irange.tpb,jrange.tpb,1);
  dim3 blocks(irange.blocks,jrange.blocks,1);
  
  forall2kernel<<<blocks,tpb>>>(irange.start,irange.end,jrange.start,jrange.end,body);
}

template<typename T1, typename T2, typename LoopBody>
void forall2(T1 &irange, T2 &jrange, LoopBody &&body){
  
  /* dim3 tpb(irange.tpb,jrange.tpb,1); */
  /* dim3 blocks(irange.blocks,jrange.blocks,1); */
  
  /* forall2kernel<<<blocks,tpb>>>(irange.start,irange.end,jrange.start,jrange.end,body); */
  forall2async(irange,jrange,body);
  cudaStreamSynchronize(0);
}
#endif
