#ifndef __POLICIES_FOR_SW4__
#define __POLICIES_FOR_SW4__
namespace Policies {
template <int T>
class Loop1 {
 public:
  const static int value = T;
};
struct Cuda {};
struct OpenMP {};
struct Default {};
template <typename T1, typename T2>
class Policy {
 public:
  typedef int type;  // to make compilation fail
};
/* template <typename T2> */
/* class Policy<Cuda, T2> { */
/*  public: */
/*   typedef RAJA::cuda_exec<T2::value> type; */
/* }; */

template <typename T1, typename... args>
class Policy2 {
 public:
  typedef int type;  // to make compilation fail
};
/* template <typename... args> */
/* class Policy2<Cuda, args...> { */
/*  public: */
/*   typedef RAJA::cuda_exec<1024> type; */
/* }; */
}  // namespace Policies
#endif




template <typename Func>
__global__ void forallkernel(int start, int N, Func f) {
  int tid = start + threadIdx.x + blockIdx.x * blockDim.x;
  if (tid < N) f(tid);
}
template <typename LoopBody>
void forall(int start, int end, LoopBody &&body) {
  int tpb = 1024;
  int blocks = (end - start) / tpb;
  blocks = ((end - start) % tpb == 0) ? blocks : blocks + 1;
  // printf("Launching the kernel blocks= %d tpb= %d \n",blocks,tpb);
  forallkernel<<<blocks, tpb>>>(start, end, body);
  cudaDeviceSynchronize();
}
template <typename LoopBody>
void forallasync(int start, int end, LoopBody &&body) {
  int tpb = 1024;
  int blocks = (end - start) / tpb;
  blocks = ((end - start) % tpb == 0) ? blocks : blocks + 1;
  // printf("Launching the kernel blocks= %d tpb= %d \n",blocks,tpb);
  forallkernel<<<blocks, tpb>>>(start, end, body);
  // cudaDeviceSynchronize();
}






template <typename T, typename Func>
__global__ void forallgskernel(T start, T N, Func f) {
  for (T i = start + threadIdx.x + blockIdx.x * blockDim.x; i < N;
       i += blockDim.x * gridDim.x)
    f(i);
}
template <int N, typename T, typename LoopBody>
void forallX(T start, T end, LoopBody &&body) {
  // On an SM7, each SM can have at most 64 warps or 2048 threads
  // V100 has 80 SMs
  // Max number of resident threads is thus 80x2048 or 80x64 warps
  int tpb = N;
  int blocks = (end - start) / tpb;
  blocks = ((end - start) % tpb == 0) ? blocks : blocks + 1;
  /* if ((blocks * tpb) <= 80 * 2048) { */
  /*   printf("Launching the kernel blocks= %d tpb= %d with N %d\n", blocks, tpb, */
  /*          N); */
  /*   forallkernel<<<blocks, tpb>>>(start, end, body); */
  /* } else { */
  /*   printf("Max resident count exceeded %d >163840 \n", blocks * tpb); */
    blocks = 80 * 2048 / tpb;
    //    blocks = 80*32 / tpb; // Really bad idea
    /* printf("Launching the GS kernel blocks= %d tpb= %dwith N %d\n", blocks, tpb, */
    /*        N); */
    forallgskernel<<<blocks, tpb>>>(start, end, body);
    //cudaDeviceSynchronize();
    // }
}


// 3D kernels

template <typename Func>
__global__ void forall3kernel(const int start0, const int N0, const int start1,
                              const int N1, const int start2, const int N2,
                              Func f) {
  int tid0 = start0 + threadIdx.x + blockIdx.x * blockDim.x;
  int tid1 = start1 + threadIdx.y + blockIdx.y * blockDim.y;
  int tid2 = start2 + threadIdx.z + blockIdx.z * blockDim.z;
  if ((tid0 < N0) && (tid1 < N1) && (tid2 < N2)) f(tid0, tid1, tid2);
}

template <int N, typename Func>
__global__ void forall3kernel(const int start0, const int N0, const int start1,
                              const int N1, const int start2, const int N2,
                              Func f) {
  int tid0 = start0 + threadIdx.x + blockIdx.x * blockDim.x;
  int tid1 = start1 + threadIdx.y + blockIdx.y * blockDim.y;
  int tid2 = start2 + threadIdx.z + blockIdx.z * blockDim.z;
  if ((tid0 < N0) && (tid1 < N1) && (tid2 < N2)) f(tid0, tid1, tid2);
}

template <typename Func>
__global__ void forall3gskernel(int start0, int N0, int start1, int N1,
                                int start2, int N2, Func f) {
  for (int i = start0 + threadIdx.x + blockIdx.x * blockDim.x; i < N0;
       i += blockDim.x * gridDim.x)
    for (int j = start1 + threadIdx.y + blockIdx.y * blockDim.y; j < N1;
         j += blockDim.y * gridDim.y)
      for (int k = start2 + threadIdx.z + blockIdx.z * blockDim.z; k < N2;
           k += blockDim.z * gridDim.z)
        f(i, j, k);
}

template <typename LoopBody>
void forall3(int start0, int end0, int start1, int end1, int start2, int end2,
             LoopBody &&body) {
  int tpb0 = 16;
  int tpb1 = 16;
  int tpb2 = 1024 / (tpb0 * tpb1);

  int block0 = (end0 - start0) / tpb0;
  block0 = ((end0 - start0) % tpb0 == 0) ? block0 : block0 + 1;
  int block1 = (end1 - start1) / tpb1;
  block1 = ((end1 - start1) % tpb1 == 0) ? block1 : block1 + 1;
  int block2 = (end2 - start2) / tpb2;
  block2 = ((end2 - start2) % tpb2 == 0) ? block2 : block2 + 1;

  // std::cout << " BLOCKS " << block0 << " " << block1 << " " << block2 <<
  // "\n";
  dim3 tpb(tpb0, tpb1, tpb2);
  dim3 blocks(block0, block1, block2);

  //printf("Launching the kernel 3d \n");
  forall3kernel<<<blocks, tpb>>>(start0, end0, start1, end1, start2, end2,
                                 body);
  //cudaDeviceSynchronize();
}

template <int N, typename LoopBody>
void forall3X(int start0, int end0, int start1, int end1, int start2, int end2,
             LoopBody &&body) {

  int tpbb = N;
  int tpb0 = 16;
  int tpb1 = 16;
  int tpb2 = tpbb / (tpb0 * tpb1);

  int blockss = 80*2048/ tpbb;
  int block0 = 4;
  int block1 = 4;
  int block2 = blockss/(block0*block1);
  

  // std::cout << " BLOCKS " << block0 << " " << block1 << " " << block2 <<
  // "\n";
  dim3 tpb(tpb0, tpb1, tpb2);
  dim3 blocks(block0, block1, block2);

  //printf("Launching the kernel 3d \n");
  forall3gskernel<<<blocks, tpb>>>(start0, end0, start1, end1, start2, end2,
                                 body);
  // cudaDeviceSynchronize();
}

// template <typename LoopBody>
// void forall3X<64>(int start0, int end0, int start1, int end1, int start2, int end2,
//              LoopBody &&body) {

//   int tpbb =64;
//   int tpb0 = 16;
//   int tpb1 = 16;
//   int tpb2 = tpbb / (tpb0 * tpb1);

//   int blockss = 80*2048/ tpbb;
//   int block0 = 4;
//   int block1 = 4;
//   int block2 = blockss/(block0*block1);
  

//   // std::cout << " BLOCKS " << block0 << " " << block1 << " " << block2 <<
//   // "\n";
//   dim3 tpb(tpb0, tpb1, tpb2);
//   dim3 blocks(block0, block1, block2);

//   //printf("Launching the kernel 3d \n");
//   forall3gskernel<<<blocks, tpb>>>(start0, end0, start1, end1, start2, end2,
//                                  body);
//   // cudaDeviceSynchronize();
// }

class Parray{
public:
  int ni,nj,nk,nc;
  double *data;
  Parray(int nii, int nji, int nki,int nci):ni(nii),nj(nji),nk(nki),nc(nci){
    cudaMallocManaged(&data,ni*nj*nk*nc*sizeof(double));
    std::cout<<"Created Parray ("<<ni<<" "<<nj<<" "<<nk<<" "<<nk<<" "<<nc<<")\n";
  }
  void set(double value){
    double *l_data = data;
    forall(0,ni*nj*nk*nc,[=]__device__(int i){
	l_data[i]=value;
      });
  }
  __device__ inline double& operator()(int i, int j, int k,int c) const {
    return data[ i+ni*j+ni*nj*k+ni*nj*nk*c];
  }
  __device__ inline double& operator()(unsigned int I, unsigned int J) const {
    int i = (I>>8) & 0x00ff;
    int j = I & 0x00ff;
    int k = (J>>8) & 0x00ff;
    int c =  J& 0x00ff;
    //printf("%d %d %d %d \n",i,j,k,c);
    return data[ i+ni*j+ni*nj*k+ni*nj*nk*c];
  }
  int offset(int i, int j, int k, int c){
    return i+ni*j+ni*nj*k+ni*nj*nk*c;
  }
  int offset2(unsigned int I, unsigned int J){
    int i = (I>>8) & 0x00ff;
    int j = I & 0x00ff;
    int k = (J>>8) & 0x00ff;
    int c =  J& 0x00ff;
    //printf("%d %d %d %d \n",i,j,k,c);
    return i+ni*j+ni*nj*k+ni*nj*nk*c;
  }
    
};
