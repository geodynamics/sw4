#ifndef __POLICIES_FOR_SW4__
#define __POLICIES_FOR_SW4__
namespace Policies {
template <sw4_type T>
class Loop1 {
 public:
  const static sw4_type value = T;
};
struct Cuda {};
struct OpenMP {};
struct Default {};
template <typename T1, typename T2>
class Policy {
 public:
  typedef sw4_type type;  // to make compilation fail
};
/* template <typename T2> */
/* class Policy<Cuda, T2> { */
/*  public: */
/*   typedef RAJA::cuda_exec<T2::value> type; */
/* }; */

template <typename T1, typename... args>
class Policy2 {
 public:
  typedef sw4_type type;  // to make compilation fail
};
/* template <typename... args> */
/* class Policy2<Cuda, args...> { */
/*  public: */
/*   typedef RAJA::cuda_exec<1024> type; */
/* }; */
}  // namespace Policies
#endif

template <typename Func>
__global__ void forallkernel(sw4_type start, sw4_type N, Func f) {
  sw4_type tid = start + threadIdx.x + blockIdx.x * blockDim.x;
  if (tid < N) f(tid);
}
template <typename LoopBody>
void forall(sw4_type start, sw4_type end, LoopBody &&body) {
  sw4_type tpb = 1024;
  sw4_type blocks = (end - start) / tpb;
  blocks = ((end - start) % tpb == 0) ? blocks : blocks + 1;
  // printf("Launching the kernel blocks= %d tpb= %d \n",blocks,tpb);
  forallkernel<<<blocks, tpb>>>(start, end, body);
  cudaDeviceSynchronize();
}
template <typename LoopBody>
void forallasync(sw4_type start, sw4_type end, LoopBody &&body) {
  sw4_type tpb = 1024;
  sw4_type blocks = (end - start) / tpb;
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
template <sw4_type N, typename T, typename LoopBody>
void forallX(T start, T end, LoopBody &&body) {
  // On an SM7, each SM can have at most 64 warps or 2048 threads
  // V100 has 80 SMs
  // Max number of resident threads is thus 80x2048 or 80x64 warps
  sw4_type tpb = N;
  sw4_type blocks = (end - start) / tpb;
  blocks = ((end - start) % tpb == 0) ? blocks : blocks + 1;
  /* if ((blocks * tpb) <= 80 * 2048) { */
  /*   printf("Launching the kernel blocks= %d tpb= %d with N %d\n", blocks,
   * tpb, */
  /*          N); */
  /*   forallkernel<<<blocks, tpb>>>(start, end, body); */
  /* } else { */
  /*   printf("Max resident count exceeded %d >163840 \n", blocks * tpb); */
  blocks = 80 * 2048 / tpb;
  //    blocks = 80*32 / tpb; // Really bad idea
  /* printf("Launching the GS kernel blocks= %d tpb= %dwith N %d\n", blocks,
   * tpb, */
  /*        N); */
  forallgskernel<<<blocks, tpb>>>(start, end, body);
  // cudaDeviceSynchronize();
  // }
}

// 3D kernels

template <typename Func>
__global__ void forall3kernel(const sw4_type start0, const sw4_type N0, const sw4_type start1,
                              const sw4_type N1, const sw4_type start2, const sw4_type N2,
                              Func f) {
  sw4_type tid0 = start0 + threadIdx.x + blockIdx.x * blockDim.x;
  sw4_type tid1 = start1 + threadIdx.y + blockIdx.y * blockDim.y;
  sw4_type tid2 = start2 + threadIdx.z + blockIdx.z * blockDim.z;
  if ((tid0 < N0) && (tid1 < N1) && (tid2 < N2)) f(tid0, tid1, tid2);
}

template <sw4_type N, typename Func>
__global__ void forall3kernel(const sw4_type start0, const sw4_type N0, const sw4_type start1,
                              const sw4_type N1, const sw4_type start2, const sw4_type N2,
                              Func f) {
  sw4_type tid0 = start0 + threadIdx.x + blockIdx.x * blockDim.x;
  sw4_type tid1 = start1 + threadIdx.y + blockIdx.y * blockDim.y;
  sw4_type tid2 = start2 + threadIdx.z + blockIdx.z * blockDim.z;
  if ((tid0 < N0) && (tid1 < N1) && (tid2 < N2)) f(tid0, tid1, tid2);
}

template <typename Func>
__global__ void forall3gskernel(sw4_type start0, sw4_type N0, sw4_type start1, sw4_type N1,
                                sw4_type start2, sw4_type N2, Func f) {
  for (sw4_type i = start0 + threadIdx.x + blockIdx.x * blockDim.x; i < N0;
       i += blockDim.x * gridDim.x)
    for (sw4_type j = start1 + threadIdx.y + blockIdx.y * blockDim.y; j < N1;
         j += blockDim.y * gridDim.y)
      for (sw4_type k = start2 + threadIdx.z + blockIdx.z * blockDim.z; k < N2;
           k += blockDim.z * gridDim.z)
        f(i, j, k);
}

template <typename LoopBody>
void forall3(sw4_type start0, sw4_type end0, sw4_type start1, sw4_type end1, sw4_type start2, sw4_type end2,
             LoopBody &&body) {
  sw4_type tpb0 = 16;
  sw4_type tpb1 = 16;
  sw4_type tpb2 = 1024 / (tpb0 * tpb1);

  sw4_type block0 = (end0 - start0) / tpb0;
  block0 = ((end0 - start0) % tpb0 == 0) ? block0 : block0 + 1;
  sw4_type block1 = (end1 - start1) / tpb1;
  block1 = ((end1 - start1) % tpb1 == 0) ? block1 : block1 + 1;
  sw4_type block2 = (end2 - start2) / tpb2;
  block2 = ((end2 - start2) % tpb2 == 0) ? block2 : block2 + 1;

  // std::cout << " BLOCKS " << block0 << " " << block1 << " " << block2 <<
  // "\n";
  dim3 tpb(tpb0, tpb1, tpb2);
  dim3 blocks(block0, block1, block2);

  // printf("Launching the kernel 3d \n");
  forall3kernel<<<blocks, tpb>>>(start0, end0, start1, end1, start2, end2,
                                 body);
  // cudaDeviceSynchronize();
}

template <sw4_type N, typename LoopBody>
void forall3X(sw4_type start0, sw4_type end0, sw4_type start1, sw4_type end1, sw4_type start2, sw4_type end2,
              LoopBody &&body) {
  sw4_type tpbb = N;
  sw4_type tpb0 = 16;
  sw4_type tpb1 = 16;
  sw4_type tpb2 = tpbb / (tpb0 * tpb1);

  sw4_type blockss = 80 * 2048 / tpbb;
  sw4_type block0 = 4;
  sw4_type block1 = 4;
  sw4_type block2 = blockss / (block0 * block1);

  // std::cout << " BLOCKS " << block0 << " " << block1 << " " << block2 <<
  // "\n";
  dim3 tpb(tpb0, tpb1, tpb2);
  dim3 blocks(block0, block1, block2);

  // printf("Launching the kernel 3d \n");
  forall3gskernel<<<blocks, tpb>>>(start0, end0, start1, end1, start2, end2,
                                   body);
  // cudaDeviceSynchronize();
}

// template <typename LoopBody>
// void forall3X<64>(sw4_type start0, sw4_type end0, sw4_type start1, sw4_type end1, sw4_type start2, sw4_type
// end2,
//              LoopBody &&body) {

//   sw4_type tpbb =64;
//   sw4_type tpb0 = 16;
//   sw4_type tpb1 = 16;
//   sw4_type tpb2 = tpbb / (tpb0 * tpb1);

//   sw4_type blockss = 80*2048/ tpbb;
//   sw4_type block0 = 4;
//   sw4_type block1 = 4;
//   sw4_type block2 = blockss/(block0*block1);

//   // std::cout << " BLOCKS " << block0 << " " << block1 << " " << block2 <<
//   // "\n";
//   dim3 tpb(tpb0, tpb1, tpb2);
//   dim3 blocks(block0, block1, block2);

//   //printf("Launching the kernel 3d \n");
//   forall3gskernel<<<blocks, tpb>>>(start0, end0, start1, end1, start2, end2,
//                                  body);
//   // cudaDeviceSynchronize();
// }

class Parray {
 public:
  sw4_type ni, nj, nk, nc;
  double *data;
  Parray(sw4_type nii, sw4_type nji, sw4_type nki, sw4_type nci)
      : ni(nii), nj(nji), nk(nki), nc(nci) {
    cudaMallocManaged(&data, ni * nj * nk * nc * sizeof(double));
    std::cout << "Created Parray (" << ni << " " << nj << " " << nk << " " << nk
              << " " << nc << ")\n";
  }
  void set(double value) {
    double *l_data = data;
    forall(0, ni * nj * nk * nc, [=] __device__(sw4_type i) { l_data[i] = value; });
  }
  __device__ inline double &operator()(sw4_type i, sw4_type j, sw4_type k, sw4_type c) const {
    return data[i + ni * j + ni * nj * k + ni * nj * nk * c];
  }
  __device__ inline double &operator()(unsigned sw4_type I, unsigned sw4_type J) const {
    sw4_type i = (I >> 8) & 0x00ff;
    sw4_type j = I & 0x00ff;
    sw4_type k = (J >> 8) & 0x00ff;
    sw4_type c = J & 0x00ff;
    // printf("%d %d %d %d \n",i,j,k,c);
    return data[i + ni * j + ni * nj * k + ni * nj * nk * c];
  }
  sw4_type offset(sw4_type i, sw4_type j, sw4_type k, sw4_type c) {
    return i + ni * j + ni * nj * k + ni * nj * nk * c;
  }
  sw4_type offset2(unsigned sw4_type I, unsigned sw4_type J) {
    sw4_type i = (I >> 8) & 0x00ff;
    sw4_type j = I & 0x00ff;
    sw4_type k = (J >> 8) & 0x00ff;
    sw4_type c = J & 0x00ff;
    // printf("%d %d %d %d \n",i,j,k,c);
    return i + ni * j + ni * nj * k + ni * nj * nk * c;
  }
};
