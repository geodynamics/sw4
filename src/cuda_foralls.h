#ifndef __CUDA_FORALLS_H__
#define __CUDA_FORALLS_H__
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

#include "mpi.h"
#include "policies.h"
std::vector<int> factors(int N);
std::vector<int> factors(int N, int start);
#ifdef ENABLE_CUDA
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

template <int N, typename LoopBody>
void forall(int start, int end, LoopBody &&body) {
  int tpb = 1024;
  int blocks = (end - start) / tpb;
  blocks = ((end - start) % tpb == 0) ? blocks : blocks + 1;
  printf("Launching the kernel blocks= %d tpb= %d on line %d\n", blocks, tpb,
         N);
  forallkernel<<<blocks, tpb>>>(start, end, body);
  cudaDeviceSynchronize();
}

template <typename Func>
__global__ void forallkernelB(int start, int N, Func f) {
  int tid = start + threadIdx.x + blockIdx.x * blockDim.x;
  const int B = 8;
  for (int i = tid; i < N; i += B * blockDim.x * gridDim.x) {
    f(i);
    // int ii=i+blockDim.x * gridDim.x;
#pragma unroll(B - 1)
    for (int ii = 1; ii < B; ii++) {
      int iii = i + ii * blockDim.x * gridDim.x;
      if (iii < N) f(iii);
    }
  }
}
template <typename LoopBody>
void forallB(int start, int end, LoopBody &&body) {
  int tpb = 1024;
  int blocks = 52;
  // blocks=((end-start)%tpb==0)?blocks:blocks+1;
  printf("Launching the kernel\n");
  forallkernelB<<<blocks, tpb>>>(start, end, body);
  cudaDeviceSynchronize();
}

template <int N>
class Range {
 public:
  Range(int istart, int iend) : start(istart), end(iend), tpb(N) {
    blocks = (end - start) / N;
    blocks = ((end - start) % N == 0) ? blocks : blocks + 1;
    invalid = false;
    if (blocks <= 0) invalid = true;
  };
  int start;
  int end;
  int blocks;
  int tpb;
  bool invalid;
};

template <int N, int M>
class RangeGS {
 public:
  RangeGS(int istart, int iend) : start(istart), end(iend), tpb(N), blocks(M){};
  int start;
  int end;
  int blocks;
  int tpb;
};

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
template <int N, typename Tag, typename Func>
__global__ void forall3kernel(Tag t, const int start0, const int N0,
                              const int start1, const int N1, const int start2,
                              const int N2, Func f) {
  int tid0 = start0 + threadIdx.x + blockIdx.x * blockDim.x;
  int tid1 = start1 + threadIdx.y + blockIdx.y * blockDim.y;
  int tid2 = start2 + threadIdx.z + blockIdx.z * blockDim.z;
  if ((tid0 < N0) && (tid1 < N1) && (tid2 < N2)) f(t, tid0, tid1, tid2);
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

  printf("Launching the kernel 3d \n");
  forall3kernel<<<blocks, tpb>>>(start0, end0, start1, end1, start2, end2,
                                 body);
  cudaDeviceSynchronize();
}

template <typename T1, typename T2, typename T3, typename LoopBody>
void forall3async(T1 &irange, T2 &jrange, T3 &krange, LoopBody &&body) {
  if (irange.invalid || jrange.invalid || krange.invalid) return;
  dim3 tpb(irange.tpb, jrange.tpb, krange.tpb);
  dim3 blocks(irange.blocks, jrange.blocks, krange.blocks);
  // std::cout<<"forall launch tpb"<<irange.tpb<<" "<<jrange.tpb<<"
  // "<<krange.tpb<<"\n"; std::cout<<"forall launch blocks"<<irange.blocks<<"
  // "<<jrange.blocks<<" "<<krange.blocks<<"\n";
  static int firstcall = 0;
  if (!firstcall) {
    firstcall = 1;
    int minGridSize, maxBlockSize;
    if (cudaOccupancyMaxPotentialBlockSize(&minGridSize, &maxBlockSize,
                                           forall3kernel<LoopBody>, 0,
                                           0) != cudaSuccess) {
      std::cerr << "cudaOccupancyMaxPotentialBlockSize Failed\n";
      abort();
    } else {
      // std::cerr<<"Min grid size "<<minGridSize<<" maxblock
      // size"<<maxBlockSize<<" Actual grid =
      // "<<irange.blocks*jrange.blocks*krange.blocks<<"\n";
    }
    if ((irange.tpb * jrange.tpb * krange.tpb) > maxBlockSize) {
      std::cerr << " Block size too large in forall3async"
                << (irange.tpb * jrange.tpb * krange.tpb) << " > "
                << maxBlockSize << "n" << std::flush;
      abort();
    }
  }

  forall3kernel<<<blocks, tpb>>>(irange.start, irange.end, jrange.start,
                                 jrange.end, krange.start, krange.end, body);
}
template <typename T1, typename T2, typename T3, typename LoopBody>
void forall3(T1 &irange, T2 &jrange, T3 &krange, LoopBody &&body) {
  // dim3 tpb(irange.tpb,jrange.tpb,krange.tpb);
  // dim3 blocks(irange.blocks,jrange.blocks,krange.blocks);

  // forall3kernel<<<blocks,tpb>>>(irange.start,irange.end,jrange.start,jrange.end,krange.start,krange.end,body);
  forall3async(irange, jrange, krange, body);
  cudaStreamSynchronize(0);
}

template <typename T1, typename T2, typename T3, typename LoopBody>
void forall3GSasync(T1 &irange, T2 &jrange, T3 &krange, LoopBody &&body) {
  dim3 tpb(irange.tpb, jrange.tpb, krange.tpb);
  dim3 blocks(irange.blocks, jrange.blocks, krange.blocks);

  forall3gskernel<<<blocks, tpb>>>(irange.start, irange.end, jrange.start,
                                   jrange.end, krange.start, krange.end, body);
}

template <typename T1, typename T2, typename T3, typename LoopBody>
void forall3GS(T1 &irange, T2 &jrange, T3 &krange, LoopBody &&body) {
  // dim3 tpb(irange.tpb,jrange.tpb,krange.tpb);
  // dim3 blocks(irange.blocks,jrange.blocks,krange.blocks);

  // forall3gskernel<<<blocks,tpb>>>(irange.start,irange.end,jrange.start,jrange.end,krange.start,krange.end,body);
  forall3async(irange, jrange, krange, body);
  cudaStreamSynchronize(0);
}

// Forall2

template <typename Func>
__global__ void forall2kernel(const int start0, const int N0, const int start1,
                              const int N1, Func f) {
  int tid0 = start0 + threadIdx.x + blockIdx.x * blockDim.x;
  int tid1 = start1 + threadIdx.y + blockIdx.y * blockDim.y;
  if ((tid0 < N0) && (tid1 < N1)) f(tid0, tid1);
}

template <typename T1, typename T2, typename LoopBody>
void forall2async(T1 &irange, T2 &jrange, LoopBody &&body) {
  dim3 tpb(irange.tpb, jrange.tpb, 1);
  dim3 blocks(irange.blocks, jrange.blocks, 1);

  forall2kernel<<<blocks, tpb>>>(irange.start, irange.end, jrange.start,
                                 jrange.end, body);
}

template <typename T1, typename T2, typename LoopBody>
void forall2(T1 &irange, T2 &jrange, LoopBody &&body) {
  /* dim3 tpb(irange.tpb,jrange.tpb,1); */
  /* dim3 blocks(irange.blocks,jrange.blocks,1); */

  /* forall2kernel<<<blocks,tpb>>>(irange.start,irange.end,jrange.start,jrange.end,body);
   */
  forall2async(irange, jrange, body);
  cudaStreamSynchronize(0);
}

// AutoTuning Code
template <int N, int M, int L>
class RangeAT {
 public:
  RangeAT(int istart, int iend, int jstart, int jend, int kstart, int kend)
      : istart(istart),
        iend(iend),
        jstart(jstart),
        jend(jend),
        kstart(kstart),
        kend(kend) {
    auto curr =
        std::make_tuple((iend - istart), (jend - jstart), (kend - kstart));
    if (files.find(curr) == files.end()) {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
      // std::cout << "Opening the file \n";
      std::stringstream s;
      s << "Times" << myRank << "_" << M << "_" << L << "_" << (iend - istart)
        << "_" << (jend - jstart) << "_" << (kend - kstart) << ".dat";
      auto &ofile = files[curr];
      ofile.open(s.str());
      // iles[curr]=std::move(ofile);
    }
    if (mintime.find(std::make_tuple((iend - istart), (jend - jstart),
                                     (kend - kstart))) == mintime.end())
      mintime[std::make_tuple((iend - istart), (jend - jstart),
                              (kend - kstart))] = 1.0e89;
    if (vloca.find(std::make_tuple((iend - istart), (jend - jstart),
                                   (kend - kstart))) == vloca.end())
      vloca[std::make_tuple((iend - istart), (jend - jstart),
                            (kend - kstart))] = 0;
    if (map.find(std::make_tuple((iend - istart), (jend - jstart),
                                 (kend - kstart))) == map.end()) {
      int max_block_size = floor(65536.0 / (N + 2)) - 1;
      max_block_size = N;
      // max_block_size = floor((256*255.0)/(N+1));
      files[curr] << "#Largest allowable block size is " << max_block_size
                  << " for Regs Per thread of " << N << " on line " << M
                  << " \n";
      files[curr] << "#LOOP SIZES " << (iend - istart) << " " << (jend - jstart)
                  << " " << (kend - kstart) << "\n";
      int ii, jj, kk;
      std::vector<std::tuple<int, int, int, int, int, int>> mins;
      if (0) {
        // Insert arbitrary configuration into vector of launch configs
        int ii = ceil(float(iend - istart) / 16);
        int jj = ceil(float(jend - jstart) / 4);
        int kk = ceil(float(kend - kstart) / 6);
        confs[curr].clear();
        dim3 tpb(16, 4, 6);
        dim3 blks(ii, jj, kk);
        confs[curr].push_back(std::make_tuple(tpb, blks));
      }
      const int offset = 0;  // Controlling the search space
      for (int block_size = max_block_size - offset;
           block_size <= max_block_size; block_size++) {
        auto f = factors(block_size, 32);
        for (const int &i : f) {
          int rest = block_size / i;
          auto ff = factors(rest);
          for (const int &j : ff) {
            int k = rest / j;
            assert(i * j * k == N);
            ii = ceil(float(iend - istart) / i);
            jj = ceil(float(jend - jstart) / j);
            kk = ceil(float(kend - kstart) / k);

            dim3 tpb(i, j, k);
            dim3 blks(ii, jj, kk);
            if (k <= 64) confs[curr].push_back(std::make_tuple(tpb, blks));

            assert(ii * i >= (iend - istart));
            assert(jj * j >= (jend - jstart));
            assert(kk * k >= (kend - kstart));
          }
        }
      }
    }
  }
  static std::map<std::tuple<int, int, int>, std::tuple<dim3, dim3>> map;
  static std::map<std::tuple<int, int, int>,
                  std::vector<std::tuple<dim3, dim3>>>
      confs;
  static std::map<std::tuple<int, int, int>, float> mintime;
  static std::map<std::tuple<int, int, int>, int> vloca;
  static std::map<std::tuple<int, int, int>, std::ofstream> files;
  static dim3 mint, minb;
  int istart, jstart, kstart;
  int iend, jend, kend;
  int ib, jb, kb;
  int ic, jc, kc;
  int blocks;
  int tpb;
};

template <int N, typename T1, typename LoopBody>
void forall3asyncAT(T1 &range3, LoopBody &&body) {
  dim3 tpb(range3.ic, range3.jc, range3.kc);
  dim3 blocks(range3.ib, range3.jb, range3.kb);

  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  auto curr = std::make_tuple((range3.iend - range3.istart),
                              (range3.jend - range3.jstart),
                              (range3.kend - range3.kstart));
#define SW4_VERBOSE
#ifdef SW4_VERBOSE
  // Calculate max block size and compare it to early calculations
  int minGridSize;
  int blockSize;

  if (range3.vloca[curr] == 0) {
    SW4_CheckDeviceError(cudaOccupancyMaxPotentialBlockSize(
        &minGridSize, &blockSize, forall3kernel<N, LoopBody>, 0, 0));
    range3.files[curr] << "#LINE " << N << " MAX BLOCK SIZE FOR OCC CALC "
                       << blockSize << " minGridSize " << minGridSize << " \n";
    int numBlocks;
    SW4_CheckDeviceError(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
        &numBlocks, forall3kernel<N, LoopBody>, blockSize, 0));
    range3.files[curr] << "#LINE " << N
                       << " maximum number of active blocks per SM "
                       << numBlocks << " \n";
  }
// end max block size
#endif

  if (range3.vloca[curr] < range3.confs[curr].size()) {
    auto &i = range3.confs[curr][range3.vloca[curr]];
    range3.vloca[curr]++;

    int n = std::get<1>(i).x * std::get<1>(i).y * std::get<1>(i).z;

    cudaEventRecord(start);
    forall3kernel<N><<<std::get<1>(i), std::get<0>(i)>>>(
        range3.istart, range3.iend, range3.jstart, range3.jend, range3.kstart,
        range3.kend, body);
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float ms = 0;
    cudaEventElapsedTime(&ms, start, stop);
    if (cudaGetLastError() == cudaSuccess) {
      range3.files[curr] << n << " " << ms << " " << std::get<0>(i).x << " "
                         << std::get<0>(i).y << " " << std::get<0>(i).z << " "
                         << range3.vloca[curr] << "\n";
      if (range3.mintime[curr] > ms) {
        range3.mintime[curr] = ms;
        range3.map[curr] = std::make_tuple(std::get<0>(i), std::get<1>(i));
      }
    } else {
      std::cout << range3.vloca[curr] << " " << N << " " << std::get<0>(i).x
                << " " << std::get<0>(i).y << " " << std::get<0>(i).z
                << " TPB FAILED \n";
      std::cout << range3.vloca[curr] << " " << N << " " << std::get<1>(i).x
                << " " << std::get<1>(i).y << " " << std::get<1>(i).z
                << " BLKS FAILED \n";
    }

    if (range3.vloca[curr] == range3.confs[curr].size()) {
      auto &mint = std::get<0>(range3.map[curr]);
      range3.files[curr] << "&\n#BEST CONFIG " << range3.mintime[curr] << " "
                         << mint.x << " " << mint.y << " " << mint.z << " "
                         << range3.confs[curr].size() << "\n";
      auto &minb = std::get<1>(range3.map[curr]);
      int minblocks = minb.x * minb.y * minb.z;

      range3.files[curr] << minblocks << " " << range3.mintime[curr] << " "
                         << mint.x << " " << mint.y << " " << mint.z << " "
                         << range3.confs[curr].size() << "\n";
    }
  } else {
    auto &v = range3.map[std::make_tuple((range3.iend - range3.istart),
                                         (range3.jend - range3.jstart),
                                         (range3.kend - range3.kstart))];
    forall3kernel<N><<<std::get<1>(v), std::get<0>(v)>>>(
        range3.istart, range3.iend, range3.jstart, range3.jend, range3.kstart,
        range3.kend, body);
  }
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
  if ((blocks * tpb) <= 80 * 2048) {
    printf("Launching the kernel blocks= %d tpb= %d with N %d\n", blocks, tpb,
           N);
    forallkernel<<<blocks, tpb>>>(start, end, body);
  } else {
    printf("Max resident count exceeded %d >163840 \n", blocks * tpb);
    blocks = 80 * 2048 / tpb;
    printf("Launching the GS kernel blocks= %d tpb= %dwith N %d\n", blocks, tpb,
           N);
    forallgskernel<<<blocks, tpb>>>(start, end, body);
    cudaDeviceSynchronize();
  }
}

#endif

template <int N, typename LoopBody>
void forall3X(int start0, int end0, int start1, int end1, int start2, int end2,
              LoopBody &&body) {
  int tpbb = N;
  int tpb0 = 8;
  int tpb1 = 8;
  int tpb2 = tpbb / (tpb0 * tpb1);

  int blockss = 80 * 2048 / tpbb;
  int block0 = 20;
  int block1 = 8;
  int block2 = blockss / (block0 * block1);

  // std::cout << " BLOCKS " << block0 << " " << block1 << " " << block2 <<
  // "\n";
  dim3 tpb(tpb0, tpb1, tpb2);
  dim3 blocks(block0, block1, block2);

  // printf("Launching the kernel 3d \n");
  forall3gskernel<<<blocks, tpb>>>(start0, end0, start1, end1, start2, end2,
                                   body);
  // cudaDeviceSynchronize();
}
template <int N, typename T1, typename T2, typename T3, typename LoopBody>
void forall3async(T1 &irange, T2 &jrange, T3 &krange, LoopBody &&body) {
  if (irange.invalid || jrange.invalid || krange.invalid) return;
  dim3 tpb(irange.tpb, jrange.tpb, krange.tpb);
  dim3 blocks(irange.blocks, jrange.blocks, krange.blocks);
  // std::cout<<"forall launch tpb"<<irange.tpb<<" "<<jrange.tpb<<"
  // "<<krange.tpb<<"\n"; std::cout<<"forall launch blocks"<<irange.blocks<<"
  // "<<jrange.blocks<<" "<<krange.blocks<<"\n";
  forall3kernel<N><<<blocks, tpb>>>(irange.start, irange.end, jrange.start,
                                    jrange.end, krange.start, krange.end, body);
}

template <int N>
class Tclass {};

template <int N, typename Tag, typename T1, typename T2, typename T3,
          typename LoopBody>
void forall3async(Tag &t, T1 &irange, T2 &jrange, T3 &krange, LoopBody &&body) {
  if (irange.invalid || jrange.invalid || krange.invalid) return;
  dim3 tpb(irange.tpb, jrange.tpb, krange.tpb);
  dim3 blocks(irange.blocks, jrange.blocks, krange.blocks);
  // std::cout<<"forall launch tpb"<<irange.tpb<<" "<<jrange.tpb<<"
  // "<<krange.tpb<<"\n"; std::cout<<"forall launch blocks"<<irange.blocks<<"
  // "<<jrange.blocks<<" "<<krange.blocks<<"\n";
  forall3kernel<N><<<blocks, tpb>>>(t, irange.start, irange.end, jrange.start,
                                    jrange.end, krange.start, krange.end, body);
}

// The multiforall for kernel fusion
template <typename T, typename F0, typename F1, typename F2, typename F3>
__global__ void multiforallkernel(T start0, T end0, F0 f0, T start1, T end1,
                                  F1 f1, T start2, T end2, F2 f2, T start3,
                                  T end3, F3 f3) {
  for (T i = start0 + threadIdx.x + blockIdx.x * blockDim.x; i < end0;
       i += blockDim.x * gridDim.x)
    f0(i);
  for (T i = start1 + threadIdx.x + blockIdx.x * blockDim.x; i < end1;
       i += blockDim.x * gridDim.x)
    f1(i);
  for (T i = start2 + threadIdx.x + blockIdx.x * blockDim.x; i < end2;
       i += blockDim.x * gridDim.x)
    f2(i);
  for (T i = start3 + threadIdx.x + blockIdx.x * blockDim.x; i < end3;
       i += blockDim.x * gridDim.x)
    f3(i);
}
template <int N, typename T, typename LB0, typename LB1, typename LB2,
          typename LB3>
void multiforall(T start0, T end0, LB0 &&body0, T start1, T end1, LB1 &&body1,
                 T start2, T end2, LB2 &&body2, T start3, T end3, LB3 &&body3) {
  int blocks = 80 * 2048 / N;  // WARNING HARDWIRED FOR V100
  multiforallkernel<<<blocks, N>>>(start0, end0, body0, start1, end1, body1,
                                   start2, end2, body2, start3, end3, body3);
}

// Generalized mforall kernel
template <typename T, typename LoopBody>
__device__ inline void runner(T start, T end, LoopBody f) {
  for (T i = start + threadIdx.x + blockIdx.x * blockDim.x; i < end;
       i += blockDim.x * gridDim.x)
    f(i);
}
template <typename T, typename LoopBody, class... Args>
__device__ inline void runner(T start, T end, LoopBody f, Args... args) {
  for (T i = start + threadIdx.x + blockIdx.x * blockDim.x; i < end;
       i += blockDim.x * gridDim.x)
    f(i);
  runner(args...);
}

template <typename T, typename LoopBody, class... Args>
__global__ void gmforallkernel(T start, T end, LoopBody body, Args... args) {
  runner(start, end, body, args...);
}

// Generalized mforall

template <int N, typename T, typename LoopBody, class... Args>
void gmforall(T start, T end, LoopBody &&body, Args... args) {
  // std::cout<<"Start is "<<start<<"\n";
  int blocks = 80 * 2048 / N;  // WARNING HARDWIRED FOR V100
  // multiforallkernel<<<blocks,N>>>(start,end, body, args...);
  gmforallkernel<<<blocks, N>>>(start, end, body, args...);
}

// Generalized gmforall3async

class MRange {
 public:
  int i, j, k;
};

template <typename T, typename LoopBody>
__device__ inline void runner3(T start, T end, LoopBody f) {
  for (decltype(start.i) i = start.i + threadIdx.x + blockIdx.x * blockDim.x;
       i < end.i; i += blockDim.x * gridDim.x)
    for (decltype(start.j) j = start.j + threadIdx.y + blockIdx.y * blockDim.y;
         j < end.j; j += blockDim.y * gridDim.y)
      for (decltype(start.k) k =
               start.k + threadIdx.z + blockIdx.z * blockDim.z;
           k < end.k; k += blockDim.z * gridDim.z)
        f(i, j, k);
}
template <typename T, typename LoopBody, class... Args>
__device__ inline void runner3(T start, T end, LoopBody f, Args... args) {
  for (decltype(start.i) i = start.i + threadIdx.x + blockIdx.x * blockDim.x;
       i < end.i; i += blockDim.x * gridDim.x)
    for (decltype(start.j) j = start.j + threadIdx.y + blockIdx.y * blockDim.y;
         j < end.j; j += blockDim.y * gridDim.y)
      for (decltype(start.k) k =
               start.k + threadIdx.z + blockIdx.z * blockDim.z;
           k < end.k; k += blockDim.z * gridDim.z)
        f(i, j, k);
  runner3(args...);
}

template <typename T, typename LoopBody, class... Args>
__global__ void gmforallkernel3(T start, T end, LoopBody body, Args... args) {
  runner3(start, end, body, args...);
}

template <int I, int J, int K, typename T, typename LoopBody, class... Args>
void gmforall3async(T &start, T &end, LoopBody &&body, Args... args) {
  //  int blocks=80 * 2048/N; // WARNING HARDWIRED FOR V100
  dim3 tpb(I, J, K);
  dim3 blocks(64 / I, 64 / J, 40 / K);  // WARNING HARDWIRED FOR V100 PBUGS

  gmforallkernel3<<<blocks, tpb>>>(start, end, body, args...);
}


// The Split Fuse ( SF) loop functions
#ifdef SW4_USE_SFK
template <int N, typename Tag, typename... Func>
__global__ void forall3kernelSF(Tag t, const int start0, const int N0,
                              const int start1, const int N1, const int start2,
                              const int N2, Func... f) {

const int STORE=5;
// NOTE: Shared memory is slightly slower, probably due to cache reduction
//#define USE_SHARED_MEMORY 1
#ifdef USE_SHARED_MEMORY
 auto off = (threadIdx.x+blockDim.x*threadIdx.y+blockDim.x*blockDim.y*threadIdx.z)*STORE;
  __shared__ double sma[512*STORE];
  double *carray = sma+off;
#else
  double carray[STORE];

#endif


  //printf("%d \n",off);
  int tid0 = start0 + threadIdx.x + blockIdx.x * blockDim.x;
  int tid1 = start1 + threadIdx.y + blockIdx.y * blockDim.y;
  int tid2 = start2 + threadIdx.z + blockIdx.z * blockDim.z;
  if ((tid0 < N0) && (tid1 < N1) && (tid2 < N2)) {
    (f(t, carray,tid0, tid1, tid2),...);
  }
}

template <int N, typename Tag, typename T1, typename T2, typename T3,
  typename... LoopBodies>
  void forall3asyncSF(Tag &t, T1 &irange, T2 &jrange, T3 &krange, LoopBodies &&... bodies) {
  if (irange.invalid || jrange.invalid || krange.invalid) return;
  dim3 tpb(irange.tpb, jrange.tpb, krange.tpb);
  dim3 blocks(irange.blocks, jrange.blocks, krange.blocks);


  // cudaEvent_t start, stop1,stop2,stop3;
  // cudaEventCreate(&start);
  // cudaEventCreate(&stop1);
  // cudaEventCreate(&stop2);
  // cudaEventCreate(&stop3);
  // bool TimeKernels=false;
  // if (TimeKernels){
  // insertEvent(start);
  // // Launch each body i a separate kernel. Gives incorrect results due to shared memory not persisting between kernels
  // (forall3kernelSM<N><<<blocks, tpb>>>(t, irange.start, irange.end, jrange.start,
  // 				     jrange.end, krange.start, krange.end, bodies),...);
  // }
  // insertEvent(stop1);
  // Launch all bodies inside 1 kernel. Regs use is 6 more than biggest kernels
  forall3kernelSF<N><<<blocks, tpb>>>(t, irange.start, irange.end, jrange.start,
				     jrange.end, krange.start, krange.end, bodies...);
  // insertEvent(stop2);

  // float ms;
  // if (TimeKernels) {
  //   ms = timeEvent(start,stop1);
  //   std::cout << "Kernel (Multiple)"<<t.value<<" runtime " << ms << " ms Best = "<<t.best<<" factor = "<<int(round(ms/t.best))<<" \n";
  // }
  // ms = timeEvent(stop1,stop2);
  // std::cout << "Kernel (Single)"<<t.value<<" runtime " << ms << " ms Best = "<<t.best<<" factor = "<<int(round(ms/t.best))<<" \n";

}
#endif // #ifdef SW4_USE_SFK

#endif  // Guards
