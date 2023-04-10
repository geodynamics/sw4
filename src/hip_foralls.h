#ifndef __HIP_FORALLS_H__
#define __HIP_FORALLS_H__

#define HSYNC_DEVICE SW4_CheckDeviceError(hipDeviceSynchronize())
#define HSYNC_STREAM SW4_CheckDeviceError(hipStreamSynchronize(0))

#include "Mspace.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

#include "mpi.h"
#include "policies.h"
std::vector<sw4_type> factors(sw4_type N);
std::vector<sw4_type> factors(sw4_type N, sw4_type start);
#ifdef ENABLE_HIP
template <typename Func>
__global__ void forallkernel(sw4_type start, sw4_type N, Func f) {
  sw4_type tid = start + threadIdx.x + blockIdx.x * blockDim.x;
  if (tid < N) f(tid);
}
template <typename LoopBody>
void forall(sw4_type start, sw4_type end, LoopBody &&body) {
  sw4_type tpb = 32;
  sw4_type blocks = (end - start) / tpb;
  blocks = ((end - start) % tpb == 0) ? blocks : blocks + 1;
  // printf("Launching the kernel blocks= %d tpb= %d \n",blocks,tpb);
  // forallkernel<<<blocks, tpb>>>(start, end, body);
  hipLaunchKernelGGL(forallkernel<LoopBody>, blocks, tpb, 0, 0, start, end,
                     body);
  HSYNC_STREAM;
}
template <typename LoopBody>
void forallasync(sw4_type start, sw4_type end, LoopBody &&body) {
  sw4_type tpb = 1024;
  sw4_type blocks = (end - start) / tpb;
  blocks = ((end - start) % tpb == 0) ? blocks : blocks + 1;
  // printf("Launching the kernel blocks= %d tpb= %d \n",blocks,tpb);
  // forallkernel<<<blocks, tpb>>>(start, end, body);
  hipLaunchKernelGGL(forallkernel<LoopBody>, blocks, tpb, 0, 0, start, end,
                     body);
}

template <sw4_type N, typename LoopBody>
void forall(sw4_type start, sw4_type end, LoopBody &&body) {
  sw4_type tpb = 1024;
  sw4_type blocks = (end - start) / tpb;
  blocks = ((end - start) % tpb == 0) ? blocks : blocks + 1;
  // printf("Launching the kernel blocks= %d tpb= %d on line %d\n", blocks, tpb,
  //    N);
  // forallkernel<<<blocks, tpb>>>(start, end, body);
  hipLaunchKernelGGL(forallkernel<LoopBody>, blocks, tpb, 0, 0, start, end,
                     body);
  HSYNC_STREAM;
}

template <typename Func>
__global__ void forallkernelB(sw4_type start, sw4_type N, Func f) {
  sw4_type tid = start + threadIdx.x + blockIdx.x * blockDim.x;
  const sw4_type B = 8;
  for (sw4_type i = tid; i < N; i += B * blockDim.x * gridDim.x) {
    f(i);
    // sw4_type ii=i+blockDim.x * gridDim.x;
#pragma unroll B - 1
    for (sw4_type ii = 1; ii < B; ii++) {
      sw4_type iii = i + ii * blockDim.x * gridDim.x;
      if (iii < N) f(iii);
    }
  }
}
template <typename LoopBody>
void forallB(sw4_type start, sw4_type end, LoopBody &&body) {
  sw4_type tpb = 1024;
  sw4_type blocks = 52;
  // blocks=((end-start)%tpb==0)?blocks:blocks+1;
  // printf("Launching the kernel\n");
  // forallkernelB<<<blocks, tpb>>>(start, end, body);
  hipLaunchKernelGGL(forallkernelB<LoopBody>, blocks, tpb, 0, 0, start, end,
                     body);
  HSYNC_STREAM;
}

template <sw4_type N>
class Range {
 public:
  Range(sw4_type istart, sw4_type iend) : start(istart), end(iend), tpb(N) {
    blocks = (end - start) / N;
    blocks = ((end - start) % N == 0) ? blocks : blocks + 1;
    invalid = false;
    if (blocks <= 0) invalid = true;
  };
  static const sw4_type value = N;
  sw4_type start;
  sw4_type end;
  sw4_type blocks;
  sw4_type tpb;
  bool invalid;
};

template <sw4_type N, sw4_type M>
class RangeGS {
 public:
  RangeGS(sw4_type istart, sw4_type iend) : start(istart), end(iend), tpb(N), blocks(M){};
  sw4_type start;
  sw4_type end;
  sw4_type blocks;
  sw4_type tpb;
};

template <sw4_type WGS, sw4_type OCC, typename Func>
__launch_bounds__(WGS, OCC) __global__
    void forall3kernel(const sw4_type start0, const sw4_type N0, const sw4_type start1,
                       const sw4_type N1, const sw4_type start2, const sw4_type N2, Func f) {
  sw4_type tid0 = start0 + threadIdx.x + blockIdx.x * blockDim.x;
  sw4_type tid1 = start1 + threadIdx.y + blockIdx.y * blockDim.y;
  sw4_type tid2 = start2 + threadIdx.z + blockIdx.z * blockDim.z;
  if ((tid0 < N0) && (tid1 < N1) && (tid2 < N2)) f(tid0, tid1, tid2);
}
template <typename Func>
__launch_bounds__(256, 2) __global__
    void forall3kernel(const sw4_type start0, const sw4_type N0, const sw4_type start1,
                       const sw4_type N1, const sw4_type start2, const sw4_type N2, Func f) {
  sw4_type tid0 = start0 + threadIdx.x + blockIdx.x * blockDim.x;
  sw4_type tid1 = start1 + threadIdx.y + blockIdx.y * blockDim.y;
  sw4_type tid2 = start2 + threadIdx.z + blockIdx.z * blockDim.z;
  if ((tid0 < N0) && (tid1 < N1) && (tid2 < N2)) f(tid0, tid1, tid2);
}

template <sw4_type N, typename Func>
__launch_bounds__(256, 2) __global__
    void forall3kernel(const sw4_type start0, const sw4_type N0, const sw4_type start1,
                       const sw4_type N1, const sw4_type start2, const sw4_type N2, Func f) {
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

  //  std::cout << " BLOCKS " << block0 << " " << block1 << " " << block2 <<
  //  "\n";
  dim3 tpb(tpb0, tpb1, tpb2);
  dim3 blocks(block0, block1, block2);

  // printf("Launching the kernel 3d \n");
  // forall3kernel<<<blocks, tpb>>>(start0, end0, start1, end1, start2, end2,
  // body);
#ifdef SW4_CHECK_LAUNCH_BOUNDS
#endif
  hipLaunchKernelGGL((forall3kernel<1024, 2>), blocks, tpb, 0, 0, start0, end0,
                     start1, end1, start2, end2, body);
  HSYNC_STREAM;
}

template <template <sw4_type> class T, sw4_type N>
constexpr sw4_type extractN(const T<N> &) {
  return N;
}

template <typename T1, typename T2, typename T3, typename LoopBody>
void forall3async(T1 &irange, T2 &jrange, T3 &krange, LoopBody &&body) {
  dim3 tpb(irange.tpb, jrange.tpb, krange.tpb);
  dim3 blocks(irange.blocks, jrange.blocks, krange.blocks);

  // forall3kernel<<<blocks, tpb>>>(irange.start, irange.end, jrange.start,
  //                              jrange.end, krange.start, krange.end, body);
#ifdef SW4_CHECK_LAUNCH_BOUNDS
  // if ((irange.tpb*jrange.tpb*krange.tpb)!=256){
  //   std::cerr<<"Check launch bounds failed !!
  //   "<<(irange.tpb*jrange.tpb*krange.tpb)<<"!=256\n"; abort();
  // }
#endif
  hipLaunchKernelGGL((forall3kernel<T1::value * T2::value * T3::value, 2>),
                     blocks, tpb, 0, 0, irange.start, irange.end, jrange.start,
                     jrange.end, krange.start, krange.end, body);
}
template <typename T1, typename T2, typename T3, typename LoopBody>
void forall3(T1 &irange, T2 &jrange, T3 &krange, LoopBody &&body) {
  // dim3 tpb(irange.tpb,jrange.tpb,krange.tpb);
  // dim3 blocks(irange.blocks,jrange.blocks,krange.blocks);

  // forall3kernel<<<blocks,tpb>>>(irange.start,irange.end,jrange.start,jrange.end,krange.start,krange.end,body);
  forall3async(irange, jrange, krange, body);
  HSYNC_STREAM;
}

template <typename T1, typename T2, typename T3, typename LoopBody>
void forall3GSasync(T1 &irange, T2 &jrange, T3 &krange, LoopBody &&body) {
  dim3 tpb(irange.tpb, jrange.tpb, krange.tpb);
  dim3 blocks(irange.blocks, jrange.blocks, krange.blocks);

  hipLaunchKernelGGL(forall3gskernel, blocks, tpb, 0, 0, irange.start,
                     irange.end, jrange.start, jrange.end, krange.start,
                     krange.end, body);
}

template <typename T1, typename T2, typename T3, typename LoopBody>
void forall3GS(T1 &irange, T2 &jrange, T3 &krange, LoopBody &&body) {
  // dim3 tpb(irange.tpb,jrange.tpb,krange.tpb);
  // dim3 blocks(irange.blocks,jrange.blocks,krange.blocks);

  // forall3gskernel<<<blocks,tpb>>>(irange.start,irange.end,jrange.start,jrange.end,krange.start,krange.end,body);
  forall3async(irange, jrange, krange, body);
  HSYNC_STREAM;
}

// Forall2

template <sw4_type WGS, sw4_type OCC, typename Func>
__launch_bounds__(WGS, OCC) __global__
    void forall2kernel(const sw4_type start0, const sw4_type N0, const sw4_type start1,
                       const sw4_type N1, Func f) {
  sw4_type tid0 = start0 + threadIdx.x + blockIdx.x * blockDim.x;
  sw4_type tid1 = start1 + threadIdx.y + blockIdx.y * blockDim.y;
  if ((tid0 < N0) && (tid1 < N1)) f(tid0, tid1);
}

template <typename Func>
__global__ void forall2kernel(const sw4_type start0, const sw4_type N0, const sw4_type start1,
                              const sw4_type N1, Func f) {
  sw4_type tid0 = start0 + threadIdx.x + blockIdx.x * blockDim.x;
  sw4_type tid1 = start1 + threadIdx.y + blockIdx.y * blockDim.y;
  if ((tid0 < N0) && (tid1 < N1)) f(tid0, tid1);
}

template <typename T1, typename T2, typename LoopBody>
void forall2async(T1 &irange, T2 &jrange, LoopBody &&body) {
  dim3 tpb(irange.tpb, jrange.tpb, 1);
  dim3 blocks(irange.blocks, jrange.blocks, 1);

  // forall2kernel<<<blocks, tpb>>>(irange.start, irange.end, jrange.start,
  // jrange.end, body);
  hipLaunchKernelGGL((forall2kernel<T1::value * T2::value, 2>), blocks, tpb, 0,
                     0, irange.start, irange.end, jrange.start, jrange.end,
                     body);
}

template <typename T1, typename T2, typename LoopBody>
void forall2(T1 &irange, T2 &jrange, LoopBody &&body) {
  /* dim3 tpb(irange.tpb,jrange.tpb,1); */
  /* dim3 blocks(irange.blocks,jrange.blocks,1); */

  /* forall2kernel<<<blocks,tpb>>>(irange.start,irange.end,jrange.start,jrange.end,body);
   */
  forall2async(irange, jrange, body);
  HSYNC_STREAM;
}

// AutoTuning Code
template <sw4_type N, sw4_type M, sw4_type L>
class RangeAT {
 public:
  RangeAT(sw4_type istart, sw4_type iend, sw4_type jstart, sw4_type jend, sw4_type kstart, sw4_type kend)
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
    if (msw4_typeime.find(std::make_tuple((iend - istart), (jend - jstart),
                                     (kend - kstart))) == msw4_typeime.end())
      msw4_typeime[std::make_tuple((iend - istart), (jend - jstart),
                              (kend - kstart))] = 1.0e89;
    if (vloca.find(std::make_tuple((iend - istart), (jend - jstart),
                                   (kend - kstart))) == vloca.end())
      vloca[std::make_tuple((iend - istart), (jend - jstart),
                            (kend - kstart))] = 0;
    if (map.find(std::make_tuple((iend - istart), (jend - jstart),
                                 (kend - kstart))) == map.end()) {
      sw4_type max_block_size = floor(65536.0 / (N + 2)) - 1;
      max_block_size = N;
      // max_block_size = floor((256*255.0)/(N+1));
      files[curr] << "#Largest allowable block size is " << max_block_size
                  << " for Regs Per thread of " << N << " on line " << M
                  << " \n";
      files[curr] << "#LOOP SIZES " << (iend - istart) << " " << (jend - jstart)
                  << " " << (kend - kstart) << "\n";
      sw4_type ii, jj, kk;
      std::vector<std::tuple<sw4_type, sw4_type, sw4_type, sw4_type, sw4_type, sw4_type>> mins;
      if (0) {
        // Insert arbitrary configuration sw4_typeo vector of launch configs
        sw4_type ii = ceil(float(iend - istart) / 16);
        sw4_type jj = ceil(float(jend - jstart) / 4);
        sw4_type kk = ceil(float(kend - kstart) / 6);
        confs[curr].clear();
        dim3 tpb(16, 4, 6);
        dim3 blks(ii, jj, kk);
        confs[curr].push_back(std::make_tuple(tpb, blks));
      }
      const sw4_type offset = 0;  // Controlling the search space
      for (sw4_type block_size = max_block_size - offset;
           block_size <= max_block_size; block_size++) {
        auto f = factors(block_size, 32);
        for (const sw4_type &i : f) {
          sw4_type rest = block_size / i;
          auto ff = factors(rest);
          for (const sw4_type &j : ff) {
            sw4_type k = rest / j;
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
  static std::map<std::tuple<sw4_type, sw4_type, sw4_type>, std::tuple<dim3, dim3>> map;
  static std::map<std::tuple<sw4_type, sw4_type, sw4_type>,
                  std::vector<std::tuple<dim3, dim3>>>
      confs;
  static std::map<std::tuple<sw4_type, sw4_type, sw4_type>, float> msw4_typeime;
  static std::map<std::tuple<sw4_type, sw4_type, sw4_type>, sw4_type> vloca;
  static std::map<std::tuple<sw4_type, sw4_type, sw4_type>, std::ofstream> files;
  static dim3 msw4_type, minb;
  sw4_type istart, jstart, kstart;
  sw4_type iend, jend, kend;
  sw4_type ib, jb, kb;
  sw4_type ic, jc, kc;
  sw4_type blocks;
  sw4_type tpb;
};

template <sw4_type N, typename T1, typename LoopBody>
void forall3asyncAT(T1 &range3, LoopBody &&body) {
  dim3 tpb(range3.ic, range3.jc, range3.kc);
  dim3 blocks(range3.ib, range3.jb, range3.kb);

  hipEvent_t start, stop;
  hipEventCreate(&start);
  hipEventCreate(&stop);
  auto curr = std::make_tuple((range3.iend - range3.istart),
                              (range3.jend - range3.jstart),
                              (range3.kend - range3.kstart));
#define SW4_VERBOSE
#ifdef SW4_VERBOSE
  // Calculate max block size and compare it to early calculations
  sw4_type minGridSize;
  sw4_type blockSize;

  if (range3.vloca[curr] == 0) {
    SW4_CheckDeviceError(hipOccupancyMaxPotentialBlockSize(
        &minGridSize, &blockSize, forall3kernel<N, LoopBody>, 0, 0));
    range3.files[curr] << "#LINE " << N << " MAX BLOCK SIZE FOR OCC CALC "
                       << blockSize << " minGridSize " << minGridSize << " \n";
    sw4_type numBlocks;
    SW4_CheckDeviceError(hipOccupancyMaxActiveBlocksPerMultiprocessor(
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

    sw4_type n = std::get<1>(i).x * std::get<1>(i).y * std::get<1>(i).z;

    hipEventRecord(start);
    hipLaunchKernelGGL(forall3kernel<N>, std::get<1>(i), std::get<0>(i), 0, 0,
                       range3.istart, range3.iend, range3.jstart, range3.jend,
                       range3.kstart, range3.kend, body);
    hipEventRecord(stop);
    hipEventSynchronize(stop);
    float ms = 0;
    hipEventElapsedTime(&ms, start, stop);
    if (hipGetLastError() == hipSuccess) {
      range3.files[curr] << n << " " << ms << " " << std::get<0>(i).x << " "
                         << std::get<0>(i).y << " " << std::get<0>(i).z << " "
                         << range3.vloca[curr] << "\n";
      if (range3.msw4_typeime[curr] > ms) {
        range3.msw4_typeime[curr] = ms;
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
      auto &msw4_type = std::get<0>(range3.map[curr]);
      range3.files[curr] << "&\n#BEST CONFIG " << range3.msw4_typeime[curr] << " "
                         << msw4_type.x << " " << msw4_type.y << " " << msw4_type.z << " "
                         << range3.confs[curr].size() << "\n";
      auto &minb = std::get<1>(range3.map[curr]);
      sw4_type minblocks = minb.x * minb.y * minb.z;

      range3.files[curr] << minblocks << " " << range3.msw4_typeime[curr] << " "
                         << msw4_type.x << " " << msw4_type.y << " " << msw4_type.z << " "
                         << range3.confs[curr].size() << "\n";
    }
  } else {
    auto &v = range3.map[std::make_tuple((range3.iend - range3.istart),
                                         (range3.jend - range3.jstart),
                                         (range3.kend - range3.kstart))];
    hipLaunchKernelGGL(forall3kernel<N>, std::get<1>(v), std::get<0>(v), 0, 0,
                       range3.istart, range3.iend, range3.jstart, range3.jend,
                       range3.kstart, range3.kend, body);
  }
}

#endif

template <sw4_type N>
class Tclass {};

template <sw4_type N, typename Tag, typename Func>
__launch_bounds__(256, 2) __global__
    void forall3kernel(Tag t, const sw4_type start0, const sw4_type N0, const sw4_type start1,
                       const sw4_type N1, const sw4_type start2, const sw4_type N2, Func f) {
  sw4_type tid0 = start0 + threadIdx.x + blockIdx.x * blockDim.x;
  sw4_type tid1 = start1 + threadIdx.y + blockIdx.y * blockDim.y;
  sw4_type tid2 = start2 + threadIdx.z + blockIdx.z * blockDim.z;
  if ((tid0 < N0) && (tid1 < N1) && (tid2 < N2)) f(t, tid0, tid1, tid2);
}
template <sw4_type WGS, sw4_type OCC, typename Tag, typename Func>
__launch_bounds__(WGS, OCC) __global__
    void forall3kernelV(Tag t, const sw4_type start0, const sw4_type N0, const sw4_type start1,
                        const sw4_type N1, const sw4_type start2, const sw4_type N2, Func f) {
  sw4_type tid0 = start0 + threadIdx.x + blockIdx.x * blockDim.x;
  sw4_type tid1 = start1 + threadIdx.y + blockIdx.y * blockDim.y;
  sw4_type tid2 = start2 + threadIdx.z + blockIdx.z * blockDim.z;
  if ((tid0 < N0) && (tid1 < N1) && (tid2 < N2)) f(t, tid0, tid1, tid2);
}

template <sw4_type WGS, sw4_type OCC, typename Func>
__launch_bounds__(WGS, OCC) __global__
    void forall3kernelV2(const sw4_type start0, const sw4_type N0, const sw4_type start1,
                         const sw4_type N1, const sw4_type start2, const sw4_type N2, Func f) {
  sw4_type tid0 = start0 + threadIdx.x + blockIdx.x * blockDim.x;
  sw4_type tid1 = start1 + threadIdx.y + blockIdx.y * blockDim.y;
  sw4_type tid2 = start2 + threadIdx.z + blockIdx.z * blockDim.z;
  if ((tid0 < N0) && (tid1 < N1) && (tid2 < N2)) f(tid0, tid1, tid2);
}
template <sw4_type N, typename Tag, typename T1, typename T2, typename T3,
          typename LoopBody>
void forall3async(Tag &t, T1 &irange, T2 &jrange, T3 &krange, LoopBody &&body) {
  if (irange.invalid || jrange.invalid || krange.invalid) return;
  dim3 tpb(irange.tpb, jrange.tpb, krange.tpb);
  dim3 blocks(irange.blocks, jrange.blocks, krange.blocks);
  // std::cout<<"forall launch tpb"<<irange.tpb<<" "<<jrange.tpb<<"
  // "<<krange.tpb<<"\n"; std::cout<<"forall launch blocks"<<irange.blocks<<"
  // "<<jrange.blocks<<" "<<krange.blocks<<"\n";
  // forall3kernel<N><<<blocks, tpb>>>(t, irange.start, irange.end,
  // jrange.start,
  //                                   jrange.end, krange.start, krange.end,
  //                                   body);
  // std::cout<<hipGetErrorString(hipPeekAtLastError())<<"\n"<<std::flush;
  // std::cout<<"Launching kernel..."<<std::flush;
#ifdef SW4_CHECK_LAUNCH_BOUNDS
  if ((irange.tpb * jrange.tpb * krange.tpb) != 256) {
    std::cerr << "Check launch bounds failed !! "
              << (irange.tpb * jrange.tpb * krange.tpb) << "!=256\n";
    abort();
  }
#endif
  hipLaunchKernelGGL(forall3kernel<N>, blocks, tpb, 0, 0, t, irange.start,
                     irange.end, jrange.start, jrange.end, krange.start,
                     krange.end, body);
  // hipStreamSynchronize(0);
  ////std::cout<<"Done\n"<<std::flush;
}
template <sw4_type N, typename Tag, typename T1, typename T2, typename T3,
          typename LoopBody>
void forall3(Tag &t, T1 &irange, T2 &jrange, T3 &krange, LoopBody &&body) {
forall3async<N, Tag>(t, irange, jrange, krange, body);
  HSYNC_STREAM;
}

template <sw4_type WGS, sw4_type OCC, typename Tag, typename T1, typename T2, typename T3,
          typename LoopBody>
void forall3asyncV(Tag &t, T1 &irange, T2 &jrange, T3 &krange,
                   LoopBody &&body) {
  if (irange.invalid || jrange.invalid || krange.invalid) return;
  dim3 tpb(irange.tpb, jrange.tpb, krange.tpb);
  dim3 blocks(irange.blocks, jrange.blocks, krange.blocks);
  // std::cout<<"forall launch tpb"<<irange.tpb<<" "<<jrange.tpb<<"
  // "<<krange.tpb<<"\n"; std::cout<<"forall launch blocks"<<irange.blocks<<"
  // "<<jrange.blocks<<" "<<krange.blocks<<"\n";
  // forall3kernel<N><<<blocks, tpb>>>(t, irange.start, irange.end,
  // jrange.start,
  //                                   jrange.end, krange.start, krange.end,
  //                                   body);
  // std::cout<<hipGetErrorString(hipPeekAtLastError())<<"\n"<<std::flush;
  // std::cout<<"Launching kernel..."<<std::flush;
#ifdef SW4_CHECK_LAUNCH_BOUNDS
  if ((irange.tpb * jrange.tpb * krange.tpb) != 256) {
    std::cerr << "Check launch bounds failed !! "
              << (irange.tpb * jrange.tpb * krange.tpb) << "!=256\n";
    abort();
  }
#endif
  hipLaunchKernelGGL((forall3kernelV<WGS, OCC>), blocks, tpb, 0, 0, t,
                     irange.start, irange.end, jrange.start, jrange.end,
                     krange.start, krange.end, body);
  // hipStreamSynchronize(0);
  ////std::cout<<"Done\n"<<std::flush;
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
template <sw4_type N, typename T, typename LB0, typename LB1, typename LB2,
          typename LB3>
void multiforall(T start0, T end0, LB0 &&body0, T start1, T end1, LB1 &&body1,
                 T start2, T end2, LB2 &&body2, T start3, T end3, LB3 &&body3) {
  // sw4_type blocks = 80 * 2048 / N;  // WARNING HARDWIRED FOR V100

  dim3 tpb(N, 1, 1);
  dim3 blocks(80 * 2048 / N, 1, 1);
  hipLaunchKernelGGL(multiforallkernel, blocks, tpb, 0, 0, start0, end0, body0,
                     start1, end1, body1, start2, end2, body2, start3, end3,
                     body3);
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

template <sw4_type N, typename T, typename LoopBody, class... Args>
void gmforall(T start, T end, LoopBody &&body, Args... args) {
  // std::cout<<"Start is "<<start<<"\n";
  dim3 tpb(N, 1, 1);
  dim3 blocks(80 * 2048 / N, 1, 1);
  // sw4_type blocks = 80 * 2048 / N;  // WARNING HARDWIRED FOR V100
  // multiforallkernel<<<blocks,N>>>(start,end, body, args...);
  // gmforallkernel<<<blocks, N>>>(start, end, body, args...);
  hipLaunchKernelGGL(gmforallkernel, blocks, tpb, 0, 0, start, end, body,
                     args...);
}

// Generalized gmforall3async

class MRange {
 public:
  sw4_type i, j, k;
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

template <sw4_type I, sw4_type J, sw4_type K, typename T, typename LoopBody, class... Args>
void gmforall3async(T &start, T &end, LoopBody &&body, Args... args) {
  //  sw4_type blocks=80 * 2048/N; // WARNING HARDWIRED FOR V100
  dim3 tpb(I, J, K);
  dim3 blocks(64 / I, 64 / J, 40 / K);  // WARNING HARDWIRED FOR V100 PBUGS

  gmforallkernel3<<<blocks, tpb>>>(start, end, body, args...);
}

// Split Fusion kernels
template <sw4_type N, typename Tag, typename... Func>
__launch_bounds__(256, 1) __global__
    void forall3kernelSF(Tag t, const sw4_type start0, const sw4_type N0,
                         const sw4_type start1, const sw4_type N1, const sw4_type start2,
                         const sw4_type N2, Func... f) {
  const sw4_type STORE = 5;
#define USE_SHARED_MEMORY 1
#ifdef USE_SHARED_MEMORY
  sw4_type off = (threadIdx.x + blockDim.x * threadIdx.y +
             blockDim.x * blockDim.y * threadIdx.z) *
            STORE;
  __shared__ double sma[512 * STORE];

  double *carray = sma + off;
#else
  double carray[STORE];

#endif
  // double sma2[3];

  sw4_type tid0 = start0 + threadIdx.x + blockIdx.x * blockDim.x;
  sw4_type tid1 = start1 + threadIdx.y + blockIdx.y * blockDim.y;
  sw4_type tid2 = start2 + threadIdx.z + blockIdx.z * blockDim.z;
  if ((tid0 < N0) && (tid1 < N1) && (tid2 < N2)) {
    (f(t, carray, tid0, tid1, tid2), ...);
  }
}

template <sw4_type N, typename Tag, typename T1, typename T2, typename T3,
          typename... LoopBodies>
void forall3asyncSF(Tag &t, T1 &irange, T2 &jrange, T3 &krange,
                    LoopBodies &&... bodies) {
  if (irange.invalid || jrange.invalid || krange.invalid) {
    std::cerr << "Invalid ranges in forall3asyncSF \n";
    return;
  }
  dim3 tpb(irange.tpb, jrange.tpb, krange.tpb);
  dim3 blocks(irange.blocks, jrange.blocks, krange.blocks);

  hipLaunchKernelGGL(forall3kernelSF<N>, blocks, tpb, 0, 0, t, irange.start,
                     irange.end, jrange.start, jrange.end, krange.start,
                     krange.end, bodies...);
}

#endif  // Guards
