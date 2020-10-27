#include "RAJA/RAJA.hpp"
#ifndef __SW4_POLICIES_H__
#define __SW4_POLICIES_H__

#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
#define ENABLE_GPU 1
#endif

#if defined(ENABLE_CUDA) && defined(ENABLE_HIP)
#define ENABLE_GPU_ERROR 1
#endif

#ifdef ENABLE_CUDA
#if RAJA_VERSION_MINOR>11
#include "cuda_policies12.h"
#else
#include "cuda_policies.h"
#endif
#endif

#ifdef ENABLE_HIP
#include "hip_policies.h"
#endif

#ifndef ENABLE_GPU
#include "omp_policies.h"
#endif

#endif  // Guards __SW4_POLICIES_H__
