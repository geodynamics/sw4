#ifndef __SW4_FORALL_H__
#define __SW4_FORALL_H__

#if defined(ENABLE_CUDA)
#include "cuda_foralls.h"
#endif

#if defined(ENABLE_HIP)
#include "hip_foralls.h"
#endif

#endif
