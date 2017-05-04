#ifdef CUDA_CODE
typedef NestedPolicy<ExecList<cuda_threadblock_x_exec<16>,cuda_threadblock_y_exec<4>,
  cuda_threadblock_z_exec<4>>>
  EXEC_CARTBC;

#define REDUCE_BLOCK_SIZE 1024
typedef RAJA::cuda_exec<REDUCE_BLOCK_SIZE> EXEC;

typedef cuda_reduce<REDUCE_BLOCK_SIZE> REDUCE_POLICY;
#define SYNC_DEVICE cudaDeviceSynchronize();
#else
typedef NestedPolicy<ExecList<omp_parallel_for_exec,omp_parallel_for_exec,
  omp_parallel_for_exec>>
  EXEC_CARTBC;

typedef RAJA::omp_parallel_for_exec EXEC;
typedef omp_reduce REDUCE_POLICY;
#define  SYNC_DEVICE
#endif

