#include <mpi.h>
#if defined(ENABLE_CUDA) && defined(OPEN_MPI)
#include <mpi-ext.h>
#endif

#include <limits.h>
#include <unistd.h>

#include <iomanip>
#include <unordered_map>

#include "GridPointSource.h"
#include "Mspace.h"
#include "caliper.h"
#include "policies.h"
#include "sys/sysinfo.h"
#include "sys/types.h"

long long node_mem();
std::string hostname();
struct global_variable_holder_struct global_variables = {0, 0, 0, 0, 0, 0,
                                                         0, 1, 0, 0, 0, 0};
using namespace std;

int presetGPUID(int mpi_rank, int local_rank, int local_size) {
  int dev_counts_local[256] = {0};
  int dev_counts_global[256] = {0};
  int device = 0;
  // SW4_CheckDeviceError(cudaErrorStartupFailure); // For testing code in
  // CheckError
#if defined(ENABLE_GPU_ERROR)
  std::cerr
      << " Compilation error. Both ENABLE_CUDA and ENABLE_HIP are defined\n";
  abort();
#endif  // ENABLE_GPU_ERROR
#if defined(ENABLE_GPU)

#ifdef USE_MAGMA
  if (magma_init() != MAGMA_SUCCESS) {
    std::cerr << "ERROR MAGMA INIT FAILED\n";
    abort();
  }
#endif  // USE_MAGMA

#ifdef ENABLE_CUDA
  int devices_per_node = 4;
  SW4_CheckDeviceError(cudaGetDeviceCount(&devices_per_node));
  global_variables.num_devices = devices_per_node;
  if (devices_per_node > 1) {
    // char *crank = getenv("SLURM_LOCALID");
    // int device = atoi(crank) % devices_per_node;
    if (local_size == devices_per_node) {
      device = local_rank;
    } else {
      device = local_rank % devices_per_node;
      std::cerr << "WARNING :: There are " << devices_per_node
                << " devices per node and " << local_size
                << " ranks per node\n";
    }
    global_variables.device = device;
    printf(" CUDA presetGPU Called ::  LOCAL RANK %d \n", device);
    SW4_CheckDeviceError(cudaSetDevice(device));
    SW4_CheckDeviceError(cudaFree(NULL));
    char uuid[80];
    if (nvmlInit() != NVML_SUCCESS) printf("NVML INIT CALL FAILED\n");
    nvmlDevice_t nvdev;
    if (nvmlDeviceGetHandleByIndex(device, &nvdev) != NVML_SUCCESS)
      printf("NVML GetHandleByIndex CALL FAILED\n");
    if (nvmlDeviceGetUUID(nvdev, uuid, 80) != NVML_SUCCESS)
      printf("UUID CALL FAILED\n");
    if (nvmlDeviceSetCpuAffinity(nvdev) != NVML_SUCCESS)
      printf("NVML SET CPU AFFINITY FAILED \n");
    else
      printf("NVML SET CPU AFFINITY CALLED SUCCESFULLY\n");
  }
  global_variables.device = device;

#endif  // ENDIF ENABLE_CUDA

#ifdef ENABLE_HIP
  int devices_per_node = 4;
  SW4_CheckDeviceError(hipGetDeviceCount(&devices_per_node));
  // printf("Number of devices is %d\n", devices_per_node);
  fflush(stdout);
  global_variables.num_devices = devices_per_node;
  if (devices_per_node > 1) {
    // char *crank = getenv("SLURM_LOCALID");
    // printf("Return from GETENV IS %s\n", crank);
    // fflush(stdout);
    // int device = atoi(crank) % devices_per_node;
    // device = mpi_rank % devices_per_node;
    if (local_size == devices_per_node) {
      device = local_rank;
    } else {
      device = local_rank % devices_per_node;
      std::cerr << "WARNING :: There are " << devices_per_node
                << " devices per node and " << local_size
                << " ranks per node\n";
    }
    global_variables.device = device;
    printf(" HIP presetGPU Called ::  LOCAL RANK %d \n", device);
    fflush(stdout);
    SW4_CheckDeviceError(hipSetDevice(device));
  }
#endif  // ENDIF ENABLE_HIP
  dev_counts_local[device] = 1;
  MPI_Reduce(dev_counts_local, dev_counts_global, devices_per_node, MPI_INT,
             MPI_SUM, 0, MPI_COMM_WORLD);
  if (!mpi_rank)
    for (int i = 0; i < devices_per_node; i++)
      if (dev_counts_global[i] != 0)
        std::cout << "Device " << i << " used by " << dev_counts_global[i]
                  << " ranks \n";
        // std::cerr << "Device set to " << global_variables.device << "\n";
#endif  // ENDIF ENABLE_GPU
  return device;
}

typedef struct {
  const char *file;
  int line;
  size_t size;
  Space type;
} pattr_t;

pattr_t *patpush(void *ptr, pattr_t *ss) {
  static std::unordered_map<void *, pattr_t *> map;
  if (ss != NULL) {
    map[ptr] = ss;
  } else {
    std::unordered_map<void *, pattr_t *>::const_iterator got = map.find(ptr);
    if (got == map.end()) {
      // std:cerr<<"ELEMENT NOT FOUND IN MAP\n";
      return NULL;
    } else
      return got->second;
  }
  return NULL;
}

void check_mem() {
  size_t mfree, mtotal;
#if defined(ENABLE_CUDA)
  SW4_CheckDeviceError(cudaMemGetInfo(&mfree, &mtotal));
#endif
#if defined(ENABLE_HIP)
  SW4_CheckDeviceError(hipMemGetInfo(&mfree, &mtotal));
#endif
  global_variables.gpu_memory_hwm =
      std::max((mtotal - mfree), global_variables.gpu_memory_hwm);
}

void print_hwm(int rank) {
  const int allocator_count = 3;
  long long host_mem_used = 0;
#ifdef SW4_CPU_HWM
  host_mem_used = node_mem();
#endif
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
  // std::cout<<"PRINT_HWM"<<std::flush;
  float hwm_local[allocator_count + 2], hwm_global[allocator_count + 2],
      hwm_global_min[allocator_count + 2];
#ifdef SW4_USE_UMPIRE
  hwm_local[0] = umpire::ResourceManager::getInstance()
                     .getAllocator("UM_pool")
                     .getHighWatermark() /
                 1024.0 / 1024.0 / 1024.0;
  hwm_local[1] = umpire::ResourceManager::getInstance()
                     .getAllocator("UM_pool_temps")
                     .getHighWatermark() /
                 1024.0 / 1024.0 / 1024.0;
  hwm_local[2] = umpire::ResourceManager::getInstance()
                     .getAllocator("UM_object_pool")
                     .getHighWatermark() /
                 1024.0 / 1024.0 / 1024.0;
  hwm_local[allocator_count] = 0;
  for (int i = 0; i < allocator_count; i++)
    hwm_local[allocator_count] += hwm_local[i];

  hwm_local[allocator_count + 1] =
      (float)host_mem_used / 1024.0 / 1024.0 / 1024.0;
  // std::cout<<getRank()<<" Umpire HWM
  // "<<umpire::ResourceManager::getInstance().getAllocator("UM_pool").getHighWatermark()/1024/1024<<"
  // MB\n";
#else
  hwm_local[0] = global_variables.gpu_memory_hwm / 1024.0 / 1024.0 / 1024.0;
  // std::cout<<getRank()<<" GPU Memory HWM =
  // "<<global_variables.gpu_memory_hwm/1024/1024<<" MB \n";
  // std::cout<<getRank()<<" GPU Memory Max =
  // "<<global_variables.max_mem/1024/1024<<" MB \n";
#endif
  float hwm_total;
  MPI_Reduce(hwm_local, hwm_global, allocator_count + 2, MPI_FLOAT, MPI_MAX, 0,
             MPI_COMM_WORLD);
  MPI_Reduce(hwm_local, hwm_global_min, allocator_count + 2, MPI_FLOAT, MPI_MIN,
             0, MPI_COMM_WORLD);
  MPI_Reduce(&hwm_local[allocator_count], &hwm_total, 1, MPI_FLOAT, MPI_SUM, 0,
             MPI_COMM_WORLD);
  if (!rank) {
    for (int i = 0; i < allocator_count; i++) {
      std::cerr << std::setprecision(2) << i << " MIN Global Device HWM is "
                << hwm_global_min[i] << " GB\n";
      std::cerr << i << " MAX Global Device HWM is " << hwm_global[i]
                << " GB\n";
    }

    std::cerr << std::setprecision(2) << " MIN Host HWM is "
              << hwm_global_min[allocator_count + 1] << " GB, MAX Host HWM is "
              << hwm_global[allocator_count + 1] << " GB\n";

    std::cout << std::fixed << std::setprecision(2)
              << " Min & Max per-device memory HWM  "
              << hwm_global_min[allocator_count] << ","
              << hwm_global[allocator_count] << " GB\n"
              << std::flush;

    std::cout << " Total device memory used " << hwm_total << " GB\n";
  }

  // umpire::util::StatisticsDatabase::getDatabase()->printStatistics(std::cout);
  if (Managed::hwm > 0) {
    std::cerr << "Space::Managed object count & HWM are " << Managed::ocount
              << " & " << Managed::hwm << " Size = "
              << Managed::hwm * sizeof(GridPointSource) / 1024.0 / 1024.0
              << " MB Each object = " << sizeof(GridPointSource) << "\n";
  }
  if (Managed::ocount != 0)
    std::cerr << "WARNING :: Managed object count should be zero at the end of "
                 "the simulation\n";
    // std::cout<<" ~HOST MEM MAX
    // "<<global_variables.host_mem_hwm/1024.0/1024.0<<" MB\n";
    // std::cout<<"PRINT_HWM DONE"<<std::flush;
#endif  // ENABLE_CUDA
}
void *operator new(std::size_t size, Space loc) throw() {
  SW4_MARK_FUNCTION;
#ifdef ENABLE_GPU
  if (loc == Space::Managed) {
    // std::cout<<"Space::Managed allocation \n";
    if (size == 0) size = 1;  // new has to return an valid pointer for 0 size.
    void *ptr;
#ifndef SW4_USE_UMPIRE
    if (SW4_MALLOC_MANAGED(&ptr, size) != SW4_DEVICE_SUCCESS) {
      std::cerr << "Managed memory allocation failed " << size << "\n";
      abort();
      // throw();
    } else {
      check_mem();
      global_variables.curr_mem += size;
      global_variables.max_mem =
          std::max(global_variables.max_mem, global_variables.curr_mem);
#if defined(ENABLE_CUDA)
      SW4_CheckDeviceError(cudaMemAdvise(ptr, size,
                                         cudaMemAdviseSetPreferredLocation,
                                         global_variables.device));
#endif
      if (ptr == NULL) {
        std::cerr << "NULL POINTER \n" << std::flush;
        abort();
      }
      return ptr;
    }
#else  // SW4_USE_UMPIRE
    umpire::ResourceManager &rma = umpire::ResourceManager::getInstance();
    auto allocator = rma.getAllocator("UM_pool");
    ptr = static_cast<void *>(allocator.allocate(size));
#if defined(ENABLE_CUDA)
    //    SW4_CheckDeviceError(cudaMemAdvise(
    //     ptr, size, cudaMemAdviseSetPreferredLocation,
    //     global_variables.device));
#endif
    // std::cout<<"PTR 1 "<<ptr<<"\n";
    // SW4_CheckDeviceError(cudaMemset(ptr,0,size));
    return ptr;
#endif  // SW4_USE_UMPIRE
  } else if (loc == Space::Host) {
    // std::cout<<"Calling my placement new \n";
    // global_variables.host_curr_mem+=size;
    // global_variables.host_max_mem=std::max(global_variables.host_max_mem,global_variables.host_curr_mem);
    return ::operator new(size);
  } else if (loc == Space::Device) {
    // std::cout<<"Managed allocation \n";
    if (size == 0) size = 1;  // new has to return an valid pointer for 0 size.
    void *ptr;
    if (SW4_MALLOC_DEVICE(&ptr, size) != SW4_DEVICE_SUCCESS) {
      std::cerr << "Device memory allocation failed " << size << "\n";
      abort();
      // throw();
    } else
      return ptr;
  } else if (loc == Space::Pinned) {
    if (size == 0) size = 1;  // new has to return an valid pointer for 0 size.
    void *ptr;
    SW4_CheckDeviceError(SW4_MALLOC_PINNED(&ptr, size));
    return ptr;
  } else if (loc == Space::Managed_temps) {
#ifdef SW4_USE_UMPIRE
    umpire::ResourceManager &rma = umpire::ResourceManager::getInstance();
    auto allocator = rma.getAllocator("UM_pool_temps");
    void *ptr = static_cast<void *>(allocator.allocate(size));
    // SW4_CheckDeviceError(cudaMemAdvise(ptr,size,cudaMemAdviseSetPreferredLocation,0));
    // std::cout<<"PTR 1 "<<ptr<<"\n";
    // SW4_CheckDeviceError(cudaMemset(ptr,0,size));
    return ptr;
#else
    // std::cerr << "Managed_temp location no defined\n";
    return ::operator new(size, Space::Managed);
#endif
  } else {
    std::cerr << "Unknown memory space for allocation request 0 " << as_int(loc)
              << "\n";
    // throw std::bad_alloc();
    abort();
  }
#else   // NOT ENABLE_GPU
  if ((loc == Space::Managed) || (loc == Space::Device) ||
      (loc == Space::Pinned)) {
    // std::cout<<"Managed location not available yet \n";
    return ::operator new(size);
  } else if (loc == Space::Host) {
    // std::cout<<"Calling my placement new \n";
    return ::operator new(size);
  } else {
    std::cerr << "Unknown memory space for allocation request 1\n";
    abort();
    // throw std::bad_alloc();
  }
#endif  // ENABE_GPU
}
void *operator new(std::size_t size, Space loc, char *file, int line) throw() {
  SW4_MARK_FUNCTION;
  std::cout << "Calling tracking new from " << line << " of " << file << "\n";
  pattr_t *ss = new pattr_t;
  ss->file = file;
  ss->line = line;
  ss->type = loc;
  ss->size = size;
  void *ret = ::operator new(size, loc);
  patpush(ret, ss);
  return ret;
}

void *operator new[](std::size_t size, Space loc) throw() {
  SW4_MARK_FUNCTION;
#ifdef ENABLE_GPU
  if (loc == Space::Managed) {
    // std::cout<<"Managed [] allocation \n";
    if (size == 0) size = 1;  // new has to return an valid pointer for 0 size.
    void *ptr;
#ifndef SW4_USE_UMPIRE
    if (SW4_MALLOC_MANAGED(&ptr, size) != SW4_DEVICE_SUCCESS) {
      std::cerr << "Managed memory allocation failed " << size << "\n";
      abort();
      // throw std::bad_alloc();
    } else {
      check_mem();
      global_variables.curr_mem += size;
      global_variables.max_mem =
          std::max(global_variables.max_mem, global_variables.curr_mem);
#if defined(ENABLE_CUDA)
      SW4_CheckDeviceError(cudaMemAdvise(ptr, size,
                                         cudaMemAdviseSetPreferredLocation,
                                         global_variables.device));
#endif
      return ptr;
    }
#else
    umpire::ResourceManager &rma = umpire::ResourceManager::getInstance();
    auto allocator = rma.getAllocator("UM_pool");
    ptr = static_cast<void *>(allocator.allocate(size));
#if defined(ENABLE_CUDA)
    SW4_CheckDeviceError(cudaMemAdvise(
        ptr, size, cudaMemAdviseSetPreferredLocation, global_variables.device));
#endif
    // std::cout<<"PTR 2 "<<ptr<<"\n";
    return ptr;
#endif
  } else if (loc == Space::Host) {
    // std::cout<<"Calling my placement new \n";
    // global_variables.host_curr_mem+=size;
    // global_variables.host_max_mem=std::max(global_variables.host_max_mem,global_variables.host_curr_mem);
    return ::operator new(size);
  } else if (loc == Space::Device) {
    // std::cout<<"Managed allocation \n";
    if (size == 0) size = 1;  // new has to return an valid pointer for 0 size.
    void *ptr;
    if (SW4_MALLOC_DEVICE(&ptr, size) != SW4_DEVICE_SUCCESS) {
      std::cerr << "Device memory allocation failed " << size << "\n";
      abort();
      // throw std::bad_alloc();
    } else
      return ptr;
  } else if (loc == Space::Pinned) {
    if (size == 0) size = 1;  // new has to return an valid pointer for 0 size.
    void *ptr;
    SW4_CheckDeviceError(SW4_MALLOC_PINNED(&ptr, size));
    return ptr;
  } else if (loc == Space::Managed_temps) {
#if defined(SW4_USE_UMPIRE)
    umpire::ResourceManager &rma = umpire::ResourceManager::getInstance();
    auto allocator = rma.getAllocator("UM_pool_temps");
    void *ptr = static_cast<void *>(allocator.allocate(size));
    // SW4_CheckDeviceError(cudaMemAdvise(ptr,size,cudaMemAdviseSetPreferredLocation,0));
    // std::cout<<"PTR 1 "<<ptr<<"\n";
    // SW4_CheckDeviceError(cudaMemset(ptr,0,size));
    return ptr;
#else
    // std::cerr << " Memory location Managed_temps is not defined\n";
    return ::operator new(size, Space::Managed);
#endif
  } else {
    // cudaHostAlloc(&ptr,size+sizeof(size_t)*MEM_PAD_LEN,cudaHostAllocMapped));
    std::cerr << "Unknown memory space for allocation request  2" << as_int(loc)
              << "\n";
    abort();
    // throw std::bad_alloc();
  }

#else   // !ENABLE_GPU
  if ((loc == Space::Managed) || (loc == Space::Device) ||
      (loc == Space::Pinned) || (loc == Space::Managed_temps)) {
    // std::cout<<"Managed location not available yet \n";
    return ::operator new(size);
  } else if (loc == Space::Host) {
    // std::cout<<"Calling my placement new \n";
    return ::operator new(size);
  } else {
    std::cerr << "Unknown memory space for allocation request 3" << as_int(loc)
              << "\n";
    throw std::bad_alloc();
  }
#endif  // ENABLE_GPU
}
void *operator new[](std::size_t size, Space loc, const char *file, int line) {
  SW4_MARK_FUNCTION;
  //  std::cout<<"Calling tracking new from "<<line<<" of "<<file<<"\n";
  pattr_t *ss = new pattr_t;
  ss->file = file;
  ss->line = line;
  ss->type = loc;
  ss->size = size;
  void *ret = ::operator new(size, loc);
  patpush(ret, ss);
  return ret;
}

void operator delete(void *ptr, Space loc) throw() {
#ifdef ENABLE_GPU
  if (loc == Space::Managed) {
    // std::cout<<"Managed delete\n";
#ifndef SW4_USE_UMPIRE
    pattr_t *ss = patpush(ptr, NULL);
    if (ss != NULL) {
      global_variables.curr_mem -= ss->size;
      // global_variables.max_mem=std::max(global_variables.max_mem,global_variables.curr_mem);
    }
    SW4_CheckDeviceError(SW4_FREE_MANAGED(ptr));
#else
    umpire::ResourceManager &rma = umpire::ResourceManager::getInstance();
    auto allocator = rma.getAllocator("UM_pool");
    allocator.deallocate(ptr);
#endif
  } else if (loc == Space::Device) {
    pattr_t *ss = patpush(ptr, NULL);
    if (ss != NULL) {
      global_variables.curr_mem -= ss->size;
      // global_variables.max_mem=std::max(global_variables.max_mem,global_variables.curr_mem);
    }
    SW4_CheckDeviceError(SW4_FREE_MANAGED(ptr));
  } else if (loc == Space::Pinned)
    SW4_CheckDeviceError(SW4_FREE_PINNED(ptr));
  else if (loc == Space::Host) {
    // std::cout<<"Calling my placement delete\n";
    ::operator delete(ptr);
  } else if (loc == Space::Managed_temps) {
#if defined(SW4_USE_UMPIRE)
    umpire::ResourceManager &rma = umpire::ResourceManager::getInstance();
    auto allocator = rma.getAllocator("UM_pool_temps");
    allocator.deallocate(ptr);
#else
    // std::cerr << "Memory location Managed_temps not defined\n";
    SW4_CheckDeviceError(SW4_FREE_MANAGED(ptr));
#endif
  } else {
    std::cerr << "Unknown memory space for de-allocation request 4 \n";
    abort();
  }
#else
  if ((loc == Space::Managed) || (loc == Space::Device)) {
    // std::cout<<"Managed delete not available yet \n";
    ::operator delete(ptr);
  } else if (loc == Space::Host) {
    // std:cout<<"Calling my placement delete\n";
    ::operator delete(ptr);
  } else {
    std::cerr << "Unknown memory space for de-allocation request "
              << as_int(loc) << "\n";
  }
#endif
}

void operator delete[](void *ptr, Space loc) throw() {
#ifdef ENABLE_GPU
  if (loc == Space::Managed) {
#ifndef SW4_USE_UMPIRE
    // std::cout<<"Managed [] delete\n";
    pattr_t *ss = patpush(ptr, NULL);
    if (ss != NULL) {
      global_variables.curr_mem -= ss->size;
      // global_variables.max_mem=std::max(global_variables.max_mem,curr_mem);
    }
    if (ptr == NULL) {
      std::cerr << "Null pointer passed to freee \n";
    } else {
      if (GML(ptr) != Space::Managed) {
        std::cerr << " Wrong space " << as_int(GML(ptr)) << " " << ptr << "\n";
        abort();
      } else {
        SW4_CheckDeviceError(SW4_FREE_MANAGED(ptr));
      }
    }
#else
    umpire::ResourceManager &rma = umpire::ResourceManager::getInstance();
    auto allocator = rma.getAllocator("UM_pool");
    allocator.deallocate(ptr);
#endif
  } else if (loc == Space::Device) {
    // std::cout<<"Managed [] delete\n";
    pattr_t *ss = patpush(ptr, NULL);
    if (ss != NULL) {
      global_variables.curr_mem -= ss->size;
      // global_variables.max_mem=std::max(global_variables.max_mem,curr_mem);
    }
    SW4_CheckDeviceError(SW4_FREE_DEVICE(ptr));

  } else if (loc == Space::Pinned)
    SW4_CheckDeviceError(SW4_FREE_PINNED(ptr));
  else if (loc == Space::Host) {
    // std:cout<<"Calling my placement delete\n";
    ::operator delete(ptr);  // No [] because it is allocated using new(size)
  } else if (loc == Space::Managed_temps) {
#ifdef SW4_USE_UMPIRE
    umpire::ResourceManager &rma = umpire::ResourceManager::getInstance();
    auto allocator = rma.getAllocator("UM_pool_temps");
    allocator.deallocate(ptr);
#else
    SW4_CheckDeviceError(SW4_FREE_MANAGED(ptr));
#endif
  } else {
    std::cerr << "Unknown memory space for de-allocation request "
              << as_int(loc) << "\n";
    abort();
  }
#else
  if ((loc == Space::Managed) || (loc == Space::Device) ||
      (loc == Space::Pinned) || (loc == Space::Managed_temps)) {
    // std::cout<<"Managed delete not available yet \n";
    ::operator delete(ptr);
  } else if (loc == Space::Host) {
    // std:cout<<"Calling my placement delete\n";
    ::operator delete(ptr);
  } else {
    std::cerr << "Unknown memory space for de-allocation request "
              << as_int(loc) << "\n";
    abort();
  }
#endif
}

void operator delete(void *ptr, Space loc, const char *file, int line) throw() {
  ::operator delete(ptr, loc);
}
void operator delete[](void *ptr, Space loc, const char *file,
                       int line) throw() {
  ::operator delete[](ptr, loc);
}

#if defined(SW4_TRACK_MEMORY_ALLOCATIONS)
void assert_check_managed(void *ptr, const char *file, int line) {
  if (ptr == NULL) return;
  pattr_t *ss = patpush(ptr, NULL);
  if (ss != NULL) {
    if (ss->type != Space::Managed) {
      std::cerr << "ASSERT_MANAGED FAILURE in line" << line << " of file "
                << file << "type = " << as_int(ss->type) << "pointer  =" << ptr
                << "\n";
      std::cerr << "ASSERT_MANAGED failed on allocation from line" << ss->line
                << "of " << ss->file << "\n";
    }
  } else {
    std::cerr << "ASSERT_MANAGED FAILURE in line" << line << " of file " << file
              << "type = " << as_int(ss->type) << "pointer  =" << ptr << "\n";
    std::cerr << "No info in map\n Use PointerAttributes\n";
    // printf("Address not in map\n Calling PrintPointerAttributes\n");
    // if ( PointerAttributes(ptr)!=HYPRE_MANAGED_POINTER){
    // fprintf(stderr,"ASSERT_MANAGED FAILURE in line %d of file %s \n NO
    // ALLOCATION INFO\n",line,file);
    //}
  }
}
void assert_check_host(void *ptr, const char *file, int line) {
  if (ptr == NULL) return;
  pattr_t *ss = patpush(ptr, NULL);
  if (ss != NULL) {
    if (ss->type != Space::Host) {
      std::cerr << "ASSERT_HOST FAILURE in line" << line << " of file " << file
                << "type = " << as_int(ss->type) << "pointer  =" << ptr << "\n";
      std::cerr << "ASSERT_HOST failed on allocation from line" << ss->line
                << "of " << ss->file << "\n";
    }
  } else {
    std::cerr << "ASSERT_HOST FAILURE in line" << line << " of file " << file
              << "type = " << as_int(ss->type) << "pointer  =" << ptr << "\n";
    std::cerr << "Ptr not in map\n Use PointerAttribs or somesuch\n";
    // printf("Address not in map\n Calling PrintPointerAttributes\n");
    // if ( PointerAttributes(ptr)!=HYPRE_HOST_POINTER){
    // 	fprintf(stderr,"ASSERT_HOST FAILURE in line %d of file %s \n NO
    // ALLOCATION INFO\n",line,file);
    // }
  }
}
void ptr_push(void *ptr, Space type, size_t size, const char *file, int line) {
  pattr_t *ss = new pattr_t;
  ss->file = file;
  ss->line = line;
  ss->type = type;
  ss->size = size;
  patpush(ptr, ss);
  return;
}

ssize_t getsize(void *ptr) {
  if (ptr == NULL) return -1;
  pattr_t *ss = patpush(ptr, NULL);
  if (ss != NULL)
    return ss->size;
  else
    return -1;
}
#endif

#if defined(ENABLE_CUDA)

void prefetch_to_device(const float_sw4 *ptr) {
#if defined(DISABLE_PREFETCH)
  return;
#else
  if (ptr == NULL) return;
  pattr_t *ss = patpush((void *)ptr, NULL);
  if (ss != NULL) {
    if (ss->size > 0) {
      SW4_MARK_BEGIN(" prefetch_to_device");
      SW4_CheckDeviceError(cudaMemPrefetchAsync(ptr, ss->size, 0, 0));
      SW4_MARK_END(" prefetch_to_device");
    }  // else std::cerr<<"Zero size prefetch \n";
  } else
    std::cerr << "NO prefetch due to unknown address\n";
#endif  // DISABLE_PREFETCH
}

void CheckError(cudaError_t const err, const char *file, char const *const fun,
                const int line) {
  if (err) {
    std::cerr << "CUDA Error Code[" << err << "]: " << cudaGetErrorString(err)
              << " " << file << " " << fun << " Line number:  " << line << "\n"
              << std::flush;
    std::cerr << "ERROR on node " << hostname() << "\n" << std::flush;
    abort();
  }
}
#endif

#if defined(ENABLE_HIP)

void prefetch_to_device(const float_sw4 *ptr) {
  return;
#if defined(DISABLE_PREFETCH)
  return;
#else
  // if (ptr == NULL) return;
  // pattr_t *ss = patpush((void *)ptr, NULL);
  // if (ss != NULL) {
  //   if (ss->size > 0) {
  //     SW4_MARK_BEGIN(" prefetch_to_device");
  //     SW4_CheckDeviceError(cudaMemPrefetchAsync(ptr, ss->size, 0, 0));
  //     SW4_MARK_END(" prefetch_to_device");
  //   }  // else std::cerr<<"Zero size prefetch \n";
  // } else
  //   std::cerr << "NO prefetch due to unknown address\n";
#endif
}

void CheckError(hipError_t const err, const char *file, char const *const fun,
                const int line) {
  if (err) {
    std::cerr << "HIP Error Code[" << err << "]: " << hipGetErrorString(err)
              << " " << file << " " << fun << " Line number:  " << line << "\n"
              << std::flush;
    std::cerr << "ERROR on node " << hostname() << "\n" << std::flush;
    abort();
  }
}
#endif  // ENABLE_HIP

#if defined(ENABLE_GPU)
void *Managed::operator new(size_t len) {
  SW4_MARK_FUNCTION;
  void *ptr;
  ocount++;
  hwm = std::max(hwm, ocount);
  // mem_total+=1;
  // std::cout<<"Total mem is now "<<mem_total<<" MB "<<len<<"\n";
  // std::cout<<"Call to Managed class "<<len<<"\n";
#if defined(SW4_USE_UMPIRE)
  // ptr=SW4_NEW(Space::Managed,char[len]);
  umpire::ResourceManager &rma = umpire::ResourceManager::getInstance();
  auto allocator = rma.getAllocator("UM_object_pool");
  ptr = static_cast<void *>(allocator.allocate(len));
#else
  SW4_CheckDeviceError(SW4_MALLOC_MANAGED(&ptr, len));
#endif
  // SW4_CheckDeviceError(cudaDeviceSynchronize());
  return ptr;
}

void *Managed::operator new[](size_t len) {
  void *ptr;
  ocount++;
  hwm = std::max(hwm, ocount);
#if defined(SW4_USE_UMPIRE)
  umpire::ResourceManager &rma = umpire::ResourceManager::getInstance();
  auto allocator = rma.getAllocator("UM_object_pool");
  ptr = static_cast<void *>(allocator.allocate(len));
  // ptr=SW4_NEW(Space::Managed,char[len]);
#else

  SW4_CheckDeviceError(SW4_MALLOC_MANAGED(&ptr, len));
#endif
  // SW4_CheckDeviceError(cudaDeviceSynchronize());
  return ptr;
}

void Managed::operator delete(void *ptr) {
  // SW4_CheckDeviceError(cudaDeviceSynchronize());
#if defined(SW4_USE_UMPIRE)
  //::operator delete(ptr,Space::Managed);
  ocount--;

  // WARNING DELETES ARE NOOPS. WORKAROUND FOR SUPER SLOW DEALLOCS IN UMPIRE
  return;
  umpire::ResourceManager &rma = umpire::ResourceManager::getInstance();
  auto allocator = rma.getAllocator("UM_object_pool");
  // mem_total+=1;
  allocator.deallocate(ptr);
  // std::cout<<"DTOR 1 "<<mem_total<<"\n";
#else
  SW4_CheckDeviceError(SW4_FREE_MANAGED(ptr));
#endif
}

void Managed::operator delete[](void *ptr) {
  // SW4_CheckDeviceError(cudaDeviceSynchronize());
#if defined(SW4_USE_UMPIRE)
  //::operator delete[](ptr,Space::Managed);
  ocount--;
  // WARNING DELETES ARE NOOPS. WORAROUND FOR SUPER SLOW DEALLOCS IN UMPIRE
  return;
  umpire::ResourceManager &rma = umpire::ResourceManager::getInstance();
  auto allocator = rma.getAllocator("UM_object_pool");
  allocator.deallocate(ptr);
  // mem_total+=1;
  // std::cout<<"DTOR 2 "<<mem_total<<"\n";
#else
  SW4_CheckDeviceError(SW4_FREE_MANAGED(ptr));
#endif
}
#endif

// size_t Managed::mem_total=0;
size_t Managed::ocount = 0;
size_t Managed::hwm = 0;

// AUTOPEEL CODE
#if defined(SW4_TRACK_MEMORY_ALLOCATIONS)
// AUTOPEEL uses getsize which only works when SW4_TRACK_MEMORY_ALLOCATIONS is
// defined
std::string line(int n, int C) {
  std::ostringstream buf;
  buf << "int arg" << C << " = " << n << ";\n";
  return buf.str();
}
std::string line(int *n, int C) {
  std::ostringstream buf;
  buf << "int arg" << C << "[6]={";
  for (int i = 0; i < 5; i++) buf << n[i] << ",";
  buf << n[5] << "};\n";
  return buf.str();
}
std::string line(double n, int C) {
  std::ostringstream buf;
  buf << "double arg" << C << " = " << n << ";\n";
  return buf.str();
}
std::string line(double *n, int C) {
  std::ostringstream buf;
  buf << "double *arg" << C << ";\n cudaMallocManaged((void*)&arg" << C << ","
      << getsize((void **)n) << ");\n";
  return buf.str();
}
std::string line(char n, int C) {
  std::ostringstream buf;
  buf << "char arg" << C << " = \"" << n << "\";\n";
  return buf.str();
}

Apc::Apc(char *s) {
  counter = 0;
  ofile.open(s, std::ios::out);
  ofile << "int main(int argc, char *argv[]){\n";
}
Apc::~Apc() {
  ofile << "FUNCTION(";
  for (int i = 0; i < counter; i++) ofile << "arg" << i << ",";
  ofile << "arg" << counter << ");\n";
  ofile << "}\n";
  ofile.close();
}
#endif

// END AUTOPEEL CODE
void global_prefetch() {
#if defined(ENABLE_CUDA)
#ifdef SW4_MASS_PREFETCH
  int count = 0;
  std::vector<std::string> allocators = {"UM_pool", "UM_pool_temps",
                                         "UM_object_pool"};
  for (auto v : global_variables.massprefetch) {
    // std::cout<<"global_prefetch "<<std::get<0>(v)<<" of
    // "<<std::get<1>(v)<<" bytes\n";
    auto size = umpire::ResourceManager::getInstance()
                    .getAllocator(allocators[count])
                    .getHighWatermark();
    std::cout << "GLOBAL PREFETCH SIZES FOR " << allocators[count] << " "
              << size << " , " << std::get<1>(v) << "\n";
#define PREFETCH_ALL 1
#ifdef PREFETCH_ALL
    SW4_CheckDeviceError(cudaMemPrefetchAsync(std::get<0>(v), std::get<1>(v),
                                              global_variables.device, 0));
#else
    SW4_CheckDeviceError(
        cudaMemPrefetchAsync(std::get<0>(v), size, global_variables.device, 0));
#endif
    SW4_CheckDeviceError(cudaStreamSynchronize(0));
    count++;
  }
#endif  // SW4_MASS_PREFETCH

#endif  // ENABLE_CUDA
}
std::vector<int> factors(int N) {
  std::vector<int> v;
  for (int i = 1; i <= N; i++)
    if (N % i == 0) v.push_back(i);
  return v;
}
std::vector<int> factors(int N, int start) {
  std::vector<int> v;
  for (int i = start; i <= N; i++)
    if (N % i == 0) v.push_back(i);
  return v;
}
// GML needs to be hipified
#ifdef ENABLE_CUDA
Space GML(const void *ptr) {
  struct cudaPointerAttributes attr;
  if (cudaPointerGetAttributes(&attr, ptr) == cudaErrorInvalidValue) {
    // This shuld go away with Cuda 11
    std::cerr << "Invalid value in GML \n";
    return Space::Host;
  }
  if (attr.type == cudaMemoryTypeUnregistered) {
    return Space::Host;
  } else if (attr.type == cudaMemoryTypeHost) {
    return Space::Pinned;
  } else if (attr.type == cudaMemoryTypeDevice) {
    return Space::Device;
  } else if (attr.type == cudaMemoryTypeManaged) {
    return Space::Managed;
  } else
    return Space::Space_Error;
}
#endif
#ifdef ENABLE_HIP
Space GML(const void *ptr) {
  return Space::Managed;
  struct hipPointerAttribute_t attr;
  if (hipPointerGetAttributes(&attr, ptr) == hipErrorInvalidValue) {
    // This shuld go away with Cuda 11
    std::cerr << "Invalid value in GML \n";
    return Space::Host;
  }
  if (attr.memoryType == hipMemoryTypeHost) {
    return Space::Pinned;
  } else if (attr.memoryType == hipMemoryTypeDevice) {
    return Space::Managed;
  } else if (attr.memoryType == hipMemoryTypeUnified) {
    return Space::Managed;
  } else
    return Space::Space_Error;
}
#endif

void invert(float_sw4 *A, int msize) {
  float_sw4 B[9];

#define M(i, j, d) A[(d * 9) + (i - 1) + (j - 1) * 3]
#define I(i, j) B[(j - 1) + (i - 1) * 3]

  for (int c = 0; c < msize; c++) {
    float_sw4 dA = M(1, 1, c) * M(2, 2, c) * M(3, 3, c) +
                   M(1, 2, c) * M(2, 3, c) * M(3, 1, c) +
                   M(1, 3, c) * M(2, 1, c) * M(3, 2, c) -
                   M(1, 3, c) * M(2, 2, c) * M(3, 1, c) -
                   M(1, 2, c) * M(2, 1, c) * M(3, 3, c) -
                   M(1, 1, c) * M(2, 3, c) * M(3, 2, c);

    I(1, 1) = M(2, 2, c) * M(3, 3, c) - M(2, 3, c) * M(3, 2, c);
    I(1, 2) = M(1, 3, c) * M(3, 2, c) - M(1, 2, c) * M(3, 3, c);
    I(1, 3) = M(1, 2, c) * M(2, 3, c) - M(2, 2, c) * M(1, 3, c);

    I(2, 1) = M(2, 3, c) * M(3, 1, c) - M(3, 3, c) * M(2, 1, c);
    I(2, 2) = M(1, 1, c) * M(3, 3, c) - M(3, 1, c) * M(1, 3, c);
    I(2, 3) = M(1, 3, c) * M(2, 1, c) - M(2, 3, c) * M(1, 1, c);

    I(3, 1) = M(2, 1, c) * M(3, 2, c) - M(3, 1, c) * M(2, 2, c);
    I(3, 2) = M(1, 2, c) * M(3, 1, c) - M(3, 2, c) * M(1, 1, c);
    I(3, 3) = M(1, 1, c) * M(2, 2, c) - M(2, 1, c) * M(1, 2, c);

#ifdef DEBUG
    float_sw4 Prod[9];
#define P(i, j) Prod[(j - 1) + (i - 1) * 3]
    for (int i = 1; i < 4; i++) {
      for (int j = 1; j < 4; j++) {
        float_sw4 sum = 0.0;
        for (int k = 1; k < 4; k++) {
          sum += M(i, k, c) * I(k, j);
        }
        P(i, j) = sum / dA;
      }
    }

    for (int i = 1; i < 4; i++) {
      for (int j = 1; j < 4; j++) {
        std::cout << P(i, j) << ",";
      }
      std::cout << "\n";
    }
    std::cout << "_______________________\n";
#endif

    for (int i = 1; i < 4; i++) {
      for (int j = 1; j < 4; j++) {
        M(j, i, c) = I(i, j) / dA;
      }
    }
  }
#undef M
#undef P
#undef I
}

bool mpi_supports_device_buffers() {
#if defined(ENABLE_CUDA)
#if defined(SMPI_VERSION)
  return MPIX_Query_cuda_support() == 1;
//   printf("Running using Spectrum MPI SMPI_VERSION \n");
//   if (MPIX_Query_cuda_support())
//     std::cout<<"Running with support for device buffer\n";
//   else
//     std::cout<<"Running WITHOUT support for device buffer\n";
#endif
#endif
  return false;
}
long long node_mem() {
  struct sysinfo memInfo;

  sysinfo(&memInfo);
  long long totalVirtualMem = memInfo.totalram;
  totalVirtualMem += memInfo.totalswap;
  totalVirtualMem *= memInfo.mem_unit;

  long long virtualMemUsed = memInfo.totalram - memInfo.freeram;
  virtualMemUsed += memInfo.totalswap - memInfo.freeswap;
  virtualMemUsed *= memInfo.mem_unit;

  long long totalPhysMem = memInfo.totalram;
  totalPhysMem *= memInfo.mem_unit;

  long long physMemUsed = memInfo.totalram - memInfo.freeram;
  physMemUsed *= memInfo.mem_unit;
  // std::cerr<<"VM Total = "<<totalVirtualMem/1024/1024.0<<" "<<"VM Used =
  // "<<virtualMemUsed/1024/1024.0<<" Phys Total =
  // "<<totalPhysMem/1024/1024.0<<" Phys Used "<<physMemUsed/1024/1024.0<<"\n";
  return physMemUsed;
}

std::string hostname() {
  char hostname[HOST_NAME_MAX];
  gethostname(hostname, HOST_NAME_MAX);
  return std::string(hostname);
}
