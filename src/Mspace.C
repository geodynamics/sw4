#include <unordered_map>

#include "Mspace.h"
#include "caliper.h"
struct global_variable_holder_struct global_variables = {0, 0, 0, 0,
                                                         0, 0, 0, 1};
using namespace std;

void presetGPUID() {
#ifdef ENABLE_CUDA
  int devices_per_node = 4;
  SW4_CheckDeviceError(cudaGetDeviceCount(&devices_per_node));
  global_variables.num_devices = devices_per_node;
  if (devices_per_node > 1) {
    char *crank = getenv("OMPI_COMM_WORLD_LOCAL_RANK");
    int device = atoi(crank) % devices_per_node;
    global_variables.device = device;
    printf(" presetGPU Called ::  LOCAL RANK %d \n", device);
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
  printf("Device set to %d \n", global_variables.device);
#endif
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
  SW4_CheckDeviceError(cudaMemGetInfo(&mfree, &mtotal));
  global_variables.gpu_memory_hwm =
      std::max((mtotal - mfree), global_variables.gpu_memory_hwm);
}

void print_hwm() {
  const int allocator_count = 3;
#if defined(ENABLE_CUDA)

  float hwm_local[allocator_count], hwm_global[allocator_count];
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
  MPI_Allreduce(&hwm_local, &hwm_global, allocator_count, MPI_FLOAT, MPI_MAX,
                MPI_COMM_WORLD);
  for (int i = 0; i < allocator_count; i++)
    if (hwm_local[i] == hwm_global[i]) {
      std::cout << i << " Global Device HWM is " << hwm_global[i] << " GB\n";
      // umpire::util::StatisticsDatabase::getDatabase()->printStatistics(std::cout);
    }
  if (Managed::hwm > 0)
    std::cout << "Managed object count & HWM are " << Managed::ocount << " & "
              << Managed::hwm << "\n";
  if (Managed::ocount != 0)
    std::cerr << "WARNING :: Managed object count should be zero at the end of "
                 "the simulation\n";
    // std::cout<<" ~HOST MEM MAX
    // "<<global_variables.host_mem_hwm/1024.0/1024.0<<" MB\n";
#endif  // ENABLE_CUDA
}
void *operator new(std::size_t size, Space loc) throw(std::bad_alloc) {
#ifdef ENABLE_CUDA
  if (loc == Managed) {
    // std::cout<<"Managed allocation \n";
    if (size == 0) size = 1;  // new has to return an valid pointer for 0 size.
    void *ptr;
#ifndef SW4_USE_UMPIRE
    if (cudaMallocManaged(&ptr, size) != cudaSuccess) {
      std::cerr << "Mananged memory allocation failed " << size << "\n";
      throw std::bad_alloc();
    } else {
      check_mem();
      global_variables.curr_mem += size;
      global_variables.max_mem =
          std::max(global_variables.max_mem, global_variables.curr_mem);
      SW4_CheckDeviceError(cudaMemAdvise(ptr, size,
                                         cudaMemAdviseSetPreferredLocation,
                                         global_variables.device));
      return ptr;
    }
#else
    umpire::ResourceManager &rma = umpire::ResourceManager::getInstance();
    auto allocator = rma.getAllocator("UM_pool");
    ptr = static_cast<void *>(allocator.allocate(size));
    SW4_CheckDeviceError(cudaMemAdvise(
        ptr, size, cudaMemAdviseSetPreferredLocation, global_variables.device));
    // std::cout<<"PTR 1 "<<ptr<<"\n";
    // SW4_CheckDeviceError(cudaMemset(ptr,0,size));
    return ptr;
#endif
  } else if (loc == Host) {
    // std::cout<<"Calling my placement new \n";
    // global_variables.host_curr_mem+=size;
    // global_variables.host_max_mem=std::max(global_variables.host_max_mem,global_variables.host_curr_mem);
    return ::operator new(size);
  } else if (loc == Device) {
    // std::cout<<"Managed allocation \n";
    if (size == 0) size = 1;  // new has to return an valid pointer for 0 size.
    void *ptr;
    if (cudaMalloc(&ptr, size) != cudaSuccess) {
      std::cerr << "Device memory allocation failed " << size << "\n";
      throw std::bad_alloc();
    } else
      return ptr;
  } else if (loc == Pinned) {
    if (size == 0) size = 1;  // new has to return an valid pointer for 0 size.
    void *ptr;
    SW4_CheckDeviceError(cudaHostAlloc(&ptr, size, cudaHostAllocMapped));
    return ptr;
  } else if (loc == Managed_temps) {
#ifdef SW4_USE_UMPIRE
    umpire::ResourceManager &rma = umpire::ResourceManager::getInstance();
    auto allocator = rma.getAllocator("UM_pool_temps");
    void *ptr = static_cast<void *>(allocator.allocate(size));
    // SW4_CheckDeviceError(cudaMemAdvise(ptr,size,cudaMemAdviseSetPreferredLocation,0));
    // std::cout<<"PTR 1 "<<ptr<<"\n";
    // SW4_CheckDeviceError(cudaMemset(ptr,0,size));
    return ptr;
#else
    std::cerr << "Managed_temp location no defined\n";
    return ::operator new(size, Managed);
#endif
  } else {
    std::cerr << "Unknown memory space for allocation request " << loc << "\n";
    throw std::bad_alloc();
  }
#else
  if ((loc == Managed) || (loc == Device) || (loc == Pinned)) {
    // std::cout<<"Managed location not available yet \n";
    return ::operator new(size);
  } else if (loc == Host) {
    // std::cout<<"Calling my placement new \n";
    return ::operator new(size);
  } else {
    std::cerr << "Unknown memory space for allocation request\n";
    throw std::bad_alloc();
  }
#endif
}
void *operator new(std::size_t size, Space loc, char *file,
                   int line) throw(std::bad_alloc) {
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

void *operator new[](std::size_t size, Space loc) throw(std::bad_alloc) {
#ifdef ENABLE_CUDA
  if (loc == Managed) {
    // std::cout<<"Managed [] allocation \n";
    if (size == 0) size = 1;  // new has to return an valid pointer for 0 size.
    void *ptr;
#ifndef SW4_USE_UMPIRE
    if (cudaMallocManaged(&ptr, size) != cudaSuccess) {
      std::cerr << "Managed memory allocation failed " << size << "\n";
      throw std::bad_alloc();
    } else {
      check_mem();
      global_variables.curr_mem += size;
      global_variables.max_mem =
          std::max(global_variables.max_mem, global_variables.curr_mem);
      SW4_CheckDeviceError(cudaMemAdvise(ptr, size,
                                         cudaMemAdviseSetPreferredLocation,
                                         global_variables.device));
      return ptr;
    }
#else
    umpire::ResourceManager &rma = umpire::ResourceManager::getInstance();
    auto allocator = rma.getAllocator("UM_pool");
    ptr = static_cast<void *>(allocator.allocate(size));
    SW4_CheckDeviceError(cudaMemAdvise(
        ptr, size, cudaMemAdviseSetPreferredLocation, global_variables.device));
    // std::cout<<"PTR 2 "<<ptr<<"\n";
    return ptr;
#endif
  } else if (loc == Host) {
    // std::cout<<"Calling my placement new \n";
    // global_variables.host_curr_mem+=size;
    // global_variables.host_max_mem=std::max(global_variables.host_max_mem,global_variables.host_curr_mem);
    return ::operator new(size);
  } else if (loc == Device) {
    // std::cout<<"Managed allocation \n";
    if (size == 0) size = 1;  // new has to return an valid pointer for 0 size.
    void *ptr;
    if (cudaMalloc(&ptr, size) != cudaSuccess) {
      std::cerr << "Device memory allocation failed " << size << "\n";
      throw std::bad_alloc();
    } else
      return ptr;
  } else if (loc == Pinned) {
    if (size == 0) size = 1;  // new has to return an valid pointer for 0 size.
    void *ptr;
    SW4_CheckDeviceError(cudaHostAlloc(&ptr, size, cudaHostAllocMapped));
    return ptr;
  } else if (loc == Managed_temps) {
#if defined(SW4_USE_UMPIRE)
    umpire::ResourceManager &rma = umpire::ResourceManager::getInstance();
    auto allocator = rma.getAllocator("UM_pool_temps");
    void *ptr = static_cast<void *>(allocator.allocate(size));
    // SW4_CheckDeviceError(cudaMemAdvise(ptr,size,cudaMemAdviseSetPreferredLocation,0));
    // std::cout<<"PTR 1 "<<ptr<<"\n";
    // SW4_CheckDeviceError(cudaMemset(ptr,0,size));
    return ptr;
#else
    std::cerr << " Memory location Managed_temps is not defined\n";
    return ::operator new(size, Managed);
#endif
  } else {
    // cudaHostAlloc(&ptr,size+sizeof(size_t)*MEM_PAD_LEN,cudaHostAllocMapped));
    std::cerr << "Unknown memory space for allocation request " << loc << "\n";
    throw std::bad_alloc();
  }

#else  // !ENABLE_CUDA
  if ((loc == Managed) || (loc == Device) || (loc == Pinned) ||
      (loc == Managed_temps)) {
    // std::cout<<"Managed location not available yet \n";
    return ::operator new(size);
  } else if (loc == Host) {
    // std::cout<<"Calling my placement new \n";
    return ::operator new(size);
  } else {
    std::cerr << "Unknown memory space for allocation request " << loc << "\n";
    throw std::bad_alloc();
  }
#endif
}
void *operator new[](std::size_t size, Space loc, const char *file, int line) {
  // std::cout<<"Calling tracking new from "<<line<<" of "<<file<<"\n";
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
#ifdef ENABLE_CUDA
  if ((loc == Managed)) {
    // std::cout<<"Managed delete\n";
#ifndef SW4_USE_UMPIRE
    pattr_t *ss = patpush(ptr, NULL);
    if (ss != NULL) {
      global_variables.curr_mem -= ss->size;
      // global_variables.max_mem=std::max(global_variables.max_mem,global_variables.curr_mem);
    }
    SW4_CheckDeviceError(cudaFree(ptr));
#else
    umpire::ResourceManager &rma = umpire::ResourceManager::getInstance();
    auto allocator = rma.getAllocator("UM_pool");
    allocator.deallocate(ptr);
#endif
  } else if (loc == Device) {
    pattr_t *ss = patpush(ptr, NULL);
    if (ss != NULL) {
      global_variables.curr_mem -= ss->size;
      // global_variables.max_mem=std::max(global_variables.max_mem,global_variables.curr_mem);
    }
    SW4_CheckDeviceError(cudaFree(ptr));
  } else if (loc == Pinned)
    SW4_CheckDeviceError(cudaFreeHost(ptr));
  else if (loc == Host) {
    // std:cout<<"Calling my placement delete\n";
    ::operator delete(ptr);
  } else if (loc == Managed_temps) {
#if defined(SW4_USE_UMPIRE)
    umpire::ResourceManager &rma = umpire::ResourceManager::getInstance();
    auto allocator = rma.getAllocator("UM_pool_temps");
    allocator.deallocate(ptr);
#else
    std::cerr << "Memory location Managed_temps not defined\n";
#endif
  } else {
    std::cerr << "Unknown memory space for de-allocation request\n";
  }
#else
  if ((loc == Managed) || (loc == Device)) {
    // std::cout<<"Managed delete not available yet \n";
    ::operator delete(ptr);
  } else if (loc == Host) {
    // std:cout<<"Calling my placement delete\n";
    ::operator delete(ptr);
  } else {
    std::cerr << "Unknown memory space for de-allocation request " << loc
              << "\n";
  }
#endif
}

void operator delete[](void *ptr, Space loc) throw() {
#ifdef ENABLE_CUDA
  if (loc == Managed) {
#ifndef SW4_USE_UMPIRE
    // std::cout<<"Managed [] delete\n";
    pattr_t *ss = patpush(ptr, NULL);
    if (ss != NULL) {
      global_variables.curr_mem -= ss->size;
      // global_variables.max_mem=std::max(global_variables.max_mem,curr_mem);
    }
    SW4_CheckDeviceError(cudaFree(ptr));
#else
    umpire::ResourceManager &rma = umpire::ResourceManager::getInstance();
    auto allocator = rma.getAllocator("UM_pool");
    allocator.deallocate(ptr);
#endif
  } else if (loc == Device) {
    // std::cout<<"Managed [] delete\n";
    pattr_t *ss = patpush(ptr, NULL);
    if (ss != NULL) {
      global_variables.curr_mem -= ss->size;
      // global_variables.max_mem=std::max(global_variables.max_mem,curr_mem);
    }
    SW4_CheckDeviceError(cudaFree(ptr));

  } else if (loc == Pinned)
    SW4_CheckDeviceError(cudaFreeHost(ptr));
  else if (loc == Host) {
    // std:cout<<"Calling my placement delete\n";
    ::operator delete(ptr);
  } else if (loc == Managed_temps) {
#ifdef SW4_USE_UMPIRE
    umpire::ResourceManager &rma = umpire::ResourceManager::getInstance();
    auto allocator = rma.getAllocator("UM_pool_temps");
    allocator.deallocate(ptr);
#endif
  } else {
    std::cerr << "Unknown memory space for de-allocation request " << loc
              << "\n";
  }
#else
  if ((loc == Managed) || (loc == Device) || (loc == Pinned)) {
    // std::cout<<"Managed delete not available yet \n";
    ::operator delete(ptr);
  } else if (loc == Host) {
    // std:cout<<"Calling my placement delete\n";
    ::operator delete(ptr);
  } else {
    std::cerr << "Unknown memory space for de-allocation request " << loc
              << "\n";
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
    if (ss->type != Managed) {
      std::cerr << "ASSERT_MANAGED FAILURE in line" << line << " of file "
                << file << "type = " << ss->type << "pointer  =" << ptr << "\n";
      std::cerr << "ASSERT_MANAGED failed on allocation from line" << ss->line
                << "of " << ss->file << "\n";
    }
  } else {
    std::cerr << "ASSERT_MANAGED FAILURE in line" << line << " of file " << file
              << "type = " << ss->type << "pointer  =" << ptr << "\n";
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
    if (ss->type != Host) {
      std::cerr << "ASSERT_HOST FAILURE in line" << line << " of file " << file
                << "type = " << ss->type << "pointer  =" << ptr << "\n";
      std::cerr << "ASSERT_HOST failed on allocation from line" << ss->line
                << "of " << ss->file << "\n";
    }
  } else {
    std::cerr << "ASSERT_HOST FAILURE in line" << line << " of file " << file
              << "type = " << ss->type << "pointer  =" << ptr << "\n";
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
#endif
}

void CheckError(cudaError_t const err, const char *file, char const *const fun,
                const int line) {
  if (err) {
    std::cerr << "CUDA Error Code[" << err << "]: " << cudaGetErrorString(err)
              << " " << file << " " << fun << " Line number:  " << line << "\n";
  }
}
#endif

#if defined(ENABLE_CUDA)
void *Managed::operator new(size_t len) {
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
  SW4_CheckDeviceError(cudaMallocManaged(&ptr, len));
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

  SW4_CheckDeviceError(cudaMallocManaged(&ptr, len));
#endif
  // SW4_CheckDeviceError(cudaDeviceSynchronize());
  return ptr;
}

void Managed::operator delete(void *ptr) {
  // SW4_CheckDeviceError(cudaDeviceSynchronize());
#if defined(SW4_USE_UMPIRE)
  //::operator delete(ptr,Space::Managed);
  ocount--;

  // WARNING DELETES ARE NOOPS. WORAROUND FOR SUPER SLOW DEALLOCS IN UMPIRE
  return;
  umpire::ResourceManager &rma = umpire::ResourceManager::getInstance();
  auto allocator = rma.getAllocator("UM_object_pool");
  // mem_total+=1;
  allocator.deallocate(ptr);
  // std::cout<<"DTOR 1 "<<mem_total<<"\n";
#else
  SW4_CheckDeviceError(cudaFree(ptr));
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
  SW4_CheckDeviceError(cudaFree(ptr));
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
#ifdef SW4_MASS_PREFETCH
  int count = 0;
  std::vector<std::string> allocators = {"UM_pool", "UM_pool_temps",
                                         "UM_object_pool"};
  for (auto v : global_variables.massprefetch) {
    // std::cout<<"global_prefetch "<<std::get<0>(v)<<" of "<<std::get<1>(v)<<"
    // bytes\n";
    auto size = umpire::ResourceManager::getInstance()
                    .getAllocator(allocators[count])
                    .getHighWatermark();
    std::cout << "GLOBAL PREFETCH SIZES " << size << " , " << std::get<1>(v) << "\n";
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
#endif
}
