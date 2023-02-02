#ifdef SW4_USE_UMPIRE
#include "umpire/ResourceManager.hpp"
//#include "umpire/Umpire.hpp"
//#include "umpire/Allocator.hpp"
//#include "umpire/strategy/DynamicPool.hpp"
#include "umpire/strategy/QuickPool.hpp"
//#include "umpire/strategy/MixedPool.hpp"
//#include "umpire/util/StatisticsDatabase.hpp"
#include "umpire/strategy/AlignedAllocator.hpp"
#include "umpire/strategy/AllocationAdvisor.hpp"
#include "umpire/strategy/MonotonicAllocationStrategy.hpp"
#include "umpire/util/Macros.hpp"
#endif
#include <string>

void init_umpire(int device){

  
 #ifdef SW4_USE_UMPIRE
  umpire::ResourceManager &rma = umpire::ResourceManager::getInstance();
#ifdef ENABLE_HIP
  auto allocator = rma.getAllocator("DEVICE::" + std::to_string(device));
#else
  auto allocator = rma.getAllocator("UM");
#endif
  // auto device_allocator = rma.getAllocator("DEVICE");
#ifdef ENABLE_HIP
  const size_t pool_size =
      static_cast<size_t>(20) * 1024 * 1024 * 1024;  //+102*1024*1024;
#else
  const size_t pool_size =
      static_cast<size_t>(15) * 1024 * 1024 * 1024;  //+102*1024*1024;
#endif

#ifdef ENABLE_HIP
  auto pref_allocator = allocator;
#else
  // auto aligned_allocator =
  // rma.makeAllocator<umpire::strategy::AlignedAllocator>(
  //   "aligned_allocator", allocator, 256);
  auto pref_allocator = rma.makeAllocator<umpire::strategy::AllocationAdvisor>(
      "preferred_location_device", allocator, "SET_PREFERRED_LOCATION",
      global_variables.device);
#endif

  const int alignment = 512; // 1024 may be 1% faster on Crusher
  auto pooled_allocator =
      rma.makeAllocator<umpire::strategy::QuickPool, true>(
							   std::string("UM_pool"), pref_allocator, pool_size, 1024 * 1024, alignment);
#ifdef ENABLE_HIP
  const size_t pool_size_small = static_cast<size_t>(1024) * 1024 * 1024;
#else
  const size_t pool_size_small = static_cast<size_t>(250) * 1024 * 1024;
#endif

  // This is a temporary workaround to the issue of Umpire always using device 0
  // for cudaMemAdvises using AllocationAdvisor.
  // if (global_variables.num_devices==1){

  // auto pooled_allocator_small =static_cast<size_t>(250)*1024*1024;
  auto pooled_allocator_small =
      rma.makeAllocator<umpire::strategy::QuickPool, true>(
							   std::string("UM_pool_temps"), pref_allocator, pool_size_small, 1024 * 1024,
          alignment);

#ifdef ENABLE_HIP
  const size_t object_pool_size = static_cast<size_t>(3) *1024* 1024 * 1024;
#else
  const size_t object_pool_size = static_cast<size_t>(500) * 1024 * 1024;
#endif

  // rma.makeAllocator<umpire::strategy::MonotonicAllocationStrategy,false>(string("UM_object_pool"),
  //					   object_pool_size,allocator);

  auto pooled_allocator_objects =
      rma.makeAllocator<umpire::strategy::QuickPool, false>(
							    std::string("UM_object_pool"), allocator, object_pool_size);

#ifdef SW4_MASS_PREFETCH
  std::cout << "Mass prefetch operational\n";
  global_variables.massprefetch.push_back(std::make_tuple(
      static_cast<char *>(pooled_allocator.allocate(1)), pool_size));
  global_variables.massprefetch.push_back(
      std::make_tuple(static_cast<char *>(pooled_allocator_small.allocate(1)),
                      pool_size_small));
  global_variables.massprefetch.push_back(
      std::make_tuple(static_cast<char *>(pooled_allocator_objects.allocate(1)),
                      object_pool_size));
#endif
  // rma.makeAllocator<umpire::strategy::MixedPool,false>(string("UM_object_pool"),
  //							   allocator,object_pool_size);

  // rma.makeAllocator<umpire::strategy::FixedPool,false>(string("UM_object_pool"),
  // 						     allocator,object_pool_size);

  // } else {
  //   auto pooled_allocator_small =
  //     rma.makeAllocator<umpire::strategy::QuickPool,true>(string("UM_pool_temps"),
  //  							   pref_allocator,pool_size_small);
  // }

  // auto pooled_allocator2 =
  //   rma.makeAllocator<umpire::strategy::QuickPool,false>(string("UM_pool_temps"),
  //                                                   allocator);
#endif
  
}
