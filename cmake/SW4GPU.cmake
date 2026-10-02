# This module is entered only for an explicitly selected GPU backend.
if(SW4_PRECISION STREQUAL single)
  message(FATAL_ERROR "The existing RAJA implementation supports double precision; select SW4_PRECISION=double")
endif()
if(SW4_BACKEND STREQUAL CUDA)
  enable_language(CUDA)
  find_package(CUDAToolkit REQUIRED)
  set(SW4_GPU_LANGUAGE CUDA)
  set(SW4_GPU_SUFFIX cu)
else()
  if(CMAKE_VERSION VERSION_LESS 3.21)
    message(FATAL_ERROR "The HIP backend requires CMake 3.21 or newer")
  endif()
  enable_language(HIP)
  find_package(hip REQUIRED CONFIG)
  set(SW4_GPU_LANGUAGE HIP)
  set(SW4_GPU_SUFFIX hip)
endif()
find_package(RAJA REQUIRED CONFIG)
if(SW4_USE_UMPIRE)
  find_package(umpire REQUIRED CONFIG)
endif()

# Use wrappers rather than changing source-file language globally: the native
# inversion target can compile the same files as C++ in the same build tree.
set(SW4_GPU_SOURCES)
foreach(source IN LISTS SW4_RAJA_SOURCES)
  if(source MATCHES "\\.C$")
    get_filename_component(name "${source}" NAME_WE)
    set(wrapper "${PROJECT_BINARY_DIR}/gpu_sources/${name}.${SW4_GPU_SUFFIX}")
    file(GENERATE OUTPUT "${wrapper}" CONTENT "#include \"${PROJECT_SOURCE_DIR}/${source}\"\n")
    set_source_files_properties("${wrapper}" PROPERTIES LANGUAGE ${SW4_GPU_LANGUAGE})
    list(APPEND SW4_GPU_SOURCES "${wrapper}")
  else()
    list(APPEND SW4_GPU_SOURCES "${source}")
  endif()
endforeach()
file(GENERATE OUTPUT "${PROJECT_BINARY_DIR}/gpu_sources/main.${SW4_GPU_SUFFIX}"
  CONTENT "#include \"${PROJECT_SOURCE_DIR}/src/main.C\"\n")
set_source_files_properties("${PROJECT_BINARY_DIR}/gpu_sources/main.${SW4_GPU_SUFFIX}"
  PROPERTIES LANGUAGE ${SW4_GPU_LANGUAGE})
add_executable(sw4 "${PROJECT_BINARY_DIR}/gpu_sources/main.${SW4_GPU_SUFFIX}"
  ${SW4_GPU_SOURCES} ${SW4_QUADPACK_SOURCES})
target_link_libraries(sw4 PRIVATE sw4_dependencies RAJA OpenMP::OpenMP_CXX)
target_compile_definitions(sw4 PRIVATE SW4_USE_RAJA=1 ENABLE_MPI_TIMING_BARRIER=1
  SW4_STAGED_MPI_BUFFERS=1 USE_DIRECT_INVERSE=1 SW4_GHCOF_NO_GP_IS_ZERO=1
  SW4_CROUTINES RAJA_USE_RESTRICT_PTR)
# Umpire already contains its compiled fmt implementation. Its exported
# header-only fmt definition triggers an NVCC 13 host-code generation bug;
# use the same compiled-library mode as the existing GPU Make workflow.
target_compile_options(sw4 PRIVATE "$<$<COMPILE_LANGUAGE:CUDA>:-UFMT_HEADER_ONLY>")
if(SW4_USE_UMPIRE)
  target_compile_definitions(sw4 PRIVATE SW4_USE_UMPIRE=1)
  if(TARGET umpire::umpire)
    target_link_libraries(sw4 PRIVATE umpire::umpire)
  else()
    target_link_libraries(sw4 PRIVATE umpire)
  endif()
endif()
set_target_properties(sw4 PROPERTIES CXX_STANDARD ${SW4_GPU_CXX_STANDARD}
  CXX_STANDARD_REQUIRED ON)
if(SW4_BACKEND STREQUAL CUDA)
  set_target_properties(sw4 PROPERTIES CUDA_STANDARD ${SW4_GPU_CXX_STANDARD}
    CUDA_STANDARD_REQUIRED ON CUDA_SEPARABLE_COMPILATION ON)
  target_compile_definitions(sw4 PRIVATE ENABLE_CUDA=1 SW4_USE_CMEM=1)
  target_compile_options(sw4 PRIVATE
    "$<$<COMPILE_LANGUAGE:CUDA>:--expt-extended-lambda;--expt-relaxed-constexpr;-Xcompiler=-fopenmp>")
  target_link_libraries(sw4 PRIVATE CUDA::cudart CUDA::cuda_driver)
  # NVML provides the existing device reporting and affinity support.
  find_library(SW4_NVML_LIBRARY NAMES nvidia-ml
    HINTS "${CUDAToolkit_LIBRARY_DIR}/stubs" "${CUDAToolkit_LIBRARY_ROOT}/lib64/stubs")
  if(NOT SW4_NVML_LIBRARY)
    message(FATAL_ERROR "CUDA backend requires NVML (libnvidia-ml)")
  endif()
  target_link_libraries(sw4 PRIVATE "${SW4_NVML_LIBRARY}")
else()
  set_target_properties(sw4 PROPERTIES HIP_STANDARD ${SW4_GPU_CXX_STANDARD}
    HIP_STANDARD_REQUIRED ON)
  target_compile_definitions(sw4 PRIVATE ENABLE_HIP=1 SW4_NO_ROCTRACER=1)
  target_compile_options(sw4 PRIVATE "$<$<COMPILE_LANGUAGE:HIP>:-fopenmp>")
  target_link_libraries(sw4 PRIVATE hip::host)
endif()

# Built and executed only when explicitly requested; never added to CTest.
if(SW4_BACKEND STREQUAL CUDA)
  add_executable(sw4_cuda_compat EXCLUDE_FROM_ALL tests/backends/cuda_compat.cu)
  target_include_directories(sw4_cuda_compat PRIVATE "${PROJECT_SOURCE_DIR}/src")
  target_link_libraries(sw4_cuda_compat PRIVATE CUDA::cudart)
  set_target_properties(sw4_cuda_compat PROPERTIES CUDA_STANDARD 11)
endif()
