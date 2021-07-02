set(CMAKE_VERBOSE_MAKEFILE ON)
cmake_minimum_required(VERSION 3.18)
project(SW4 VERSION 1.0 LANGUAGES CXX Fortran)

set(BLT_CXX_STD "c++11" CACHE STRING "")
include(tpl/blt/SetupBLT.cmake)

configure_file(SW4Config.h.in SW4Config.h)


# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if (NOT BLT_CXX_STD)
    set(BLT_CXX_STD "c++11" CACHE STRING "")
endif()

option(ENABLE_RAJA "ENABLE_RAJA" OFF)
if (ENABLE_RAJA)
  set (ENABLE_TESTS Off CACHE BOOL "" FORCE)
  set (ENABLE_EXAMPLES Off CACHE BOOL "" FORCE)
  set (ENABLE_REPRODUCERS Off CACHE BOOL "" FORCE)
  set (ENABLE_EXERCISES Off CACHE BOOL "" FORCE)
  set (ENABLE_DOCUMENTATION Off CACHE BOOL "" FORCE)
  set (ENABLE_BENCHMARKS Off CACHE BOOL "" FORCE)


  if (DEFINED RAJA_DIR)
    find_package(RAJA REQUIRED)
    blt_print_target_properties(TARGET RAJA)
  else ()
    add_subdirectory(tpl/RAJA)
  endif ()
endif()
add_subdirectory(src)
