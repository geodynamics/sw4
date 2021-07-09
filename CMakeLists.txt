set(CMAKE_VERBOSE_MAKEFILE ON)
cmake_minimum_required(VERSION 3.18)
project(SW4 VERSION 1.0 LANGUAGES CXX Fortran)



configure_file(SW4Config.h.in SW4Config.h)


# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if (NOT BLT_CXX_STD)
    set(BLT_CXX_STD "c++11" CACHE STRING "")
endif()

# set options for the RAJA submodule build
set (ENABLE_TESTS Off CACHE BOOL "" FORCE)
set (ENABLE_EXAMPLES Off CACHE BOOL "" FORCE)
set (ENABLE_REPRODUCERS Off CACHE BOOL "" FORCE)
set (ENABLE_EXERCISES Off CACHE BOOL "" FORCE)
set (ENABLE_DOCUMENTATION Off CACHE BOOL "" FORCE)
set (ENABLE_BENCHMARKS Off CACHE BOOL "" FORCE)
set(BLT_CXX_STD "c++11" CACHE STRING "")
include(tpl/blt/SetupBLT.cmake)

if (DEFINED UMPIRE_DIR)
  find_package(Umpire REQUIRED)
  blt_print_target_properties(TARGET Umpire)
else ()
add_subdirectory(tpl/Umpire)
endif ()

if (DEFINED RAJA_DIR)
  find_package(RAJA REQUIRED)
  blt_print_target_properties(TARGET RAJA)
else ()
  add_subdirectory(tpl/RAJA)
endif ()


MACRO (TODAY RESULT)
    IF (WIN32)
        EXECUTE_PROCESS(COMMAND "date" "/T" OUTPUT_VARIABLE ${RESULT})
        STRING(REGEX REPLACE "(..)/(..)/..(..).*" "\\3\\2\\1" ${RESULT} ${${RESULT}})
    ELSEIF(UNIX)
        EXECUTE_PROCESS(COMMAND "date" "+%d/%m/%Y" OUTPUT_VARIABLE ${RESULT})
        STRING(REGEX REPLACE "(..)/(..)/..(..).*" "\\3\\2\\1" ${RESULT} ${${RESULT}})
    ELSE (WIN32)
        MESSAGE(SEND_ERROR "date not implemented")
        SET(${RESULT} 000000)
    ENDIF (WIN32)
ENDMACRO (TODAY)

MACRO (GIT_VERSION RESULT)
  EXECUTE_PROCESS(COMMAND "git" "describe" "--abbrev=4" "--dirty" "--always" "--tags" OUTPUT_VARIABLE ${RESULT})
  string(REGEX REPLACE "\n$" "" ${RESULT} ${${RESULT}})
ENDMACRO(GIT_VERSION)
# Use the preprocessor to handle #ifdefs
SITE_NAME(HOSTNAME)
SET(USERNAME $ENV{USER})
TODAY(CURDATETIME)
GIT_VERSION(SW4_VERSION)
ADD_DEFINITIONS(-DEW_MADEBY="${USERNAME}" -DEW_OPT_LEVEL="${CMAKE_BUILD_TYPE}" -DEW_COMPILER="${CMAKE_CXX_COMPILER}" -DEW_LIBDIR="NA" -DEW_INCDIR="NA" -DEW_HOSTNAME="${HOSTNAME}" -DEW_WHEN="${CURDATETIME}" -DSW4_GIT_VERSION="${SW4_VERSION}")


SET(SW4_BASE_COMPILE_FLAGS "-DSW4_USE_UMPIRE=1 -DENABLE_PROJ4 -DENABLE_MPI_TIMING_BARRIER=1 -DSW4_STAGED_MPI_BUFFERS=1 -DUSE_DIRECT_INVERSE=1 -I $ENV{LAPACK_INC}  -DCAMP_USE_PLATFORM_DEFAULT_STREAM=1 -DSW4_USE_CMEM=1 -DSW4_GHCOF_NO_GP_IS_ZERO=1 --expt-relaxed-constexpr")
SET(CMAKE_CUDA_FLAGS  "${CMAKE_CUDA_FLAGS} ${SW4_BASE_COMPILE_FLAGS}")

#add_subdirectory(tpl/PROJ)

add_subdirectory(src)
