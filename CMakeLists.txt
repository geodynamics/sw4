# The name of our project is "SW4". CMakeLists files in this project can 
# refer to the root source directory of the project as ${SW4_SOURCE_DIR} and 
# to the root binary directory of the project as ${SW4_BINARY_DIR}. 
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(SW4 C CXX Fortran)
ENABLE_LANGUAGE(Fortran)

# Enable testing framework
INCLUDE(CTest)

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

# If the user doesn't specify, build the release version by default
IF (NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release)
ENDIF (NOT CMAKE_BUILD_TYPE)

SET(CMAKE_Fortran_MODULE_DIRECTORY "${SW4_BINARY_DIR}/include")
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${SW4_BINARY_DIR}/bin")
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${SW4_BINARY_DIR}/lib")

# Use the preprocessor to handle #ifdefs
SITE_NAME(HOSTNAME)
SET(USERNAME $ENV{USER})
TODAY(CURDATETIME)
ADD_DEFINITIONS(-DEW_MADEBY="${USERNAME}" -DEW_OPT_LEVEL="${CMAKE_BUILD_TYPE}" -DEW_COMPILER="${CMAKE_CXX_COMPILER}" -DEW_LIBDIR="???" -DEW_INCDIR="???" -DEW_HOSTNAME="${HOSTNAME}" -DEW_WHEN="${CURDATETIME}")

# Look for MPI and add the required flags if found
FIND_PACKAGE(MPI REQUIRED)
ADD_DEFINITIONS(${MPI_C_COMPILE_FLAGS})
INCLUDE_DIRECTORIES(${MPI_C_INCLUDE_PATH})

# Find BLAS and LAPACK installation
FIND_PACKAGE(BLAS REQUIRED)
FIND_PACKAGE(LAPACK REQUIRED)

INCLUDE(FortranCInterface)
FortranCInterface_HEADER(FC.h MACRO_NAMESPACE "FC_")
INCLUDE_DIRECTORIES(${SW4_BINARY_DIR})
ADD_DEFINITIONS(-DCMAKE_FC_FUNC="1")

# Check validity of different flags
INCLUDE(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG(-O0 COMPILER_SUPPORTS_OPT0_FLAG)
CHECK_CXX_COMPILER_FLAG(-O3 COMPILER_SUPPORTS_OPT3_FLAG)
CHECK_CXX_COMPILER_FLAG(-g COMPILER_SUPPORTS_DEBUG_FLAG)

# Define debug vs release compiler flags
IF(COMPILER_SUPPORTS_OPT0_FLAG)
    SET(DEBUG_FLAGS "${DEBUG_FLAGS} -O0")
ENDIF(COMPILER_SUPPORTS_OPT0_FLAG)

IF(COMPILER_SUPPORTS_OPT3_FLAG)
    SET(RELEASE_FLAGS "${RELEASE_FLAGS} -O3")
    SET(RELWITHDEBINFO_FLAGS "${RELWITHDEBINFO_FLAGS} -O3")
ENDIF(COMPILER_SUPPORTS_OPT3_FLAG)

IF(COMPILER_SUPPORTS_DEBUG_FLAG)
    SET(DEBUG_FLAGS "${DEBUG_FLAGS} -g")
    SET(RELWITHDEBINFO_FLAGS "${RELWITHDEBINFO_FLAGS} -g")
ENDIF(COMPILER_SUPPORTS_DEBUG_FLAG)

SET(CMAKE_CXX_FLAGS_DEBUG ${DEBUG_FLAGS})
SET(CMAKE_CXX_FLAGS_RELEASE ${RELEASE_FLAGS})
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO ${RELWITHDEBINFO_FLAGS})

# Define source files
FILE(GLOB SW4_QUADPACK_SRC_FILES src/quadpack/*.f)

SET(SW4_SRC_C_FILES src/main.C src/EW.C src/Sarray.C src/version.C src/parseInputFile.C 
    src/ForcingTwilight.C src/curvilinearGrid.C src/parallelStuff.C src/Source.C 
    src/MaterialProperty.C src/MaterialData.C src/material.C src/setupRun.C 
    src/solve.C src/Parallel_IO.C src/Image.C src/GridPointSource.C src/MaterialBlock.C
    src/TimeSeries.C src/sacsubc.C src/SuperGrid.C
    src/TestRayleighWave.C src/MaterialPfile.C src/Filter.C src/Polynomial.C
    src/SecondOrderSection.C src/time_functions.C src/Qspline.C
    src/EtreeFile.C src/MaterialIfile.C src/GeographicProjection.C 
    src/Image3D.C src/MaterialVolimagefile.C
    src/MaterialRfile.C 
    src/AnisotropicMaterialBlock.C )

SET(SW4_SRC_F_FILES src/addsgd.f src/bcfort.f 
    src/randomfield3d.f src/testsrc.f src/twilightfort.f
    src/boundaryOp.f src/rayleighfort.f src/twilightsgfort.f
    src/curvilinear4.f src/curvilinear4sg.f 
    src/innerloop-ani-sgstr-vc.f src/checkanisomtrl.f src/computedtaniso.f
    src/rhs4curvilinear.f src/velsum.f src/gradients.f
    src/rhs4curvilinearsg.f src/energy4.f 
    src/rhs4th3fort.f src/lamb_exact_numquad.f src/solerr3.f src/bcfortanisg.f )

# Unused:
# bcsf.f
# exactRhs.f
# lamb_one_point.f
# facc.f
# test.f
# test3.f
# ffsurf.f
# bctwilight.f
# fg.f
# rayleigh.f
# boundaryDisc.f
# fgt.f
# fgtt.f
# consinterface.f
# erhs3.f

ADD_EXECUTABLE(sw4 ${SW4_SRC_C_FILES} ${SW4_SRC_F_FILES} ${SW4_QUADPACK_SRC_FILES})
TARGET_LINK_LIBRARIES(sw4 ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

# Show the banner once sw4 is built
ADD_CUSTOM_COMMAND(TARGET sw4 POST_BUILD COMMAND cat ${SW4_SOURCE_DIR}/wave.txt)

SET(EXAMPLE_DIR ${SW4_SOURCE_DIR}/examples)

# Run the shorter twilight tests
ADD_TEST(
    NAME twilight_flat-twi-1
    WORKING_DIRECTORY ${SW4_BINARY_DIR}
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 8 ${MPIEXEC_PREFLAGS} ${MPIEXEC_POSTFLAGS}
            ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/sw4 ${EXAMPLE_DIR}/twilight/flat-twi-1.in
    )

ADD_TEST(
    NAME twilight_flat-twi-2
    WORKING_DIRECTORY ${SW4_BINARY_DIR}
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 8 ${MPIEXEC_PREFLAGS} ${MPIEXEC_POSTFLAGS}
            ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/sw4 ${EXAMPLE_DIR}/twilight/flat-twi-2.in
    )

ADD_TEST(
    NAME twilight_gauss-twi-1
    WORKING_DIRECTORY ${SW4_BINARY_DIR}
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 8 ${MPIEXEC_PREFLAGS} ${MPIEXEC_POSTFLAGS}
            ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/sw4 ${EXAMPLE_DIR}/twilight/gauss-twi-1.in
    )

ADD_TEST(
    NAME twilight_gauss-twi-2
    WORKING_DIRECTORY ${SW4_BINARY_DIR}
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 8 ${MPIEXEC_PREFLAGS} ${MPIEXEC_POSTFLAGS}
            ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/sw4 ${EXAMPLE_DIR}/twilight/gauss-twi-2.in
    )

