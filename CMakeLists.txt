# The name of our project is "SW4". CMakeLists files in this project can 
# refer to the root source directory of the project as ${SW4_SOURCE_DIR} and 
# to the root binary directory of the project as ${SW4_BINARY_DIR}. 
CMAKE_MINIMUM_REQUIRED(VERSION 3.12)
PROJECT(SW4 C CXX Fortran)
ENABLE_LANGUAGE(Fortran)
CMAKE_POLICY(SET CMP0028 NEW) # Double colon in target name means ALIAS or IMPORTED target.

# Add extra CMake module path (is it still neeed?)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/configs/")

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
ADD_DEFINITIONS(-DEW_MADEBY="${USERNAME}" -DEW_OPT_LEVEL="${CMAKE_BUILD_TYPE}" -DEW_COMPILER="${CMAKE_CXX_COMPILER}" -DEW_LIBDIR="NA" -DEW_INCDIR="NA" -DEW_HOSTNAME="${HOSTNAME}" -DEW_WHEN="${CURDATETIME}")

# Look for MPI and add the required flags if found
FIND_PACKAGE(MPI REQUIRED)
ADD_DEFINITIONS(${MPI_C_COMPILE_FLAGS})
INCLUDE_DIRECTORIES(${MPI_C_INCLUDE_PATH})

# Allow a shorter version of the PROJ directory definition
IF (NOT PROJ_OSGEO4W_HOME)
    SET(PROJ_OSGEO4W_HOME ${PROJ_HOME})
ENDIF (NOT PROJ_OSGEO4W_HOME)

OPTION(USE_PROJ "Use PROJ." OFF)
IF(USE_PROJ)
    FIND_PACKAGE(PROJ REQUIRED)
    IF (PROJ_FOUND)
        INCLUDE_DIRECTORIES(${PROJ_INCLUDE_DIR})
        if(PROJ_VERSION VERSION_LESS "6.0.0")
            MESSAGE(SEND_ERROR "PROJ library must be version 6 or higher")
        else()
            ADD_DEFINITIONS(-DENABLE_PROJ)
        endif()
    ENDIF (PROJ_FOUND)
ENDIF (USE_PROJ)

OPTION(USE_HDF5 "Use HDF5." OFF)
IF(USE_HDF5)
    FIND_PACKAGE(HDF5 REQUIRED)
    IF (HDF5_FOUND)
        ADD_DEFINITIONS(-DUSE_HDF5)
        INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
    ENDIF (HDF5_FOUND)
ENDIF (USE_HDF5)

OPTION(USE_ZFP "Use ZFP compression." OFF)
IF(USE_ZFP)
    SET(H5Z_ZFP_USE_STATIC_LIBS ON)
    FIND_PACKAGE(H5Z_ZFP 1.1.0 CONFIG)
    ADD_DEFINITIONS(-DUSE_ZFP)
ENDIF (USE_ZFP)

OPTION(USE_SZ "Use SZ compression." OFF)
IF(USE_SZ)
    FIND_PACKAGE(SZ REQUIRED)
    IF (SZ_FOUND)
        ADD_DEFINITIONS(-DUSE_SZ)
        INCLUDE_DIRECTORIES(${SZ_INCLUDE_DIR})
    ENDIF (SZ_FOUND)
ENDIF (USE_SZ)

OPTION(USE_FFTW3 "Use FFTW3." OFF)
IF(USE_FFTW3)
    FIND_PACKAGE(FFTW3 REQUIRED)
    IF (FFTW3_FOUND)
        ADD_DEFINITIONS(-DUSE_FFTW)
        INCLUDE_DIRECTORIES(${FFTW3_INCLUDE_DIR})
    ENDIF (FFTW3_FOUND)
ENDIF (USE_FFTW3)

# Find BLAS and LAPACK installation
FIND_PACKAGE(BLAS REQUIRED)
FIND_PACKAGE(LAPACK REQUIRED)
FIND_PACKAGE(OpenMP)

# Determine how to mangle names between Fortran and C
INCLUDE(FortranCInterface)
FortranCInterface_HEADER(FC.h MACRO_NAMESPACE "FC_")
INCLUDE_DIRECTORIES(${SW4_BINARY_DIR})
ADD_DEFINITIONS(-DCMAKE_FC_FUNC="1")
INCLUDE_DIRECTORIES(${SW4_SOURCE_DIR}/src/double)
INCLUDE_DIRECTORIES(${SW4_SOURCE_DIR}/src)

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
    src/TimeSeries.C src/sacsubc.C src/SuperGrid.C src/TestRayleighWave.C src/MaterialPfile.C
    src/Filter.C src/Polynomial.C src/SecondOrderSection.C src/time_functions.C src/Qspline.C
    src/MaterialIfile.C src/GeographicProjection.C src/Image3D.C src/ESSI3D.C src/ESSI3DHDF5.C
    src/MaterialVolimagefile.C src/MaterialRfile.C src/MaterialSfile.C src/AnisotropicMaterialBlock.C
    src/sacutils.C src/DataPatches.C src/addmemvarforcing2.C src/consintp.C src/oddIoddJinterp.C
    src/evenIoddJinterp.C src/oddIevenJinterp.C src/evenIevenJinterp.C src/CheckPoint.C src/geodyn.C
    src/AllDims.C src/Patch.C src/RandomizedMaterial.C src/MaterialInvtest.C src/sw4-prof.C
    src/sachdf5.C src/readhdf5.C src/TestTwilight.C src/TestPointSource.C src/curvilinear4sgwind.C
    src/TestEcons.C src/GridGenerator.C src/GridGeneratorGeneral.C src/GridGeneratorGaussianHill.C
    src/CurvilinearInterface2.C src/SfileOutput.C src/pseudohess.C src/fastmarching.C src/solveTT.C
    src/rhs4th3point.C src/MaterialGMG.C src/addsgdc.C src/bcfortc.C src/bcfortanisgc.C 
    src/bcfreesurfcurvanic.C src/boundaryOpc.C src/energy4c.C src/checkanisomtrlc.C src/computedtanisoc.C 
    src/curvilinear4sgc.C src/gradientsc.C src/randomfield3dc.C src/innerloop-ani-sgstr-vcc.C 
    src/ilanisocurvc.C src/rhs4curvilinearc.C src/rhs4curvilinearsgc.C src/rhs4th3fortc.C src/solerr3c.C 
    src/testsrcc.C src/rhs4th3windc.C src/tw_aniso_forcec.C src/tw_aniso_force_ttc.C src/velsumc.C 
    src/twilightfortc.C src/twilightsgfortc.C src/tw_ani_stiffc.C src/anisomtrltocurvilinearc.C 
    src/scalar_prodc.C src/updatememvarc.C src/addsg4windc.C src/bndryOpNoGhostc.C src/rhs4th3windc2.C
    src/Sarray.h src/double/sw4.h)


SET(SW4_SRC_F_FILES src/rayleighfort.f src/lamb_exact_numquad.f)

ADD_EXECUTABLE(sw4 ${SW4_SRC_C_FILES} ${SW4_SRC_F_FILES} ${SW4_QUADPACK_SRC_FILES})
TARGET_LINK_LIBRARIES(sw4 PUBLIC ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES}
                        ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}
                        ${PROJ_LIBRARIES} ${HDF5_LIBRARIES})

IF(OpenMP_FOUND)
    TARGET_LINK_LIBRARIES(sw4 PRIVATE OpenMP::OpenMP_CXX)
ENDIF (OpenMP_FOUND)

IF (USE_FFTW3)
    TARGET_LINK_LIBRARIES(sw4 PRIVATE fftw3)
ENDIF (USE_FFTW3)

IF (USE_PROJ)
    TARGET_LINK_LIBRARIES(sw4 PRIVATE PROJ::proj)
ENDIF (USE_PROJ)

IF(USE_ZFP)
    TARGET_LINK_LIBRARIES(sw4 PRIVATE h5z_zfp::h5z_zfp)
ENDIF (USE_ZFP)

# Show the banner once sw4 is built
ADD_CUSTOM_COMMAND(TARGET sw4 POST_BUILD COMMAND cat ${SW4_SOURCE_DIR}/wave.txt)

# Level of tests to run, corresponding to computational time
# Level 0 tests should take just less than a minute each, level 1 tests take less than 10 minutes (includes level 0 tests)
# and level 2 consists of all tests
SET(TESTING_LEVEL 0 CACHE STRING "Level of testing thoroughness: 0=short tests (<1 minute each), 2=all tests (may take an hour)")
SET(MPI_NUM_TEST_PROCS 4 CACHE STRING "Number of MPI processes to use in tests")

SET(PYTEST_DIR ${SW4_SOURCE_DIR}/pytest)
SET(REF_DIR ${SW4_SOURCE_DIR}/pytest/reference)

# Define tests and tolerances
SET(TEST_DIRS       twilight               twilight         twilight
                    twilight               twilight         twilight
                    attenuation            attenuation      attenuation
                    attenuation            attenuation      attenuation
		    meshrefine             meshrefine       meshrefine
		    meshrefine             meshrefine       meshrefine
                    lamb                   lamb             lamb
                    pointsource            pointsource      pointsource)

SET(RESULT_SUBDIRS  flat-twi-1             flat-twi-2       flat-twi-3
                    gauss-twi-1            gauss-twi-2      gauss-twi-3
                    tw-att-1               tw-att-2         tw-att-3
                    tw-topo-att-1          tw-topo-att-2    tw-topo-att-3
                    refine-el-1            refine-att-1     refine-att-2nd-1
                    refine-el-2            refine-att-2     refine-att-2nd-2
                    lamb-1                 lamb-2           lamb-3
                    pointsource-sg-1       pointsource-sg-2 pointsource-sg-3)

SET(TEST_IN_FILES   flat-twi-1.in          flat-twi-2.in       flat-twi-3.in
                    gauss-twi-1.in         gauss-twi-2.in      gauss-twi-3.in
                    tw-att-1.in            tw-att-2.in         tw-att-3.in
                    tw-topo-att-1.in       tw-topo-att-2.in    tw-topo-att-3.in
		    refine-el-1.in         refine-att-1.in     refine-att-2nd-1.in
		    refine-el-2.in         refine-att-2.in     refine-att-2nd-2.in
                    lamb-1.in              lamb-2.in           lamb-3.in
                    pointsource-sg-1.in    pointsource-sg-2.in pointsource-sg-3.in)

SET(TEST_BASE_FILES TwilightErr            TwilightErr             TwilightErr    
                    TwilightErr            TwilightErr             TwilightErr    
                    TwilightErr            TwilightErr             TwilightErr    
                    TwilightErr            TwilightErr             TwilightErr    
                    TwilightErr            TwilightErr             TwilightErr    
                    TwilightErr            TwilightErr             TwilightErr
                    LambErr                LambErr                 LambErr
                    PointSourceErr         PointSourceErr          PointSourceErr)

SET(TEST_CHECKS     compare             compare                 compare 
                    compare             compare                 compare 
                    compare             compare                 compare 
                    compare             compare                 compare 
                    compare             compare                 compare 
		    compare             compare                 compare
                    compare             compare                 compare
                    compare             compare                 compare)

SET(TEST_LEVELS     0                       0                    1
                    0                       0                    1
                    0                       0                    1
                    0                       1                    2
		    0                       0                    0
		    1                       1                    1
                    0                       1                    2
                    0                       1                    2)

SET(TEST_ERRINF_TOL 1e-5                    1e-5                 1e-5
                    1e-5                    1e-5                 1e-5
                    1e-5                    1e-5                 1e-5
                    1e-5                    1e-5                 1e-5
		    1e-5                    1e-5                 1e-5
		    1e-5                    1e-5                 1e-5
                    1e-5                    1e-5                 1e-5
                    1e-5                    1e-5                 1e-5)

SET(TEST_ERRL2_TOL  1e-5                    1e-5                 1e-5
                    1e-5                    1e-5                 1e-5
                    1e-5                    1e-5                 1e-5
                    1e-5                    1e-5                 1e-5
		    1e-5                    1e-5                 1e-5
		    1e-5                    1e-5                 1e-5
                    1e-5                    1e-5                 1e-5
                    1e-5                    1e-5                 1e-5)

SET(TEST_SOLINF_TOL 1e-2                    1e-2                 1e-2
                    1e-2                    1e-2                 1e-2
                    1e-2                    1e-2                 1e-2
                    1e-2                    1e-2                 1e-2
		    1e-2                    1e-2                 1e-2
		    1e-2                    1e-2                 1e-2
                    0                       0                    0
                    0                       0                    0)

LIST(LENGTH TEST_DIRS N)
MATH(EXPR NUM_TESTS "${N}-1")

# Run through and register all tests within the current testing level
FOREACH(TEST_IND RANGE ${NUM_TESTS})
    LIST(GET TEST_LEVELS ${TEST_IND} TEST_LEVEL)
    IF (NOT ${TESTING_LEVEL} LESS ${TEST_LEVEL})
        LIST(GET TEST_DIRS ${TEST_IND} TEST_DIR)
        LIST(GET RESULT_SUBDIRS ${TEST_IND} RESULT_SUBDIR)
        SET(Test_Name ${TEST_DIR}/${RESULT_SUBDIR})
        LIST(GET TEST_IN_FILES ${TEST_IND} TEST_IN_FILE)
        LIST(GET TEST_BASE_FILES ${TEST_IND} TEST_BASE_FILE)
        LIST(GET TEST_CHECKS ${TEST_IND} TEST_CHECK_TYPE)
        LIST(GET TEST_ERRINF_TOL ${TEST_IND} ERRINF_TOL)
        LIST(GET TEST_ERRL2_TOL ${TEST_IND} ERRL2_TOL)
        LIST(GET TEST_SOLINF_TOL ${TEST_IND} SOLINF_TOL)
        ADD_TEST(
            NAME Run_${Test_Name}
            WORKING_DIRECTORY ${SW4_BINARY_DIR}
            COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPI_NUM_TEST_PROCS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_POSTFLAGS}
                ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/sw4 ${REF_DIR}/${TEST_DIR}/${TEST_IN_FILE}
            )

        SET(TEST_REF_FILE ${REF_DIR}/${TEST_DIR}/${RESULT_SUBDIR}/${TEST_BASE_FILE}.txt)
        SET(TEST_OUT_FILE ${SW4_BINARY_DIR}/${RESULT_SUBDIR}/${TEST_BASE_FILE}.txt)
        ADD_TEST(
            NAME Check_Result_${Test_Name}
            WORKING_DIRECTORY ${SW4_BINARY_DIR}/${TEST_OUTPUT_DIR}
            COMMAND ${PYTEST_DIR}/check_results.py ${TEST_CHECK_TYPE} ${TEST_REF_FILE} ${TEST_OUT_FILE} ${ERRINF_TOL} ${ERRL2_TOL} ${SOLINF_TOL}
            )
        SET_TESTS_PROPERTIES(
            Check_Result_${Test_Name}
            PROPERTIES DEPENDS ${Test_Name}
            )
    ENDIF (NOT ${TESTING_LEVEL} LESS ${TEST_LEVEL})
ENDFOREACH(TEST_IND RANGE ${NUM_AUTO_TESTS})

