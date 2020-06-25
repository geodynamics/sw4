# The name of our project is "SW4". CMakeLists files in this project can 
# refer to the root source directory of the project as ${SW4_SOURCE_DIR} and 
# to the root binary directory of the project as ${SW4_BINARY_DIR}. 
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(SW4 C CXX Fortran)
ENABLE_LANGUAGE(Fortran)

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

# Allow a shorter version of the PROJ4 directory definition
IF (NOT PROJ4_OSGEO4W_HOME)
    SET(PROJ4_OSGEO4W_HOME ${PROJ4_HOME})
ENDIF (NOT PROJ4_OSGEO4W_HOME)

# Find proj4 installation
FIND_PACKAGE(Proj4)
IF (PROJ4_FOUND)
    ADD_DEFINITIONS(-DENABLE_PROJ4)
    INCLUDE_DIRECTORIES(${PROJ4_INCLUDE_DIR})
    # If proj4 is installed, we can also look for cencalvm
    FIND_PACKAGE(CENCALVM)
    IF (CENCALVM_FOUND)
        ADD_DEFINITIONS(-DENABLE_ETREE)
        INCLUDE_DIRECTORIES(${CENCALVM_INCLUDE_DIR})
    ENDIF (CENCALVM_FOUND)
ENDIF (PROJ4_FOUND)

# Find BLAS and LAPACK installation
FIND_PACKAGE(BLAS REQUIRED)
FIND_PACKAGE(LAPACK REQUIRED)

# Determine how to mangle names between Fortran and C
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
    src/Image3D.C src/MaterialVolimagefile.C src/SfileOutput.C
    src/MaterialRfile.C src/MaterialSfile.C src/AnisotropicMaterialBlock.C src/sacutils.C
    src/consintp.C src/addmemvarforcing2.C src/ESSI3D.C src/ESSI3DHDF5.C  src/sachdf5.C src/readhdf5.C)

SET(SW4_SRC_F_FILES src/addsgd.f src/bcfort.f 
    src/randomfield3d.f src/testsrc.f src/twilightfort.f
    src/boundaryOp.f src/bndryOpNoGhost.f90 src/rayleighfort.f src/twilightsgfort.f
    src/curvilinear4.f src/curvilinear4sg.f 
    src/innerloop-ani-sgstr-vc.f src/checkanisomtrl.f src/computedtaniso.f
    src/rhs4curvilinear.f src/velsum.f src/gradients.f
    src/rhs4curvilinearsg.f src/energy4.f 
    src/rhs4th3fort.f src/lamb_exact_numquad.f src/solerr3.f
    src/bcfortanisg.f src/ilanisocurv.f
    src/anisomtrltocurvilinear.f src/bcfreesurfcurvani.f 
    src/tw_aniso_force_tt.f src/tw_aniso_force.f src/tw_ani_stiff.f90
    src/rhs4th3fortwind.f src/scalar_prod.f90 src/addsg4wind.f90 src/updatememvar.f90)

ADD_EXECUTABLE(sw4 ${SW4_SRC_C_FILES} ${SW4_SRC_F_FILES} ${SW4_QUADPACK_SRC_FILES})
TARGET_LINK_LIBRARIES(sw4 ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES}
                        ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}
                        ${PROJ4_LIBRARIES} ${CENCALVM_LIBRARIES})

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

