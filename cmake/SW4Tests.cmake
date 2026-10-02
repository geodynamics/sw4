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
            PROPERTIES DEPENDS Run_${Test_Name}
            )
    ENDIF (NOT ${TESTING_LEVEL} LESS ${TEST_LEVEL})
ENDFOREACH(TEST_IND RANGE ${NUM_AUTO_TESTS})

