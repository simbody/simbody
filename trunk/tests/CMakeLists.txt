# Adhoc tests
SUBDIRS(BigMatrixTest)

# Old-style regression tests
# TODO: convert these to the quiet, automated regressiosn we've come
# to know and love
SUBDIRS (ListTest BNTTest MatVecTest SpatialAlgebraTest OrientationTest 
          StateTest RandomTest SFMTTest PolynomialTest 
         EnumerationTest RotationTest)


ADD_TEST(ListTest                 ${EXECUTABLE_OUTPUT_PATH}/ListTest)
ADD_TEST(ListTestStatic           ${EXECUTABLE_OUTPUT_PATH}/ListTestStatic)
ADD_TEST(BNTTest                  ${EXECUTABLE_OUTPUT_PATH}/BNTTest)
ADD_TEST(BNTTestStatic            ${EXECUTABLE_OUTPUT_PATH}/BNTTestStatic)
ADD_TEST(MatVecTest               ${EXECUTABLE_OUTPUT_PATH}/MatVecTest)
ADD_TEST(MatVecTestStatic         ${EXECUTABLE_OUTPUT_PATH}/MatVecTestStatic)
ADD_TEST(SpatialAlgebraTest      	${EXECUTABLE_OUTPUT_PATH}/SpatialAlgebraTest)
ADD_TEST(SpatialAlgebraTestStatic 	${EXECUTABLE_OUTPUT_PATH}/SpatialAlgebraTestStatic)
ADD_TEST(OrientationTest          ${EXECUTABLE_OUTPUT_PATH}/OrientationTest)
ADD_TEST(OrientationTestStatic    	${EXECUTABLE_OUTPUT_PATH}/OrientationTestStatic)
ADD_TEST(StateTest                ${EXECUTABLE_OUTPUT_PATH}/StateTest)
ADD_TEST(StateTestStatic          ${EXECUTABLE_OUTPUT_PATH}/StateTestStatic)
ADD_TEST(RandomTest               ${EXECUTABLE_OUTPUT_PATH}/RandomTest)
ADD_TEST(RandomTestStatic         ${EXECUTABLE_OUTPUT_PATH}/RandomTestStatic)
ADD_TEST(SFMTTest                 ${EXECUTABLE_OUTPUT_PATH}/SFMTTest)
ADD_TEST(SFMTTestStatic           ${EXECUTABLE_OUTPUT_PATH}/SFMTTestStatic)
ADD_TEST(PolynomialTest           ${EXECUTABLE_OUTPUT_PATH}/PolynomialTest)
ADD_TEST(PolynomialTestStatic     	${EXECUTABLE_OUTPUT_PATH}/PolynomialTestStatic)
ADD_TEST(EnumerationTest          ${EXECUTABLE_OUTPUT_PATH}/EnumerationTest)
ADD_TEST(EnumerationTestStatic    	${EXECUTABLE_OUTPUT_PATH}/EnumerationTestStatic)
ADD_TEST(RotationTest             ${EXECUTABLE_OUTPUT_PATH}/RotationTest)
ADD_TEST(RotationTestStatic       	${EXECUTABLE_OUTPUT_PATH}/RotationTestStatic)

# Regression tests

# Thanks to Chris Bruns for this next part. Use this only for
# actual CTest tests, that is, quiet tests that know the right
# answers and either succeed or fail.

# Automatically create tests using files named "Test*.cpp"
FILE(GLOB TEST_PROGS "Test*.cpp")
FOREACH(TEST_PROG ${TEST_PROGS})
    GET_FILENAME_COMPONENT(TEST_ROOT ${TEST_PROG} NAME_WE)

    # Link with shared library
    ADD_EXECUTABLE(${TEST_ROOT} ${TEST_PROG})
    TARGET_LINK_LIBRARIES(${TEST_ROOT} ${SHARED_TARGET})
    ADD_TEST(${TEST_ROOT} ${EXECUTABLE_OUTPUT_PATH}/${TEST_ROOT})

    # Link with static library
    SET(TEST_STATIC ${TEST_ROOT}Static)
    ADD_EXECUTABLE(${TEST_STATIC} ${TEST_PROG})
    SET_TARGET_PROPERTIES(${TEST_STATIC}
		PROPERTIES
		COMPILE_FLAGS "-DSimTK_USE_STATIC_LIBRARIES"
		)
    TARGET_LINK_LIBRARIES(${TEST_STATIC} ${STATIC_TARGET})
    ADD_TEST(${TEST_STATIC} ${EXECUTABLE_OUTPUT_PATH}/${TEST_STATIC})

ENDFOREACH(TEST_PROG ${TEST_PROGS})


