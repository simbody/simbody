SUBDIRS (ListTest BNTTest MatVecTest SpatialAlgebraTest OrientationTest 
         BigMatrixTest StateTest RandomTest SFMTTest PolynomialTest 
         EnumerationTest RotationTest)

# Thanks to Chris Bruns for this next part. Use this only for
# actual CTest tests, that is, quiet tests that know the right
# answers and either succeed or fail.
# TODO: convert the above tests.

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


