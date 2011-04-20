# Adhoc tests are those test or demo programs which are not intended,
# or not ready, to be part of the regression suite.
ADD_SUBDIRECTORY(adhoc)

# Generate regression tests.
#
# This is boilerplate code for generating a set of executables, one per
# .cpp file in a "tests" directory. These are for test cases which 
# have been engineered to be suitable as part of the nightly regression
# test suite. Ideally, they know their correct answers and return
# a simple thumbs-up/thumbs-down result, although some regressions 
# may be present and useful just to determine if they compile and
# run to completion.
#
# For IDEs that can deal with PROJECT_LABEL properties (at least
# Visual Studio) the projects for building each of these adhoc
# executables will be labeled "Regr - TheTestName" if a file
# TheTestName.cpp is found in this directory.
#
# We check the BUILD_TESTING_{SHARED,STATIC} variables to determine
# whether to build dynamically linked, statically linked, or both
# versions of the executable.

FILE(GLOB REGR_TESTS "*.cpp")
FOREACH(TEST_PROG ${REGR_TESTS})
    GET_FILENAME_COMPONENT(TEST_ROOT ${TEST_PROG} NAME_WE)

    IF (BUILD_TESTING_SHARED)
        # Link with shared library
        ADD_EXECUTABLE(${TEST_ROOT} ${TEST_PROG})
        SET_TARGET_PROPERTIES(${TEST_ROOT}
		PROPERTIES
	      PROJECT_LABEL "Test_Regr - ${TEST_ROOT}")
        TARGET_LINK_LIBRARIES(${TEST_ROOT} 
					${TEST_SHARED_TARGET})
        ADD_TEST(${TEST_ROOT} ${EXECUTABLE_OUTPUT_PATH}/${TEST_ROOT})
    ENDIF (BUILD_TESTING_SHARED)

    IF (BUILD_STATIC_LIBRARIES AND BUILD_TESTING_STATIC)
        # Link with static library
        SET(TEST_STATIC ${TEST_ROOT}Static)
        ADD_EXECUTABLE(${TEST_STATIC} ${TEST_PROG})
        SET_TARGET_PROPERTIES(${TEST_STATIC}
		PROPERTIES
		COMPILE_FLAGS "-DSimTK_USE_STATIC_LIBRARIES"
		PROJECT_LABEL "Test_Regr - ${TEST_STATIC}")
        TARGET_LINK_LIBRARIES(${TEST_STATIC}
					${TEST_STATIC_TARGET})
        ADD_TEST(${TEST_STATIC} ${EXECUTABLE_OUTPUT_PATH}/${TEST_STATIC})
    ENDIF()

ENDFOREACH(TEST_PROG ${REGR_TESTS})

