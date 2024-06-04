# Adhoc tests are those test or demo programs which are not intended,
# or not ready, to be part of the regression suite.
add_subdirectory(adhoc)

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
# Tests are statically or dynamically linked to the library targets as-built,
# depending on BUILD_SHARED_LIBS

# Some tests here share headers. We don't know the actual dependencies
# so we'll just make each test depend on all the headers which is
# at least conservative.
file(GLOB ALL_HEADERS "*.h")

file(GLOB REGR_TESTS "*.cpp")
foreach(TEST_PROG ${REGR_TESTS})
    get_filename_component(TEST_ROOT ${TEST_PROG} NAME_WE)

    add_executable(${TEST_ROOT} ${TEST_PROG} ${ALL_HEADERS})
    set_target_properties(${TEST_ROOT} PROPERTIES
        PROJECT_LABEL "Test_Regr - ${TEST_ROOT}")

    # Link with library (shared or static determined by BUILD_SHARED_LIBS)
    target_link_libraries(${TEST_ROOT} SimTKmath)
    add_test(${TEST_ROOT} ${EXECUTABLE_OUTPUT_PATH}/${TEST_ROOT})

endforeach()

