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
# We check the BUILD_TESTS_AND_EXAMPLES_{SHARED,STATIC} variables to determine
# whether to build dynamically linked, statically linked, or both
# versions of the executable.

file(GLOB REGR_TESTS "*.cpp")
foreach(TEST_PROG ${REGR_TESTS})
    get_filename_component(TEST_ROOT ${TEST_PROG} NAME_WE)

    if(BUILD_TESTS_AND_EXAMPLES_SHARED)
        # Link with shared library
        add_executable(${TEST_ROOT} ${TEST_PROG})
        set_target_properties(${TEST_ROOT}
        PROPERTIES
          PROJECT_LABEL "Test_Regr - ${TEST_ROOT}")
        target_link_libraries(${TEST_ROOT}
                    ${TEST_SHARED_TARGET})
        add_test(${TEST_ROOT} ${EXECUTABLE_OUTPUT_PATH}/${TEST_ROOT})
    endif(BUILD_TESTS_AND_EXAMPLES_SHARED)

    if(BUILD_STATIC_LIBRARIES AND BUILD_TESTS_AND_EXAMPLES_STATIC)
        # Link with static library
        set(TEST_STATIC ${TEST_ROOT}Static)
        add_executable(${TEST_STATIC} ${TEST_PROG})
        set_target_properties(${TEST_STATIC}
        PROPERTIES
        COMPILE_FLAGS "-DSimTK_USE_STATIC_LIBRARIES"
        PROJECT_LABEL "Test_Regr - ${TEST_STATIC}")
        target_link_libraries(${TEST_STATIC}
                    ${TEST_STATIC_TARGET})
        add_test(${TEST_STATIC} ${EXECUTABLE_OUTPUT_PATH}/${TEST_STATIC})
    endif()

endforeach(TEST_PROG ${REGR_TESTS})




