# This set of tests was written by Radu Serban to test CPodes within
# the Sundials framework. They do not use the SimTK C++ interface at all.

# Adhoc tests are those test or demo programs which are not intended,
# or not ready, to be part of the regression suite.
#add_subdirectory(adhoc)

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

file(GLOB REGR_TESTS "*.c")
foreach(TEST_PROG ${REGR_TESTS})
    get_filename_component(TEST_ROOT ${TEST_PROG} NAME_WE)

    add_executable(${TEST_ROOT} ${TEST_PROG})
    set_target_properties(${TEST_ROOT} PROPERTIES
        FOLDER "Tests/SimTKmath/CPodes/sundials" # Organize target list in IDE
        PROJECT_LABEL "Test_Regr - ${TEST_ROOT}")
    target_link_libraries(${TEST_ROOT} 
                SimTKmath)
    target_include_directories(${TEST_ROOT} PRIVATE
        ${PROJECT_SOURCE_DIR}/Integrators/src/CPodes/sundials/include
        ${SimTK_SDK}/include)
    add_test(NAME ${TEST_ROOT}
        COMMAND $<TARGET_FILE:${TEST_ROOT}>)

endforeach()
