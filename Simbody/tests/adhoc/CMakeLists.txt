# Generate adhoc tests.
#
# This is boilerplate code for generating a set of executables, one per
# .cpp file in a "tests/adhoc" subdirectory. These are for "loose" test
# cases, not suitable as regression tests but perhaps useful to the
# developer or as demos. Unlike the similar boilerplate code in the "tests"
# directory, this does not generate CMake ADD_TEST commands so these
# will never run automatically.
#
# For IDEs that can deal with PROJECT_LABEL properties (at least
# Visual Studio) the projects for building each of these adhoc
# executables will be labeled "Adhoc - TheTestName" if a file
# TheTestName.cpp is found in this directory.
#
# We check the BUILD_TESTS_AND_EXAMPLES_{SHARED,STATIC} variables to determine
# whether to build dynamically linked, statically linked, or both
# versions of the executable.


file(GLOB ADHOC_TESTS "*.cpp")
foreach(TEST_PROG ${ADHOC_TESTS})
    get_filename_component(TEST_ROOT ${TEST_PROG} NAME_WE)

    if(BUILD_TESTS_AND_EXAMPLES_SHARED)
        # Link with shared library
        add_executable(${TEST_ROOT} ${TEST_PROG})
    if(GUI_NAME)
        add_dependencies(${TEST_ROOT} ${GUI_NAME})
    endif()
        set_target_properties(${TEST_ROOT}
        PROPERTIES
          PROJECT_LABEL "Test_Adhoc - ${TEST_ROOT}")
        target_link_libraries(${TEST_ROOT} 
                      ${TEST_SHARED_TARGET})
    endif(BUILD_TESTS_AND_EXAMPLES_SHARED)

    if(BUILD_STATIC_LIBRARIES AND BUILD_TESTS_AND_EXAMPLES_STATIC)
        # Link with static library
        set(TEST_STATIC ${TEST_ROOT}Static)
        add_executable(${TEST_STATIC} ${TEST_PROG})
    if(GUI_NAME)
        add_dependencies(${TEST_STATIC} ${GUI_NAME})
    endif()
        set_target_properties(${TEST_STATIC}
        PROPERTIES
        COMPILE_FLAGS "-DSimTK_USE_STATIC_LIBRARIES"
        PROJECT_LABEL "Test_Adhoc - ${TEST_STATIC}")
        target_link_libraries(${TEST_STATIC}
                      ${TEST_STATIC_TARGET})
    endif()

endforeach(TEST_PROG ${ADHOC_TESTS})

file(GLOB ADHOC_OBJS "*.obj" "*.vtp")
foreach(TEST_OBJ ${ADHOC_OBJS})
    get_filename_component(OBJ_ROOT ${TEST_OBJ} NAME)
    configure_file(${TEST_OBJ} 
            ${CMAKE_CURRENT_BINARY_DIR}/${OBJ_ROOT} COPYONLY)
endforeach(TEST_OBJ ${ADHOC_OBJS})
