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
# Tests are statically or dynamically linked to the library targets as-built,
# depending on BUILD_SHARED_LIBS


file(GLOB ADHOC_TESTS "*.cpp")
foreach(TEST_PROG ${ADHOC_TESTS})
    get_filename_component(TEST_ROOT ${TEST_PROG} NAME_WE)

    add_executable(${TEST_ROOT} ${TEST_PROG})
    if(TARGET simbody-visualizer)
        add_dependencies(${TEST_ROOT} simbody-visualizer)
    endif()
    set_target_properties(${TEST_ROOT} PROPERTIES
        PROJECT_LABEL "Test_Adhoc - ${TEST_ROOT}")

    # Link with library (shared or static determined by BUILD_SHARED_LIBS)
    target_link_libraries(${TEST_ROOT} SimTKsimbody)

endforeach()

file(GLOB ADHOC_OBJS "*.obj" "*.vtp")
foreach(TEST_OBJ ${ADHOC_OBJS})
    get_filename_component(OBJ_ROOT ${TEST_OBJ} NAME)
    configure_file(${TEST_OBJ} 
            ${CMAKE_CURRENT_BINARY_DIR}/${OBJ_ROOT} COPYONLY)
endforeach()
