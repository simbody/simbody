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
# We check the BUILD_TESTING_{SHARED,STATIC} variables to determine
# whether to build dynamically linked, statically linked, or both
# versions of the executable.


FILE(GLOB ADHOC_TESTS "*.cpp")
FOREACH(TEST_PROG ${ADHOC_TESTS})
    GET_FILENAME_COMPONENT(TEST_ROOT ${TEST_PROG} NAME_WE)

    IF (BUILD_TESTING_SHARED)
        # Link with shared library
        ADD_EXECUTABLE(${TEST_ROOT} ${TEST_PROG})
        SET_TARGET_PROPERTIES(${TEST_ROOT}
		PROPERTIES
	      PROJECT_LABEL "Adhoc - ${TEST_ROOT}")
        TARGET_LINK_LIBRARIES(${TEST_ROOT} 
					${SHARED_TARGET} ${AUX_SHARED_TARGET})
    ENDIF (BUILD_TESTING_SHARED)

    IF (BUILD_TESTING_STATIC)
        # Link with static library
        SET(TEST_STATIC ${TEST_ROOT}Static)
        ADD_EXECUTABLE(${TEST_STATIC} ${TEST_PROG})
        SET_TARGET_PROPERTIES(${TEST_STATIC}
		PROPERTIES
		COMPILE_FLAGS "-DSimTK_USE_STATIC_LIBRARIES"
		PROJECT_LABEL "Adhoc - ${TEST_STATIC}")
        TARGET_LINK_LIBRARIES(${TEST_STATIC}
					${STATIC_TARGET} ${AUX_STATIC_TARGET})
    ENDIF (BUILD_TESTING_STATIC)

ENDFOREACH(TEST_PROG ${ADHOC_TESTS})