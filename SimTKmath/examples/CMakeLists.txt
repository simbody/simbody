
#FILE( GLOB SOURCES ${PROJECT_SOURCE_DIR}/examples/*.cpp )

#FOREACH (src_file ${SOURCES} )
#   GET_FILENAME_COMPONENT( file_name ${src_file} NAME )
#   INSTALL_FILES( "/core/examples/simmath/" FILES ${file_name} )
#ENDFOREACH (src_file ${SOURCES})

#INSTALL_FILES( "/core/examples/simmath/" FILES
# ${PROJECT_SOURCE_DIR}/examples/README.txt)

# Generate and install examples.
#
# This is boilerplate code for generating a set of executables, one per
# .cpp file in an "examples" subdirectory. These are intended to
# be installed with the library but we don't handle installation
# here. Unlike the similar boilerplate code in the "tests"
# directory (but like the "tests/adhoc" boilerplate), this does
# not generate CMake ADD_TEST commands so these will never run
# automatically.
#
# For IDEs that can deal with PROJECT_LABEL properties (at least
# Visual Studio) the projects for building each of these adhoc
# executables will be labeled "Example - TheExampleName" if a file
# TheExampleName.cpp is found in this directory.
#
# We check the BUILD_TESTS_AND_EXAMPLES_{SHARED,STATIC} variables to determine
# whether to build dynamically linked, statically linked, or both
# versions of the executable.


SET(EXAMPLES_INSTALL ${CMAKE_INSTALL_DOCDIR}/examples/simmath/)
FILE(GLOB EXAMPLES "*.cpp")
FOREACH(EX_PROG ${EXAMPLES})
    GET_FILENAME_COMPONENT(EX_SRC  ${EX_PROG} NAME)
    GET_FILENAME_COMPONENT(EX_ROOT ${EX_PROG} NAME_WE)

    IF (BUILD_TESTS_AND_EXAMPLES_SHARED)
        # Link with shared library
        ADD_EXECUTABLE(${EX_ROOT} ${EX_PROG})
        SET_TARGET_PROPERTIES(${EX_ROOT}
		PROPERTIES
	      PROJECT_LABEL "Example - ${EX_ROOT}")
        TARGET_LINK_LIBRARIES(${EX_ROOT} ${TEST_SHARED_TARGET})
    ENDIF (BUILD_TESTS_AND_EXAMPLES_SHARED)

    IF (BUILD_STATIC_LIBRARIES AND BUILD_TESTS_AND_EXAMPLES_STATIC)
        # Link with static library
        SET(EX_STATIC ${EX_ROOT}Static)
        ADD_EXECUTABLE(${EX_STATIC} ${EX_PROG})
        SET_TARGET_PROPERTIES(${EX_STATIC}
		PROPERTIES
		COMPILE_FLAGS "-DSimTK_USE_STATIC_LIBRARIES"
		PROJECT_LABEL "Example - ${EX_STATIC}")
        TARGET_LINK_LIBRARIES(${EX_STATIC} ${TEST_STATIC_TARGET})
    ENDIF()

    INSTALL(FILES ${EX_SRC} DESTINATION ${EXAMPLES_INSTALL})

ENDFOREACH(EX_PROG ${ADHOC_TESTS})

INSTALL(FILES README.txt DESTINATION ${EXAMPLES_INSTALL})

