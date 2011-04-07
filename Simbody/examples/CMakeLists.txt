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
# We check the BUILD_TESTING_{SHARED,STATIC} variables to determine
# whether to build dynamically linked, statically linked, or both
# versions of the executable.


FILE(GLOB EXAMPLES "*.cpp")
FOREACH(EX_PROG ${EXAMPLES})
    GET_FILENAME_COMPONENT(EX_SRC  ${EX_PROG} NAME)
    GET_FILENAME_COMPONENT(EX_ROOT ${EX_PROG} NAME_WE)

    IF (BUILD_TESTING_SHARED)
        # Link with shared library
        ADD_EXECUTABLE(${EX_ROOT} ${EX_PROG})
	IF(GUI_NAME)
	    ADD_DEPENDENCIES(${EX_ROOT} ${GUI_NAME})
	ENDIF()
        SET_TARGET_PROPERTIES(${EX_ROOT}
		PROPERTIES
	      PROJECT_LABEL "Example - ${EX_ROOT}")
        TARGET_LINK_LIBRARIES(${EX_ROOT} ${TEST_SHARED_TARGET})

	INSTALL(TARGETS ${EX_ROOT} DESTINATION examples/bin)
    ENDIF (BUILD_TESTING_SHARED)

    IF (BUILD_TESTING_STATIC)
        # Link with static library
        SET(EX_STATIC ${EX_ROOT}Static)
        ADD_EXECUTABLE(${EX_STATIC} ${EX_PROG})
	IF(GUI_NAME)
	    ADD_DEPENDENCIES(${EX_STATIC} ${GUI_NAME})
	ENDIF()
        SET_TARGET_PROPERTIES(${EX_STATIC}
		PROPERTIES
		COMPILE_FLAGS "-DSimTK_USE_STATIC_LIBRARIES"
		PROJECT_LABEL "Example - ${EX_STATIC}")
        TARGET_LINK_LIBRARIES(${EX_STATIC} ${TEST_STATIC_TARGET})
	# Don't install static examples
    ENDIF (BUILD_TESTING_STATIC)

ENDFOREACH()

# Copy extra files associated with these examples
# to the binary directory that will be the working
# directory when the tests are run; and include them in
# the examples installation. ".obj" files contain geometry.
FILE(GLOB EXTRA_FILES "*.obj")
FOREACH(XTRA ${EXTRA_FILES})
    GET_FILENAME_COMPONENT(XTRA_ROOT ${XTRA} NAME)
    CONFIGURE_FILE(${XTRA} 
	           ${CMAKE_CURRENT_BINARY_DIR}/${XTRA_ROOT} COPYONLY)
    INSTALL(FILES ${XTRA} DESTINATION examples/src)
    INSTALL(FILES ${XTRA} DESTINATION examples/bin)
ENDFOREACH()

# install source for examples
FILE(GLOB SRC_FILES "*.cpp" "*.h")
INSTALL(FILES ${SRC_FILES} DESTINATION examples/src)

# install .txt files except CMakeLists.txt
FILE(GLOB TXT_FILES "*.txt")
FOREACH(TXTF ${TXT_FILES})
    IF(NOT "${TXTF}" MATCHES "CMakeLists.txt")
	INSTALL(FILES ${TXTF} DESTINATION examples/src)
    ENDIF()
ENDFOREACH()

# install the installed version of CMakeLists.txt
INSTALL(FILES InstalledCMakeLists.txt DESTINATION examples/src
	RENAME CMakeLists.txt)

# install Makefile
INSTALL(FILES "${PROJECT_SOURCE_DIR}/examples/Makefile"
    DESTINATION examples/src)

