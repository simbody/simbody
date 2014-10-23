# Build and install examples.
#
# This is for use when building Simbody from source; it is not for
# end users to build examples themselves. This guarantees that the 
# examples built with Simbody are using the libraries just built regardless
# of what is installed. We install a different CMakeLists.txt along with
# the examples for end users to use.
#
# Examples are built from an explicit list of subdirectories given below,
# and then simpler examples are built one per .cpp file in the top level
# examples directory. If your example needs more than just a single .cpp,
# put it in its own subdirectory where it can have its own CMakeLists.txt.
#
# Note that unlike code in the "tests" directory (but like the 
# "tests/adhoc" code), this does not generate CMake ADD_TEST commands so 
# these will never run automatically.
#
# For IDEs that can deal with PROJECT_LABEL properties (at least
# Visual Studio) the projects for building each of these adhoc
# executables will be labeled "Example - TheExampleName" if a file
# TheExampleName.cpp is found in this directory, or if we build from 
# a subdirectory called TheExampleName.
#
# We check the BUILD_TESTS_AND_EXAMPLES_{SHARED,STATIC} variables to 
# determine whether to build dynamically linked, statically linked, or both
# versions of the executable.
IF(WIN32)
  SET(EXAMPLES_INSTALL_BIN ${CMAKE_INSTALL_DOCDIR}/examples/bin/)
ELSE()
  SET(EXAMPLES_INSTALL_BIN ${CMAKE_INSTALL_FULL_LIBDIR}/simbody/examples/)
  SET(EXAMPLES_SYMLINK_BIN ${CMAKE_INSTALL_DOCDIR}/examples/)
ENDIF()
SET(EXAMPLES_INSTALL_SRC ${CMAKE_INSTALL_DOCDIR}/examples/src/)


# Examples can depend on any Simbody includes.
INCLUDE_DIRECTORIES(${PLATFORM_INCLUDE_DIRECTORIES}
	                ${SimTKCOMMON_INCLUDE_DIRECTORIES}
	                ${SimTKMATH_INCLUDE_DIRECTORIES}
	                ${SimTKSIMBODY_INCLUDE_DIRECTORIES})

## Build examples against the unversioned libraries if they are being 
## built; otherwise against the versioned libraries.
IF(BUILD_UNVERSIONED_LIBRARIES)
	SET(EXAMPLES_SHARED_TARGET ${SimTKSIMBODY_LIBRARY_NAME})
	SET(EXAMPLES_STATIC_TARGET ${SimTKSIMBODY_LIBRARY_NAME}_static)
ELSE()
	SET(EXAMPLES_SHARED_TARGET ${SimTKSIMBODY_LIBRARY_NAME}${VN})
	SET(EXAMPLES_STATIC_TARGET ${SimTKSIMBODY_LIBRARY_NAME}${VN}_static)
ENDIF()


# First build examples that have their own subdirectories and CMakeLists
# files.
ADD_SUBDIRECTORY(BricardMechanism)
ADD_SUBDIRECTORY(Gazebo2Simbody)
ADD_SUBDIRECTORY(TaskSpaceControl-UR10)
ADD_SUBDIRECTORY(TaskSpaceControl-Atlas)
ADD_SUBDIRECTORY(TaskSpaceControl-Humanoid)

# Next pick up stragglers that consist solely of a .cpp file in the top
# level examples directory.

FILE(GLOB EXAMPLES "*.cpp")
FOREACH(EX_PROG ${EXAMPLES})
    GET_FILENAME_COMPONENT(EX_SRC  ${EX_PROG} NAME)
    GET_FILENAME_COMPONENT(EX_ROOT ${EX_PROG} NAME_WE)

    IF (BUILD_TESTS_AND_EXAMPLES_SHARED)
        # Link with shared library
        ADD_EXECUTABLE(${EX_ROOT} ${EX_PROG})
        IF(GUI_NAME)
            ADD_DEPENDENCIES(${EX_ROOT} ${GUI_NAME})
        ENDIF()
        SET_TARGET_PROPERTIES(${EX_ROOT}
        PROPERTIES
          PROJECT_LABEL "Example - ${EX_ROOT}")
        TARGET_LINK_LIBRARIES(${EX_ROOT} ${EXAMPLES_SHARED_TARGET})
        INSTALL(TARGETS ${EX_ROOT} DESTINATION ${EXAMPLES_INSTALL_BIN})
    ENDIF()

    IF (BUILD_STATIC_LIBRARIES AND BUILD_TESTS_AND_EXAMPLES_STATIC)
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
        TARGET_LINK_LIBRARIES(${EX_STATIC} ${EXAMPLES_STATIC_TARGET})
        # Don't install static examples
    ENDIF()
ENDFOREACH()

IF(NOT WIN32)
  EXECUTE_PROCESS(COMMAND "${CMAKE_COMMAND}" -E create_symlink
                  ${EXAMPLES_INSTALL_BIN}
                  ${CMAKE_CURRENT_BINARY_DIR}/bin
  )
  INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin
          DESTINATION ${EXAMPLES_SYMLINK_BIN}
  )
ENDIF()

# install top-level source for examples (shouldn't be any headers though)
FILE(GLOB SRC_FILES "*.cpp" "*.h")
INSTALL(FILES ${SRC_FILES} DESTINATION ${EXAMPLES_INSTALL_SRC})

# install shared source and other files
FILE(GLOB SHARED_FILES "shared/*")
INSTALL(FILES ${SHARED_FILES} DESTINATION ${EXAMPLES_INSTALL_SRC}/shared)

# install .txt files except CMakeLists.txt
FILE(GLOB TXT_FILES "*.txt")
FOREACH(TXTF ${TXT_FILES})
    IF(NOT "${TXTF}" MATCHES "CMakeLists.txt")
    INSTALL(FILES ${TXTF} DESTINATION ${EXAMPLES_INSTALL_SRC})
    ENDIF()
ENDFOREACH()

# install the installed version of CMakeLists.txt
INSTALL(FILES InstalledCMakeLists.txt DESTINATION ${EXAMPLES_INSTALL_SRC}
    RENAME CMakeLists.txt)

# install Makefile
INSTALL(FILES "${PROJECT_SOURCE_DIR}/examples/Makefile"
    DESTINATION ${EXAMPLES_INSTALL_SRC})

