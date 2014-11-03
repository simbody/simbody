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
# examples directory. If your example needs more than just a single .cpp
# with a single .h of the same name, then put it in its own subdirectory 
# where it can have its own CMakeLists.txt. You must list those
# subdirectories explicitly below; the simpler ones will be discovered
# automatically.
#
# Examples get built with these macros defined:
#    SIMBODY_EXAMPLE_NAME     
#    SIMBODY_EXAMPLES_INSTALL_SRC  (has trailing "/")
# Subdirectory examples should locate their auxiliary files by searching
# for one file, say "geometry/pelvis.stl", using this search path:
#    .  (current directory first)
#    ./examples/${SIMBODY_EXAMPLE_NAME}
#    ${SIMBODY_EXAMPLES_INSTALL_SRC}${SIMBODY_EXAMPLE_NAME}
#
# There is also a "shared" directory that has shared header files for
# the examples, but also may contain shared auxiliary files like common
# geometric shapes. In that case it must also be located at run time
# using the search path:
#    ../shared   (parallel to current directory first)
#    ./examples/shared
#    ${SIMBODY_EXAMPLES_INSTALL_SRC}shared
# 
# Installed shared files can always be found in
#    ${SIMBODY_EXAMPLES_INSTALL_SRC}shared
# TODO: how to access shared files from build?
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

# All subdirectory examples must be listed here.  A subdirectory is needed 
# for any example that consists of more than just one source file and one 
# header file.
set(SUBDIR_EXAMPLES
    BricardMechanism
    Gazebo2Simbody
    TaskSpaceControl-UR10
    TaskSpaceControl-Atlas
    TaskSpaceControl-Humanoid
)

IF(WIN32)
  SET(EXAMPLES_INSTALL_BIN ${CMAKE_INSTALL_DOCDIR}/examples/bin/)
ELSE()
  SET(EXAMPLES_INSTALL_BIN ${CMAKE_INSTALL_FULL_LIBDIR}/simbody/examples/)
  SET(EXAMPLES_SYMLINK_BIN ${CMAKE_INSTALL_DOCDIR}/examples/)
ENDIF()
SET(EXAMPLES_INSTALL_SRC ${CMAKE_INSTALL_DOCDIR}/examples/src/)

ADD_DEFINITIONS(
    "-DSIMBODY_EXAMPLES_INSTALL_SRC=\"${CMAKE_INSTALL_PREFIX}/${EXAMPLES_INSTALL_SRC}\"")

# Examples can depend on any Simbody includes, and may depend on headers
# in the examples/shared directory, in which case they should be
# included as #include "shared/SharedHeaderFile.h". This also permits use 
# of headers between examples, via #include "OtherExample/SomeHeader.h",
# although it is not clear that's a good idea.
INCLUDE_DIRECTORIES(${PLATFORM_INCLUDE_DIRECTORIES}
	                ${SimTKCOMMON_INCLUDE_DIRECTORIES}
	                ${SimTKMATH_INCLUDE_DIRECTORIES}
	                ${SimTKSIMBODY_INCLUDE_DIRECTORIES})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

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
foreach(EX_NAME ${SUBDIR_EXAMPLES})
    ADD_SUBDIRECTORY(${EX_NAME})
endforeach()

# Next pick up stragglers that consist solely of a .cpp file in the top
# level examples directory, with a single .h file of the same name.
# The example name is the source file name with ".cpp" removed, and
# the generated executable will have that name.

FILE(GLOB EXAMPLES "*.cpp")
FOREACH(EX_PROG ${EXAMPLES})
    GET_FILENAME_COMPONENT(EX_SRC  ${EX_PROG} NAME)
    GET_FILENAME_COMPONENT(EX_NAME ${EX_PROG} NAME_WE)
    set(EX_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/${EX_NAME}.h")
    if(NOT EXISTS ${EX_HEADER})
        unset(EX_HEADER)
    endif()

    IF (BUILD_TESTS_AND_EXAMPLES_SHARED)
        # Link with shared library
        ADD_EXECUTABLE(${EX_NAME} ${EX_PROG} ${EX_HEADER})
        IF(GUI_NAME)
            ADD_DEPENDENCIES(${EX_NAME} ${GUI_NAME})
        ENDIF()
        SET_TARGET_PROPERTIES(${EX_NAME}
        PROPERTIES
          COMPILE_FLAGS "-DSIMBODY_EXAMPLE_NAME=\"${EX_NAME}\""
          PROJECT_LABEL "Example - ${EX_NAME}")
        TARGET_LINK_LIBRARIES(${EX_NAME} ${EXAMPLES_SHARED_TARGET})
        INSTALL(TARGETS ${EX_NAME} DESTINATION ${EXAMPLES_INSTALL_BIN})
    ENDIF()

    IF (BUILD_STATIC_LIBRARIES AND BUILD_TESTS_AND_EXAMPLES_STATIC)
        # Link with static library
        SET(EX_STATIC ${EX_NAME}Static)
        ADD_EXECUTABLE(${EX_STATIC} ${EX_PROG} ${EX_HEADER})
        IF(GUI_NAME)
            ADD_DEPENDENCIES(${EX_STATIC} ${GUI_NAME})
        ENDIF()
        SET_TARGET_PROPERTIES(${EX_STATIC}
        PROPERTIES # note that the name is not EX_STATIC
          COMPILE_FLAGS "-DSIMBODY_EXAMPLE_NAME=\"${EX_NAME}\""
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

# install top-level source for examples
FILE(GLOB SRC_FILES "*.cpp" "*.h")
INSTALL(FILES ${SRC_FILES} DESTINATION ${EXAMPLES_INSTALL_SRC})

# install shared source and other files; also copy to binary directory
# so we can use them when running from the build
FILE(GLOB_RECURSE SHARED_FILES 
    RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/shared"
    "${CMAKE_CURRENT_SOURCE_DIR}/shared/*")
FOREACH(one_shared_file ${SHARED_FILES})
    CONFIGURE_FILE(shared/${one_shared_file}
       ${CMAKE_CURRENT_BINARY_DIR}/shared/${one_shared_file} COPYONLY)
ENDFOREACH()

INSTALL(DIRECTORY shared/
        DESTINATION ${EXAMPLES_INSTALL_SRC}shared) # already has "/"

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

