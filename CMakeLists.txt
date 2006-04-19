#---------------------------------------------------
# Simbody 
#
# Creates SimTK Core library, base name=simbody.
# Default libraries are shared & optimized. Variants
# are created for static (_static) and debug (_d).
#
# Windows:
#   simbody[_d].dll
#   simbody[_d].lib
#   simbody_static[_d].lib
# Unix:
#   libsimbody[_d].so
#   libsimbody_static[_d].a
#
# All libraries are installed in 
#   C:\Program Files\SimTK\core\lib  (Windows)
#   /usr/local/SimTK/core/lib        (UNIX)
#
# TODO: should we just install the default libraries
# at the top level, and install the static & debug
# libraries in .../core/lib/simbody/lib
# or           .../simbody/lib ?
#----------------------------------------------------

PROJECT (Simbody)
SUBDIRS (staticTarget sharedTarget tests)

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
SET(SIMBODY_SOURCE_SUBDIRS .)

# Collect up information about the version of the simbody library we're building
# and make it available to the code so it can be built into the binaries.

SET(SIMBODY_LIBRARY_NAME simbody)
SET(SIMBODY_MAJOR_VERSION 0)
SET(SIMBODY_MINOR_VERSION 5)
SET(SIMBODY_BUILD_VERSION 1)

SET(SIMBODY_COPYRIGHT_YEARS "2005-6")
SET(SIMBODY_AUTHORS         "Michael Sherman") # comma separated list of authors


# Get the Subversion revision number
EXEC_PROGRAM ( svnversion
                ARGS  ${PROJECT_SOURCE_DIR}
                OUTPUT_VARIABLE SIMBODY_SVN_REVISION )

ADD_DEFINITIONS(-DSimTK_SIMBODY_LIBRARY_NAME=${SIMBODY_LIBRARY_NAME}
                -DSimTK_SIMBODY_MAJOR_VERSION=${SIMBODY_MAJOR_VERSION}
                -DSimTK_SIMBODY_MINOR_VERSION=${SIMBODY_MINOR_VERSION}
		-DSimTK_SIMBODY_BUILD_VERSION=${SIMBODY_BUILD_VERSION})

# CMake quotes automatically when building Visual Studio projects but we need
# to add them ourselves for Linux or Cygwin. Two cases to avoid duplicate quotes
# in Visual Studio which end up in the binary.
IF(UNIX)
   ADD_DEFINITIONS(-DSimTK_SIMBODY_SVN_REVISION="${SIMBODY_SVN_REVISION}"
                   -DSimTK_SIMBODY_COPYRIGHT_YEARS="${SIMBODY_COPYRIGHT_YEARS}"
                   -DSimTK_SIMBODY_AUTHORS="${SIMBODY_AUTHORS}")
ELSE(UNIX)
   ADD_DEFINITIONS(-DSimTK_SIMBODY_SVN_REVISION=${SIMBODY_SVN_REVISION}
                   -DSimTK_SIMBODY_COPYRIGHT_YEARS=${SIMBODY_COPYRIGHT_YEARS}
                   -DSimTK_SIMBODY_AUTHORS=${SIMBODY_AUTHORS})
ENDIF(UNIX)

# -DSimTK_SIMBODY_TYPE has to be defined in the target subdirectories.
# -Dsimbody_EXPORTS defined automatically when Windows DLL build is being done.

# Report the version number to the CMake UI
SET(SIMBODY_VERSION
    "${SIMBODY_MAJOR_VERSION}.${SIMBODY_MINOR_VERSION}.${SIMBODY_BUILD_VERSION}"
     CACHE STRING "This is the version of Simbody which will be built." FORCE)

SET(SHARED_TARGET ${SIMBODY_LIBRARY_NAME})
SET(STATIC_TARGET ${SIMBODY_LIBRARY_NAME}_static)

## If no one says otherwise, change the executable path to drop into the same binary
## location as the DLLs so that the test cases will use the just-build DLLs.
IF(NOT EXECUTABLE_OUTPUT_PATH)
  SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR} 
      CACHE INTERNAL "Single output directory for building all executables.")
ENDIF(NOT EXECUTABLE_OUTPUT_PATH)
IF(NOT LIBRARY_OUTPUT_PATH)
  SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR} 
      CACHE INTERNAL "Single output directory for building all libraries.")
ENDIF(NOT LIBRARY_OUTPUT_PATH)
SET(${PROJECT_NAME}_EXECUTABLE_DIR ${EXECUTABLE_OUTPUT_PATH}/${CMAKE_CFG_INTDIR})
SET(${PROJECT_NAME}_LIBRARY_DIR    ${LIBRARY_OUTPUT_PATH}/${CMAKE_CFG_INTDIR})

# Ensure that debug libraries have "_d" appended to their names.
# CMake gets this right on Windows automatically with this definition.
SET(CMAKE_DEBUG_POSTFIX "_d" CACHE INTERNAL "" FORCE)

# But on Unix or Cygwin we have to add the suffix manually
IF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)
    SET(SHARED_TARGET ${SHARED_TARGET}${CMAKE_DEBUG_POSTFIX})
    SET(STATIC_TARGET ${STATIC_TARGET}${CMAKE_DEBUG_POSTFIX})
ENDIF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)


# Choose a LAPACK to link with.
SET(SimTK_LAPACK SimTK CACHE STRING "What LAPACK to link with (ACML or SimTK)")
IF(NOT SimTK_LAPACK)
    SET(SimTK_LAPACK SimTK)
ENDIF(NOT SimTK_LAPACK)
STRING(TOLOWER ${SimTK_LAPACK} SimTK_LAPACK_DOWNSHIFT)


# Find the right place to look for dependency resolution (which
# can be the default install directory or specified SimTK_SDK directory).
#
# Also, determine which math libraries to use for this platform.
#

IF(${CMAKE_C_COMPILER} MATCHES "gcc")

    ## A gcc build
    IF (SimTK_LAPACK_DOWNSHIFT STREQUAL "acml")
        SET(MATH_LIBS_TO_USE acml)
        ADD_DEFINITIONS(-DSimTK_USE_ACML_LAPACK)
    ELSE (SimTK_LAPACK_DOWNSHIFT STREQUAL "acml")
        SET(MATH_LIBS_TO_USE SimTKlapack)
    ENDIF (SimTK_LAPACK_DOWNSHIFT STREQUAL "acml")

    SET(MATH_LIBS_TO_USE ${MATH_LIBS_TO_USE} c m g2c)
    LINK_DIRECTORIES(/usr/local/lib /lib)

ELSE(${CMAKE_C_COMPILER} MATCHES "gcc")

    ## Assume Microsoft Visual Studio
    IF (SimTK_LAPACK_DOWNSHIFT STREQUAL "acml")
        SET(MATH_LIBS_TO_USE libacml_dll)
        ADD_DEFINITIONS(-DSimTK_USE_ACML_LAPACK)
    ELSE (SimTK_LAPACK_DOWNSHIFT STREQUAL "acml")
        SET(MATH_LIBS_TO_USE SimTKlapack)
    ENDIF (SimTK_LAPACK_DOWNSHIFT STREQUAL "acml")

ENDIF(${CMAKE_C_COMPILER} MATCHES "gcc")


# These are all the places to search for header files which are
# to be part of the API.
SET(API_INCLUDE_DIRS) # start empty
FOREACH(subdir ${SIMBODY_SOURCE_SUBDIRS})
    # append
    SET(API_INCLUDE_DIRS ${API_INCLUDE_DIRS}
                         ${subdir}/include 
                         ${subdir}/include/simbody 
                         ${subdir}/include/simbody/internal)
ENDFOREACH(subdir)

# We'll need both *relative* path names, starting with their API_INCLUDE_DIRS,
# and absolute pathnames.
SET(API_REL_INCLUDE_FILES)   # start these out empty
SET(API_ABS_INCLUDE_FILES)

FOREACH(dir ${API_INCLUDE_DIRS})
    FILE(GLOB fullpaths ${dir}/*.h)	# returns full pathnames
    SET(API_ABS_INCLUDE_FILES ${API_ABS_INCLUDE_FILES} ${fullpaths})

    FOREACH(pathname ${fullpaths})
        GET_FILENAME_COMPONENT(filename ${pathname} NAME)
        SET(API_REL_INCLUDE_FILES ${API_REL_INCLUDE_FILES} ${dir}/${filename})
    ENDFOREACH(pathname)
ENDFOREACH(dir)

# collect up source files
SET(SOURCE_FILES) # empty
SET(SOURCE_INCLUDE_FILES)

FOREACH(subdir ${SIMBODY_SOURCE_SUBDIRS})
    FILE(GLOB src_files  ${subdir}/src/*.cpp)
    FILE(GLOB incl_files ${subdir}/src/*.h)
    SET(SOURCE_FILES         ${SOURCE_FILES}         ${src_files})   #append
    SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})

    INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/${subdir}/include)
ENDFOREACH(subdir)

#
# Allow automated build and dashboard.
#
INCLUDE (Dart)

IF (UNIX AND NOT CYGWIN)
  IF (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE MATCHES Debug)
    ADD_DEFINITIONS(-fprofile-arcs -ftest-coverage)
    LINK_LIBRARIES(gcov)
  ENDIF (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE MATCHES Debug)
ENDIF (UNIX AND NOT CYGWIN)


#
# Testing
#

ENABLE_TESTING()

IF (EXECUTABLE_OUTPUT_PATH)
  SET (TEST_PATH ${EXECUTABLE_OUTPUT_PATH})
ELSE (EXECUTABLE_OUTPUT_PATH)
  SET (TEST_PATH .)
ENDIF (EXECUTABLE_OUTPUT_PATH)

ADD_TEST(SBPendulum1              ${TEST_PATH}/SBPendulum1)
ADD_TEST(SBPendulum1Static        ${TEST_PATH}/SBPendulum1Static)
## add more

#
# Installation
#
# CMake installs go into CMAKE_INSTALL_PREFIX, which is always
# set to something incorrect by CMake. If SimTK_INSTALL_PREFIX is set then
# it is a sandbox installation area, otherwise we want to install
# in /usr/local/SimTK or c:\program files\SimTK. Similarly, SimTK_SDK
# tells us where to find other SimTK Core modules on which this one
# depends. If it is not set, we'll look in /usr/local/SimTK or
# c:\program files\SimTK *regardless* of SimTK_INSTALL_PREFIX.

IF(${CMAKE_C_COMPILER} MATCHES "gcc")

    IF(NOT SimTK_INSTALL_PREFIX)
        SET(SimTK_INSTALL_PREFIX "/usr/local/SimTK")
    ENDIF(NOT SimTK_INSTALL_PREFIX)

    # unix or cygwin; cygwin requires some trickery
    IF(CYGWIN)
        SET(SimTK_INSTALL_PREFIX C:/cygwin${SimTK_INSTALL_PREFIX})
    ENDIF(CYGWIN)

    IF(NOT SimTK_SDK)
        SET(SimTK_SDK "/usr/local/SimTK")
    ENDIF(NOT SimTK_SDK)

ELSE(${CMAKE_C_COMPILER} MATCHES "gcc")

    IF(NOT SimTK_INSTALL_PREFIX)
        SET(SimTK_INSTALL_PREFIX "C:/Program Files/SimTK")
    ENDIF(NOT SimTK_INSTALL_PREFIX)

    IF(NOT SimTK_SDK)
        SET(SimTK_SDK "C:/Program Files/SimTK")
    ENDIF(NOT SimTK_SDK)

ENDIF(${CMAKE_C_COMPILER} MATCHES "gcc")

SET(CMAKE_INSTALL_PREFIX ${SimTK_INSTALL_PREFIX} CACHE STRING "Install path prefix." FORCE)

INCLUDE_DIRECTORIES(${SimTK_SDK}/core/include)
LINK_DIRECTORIES(${SimTK_SDK}/core/lib)

# libraries are installed from their subdirectories; headers here

# install headers
FILE(GLOB CORE_HEADERS     include/*.h                  */include/*.h)
FILE(GLOB TOP_HEADERS      include/simbody/*.h          */include/simbody/*.h)
FILE(GLOB INTERNAL_HEADERS include/simbody/internal/*.h */include/simbody/internal/*.h)
INSTALL_FILES(/core/include/                 FILES ${CORE_HEADERS})
INSTALL_FILES(/core/include/simbody/         FILES ${TOP_HEADERS})
INSTALL_FILES(/core/include/simbody/internal FILES ${INTERNAL_HEADERS})
