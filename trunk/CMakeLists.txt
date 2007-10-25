#---------------------------------------------------
# Simmath 
#
# Creates SimTK Core library, base name=simmath.
# Default libraries are shared & optimized. Variants
# are created for static (_static) and debug (_d).
#
# Windows:
#   simmath[_d].dll
#   simmath[_d].lib
#   simmath_static[_d].lib
# Unix:
#   libsimmath[_d].so
#   libsimmath_static[_d].a
#
# All libraries are installed in 
#   C:\Program Files\SimTK\core\lib  (Windows)
#   /usr/local/SimTK/core/lib        (UNIX)
#
# TODO: should we just install the default libraries
# at the top level, and install the static & debug
# libraries in .../core/lib/simmath/lib
# or           .../simmath/lib ?
#----------------------------------------------------

PROJECT ( SimTKsimmath )
SUBDIRS (staticTarget sharedTarget Integrator tests examples)

# The build system will set ARCH64 for 64 bit builds, which require
# use of the lib64/ library directories rather than lib/.
SET( ARCH64 OFF CACHE BOOL "ON for 64bit builds, OFF for 32bit builds")
MARK_AS_ADVANCED( ARCH64  )

IF(${CMAKE_C_COMPILER} MATCHES "gcc")

   IF (ARCH64)
       SET( M32 -m64 )
       SET(LIB64 64)
   ELSE (ARCH64)
       SET( M32 -m32 )
       SET(LIB64) # nothing
   ENDIF (ARCH64)

   IF( WIN32 )
      SET( NO_CYGWIN -mno-cygwin )
   ENDIF( WIN32 )

   SET( CMAKE_CXX_FLAGS ${M32} )
   SET( CMAKE_SHARED_LINKER_FLAGS "${M32} ${NO_CYGWIN}" )
   SET( CMAKE_EXE_LINKER_FLAGS    "${M32} ${NO_CYGWIN}" )

ENDIF(${CMAKE_C_COMPILER} MATCHES "gcc")

IF(UNIX AND NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Debug CACHE STRING "Debug or Release build" FORCE)
ENDIF (UNIX AND NOT CMAKE_BUILD_TYPE)

IF (NOT CMAKE_CXX_FLAGS_DEBUG)
    SET(CMAKE_CXX_FLAGS_DEBUG "-g" CACHE STRING "To use when CMAKE_BUILD_TYPE=Debug" FORCE)
ENDIF (NOT CMAKE_CXX_FLAGS_DEBUG)

IF (NOT CMAKE_CXX_FLAGS_RELEASE)
    SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING 
          "To use when CMAKE_BUILD_TYPE=Release" FORCE)
ENDIF (NOT CMAKE_CXX_FLAGS_RELEASE)

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
SET(SIMMATH_SOURCE_SUBDIRS  src Ipopt include include/simmath include/simmath/internal 
    Integrator/src Integrator/include Integrator/include/simmath Integrator/include/simmath/internal)


# Collect up information about the version of the simmath library we're building
# and make it available to the code so it can be built into the binaries.

SET(SIMMATH_LIBRARY_NAME SimTKmath)
SET(SIMMATH_MAJOR_VERSION 0)
SET(SIMMATH_MINOR_VERSION 7)
SET(SIMMATH_BUILD_VERSION 1)

SET(SIMMATH_COPYRIGHT_YEARS "2005-6")

# underbar separated list of dotted authors, no spaces or commas
SET(SIMMATH_AUTHORS         "Jack.Middleton_Michael.Sherman")

# Get the subversion revision number if we can
# It's possible that WIN32 installs use svnversion through cygwin
# so we'll try for both svnversion.exe and svnversion. Note that
# this will result in warnings if all you have is Tortoise without
# Cygwin, and your "about" string will say "unknown" rather than
# providing the SVN version of the source.


FIND_PROGRAM (SVNVERSION_EXE svnversion.exe)
IF (SVNVERSION_EXE)
    EXEC_PROGRAM (${SVNVERSION_EXE}
                  ARGS  \"${PROJECT_SOURCE_DIR}\"
                  OUTPUT_VARIABLE SIMMATH_SVN_REVISION )
 ELSE (SVNVERSION_EXE) 
    FIND_PROGRAM (SVNVERSION svnversion) 
    IF (SVNVERSION)
        EXEC_PROGRAM (${SVNVERSION}
                      ARGS ${PROJECT_SOURCE_DIR}
                      OUTPUT_VARIABLE SIMMATH_SVN_REVISION)
    ELSE (SVNVERSION)
        MESSAGE (STATUS 
          "Could not find 'svnversion' executable; 'about' will be wrong. (Cygwin provides one on Windows.)"
        )
        SET (SIMMATH_SVN_REVISION unknown) 
    ENDIF (SVNVERSION)
ENDIF (SVNVERSION_EXE)
# Remove colon from build version, for easier placement in directory names
STRING(REPLACE ":" "_" SIMMATH_SVN_REVISION ${SIMMATH_SVN_REVISION})

ADD_DEFINITIONS(-DSimTK_SIMMATH_LIBRARY_NAME=${SIMMATH_LIBRARY_NAME}
                -DSimTK_SIMMATH_MAJOR_VERSION=${SIMMATH_MAJOR_VERSION}
                -DSimTK_SIMMATH_MINOR_VERSION=${SIMMATH_MINOR_VERSION}
		-DSimTK_SIMMATH_BUILD_VERSION=${SIMMATH_BUILD_VERSION})

# CMake quotes automatically when building Visual Studio projects but we need
# to add them ourselves for Linux or Cygwin. Two cases to avoid duplicate quotes
# in Visual Studio which end up in the binary.

IF (${CMAKE_GENERATOR} MATCHES "Visual Studio")
   SET(NEED_QUOTES FALSE)
ELSE (${CMAKE_GENERATOR} MATCHES "Visual Studio")
   SET(NEED_QUOTES TRUE)
ENDIF (${CMAKE_GENERATOR} MATCHES "Visual Studio")

##TODO: doesn't work without quotes in nightly build
SET(NEED_QUOTES TRUE)

IF(NEED_QUOTES)
   ADD_DEFINITIONS(-DSimTK_SIMMATH_SVN_REVISION="${SIMMATH_SVN_REVISION}"
                   -DSimTK_SIMMATH_COPYRIGHT_YEARS="${SIMMATH_COPYRIGHT_YEARS}"
                   -DSimTK_SIMMATH_AUTHORS="${SIMMATH_AUTHORS}")
ELSE(NEED_QUOTES)
   ADD_DEFINITIONS(-DSimTK_SIMMATH_SVN_REVISION=${SIMMATH_SVN_REVISION}
                   -DSimTK_SIMMATH_COPYRIGHT_YEARS=${SIMMATH_COPYRIGHT_YEARS}
                   -DSimTK_SIMMATH_AUTHORS=${SIMMATH_AUTHORS})
ENDIF(NEED_QUOTES)

# -DSimTK_SIMMATH_TYPE has to be defined in the target subdirectories.
# -Dsimmath_EXPORTS defined automatically when Windows DLL build is being done.

# Report the version number to the CMake UI
SET(SIMMATH_VERSION
    "${SIMMATH_MAJOR_VERSION}.${SIMMATH_MINOR_VERSION}.${SIMMATH_BUILD_VERSION}"
     CACHE STRING "This is the version of Simmath which will be built." FORCE)


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

SET(SHARED_TARGET ${SIMMATH_LIBRARY_NAME})
SET(STATIC_TARGET ${SIMMATH_LIBRARY_NAME}_static)

# Ensure that debug libraries have "_d" appended to their names.
# CMake gets this right on Windows automatically with this definition.
IF (${CMAKE_GENERATOR} MATCHES "Visual Studio")
    SET(CMAKE_DEBUG_POSTFIX "_d" CACHE INTERNAL "" FORCE)
ENDIF (${CMAKE_GENERATOR} MATCHES "Visual Studio")

# But on Unix or Cygwin we have to add the suffix manually
IF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)
    SET(SHARED_TARGET ${SHARED_TARGET}_d)
    SET(STATIC_TARGET ${STATIC_TARGET}_d)
ENDIF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)

SET(SimTKCOMMON_LIBRARY_NAME SimTKcommon)
SET(SimTKCOMMON_SHARED_LIBRARY ${SimTKCOMMON_LIBRARY_NAME})
SET(SimTKCOMMON_STATIC_LIBRARY ${SimTKCOMMON_LIBRARY_NAME}_static)

SET(SimTKCPODES_LIBRARY_NAME SimTKcpodes)
SET(SimTKCPODES_SHARED_LIBRARY ${SimTKCPODES_LIBRARY_NAME})
SET(SimTKCPODES_STATIC_LIBRARY ${SimTKCPODES_LIBRARY_NAME}_static)

# Choose a LAPACK to link with.
SET(SimTK_LAPACK SimTK CACHE STRING "What LAPACK to link with (ACML or SimTK)")
IF(NOT SimTK_LAPACK)
    SET(SimTK_LAPACK SimTK)
ENDIF(NOT SimTK_LAPACK)
STRING(TOLOWER ${SimTK_LAPACK} SimTK_LAPACK_DOWNSHIFT)

#
# Find the right place to look for dependency resolution (which
# can be the default install directory or specified SimTK_SDK directory).
#
# Also, determine which math libraries to use for this platform.
#

IF(${CMAKE_C_COMPILER} MATCHES "gcc")

    ## A gcc build
   SET(MATH_LIBS_TO_USE SimTKlapack pthread dl)

 ##   LINK_DIRECTORIES(/usr/local/lib${LIB64} /lib${LIB64})

ELSE(${CMAKE_C_COMPILER} MATCHES "gcc")

    ## Assume Microsoft Visual Studio
    SET(MATH_LIBS_TO_USE SimTKlapack)

ENDIF(${CMAKE_C_COMPILER} MATCHES "gcc")


# These are all the places to search for header files which are
# to be part of the API.

SET(API_INCLUDE_DIRS ${PROJECT_SRC_DIR}/include ${PROJECT_SRC_DIR}/Integrator/include) 

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

FOREACH(subdir ${SIMMATH_SOURCE_SUBDIRS})
    FILE(GLOB src_files  ${subdir}/*.cpp)
    FILE(GLOB incl_files ${subdir}/*.h)
    SET(SOURCE_FILES         ${SOURCE_FILES}         ${src_files})   #append
    SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})

    INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/${subdir})
ENDFOREACH(subdir)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/Integrator/include)

IF( APPLE )
  INCLUDE_DIRECTORIES( "/usr/include/malloc")
ENDIF( APPLE )

#
# Allow automated build and dashboard.
#
INCLUDE (Dart)

#IF (UNIX AND NOT CYGWIN AND NOT APPLE)
#  IF (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE MATCHES Debug)
#    ADD_DEFINITIONS(-fprofile-arcs -ftest-coverage)
#    LINK_LIBRARIES(gcov)
#  ENDIF (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE MATCHES Debug)
#ENDIF (UNIX AND NOT CYGWIN AND NOT APPLE)


#
# Testing
#

ENABLE_TESTING()

IF (EXECUTABLE_OUTPUT_PATH)
  SET (TEST_PATH ${EXECUTABLE_OUTPUT_PATH})
ELSE (EXECUTABLE_OUTPUT_PATH)
  SET (TEST_PATH .)
ENDIF (EXECUTABLE_OUTPUT_PATH)

ADD_TEST(IpoptTest                      ${TEST_PATH}/IpoptTest)
ADD_TEST(IpoptTestStatic                ${TEST_PATH}/IpoptTestStatic)
ADD_TEST(LBFGSTest                      ${TEST_PATH}/LBFGSTest)
ADD_TEST(LBFGSTestStatic                ${TEST_PATH}/LBFGSTestStatic)
ADD_TEST(LBFGSBTest                     ${TEST_PATH}/LBFGSBTest)
ADD_TEST(LBFGSBTestStatic               ${TEST_PATH}/LBFGSBTestStatic)
ADD_TEST(IpoptDiffTest                  ${TEST_PATH}/IpoptDiffTest)
ADD_TEST(IpoptDiffTestStatic            ${TEST_PATH}/IpoptDiffTestStatic)
ADD_TEST(LBFGSDiffTest                  ${TEST_PATH}/LBFGSDiffTest)
ADD_TEST(LBFGSDiffTestStatic            ${TEST_PATH}/LBFGSDiffTestStatic)
ADD_TEST(LBFGSBDiffTest                 ${TEST_PATH}/LBFGSBDiffTest)
ADD_TEST(LBFGSBDiffTestStatic           ${TEST_PATH}/LBFGSBDiffTestStatic)
ADD_TEST(DifferentiatorTest             ${TEST_PATH}/DifferentiatorTest)
ADD_TEST(DifferentiatorTestStatic       ${TEST_PATH}/DifferentiatorTestStatic)
ADD_TEST(SimpleDifferentiatorTest       ${TEST_PATH}/SimpleDifferentiatorTest)
ADD_TEST(SimpleDifferentiatorTestStatic ${TEST_PATH}/SimpleDifferentiatorTestStatic )
ADD_TEST(IntegratorTest			${TEST_PATH}/IntegratorTest)
ADD_TEST(IntegratorTestStatic		${TEST_PATH}/IntegratorTestStatic)
ADD_TEST(TimeStepperTest		${TEST_PATH}/TimeStepperTest)
ADD_TEST(TimeStepperTestStatic		${TEST_PATH}/TimeStepperTestStatic)
ADD_TEST(CPodesIntegratorTest		${TEST_PATH}/CPodesIntegratorTest)
ADD_TEST(CPodesIntegratorTestStatic	${TEST_PATH}/CPodesIntegratorTestStatic)

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
        SET(SimTK_INSTALL_PREFIX "/usr/local/SimTK" CACHE PATH "Install Directory")
    ENDIF(NOT SimTK_INSTALL_PREFIX)

    # unix or cygwin; cygwin requires some trickery
    IF(CYGWIN)
        SET(SimTK_INSTALL_PREFIX C:/cygwin${SimTK_INSTALL_PREFIX} CACHE PATH "Install Directory")
    ENDIF(CYGWIN)

    IF(NOT SimTK_SDK)
        SET(SimTK_SDK "/usr/local/SimTK" CACHE PATH "SimTK Dependency path")
    ENDIF(NOT SimTK_SDK)

ELSE(${CMAKE_C_COMPILER} MATCHES "gcc")

    IF(NOT SimTK_INSTALL_PREFIX)
        SET(SimTK_INSTALL_PREFIX "C:/Program Files/SimTK" CACHE PATH "Install Directory")
#        SET(SimTK_INSTALL_PREFIX "C:/Documents and Settings/jacklm/SimTK" CACHE PATH "Install Directory")
    ENDIF(NOT SimTK_INSTALL_PREFIX)

    IF(NOT SimTK_SDK)
        SET(SimTK_SDK "C:/Program Files/SimTK" CACHE PATH "SimTK Dependency path")
#        SET(SimTK_SDK "C:/Documents and Settings/jacklm/SimTK" CACHE PATH "SimTK Dependency path")
    ENDIF(NOT SimTK_SDK)

ENDIF(${CMAKE_C_COMPILER} MATCHES "gcc")

MARK_AS_ADVANCED( SimTK_SDK )
MARK_AS_ADVANCED( SimTK_INSTALL_PREFIX )

SET(CMAKE_INSTALL_PREFIX ${SimTK_INSTALL_PREFIX} CACHE STRING "Install path prefix." FORCE)

INCLUDE_DIRECTORIES(${SimTK_SDK}/core/include)

## Watch out -- this next line will cause CMake to find installed versions of
## simmath libraries before the local one. Be sure to delete installed ones first.

LINK_DIRECTORIES(${SimTK_SDK}/core/lib${LIB64})


# libraries are installed from their subdirectories; headers here

# install headers
FILE(GLOB CORE_HEADERS     include/*.h                  */include/*.h)
FILE(GLOB TOP_HEADERS      include/simmath/*.h          */include/simmath/*.h)
FILE(GLOB INTERNAL_HEADERS include/simmath/internal/*.h */include/simmath/internal/*.h)

INSTALL( FILES ${CORE_HEADERS}     DESTINATION core/include     )
INSTALL( FILES ${TOP_HEADERS}      DESTINATION core/include/simmath)
INSTALL( FILES ${INTERNAL_HEADERS} DESTINATION core/include/simmath/internal )

