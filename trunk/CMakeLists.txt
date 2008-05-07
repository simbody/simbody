#---------------------------------------------------
# SimTKcommon 
#
# Creates SimTK Core library, base name=SimTKcommon.
# Default libraries are shared & optimized. Variants
# are created for static (_static) and debug (_d).
#
# Windows:
#   SimTKcommon[_d].dll
#   SimTKcommon[_d].lib
#   SimTKcommon_static[_d].lib
# Unix:
#   libSimTKcommon[_d].so
#   libSimTKcommon_static[_d].a
#
# All libraries are installed in 
#   C:\Program Files\SimTK\core\lib[64]  (Windows)
#   /usr/local/SimTK/core/lib[64]        (UNIX)
#
# TODO: should we just install the default libraries
# at the top level, and install the static & debug
# libraries in .../core/lib/SimTKcommon/lib
# or           .../SimTKcommon/lib ?
#----------------------------------------------------

PROJECT (SimTKcommon)
SUBDIRS (staticTarget sharedTarget tests)

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
SET(SimTKCOMMON_SOURCE_SUBDIRS . Scalar SmallMatrix Mechanics BigMatrix Geometry Simulation Random Polynomial)


# 
# Make sure "sandbox" input & output directories are set. During
# SimTK Core build, the build system will set them. Otherwise, we'll
# set them to sensible local values.
# If SimTK_INSTALL_PREFIX is set then
# it is a sandbox installation area, otherwise we want to install
# in /usr/local/SimTK or c:\program files\SimTK. Similarly, SimTK_SDK
# tells us where to find other SimTK Core modules on which this one
# depends. If it is not set, we'll look in /usr/local/SimTK or
# c:\program files\SimTK *regardless* of SimTK_INSTALL_PREFIX.
#

IF(${CMAKE_C_COMPILER} MATCHES "gcc")

    IF(NOT SimTK_INSTALL_PREFIX)
        SET(SimTK_INSTALL_PREFIX "/usr/local/SimTK" CACHE PATH "Install directory")
    ENDIF(NOT SimTK_INSTALL_PREFIX)

    IF(NOT SimTK_SDK)
        SET(SimTK_SDK "/usr/local/SimTK" CACHE PATH "SimTK dependency path" )
    ENDIF(NOT SimTK_SDK)

ELSE(${CMAKE_C_COMPILER} MATCHES "gcc")

    IF(NOT SimTK_INSTALL_PREFIX)
        SET(SimTK_INSTALL_PREFIX "C:/Program Files/SimTK" CACHE PATH "Install directory")
    ENDIF(NOT SimTK_INSTALL_PREFIX)

    IF(NOT SimTK_SDK)
        SET(SimTK_SDK "C:/Program Files/SimTK" CACHE PATH "SimTK dependency path")
    ENDIF(NOT SimTK_SDK)

ENDIF(${CMAKE_C_COMPILER} MATCHES "gcc")


MARK_AS_ADVANCED( SimTK_SDK  )
MARK_AS_ADVANCED( SimTK_INSTALL_PREFIX  )

INCLUDE_DIRECTORIES(${SimTK_SDK}/core/include)

IF( ${CMAKE_SIZEOF_VOID_P} EQUAL 8 )
  SET( LIB64 64 )
ELSE( ${CMAKE_SIZEOF_VOID_P} EQUAL 8 )
  SET( LIB64  )
ENDIF( ${CMAKE_SIZEOF_VOID_P} EQUAL 8 )

# The build system will set ARCH64 for 64 bit builds, which require
# use of the lib64/ library directories rather than lib/.
# uncomment the folloing lines to support cross compiling 32/64 bit builds
#SET( ARCH64 OFF CACHE BOOL "ON for 64bit builds, OFF for 32bit builds")
#MARK_AS_ADVANCED( ARCH64  )
#IF (ARCH64)
#    SET(LIB64 64)
#ELSE (ARCH64)
#    SET(LIB64) # nothing
#ENDIF (ARCH64)



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




# Collect up information about the version of the SimTKcommon library we're building
# and make it available to the code so it can be built into the binaries.

SET(SimTKCOMMON_LIBRARY_NAME SimTKcommon)
SET(SimTKCOMMON_MAJOR_VERSION 0)
SET(SimTKCOMMON_MINOR_VERSION 6)
SET(SimTKCOMMON_BUILD_VERSION 3)

SET(SimTKCOMMON_COPYRIGHT_YEARS "2005-7")

# underbar separated list of dotted authors, no spaces or commas
SET(SimTKCOMMON_AUTHORS "Michael.Sherman")

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
                  OUTPUT_VARIABLE SimTKCOMMON_SVN_REVISION )
ELSE (SVNVERSION_EXE) 
    FIND_PROGRAM (SVNVERSION svnversion) 
    IF (SVNVERSION)
        EXEC_PROGRAM (${SVNVERSION}
                      ARGS \"${PROJECT_SOURCE_DIR}\"
                      OUTPUT_VARIABLE SimTKCOMMON_SVN_REVISION)
    ELSE (SVNVERSION)
        MESSAGE (STATUS 
          "Could not find 'svnversion' executable; 'about' will be wrong. (Cygwin provides one on Windows.)"
        )
        SET (SimTKCOMMON_SVN_REVISION unknown) 
    ENDIF (SVNVERSION)
ENDIF (SVNVERSION_EXE)
# Remove colon from build version, for easier placement in directory names
STRING(REPLACE ":" "_" SimTKCOMMON_SVN_REVISION ${SimTKCOMMON_SVN_REVISION})


ADD_DEFINITIONS(-DSimTK_SimTKCOMMON_LIBRARY_NAME=${SimTKCOMMON_LIBRARY_NAME}
                -DSimTK_SimTKCOMMON_MAJOR_VERSION=${SimTKCOMMON_MAJOR_VERSION}
                -DSimTK_SimTKCOMMON_MINOR_VERSION=${SimTKCOMMON_MINOR_VERSION}
		-DSimTK_SimTKCOMMON_BUILD_VERSION=${SimTKCOMMON_BUILD_VERSION})


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
   ADD_DEFINITIONS(-DSimTK_SimTKCOMMON_SVN_REVISION="${SimTKCOMMON_SVN_REVISION}"
                   -DSimTK_SimTKCOMMON_COPYRIGHT_YEARS="${SimTKCOMMON_COPYRIGHT_YEARS}"
                   -DSimTK_SimTKCOMMON_AUTHORS="${SimTKCOMMON_AUTHORS}")
ELSE(NEED_QUOTES)
   ADD_DEFINITIONS(-DSimTK_SimTKCOMMON_SVN_REVISION=${SimTKCOMMON_SVN_REVISION}
                   -DSimTK_SimTKCOMMON_COPYRIGHT_YEARS=${SimTKCOMMON_COPYRIGHT_YEARS}
                   -DSimTK_SimTKCOMMON_AUTHORS=${SimTKCOMMON_AUTHORS})
ENDIF(NEED_QUOTES)

# -DSimTK_SimTKCOMMON_LIBRARY_TYPE has to be defined in the target subdirectories.
# -DSimTKcommon_EXPORTS defined automatically when Windows DLL build is being done.



# Report the version number to the CMake UI
SET(SimTKCOMMON_VERSION
    "${SimTKCOMMON_MAJOR_VERSION}.${SimTKCOMMON_MINOR_VERSION}.${SimTKCOMMON_BUILD_VERSION}"
     CACHE STRING "This is the version which will be built." FORCE)

SET(SHARED_TARGET ${SimTKCOMMON_LIBRARY_NAME})
SET(STATIC_TARGET ${SimTKCOMMON_LIBRARY_NAME}_static)

## If no one says otherwise, change the executable path to drop into the same binary
## location as the DLLs so that the test cases will use the just-built DLLs.
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
IF (${CMAKE_GENERATOR} MATCHES "Visual Studio")
    SET(CMAKE_DEBUG_POSTFIX "_d" CACHE INTERNAL "" FORCE)
ENDIF (${CMAKE_GENERATOR} MATCHES "Visual Studio")

# But on Unix or Cygwin we have to add the suffix manually
IF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)
    SET(SHARED_TARGET ${SHARED_TARGET}_d)
    SET(STATIC_TARGET ${STATIC_TARGET}_d)
ENDIF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)

# Determine which math libraries to use for this platform.

IF(${CMAKE_C_COMPILER} MATCHES "gcc")
    SET(MATH_LIBS_TO_USE SimTKlapack pthread)
ELSE(${CMAKE_C_COMPILER} MATCHES "gcc")
    ## Assume Microsoft Visual Studio
    SET(MATH_LIBS_TO_USE SimTKlapack pthreadVC2)
ENDIF(${CMAKE_C_COMPILER} MATCHES "gcc")


# These are all the places to search for header files which are
# to be part of the API.
SET(API_INCLUDE_DIRS) # start empty
FOREACH(subdir ${SimTKCOMMON_SOURCE_SUBDIRS})
    # append
    SET(API_INCLUDE_DIRS ${API_INCLUDE_DIRS}
                         ${subdir}/include 
                         ${subdir}/include/SimTKcommon 
                         ${subdir}/include/SimTKcommon/internal)
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

FOREACH(subdir ${SimTKCOMMON_SOURCE_SUBDIRS})
    FILE(GLOB src_files  ${subdir}/src/*.cpp)
    FILE(GLOB incl_files ${subdir}/src/*.h)
    SET(SOURCE_FILES         ${SOURCE_FILES}         ${src_files})   #append
    SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})

    ## Make sure we find these locally before looking in SimTK/core/include if
    ## Simbody was previously installed there.
    INCLUDE_DIRECTORIES(BEFORE ${PROJECT_SOURCE_DIR}/${subdir}/include)
ENDFOREACH(subdir)



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

ADD_TEST(ListTest                 ${TEST_PATH}/ListTest)
ADD_TEST(ListTestStatic           ${TEST_PATH}/ListTestStatic)
ADD_TEST(BNTTest                  ${TEST_PATH}/BNTTest)
ADD_TEST(BNTTestStatic            ${TEST_PATH}/BNTTestStatic)
ADD_TEST(MatVecTest               ${TEST_PATH}/MatVecTest)
ADD_TEST(MatVecTestStatic         ${TEST_PATH}/MatVecTestStatic)
ADD_TEST(SpatialAlgebraTest       ${TEST_PATH}/SpatialAlgebraTest)
ADD_TEST(SpatialAlgebraTestStatic ${TEST_PATH}/SpatialAlgebraTestStatic)
ADD_TEST(OrientationTest          ${TEST_PATH}/OrientationTest)
ADD_TEST(OrientationTestStatic    ${TEST_PATH}/OrientationTestStatic)
#ADD_TEST(BigMatrixTest            ${TEST_PATH}/BigMatrixTest)
#ADD_TEST(BigMatrixTestStatic      ${TEST_PATH}/BigMatrixTestStatic)
ADD_TEST(StateTest                ${TEST_PATH}/StateTest)
ADD_TEST(StateTestStatic          ${TEST_PATH}/StateTestStatic)
ADD_TEST(RandomTest               ${TEST_PATH}/RandomTest)
ADD_TEST(RandomTestStatic         ${TEST_PATH}/RandomTestStatic)
ADD_TEST(SFMTTest                 ${TEST_PATH}/SFMTTest)
ADD_TEST(SFMTTestStatic           ${TEST_PATH}/SFMTTestStatic)
ADD_TEST(PolynomialTest           ${TEST_PATH}/PolynomialTest)
ADD_TEST(PolynomialTestStatic     ${TEST_PATH}/PolynomialTestStatic)
ADD_TEST(EnumerationTest          ${TEST_PATH}/EnumerationTest)
ADD_TEST(EnumerationTestStatic    ${TEST_PATH}/EnumerationTestStatic)
ADD_TEST(RotationTest             ${TEST_PATH}/RotationTest)
ADD_TEST(RotationTestStatic       ${TEST_PATH}/RotationTestStatic)

#
# Installation
#
# CMake installs go into CMAKE_INSTALL_PREFIX, which is always
# set to something incorrect by CMake. We'll use SimTK_INSTALL_PREFIX
# instead which is passed in by the build system or set above.

SET(CMAKE_INSTALL_PREFIX ${SimTK_INSTALL_PREFIX} CACHE STRING "Install path prefix." FORCE)


## Watch out -- this next line will cause CMake to find installed versions of
## SimTKcommon library before the local one. Be sure to delete installed one first.
LINK_DIRECTORIES(${SimTK_SDK}/core/lib${LIB64})

INCLUDE(ApiDoxygen.cmake)

# libraries are installed from their subdirectories; headers here

# install headers
FILE(GLOB CORE_HEADERS     include/*.h                      */include/*.h)
FILE(GLOB TOP_HEADERS      include/SimTKcommon/*.h          */include/SimTKcommon/*.h)
FILE(GLOB INTERNAL_HEADERS include/SimTKcommon/internal/*.h */include/SimTKcommon/internal/*.h)
INSTALL_FILES("/core/include"                      FILES ${CORE_HEADERS})
INSTALL_FILES("/core/include/SimTKcommon"          FILES ${TOP_HEADERS})
INSTALL_FILES("/core/include/SimTKcommon/internal" FILES ${INTERNAL_HEADERS})
INSTALL( FILES ${PROJECT_SOURCE_DIR}/doc/Simmatrix.pdf DESTINATION core/doc/simtkcommon )

