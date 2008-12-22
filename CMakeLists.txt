#---------------------------------------------------
# Simbody 
#
# Creates SimTK Core library, base name=SimTKsimbody.
# Default libraries are shared & optimized. Variants
# are created for static (_static) and debug (_d).
#
# Windows:
#   SimTKsimbody[_d].dll
#   SimTKsimbody[_d].lib
#   SimTKsimbody_static[_d].lib
# Unix:
#   libSimTKsimbody[_d].so
#   libSimTKsimbody_static[_d].a
#
# All libraries are installed in 
#   C:\Program Files\SimTK\core\lib[64]  (Windows)
#   /usr/local/SimTK/core/lib[64]        (UNIX)
#
# TODO: should we just install the default libraries
# at the top level, and install the static & debug
# libraries in .../core/lib/simbody/lib
# or           .../simbody/lib ?
#----------------------------------------------------

cmake_minimum_required(VERSION 2.4.6)

if(COMMAND cmake_policy)
        cmake_policy(SET CMP0003 NEW)
        cmake_policy(SET CMP0005 NEW)
endif(COMMAND cmake_policy)

PROJECT (Simbody)

SET( BUILD_AUX OFF CACHE BOOL
	"Control building of Simbody aux library and tests" )

# Caution: this variable is automatically created by the CMake
# ENABLE_TESTING() command, but we'll take it over here for
# our own purposes too.
SET( BUILD_TESTING ON CACHE BOOL
	"Control building of Simbody test programs." )

# Turning this off reduces the build time (and space) substantially,
# but you may miss the occasional odd bug. Also currently on Windows it
# is easier to debug the static tests than the DLL-liked ones.
SET( BUILD_TESTING_STATIC ON CACHE BOOL
	"If building test programs, include statically-linked ones?" )

SET( BUILD_TESTING_SHARED ON CACHE BOOL
 	"If building test programs, include dynamically-linked ones?" )


# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
# (Only "./src" currently.)
SET(SIMBODY_SOURCE_SUBDIRS .)

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
        SET(SimTK_INSTALL_PREFIX "/usr/local/SimTK")
    ENDIF(NOT SimTK_INSTALL_PREFIX)

    IF(NOT SimTK_SDK)
        SET(SimTK_SDK "/usr/local/SimTK")
    ENDIF(NOT SimTK_SDK)

ELSE(${CMAKE_C_COMPILER} MATCHES "gcc")

    IF(NOT SimTK_INSTALL_PREFIX)
        FIND_PATH(SimTK_INSTALL_PREFIX
            NAME core
            PATHS
                "C:/Program Files/SimTK"
    )
    ENDIF(NOT SimTK_INSTALL_PREFIX)

    IF(NOT SimTK_SDK)
        FIND_PATH(SimTK_SDK
            NAME core
            PATHS
                "C:/Program Files/SimTK"
    )
    ENDIF(NOT SimTK_SDK)

ENDIF(${CMAKE_C_COMPILER} MATCHES "gcc")

MARK_AS_ADVANCED( SimTK_SDK  )
MARK_AS_ADVANCED( SimTK_INSTALL_PREFIX  )

INCLUDE_DIRECTORIES(${SimTK_SDK}/core/include)


# The build system will set ARCH64 for 64 bit builds, which require
# use of the lib64/ library directories rather than lib/.
#SET( ARCH64 OFF CACHE BOOL "ON for 64bit builds, OFF for 32bit builds")
#MARK_AS_ADVANCED( ARCH64  )
#IF (ARCH64)
#    SET(LIB64 64)
#ELSE (ARCH64)
#    SET(LIB64) # nothing
#ENDIF (ARCH64)

IF( ${CMAKE_SIZEOF_VOID_P} EQUAL 8 )
  SET( LIB64 64 )
ELSE( ${CMAKE_SIZEOF_VOID_P} EQUAL 8 )
  SET( LIB64  )
ENDIF( ${CMAKE_SIZEOF_VOID_P} EQUAL 8 )



IF( APPLE )
   SET( CMAKE_CXX_FLAGS "-mmacosx-version-min=10.4" )
ENDIF( APPLE )

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



# Collect up information about the version of the simbody library we're building
# and make it available to the code so it can be built into the binaries.

SET(SIMBODY_LIBRARY_NAME SimTKsimbody)
SET(SIMBODY_MAJOR_VERSION 1)
SET(SIMBODY_MINOR_VERSION 6)
SET(SIMBODY_BUILD_VERSION 0)

SET(SIMBODY_COPYRIGHT_YEARS "2005-9")

# underbar separated list of dotted authors, no spaces or commas
SET(SIMBODY_AUTHORS "Michael.Sherman_Peter.Eastman")

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
                  OUTPUT_VARIABLE SIMBODY_SVN_REVISION )
ELSE (SVNVERSION_EXE) 
    FIND_PROGRAM (SVNVERSION svnversion) 
    IF (SVNVERSION)
        EXEC_PROGRAM (${SVNVERSION}
                      ARGS "${PROJECT_SOURCE_DIR}"
                      OUTPUT_VARIABLE SIMBODY_SVN_REVISION)
    ELSE (SVNVERSION)
        MESSAGE (STATUS 
          "Could not find 'svnversion' executable; 'about' will be wrong. (Cygwin provides one on Windows.)"
        )
        SET (SIMBODY_SVN_REVISION unknown) 
    ENDIF (SVNVERSION)
ENDIF (SVNVERSION_EXE)
# Remove colon from build version, for easier placement in directory names
STRING(REPLACE ":" "_" SIMBODY_SVN_REVISION ${SIMBODY_SVN_REVISION})

ADD_DEFINITIONS(-DSimTK_SIMBODY_LIBRARY_NAME=${SIMBODY_LIBRARY_NAME}
                -DSimTK_SIMBODY_MAJOR_VERSION=${SIMBODY_MAJOR_VERSION}
                -DSimTK_SIMBODY_MINOR_VERSION=${SIMBODY_MINOR_VERSION}
		-DSimTK_SIMBODY_BUILD_VERSION=${SIMBODY_BUILD_VERSION})

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
   ADD_DEFINITIONS(-DSimTK_SIMBODY_SVN_REVISION="${SIMBODY_SVN_REVISION}"
                   -DSimTK_SIMBODY_COPYRIGHT_YEARS="${SIMBODY_COPYRIGHT_YEARS}"
                   -DSimTK_SIMBODY_AUTHORS="${SIMBODY_AUTHORS}")
ELSE(NEED_QUOTES)
   ADD_DEFINITIONS(-DSimTK_SIMBODY_SVN_REVISION=${SIMBODY_SVN_REVISION}
                   -DSimTK_SIMBODY_COPYRIGHT_YEARS=${SIMBODY_COPYRIGHT_YEARS}
                   -DSimTK_SIMBODY_AUTHORS=${SIMBODY_AUTHORS})
ENDIF(NEED_QUOTES)

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

SET(SimTKMATH_LIBRARY_NAME SimTKmath)
SET(SimTKMATH_SHARED_LIBRARY ${SimTKMATH_LIBRARY_NAME})
SET(SimTKMATH_STATIC_LIBRARY ${SimTKMATH_LIBRARY_NAME}_static)


# Determine which math libraries to use for this platform.
#

IF(${CMAKE_C_COMPILER} MATCHES "gcc")

    ## A gcc build
    SET(MATH_LIBS_TO_USE SimTKlapack pthread dl)

ELSE(${CMAKE_C_COMPILER} MATCHES "gcc")

    ## Assume Microsoft Visual Studio
    SET(MATH_LIBS_TO_USE SimTKlapack pthreadVC2)

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
    FILE(GLOB src_files  ${subdir}/src/*.cpp ${subdir}/src/*/*.cpp)
    FILE(GLOB incl_files ${subdir}/src/*.h)
    SET(SOURCE_FILES         ${SOURCE_FILES}         ${src_files})   #append
    SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})

    ## Make sure we find these locally before looking in SimTK/core/include if
    ## Simbody was previously installed there.
    INCLUDE_DIRECTORIES(BEFORE ${PROJECT_SOURCE_DIR}/${subdir}/include)
ENDFOREACH(subdir)

INCLUDE_DIRECTORIES(BEFORE ${PROJECT_SOURCE_DIR}/src)

IF (BUILD_TESTING)
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

ENDIF (BUILD_TESTING)

INCLUDE(ApiDoxygen.cmake)

#
# Installation
#
# CMake installs go into CMAKE_INSTALL_PREFIX, which is always
# set to something incorrect by CMake. We'll use SimTK_INSTALL_PREFIX
# instead which is passed in by the build system or set above.


SET(CMAKE_INSTALL_PREFIX ${SimTK_INSTALL_PREFIX} CACHE STRING "Install path prefix." FORCE)

## Watch out -- this next line will cause CMake to find installed versions of
## simmatrix libraries before the local one. Be sure to delete installed ones first.
LINK_DIRECTORIES(${SimTK_SDK}/core/lib${LIB64})

# libraries are installed from their subdirectories; headers here

# install headers
FILE(GLOB CORE_HEADERS     include/*.h                  */include/*.h)
FILE(GLOB TOP_HEADERS      include/simbody/*.h          */include/simbody/*.h)
FILE(GLOB INTERNAL_HEADERS include/simbody/internal/*.h */include/simbody/internal/*.h)
INSTALL_FILES(/core/include/                 FILES ${CORE_HEADERS})
INSTALL_FILES(/core/include/simbody/         FILES ${TOP_HEADERS})
INSTALL_FILES(/core/include/simbody/internal FILES ${INTERNAL_HEADERS})

# Notice that the (deprecated) INSTALL_FILES(dir ...) command works differently than the
# newer (recommended) INSTALL(dir FILES ...) command
# when the target directory starts with a slash
INSTALL( FILES ${PROJECT_SOURCE_DIR}/doc/SimbodyUserGuide.pdf DESTINATION core/doc/simbody )


ADD_SUBDIRECTORY( staticTarget )
ADD_SUBDIRECTORY( sharedTarget )
IF( BUILD_AUX )
  ADD_SUBDIRECTORY( auxiliary )
ENDIF( BUILD_AUX )

IF( BUILD_TESTING )
  ADD_SUBDIRECTORY( tests )
ENDIF( BUILD_TESTING )

