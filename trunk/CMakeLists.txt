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
#   %ProgramFiles%\SimTK\lib  (Windows)
#   /usr/local/SimTK/lib        (UNIX)
#
#----------------------------------------------------

cmake_minimum_required(VERSION 2.4.6)

if(COMMAND cmake_policy)
        cmake_policy(SET CMP0003 NEW)
        cmake_policy(SET CMP0005 NEW)
endif(COMMAND cmake_policy)

PROJECT ( SimTKsimmath )

## When in Debug mode and running valgrind, some of the test
## cases take longer than the default 1500 seconds.
SET(DART_TESTING_TIMEOUT 7200)

# Caution: this variable is automatically created by the CMake
# ENABLE_TESTING() command, but we'll take it over here for
# our own purposes too.
SET( BUILD_TESTING ON CACHE BOOL
	"Control building of Simmath test programs." )

SET( BUILD_EXAMPLES ON CACHE BOOL
	"Control building of Simmath example programs.")

# Turning this off reduces the build time (and space) substantially,
# but you may miss the occasional odd bug. Also currently on Windows it
# is easier to debug the static tests than the DLL-liked ones.
SET( BUILD_TESTING_STATIC ON CACHE BOOL
	"If building test or example programs, include statically-linked ones?" )

SET( BUILD_TESTING_SHARED ON CACHE BOOL
 	"If building test or example programs, include dynamically-linked ones?" )

# 
# Make sure "sandbox" input & output directories are set. During
# SimTK Core build, the build system will set them. Otherwise, we'll
# set them to sensible local values.
# If SimTK_INSTALL_PREFIX is set then
# it is a sandbox installation area, otherwise we want to install
# in /usr/local/SimTK or %ProgramFiles%\SimTK. Similarly, SimTK_SDK
# tells us where to find other SimTK Core modules on which this one
# depends. If it is not set, we'll look in /usr/local/SimTK or
# %ProgramFiles%\SimTK *regardless* of SimTK_INSTALL_PREFIX.
#

IF(${CMAKE_C_COMPILER} MATCHES "gcc")

    IF(NOT SimTK_INSTALL_PREFIX)
        SET(SimTK_INSTALL_PREFIX "/usr/local/SimTK" 
            CACHE PATH "Install directory")
    ENDIF(NOT SimTK_INSTALL_PREFIX)

    IF(NOT SimTK_SDK)
        SET(SimTK_SDK "/usr/local/SimTK" 
            CACHE PATH "SimTK dependency path" )
    ENDIF(NOT SimTK_SDK)

ELSE(${CMAKE_C_COMPILER} MATCHES "gcc")

    IF(NOT SimTK_INSTALL_PREFIX)
        SET(SimTK_INSTALL_PREFIX "$ENV{ProgramFiles}/SimTK" 
            CACHE PATH "Install directory")
    ENDIF(NOT SimTK_INSTALL_PREFIX)

    IF(NOT SimTK_SDK)
        SET(SimTK_SDK "$ENV{ProgramFiles}/SimTK" 
            CACHE PATH "SimTK dependency path")
    ENDIF(NOT SimTK_SDK)

ENDIF(${CMAKE_C_COMPILER} MATCHES "gcc")

MARK_AS_ADVANCED( SimTK_SDK  )
MARK_AS_ADVANCED( SimTK_INSTALL_PREFIX  )

INCLUDE_DIRECTORIES(${SimTK_SDK}/include)

IF( ${CMAKE_SIZEOF_VOID_P} EQUAL 8 )
  SET( LIB64 64 )
ELSE( ${CMAKE_SIZEOF_VOID_P} EQUAL 8 )
  SET( LIB64  )
ENDIF( ${CMAKE_SIZEOF_VOID_P} EQUAL 8 )

IF(${CMAKE_C_COMPILER} MATCHES "gcc")

   IF( WIN32 )
      SET( NO_CYGWIN -mno-cygwin )
   ENDIF( WIN32 )

   SET( CMAKE_CXX_FLAGS ${M32} )
   SET( CMAKE_SHARED_LINKER_FLAGS "${M32} ${NO_CYGWIN}" )
   SET( CMAKE_EXE_LINKER_FLAGS    "${M32} ${NO_CYGWIN}" )

ENDIF(${CMAKE_C_COMPILER} MATCHES "gcc")

# Allow build of binaries on Leopard that work on Tiger
IF( APPLE )
   SET( CMAKE_CXX_FLAGS "-mmacosx-version-min=10.4" )
   SET( CMAKE_C_FLAGS "-mmacosx-version-min=10.4" )
ENDIF( APPLE )

IF(UNIX AND NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Debug CACHE STRING "Debug or Release build" FORCE)
ENDIF (UNIX AND NOT CMAKE_BUILD_TYPE)

## Choose the maximum level of x86 instruction set that the compiler is 
## allowed to use. SSE2 is ubiquitous enough now that we don't mind
## abandoning machines that can't handle those instructions. SSE3 migh
## also be reasonable by now (April 2009) so this default should be
## revisited soon. This can be set to a different value by the person
## running CMake.
SET(BUILD_INST_SET "sse2"  # use SSE2 instruction set by default
    CACHE STRING "CPU instruction level compiler is permitted to use.")
MARK_AS_ADVANCED( BUILD_INST_SET )

## When building in any of the Release modes, tell gcc to use full optimization and
## to generate SSE2 floating point instructions. Here we are specifying *all* of the
## Release flags, overriding CMake's defaults.

IF(${CMAKE_C_COMPILER} MATCHES "gcc")
    STRING(TOLOWER ${BUILD_INST_SET} GCC_INST_SET)
    # C++
    SET(CMAKE_CXX_FLAGS_DEBUG          "-g -m${GCC_INST_SET}" 
        CACHE STRING "g++ Debug build compile flags" FORCE)
    SET(CMAKE_CXX_FLAGS_RELEASE        "-DNDEBUG -O3 -m${GCC_INST_SET}" 
        CACHE STRING "g++ Release build compile flags" FORCE)
    SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-DNDEBUG -O3 -g -m${GCC_INST_SET}" 
        CACHE STRING "g++ RelWithDebInfo build compile flags" FORCE)
    SET(CMAKE_CXX_FLAGS_MINSIZEREL     "-DNDEBUG -Os -m${GCC_INST_SET}" 
        CACHE STRING "g++ MinSizeRel build compile flags" FORCE)

    # C
    SET(CMAKE_C_FLAGS_DEBUG            "-g -m${GCC_INST_SET}" 
        CACHE STRING "gcc Debug build compile flags" FORCE)
    SET(CMAKE_C_FLAGS_RELEASE          "-DNDEBUG -O3 -m${GCC_INST_SET}" 
        CACHE STRING "gcc Release build compile flags" FORCE)
    SET(CMAKE_C_FLAGS_RELWITHDEBINFO   "-DNDEBUG -O3 -g -m${GCC_INST_SET}" 
        CACHE STRING "gcc RelWithDebInfo build compile flags" FORCE)
    SET(CMAKE_C_FLAGS_MINSIZEREL       "-DNDEBUG -Os -m${GCC_INST_SET}" 
        CACHE STRING "gcc MinSizeRel build compile flags" FORCE)

ENDIF(${CMAKE_C_COMPILER} MATCHES "gcc")

## When building in any of the Release modes, tell VC++ cl compiler to use intrinsics 
## (i.e. sqrt instruction rather than sqrt subroutine) with flag /Oi.

IF(${CMAKE_C_COMPILER} MATCHES "cl")
    STRING(TOUPPER ${BUILD_INST_SET} CL_INST_SET)

    ## C++
    SET(CMAKE_CXX_FLAGS_DEBUG        
	"/D _DEBUG /MDd /Od /Ob0 /RTC1 /Zi /arch:${CL_INST_SET}" 
        CACHE STRING "VC++ Debug build compile flags" FORCE)
    SET(CMAKE_CXX_FLAGS_RELEASE        
	"/D NDEBUG /MD  /O2 /Ob2 /Oi /arch:${CL_INST_SET}" 
        CACHE STRING "VC++ Release build compile flags" FORCE)
    SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO 
	"/D NDEBUG /MD  /O2 /Ob2 /Oi /Zi /arch:${CL_INST_SET}" 
        CACHE STRING "VC++ RelWithDebInfo build compile flags" FORCE)
    SET(CMAKE_CXX_FLAGS_MINSIZEREL 
	"/D NDEBUG /MD  /O1 /Ob1 /Oi /arch:${CL_INST_SET}" 
        CACHE STRING "VC++ MinSizeRel build compile flags" FORCE)

    ## C
    SET(CMAKE_C_FLAGS_DEBUG        
	"/D _DEBUG /MDd /Od /Ob0 /RTC1 /Zi /arch:${CL_INST_SET}" 
        CACHE STRING "VC++ Debug build compile flags" FORCE)
    SET(CMAKE_C_FLAGS_RELEASE        
	"/D NDEBUG /MD  /O2 /Ob2 /Oi /arch:${CL_INST_SET}" 
        CACHE STRING "VC++ Release build compile flags" FORCE)
    SET(CMAKE_C_FLAGS_RELWITHDEBINFO 
	"/D NDEBUG /MD  /O2 /Ob2 /Oi /Zi /arch:${CL_INST_SET}" 
        CACHE STRING "VC++ RelWithDebInfo build compile flags" FORCE)
    SET(CMAKE_C_FLAGS_MINSIZEREL 
	"/D NDEBUG /MD  /O1 /Ob1 /Oi /arch:${CL_INST_SET}" 
        CACHE STRING "VC++ MinSizeRel build compile flags" FORCE)

ENDIF(${CMAKE_C_COMPILER} MATCHES "cl")

SET(SIMMATH_SOURCE_SUBDIRS  . Ipopt LinearAlgebra Integrator)

# Collect up information about the version of the simmath library we're building
# and make it available to the code so it can be built into the binaries.

SET(SIMMATH_LIBRARY_NAME SimTKmath)
SET(SIMMATH_MAJOR_VERSION 1)
SET(SIMMATH_MINOR_VERSION 6)
SET(SIMMATH_BUILD_VERSION 0)

SET(SIMMATH_COPYRIGHT_YEARS "2005-9")

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

SET(SHARED_TARGET ${SIMMATH_LIBRARY_NAME})
SET(STATIC_TARGET ${SIMMATH_LIBRARY_NAME}_static)

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

SET(SimTKCOMMON_LIBRARY_NAME SimTKcommon)
SET(SimTKCOMMON_SHARED_LIBRARY ${SimTKCOMMON_LIBRARY_NAME})
SET(SimTKCOMMON_STATIC_LIBRARY ${SimTKCOMMON_LIBRARY_NAME}_static)

SET(SimTKCPODES_LIBRARY_NAME SimTKcpodes)
SET(SimTKCPODES_SHARED_LIBRARY ${SimTKCPODES_LIBRARY_NAME})
SET(SimTKCPODES_STATIC_LIBRARY ${SimTKCPODES_LIBRARY_NAME}_static)


#
# Find the right place to look for dependency resolution (which
# can be the default install directory or specified SimTK_SDK directory).
#
# Also, determine which math libraries to use for this platform.
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
FOREACH(subdir ${SIMMATH_SOURCE_SUBDIRS})
    # append
    SET(API_INCLUDE_DIRS ${API_INCLUDE_DIRS}
                         ${subdir}/include 
                         ${subdir}/include/simmath 
                         ${subdir}/include/simmath/internal)
ENDFOREACH(subdir)

# We'll need both *relative* path names, starting with
# their API_INCLUDE_DIRS, and absolute pathnames.
SET(API_REL_INCLUDE_FILES)   # start these out empty
SET(API_ABS_INCLUDE_FILES)

FOREACH(dir ${API_INCLUDE_DIRS})
    FILE(GLOB fullpaths ${dir}/*.h)	# returns full pathnames
    SET(API_ABS_INCLUDE_FILES ${API_ABS_INCLUDE_FILES} ${fullpaths})

    FOREACH(pathname ${fullpaths})
        GET_FILENAME_COMPONENT(filename ${pathname} NAME)
        SET(API_REL_INCLUDE_FILES ${API_REL_INCLUDE_FILES}
		${dir}/${filename})
    ENDFOREACH(pathname)
ENDFOREACH(dir)

# collect up source files
SET(SOURCE_FILES) # empty
SET(SOURCE_INCLUDE_FILES) # for dependency tracking only, I believe (sherm)

FOREACH(subdir ${SIMMATH_SOURCE_SUBDIRS})
    FILE(GLOB src_files  ${subdir}/src/*.cpp ${subdir}/src/*.c
		             ${subdir}/src/*/*.cpp ${subdir}/src/*/*.c)
    FILE(GLOB incl_files ${subdir}/src/*.h)
    SET(SOURCE_FILES         ${SOURCE_FILES}         ${src_files})   #append
    SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})

    ## Make sure we find these locally before looking in
    ## SimTK/include if Simmath was previously installed there.
    INCLUDE_DIRECTORIES(BEFORE ${PROJECT_SOURCE_DIR}/${subdir}/include)
ENDFOREACH(subdir)

## Special cases for Ipopt since it doesn't have a src subdirectory.
FILE(GLOB src_files  Ipopt/*.cpp)
FILE(GLOB incl_files Ipopt/*.h)
SET(SOURCE_FILES         ${SOURCE_FILES}         ${src_files})   #append
SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})

## Watch out -- this next line will cause CMake to find installed versions of
## simmath libraries before the local one. Be sure to delete installed ones first.
LINK_DIRECTORIES(${SimTK_SDK}/lib${LIB64})


IF( APPLE )
  INCLUDE_DIRECTORIES( "/usr/include/malloc")
ENDIF( APPLE )

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

# libraries and examples are installed from their subdirectories; headers here

# install headers
FILE(GLOB CORE_HEADERS     include/*.h                  */include/*.h)
FILE(GLOB TOP_HEADERS      include/simmath/*.h          */include/simmath/*.h)
FILE(GLOB INTERNAL_HEADERS include/simmath/internal/*.h */include/simmath/internal/*.h)

INSTALL_FILES(/include/                 FILES ${CORE_HEADERS})
INSTALL_FILES(/include/simmath/         FILES ${TOP_HEADERS})
INSTALL_FILES(/include/simmath/internal FILES ${INTERNAL_HEADERS})
INSTALL(FILES ${CMAKE_SOURCE_DIR}/doc/SimmathUserGuide.pdf DESTINATION doc/simmath)

# These are at the end because we want them processed after
# all the various variables have been set above.

ADD_SUBDIRECTORY( staticTarget )
ADD_SUBDIRECTORY( sharedTarget )
IF( BUILD_EXAMPLES )
  ADD_SUBDIRECTORY( examples )
ENDIF( BUILD_EXAMPLES )

IF( BUILD_TESTING )
  ADD_SUBDIRECTORY( tests )
ENDIF( BUILD_TESTING )
