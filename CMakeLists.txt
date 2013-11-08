#---------------------------------------------------
# Simbody 
#
# This is the master CMake file that coordinates 
# the build of Simbody. There are four steps:
#    (1) Choose appropriate platform
#    (2) Build SimTKcommon
#    (3) Build SimTKmath
#    (4) Build SimTKsimbody
#
#----------------------------------------------------

cmake_minimum_required(VERSION 2.8)

if(COMMAND cmake_policy)
        cmake_policy(SET CMP0003 NEW)
        cmake_policy(SET CMP0005 NEW)
endif(COMMAND cmake_policy)

PROJECT(Simbody)

SET(SIMBODY_MAJOR_VERSION 3)
SET(SIMBODY_MINOR_VERSION 4)
SET(SIMBODY_PATCH_VERSION 0)

SET(SIMBODY_COPYRIGHT_YEARS "2005-13")

# underbar separated list of dotted authors, no spaces or commas
SET(SIMBODY_AUTHORS "Michael.Sherman_Peter.Eastman")

# Report the version number to the CMake UI. Don't include the 
# build version if it is zero.
SET(PATCH_VERSION_STRING)
IF(SIMBODY_PATCH_VERSION)
    SET(PATCH_VERSION_STRING ".${SIMBODY_PATCH_VERSION}")
ENDIF()

SET(SIMBODY_VERSION 
    "${SIMBODY_MAJOR_VERSION}.${SIMBODY_MINOR_VERSION}${PATCH_VERSION_STRING}" 
    CACHE STRING 
    "This is the version that will be built (can't be changed here)." 
    FORCE)

SET(SIMBODY_SONAME_VERSION
    "${SIMBODY_MAJOR_VERSION}.${SIMBODY_MINOR_VERSION}"
    CACHE STRING
    "Soname version to use when generating shared libs")

# This is the suffix if we're building and depending on versioned libraries.
SET(VN "_${SIMBODY_VERSION}")

SET(BUILD_BINARY_DIR ${CMAKE_BINARY_DIR}
    CACHE PATH 
    "The Simbody build (not the install) puts all the libraries and executables together here (with /Release, etc. appended on some platforms).")

# Make everything go in the same binary directory. (These are CMake-defined
# variables.)
SET(EXECUTABLE_OUTPUT_PATH ${BUILD_BINARY_DIR})
SET(LIBRARY_OUTPUT_PATH ${BUILD_BINARY_DIR})

# Static libraries, tests, and examples won't be built unless this 
# is set.
SET(BUILD_STATIC_LIBRARIES FALSE CACHE BOOL
"Build '_static' versions of libraries in addition to dynamic libraries?")

# Use this to generate a private set of libraries whose names
# won't conflict with installed versions.
SET(BUILD_USING_NAMESPACE "" CACHE STRING
	"All library names will be prefixed with 'xxx_' if this is set to xxx.")

SET(BUILD_UNVERSIONED_LIBRARIES TRUE CACHE BOOL
 "Build library names, and assume dependency names, with no version numbers?")

SET(BUILD_VERSIONED_LIBRARIES FALSE CACHE BOOL
 "Build library names, and assume dependency names, with version numbers?")

SET(NS)
IF(BUILD_USING_NAMESPACE)
    SET(NS "${BUILD_USING_NAMESPACE}_")
ENDIF(BUILD_USING_NAMESPACE)


#
# These are the names of all the libraries we may generate. These are
# target names so can be used to specify dependencies of one library
# on another. (In Debug mode the actual targets will have "_d" appended.)
#

SET(SimTKSIMBODY_LIBRARY_NAME ${NS}SimTKsimbody CACHE STRING
"Base name of the library being built; can't be changed here; see BUILD_USING_NAMESPACE variable."
FORCE)
SET(SimTKMATH_LIBRARY_NAME ${NS}SimTKmath CACHE STRING
"Base name of the library being built; can't be changed here; see BUILD_USING_NAMESPACE variable."
FORCE)
SET(SimTKCOMMON_LIBRARY_NAME ${NS}SimTKcommon CACHE STRING
"Base name of the library being built; can't be changed here; see BUILD_USING_NAMESPACE variable."
FORCE)


SET(SimTKCOMMON_SHARED_LIBRARY ${SimTKCOMMON_LIBRARY_NAME})
SET(SimTKCOMMON_STATIC_LIBRARY ${SimTKCOMMON_LIBRARY_NAME}_static)

SET(SimTKCOMMON_LIBRARY_NAME_VN ${NS}SimTKcommon${VN})
SET(SimTKCOMMON_SHARED_LIBRARY_VN ${SimTKCOMMON_LIBRARY_NAME_VN})
SET(SimTKCOMMON_STATIC_LIBRARY_VN ${SimTKCOMMON_LIBRARY_NAME_VN}_static)

SET(SimTKMATH_SHARED_LIBRARY ${SimTKMATH_LIBRARY_NAME})
SET(SimTKMATH_STATIC_LIBRARY ${SimTKMATH_LIBRARY_NAME}_static)

SET(SimTKMATH_LIBRARY_NAME_VN ${NS}SimTKmath${VN})
SET(SimTKMATH_SHARED_LIBRARY_VN ${SimTKMATH_LIBRARY_NAME_VN})
SET(SimTKMATH_STATIC_LIBRARY_VN ${SimTKMATH_LIBRARY_NAME_VN}_static)

SET(SimTKSIMBODY_SHARED_LIBRARY ${SimTKSIMBODY_LIBRARY_NAME})
SET(SimTKSIMBODY_STATIC_LIBRARY ${SimTKSIMBODY_LIBRARY_NAME}_static)

SET(SimTKSIMBODY_LIBRARY_NAME_VN ${NS}SimTKsimbody${VN})
SET(SimTKSIMBODY_SHARED_LIBRARY_VN ${SimTKSIMBODY_LIBRARY_NAME_VN})
SET(SimTKSIMBODY_STATIC_LIBRARY_VN ${SimTKSIMBODY_LIBRARY_NAME_VN}_static)


# Caution: this variable is automatically created by the CMake
# ENABLE_TESTING() command, but we'll take it over here for
# our own purposes too.
SET( BUILD_TESTING ON CACHE BOOL
	"Control building of Simbody test programs." )

SET( BUILD_EXAMPLES ON CACHE BOOL
	"Control building of Simbody example programs.")

# Set whether to build the Visualizer code.
SET(BUILD_VISUALIZER ON CACHE BOOL 
	"Control building of the visualizer component")

# Turning this off reduces the build time (and space) substantially,
# but you may miss the occasional odd bug. Also currently on Windows it
# is easier to debug the static tests than the DLL-liked ones.
SET( BUILD_TESTING_STATIC ON CACHE BOOL
    "If building static libraries, build static test and example programs too?" )

SET( BUILD_TESTING_SHARED ON CACHE BOOL
    "If building test or example programs, include dynamically-linked ones?" )

#
# Create a platform name useful for finding things in the Platform
# directory.
IF(WIN32)
    SET(PLATFORM_NAME Windows)
    SET(NATIVE_COPY_CMD copy)
ELSEIF(APPLE)
    SET(PLATFORM_NAME Mac)
    SET(NATIVE_COPY_CMD cp)
ELSE()
    SET(PLATFORM_NAME Linux)
    SET(NATIVE_COPY_CMD cp)
ENDIF()

# In addition to the platform name we need to know the Application Binary
# Interface (ABI) we're building for. Currently that is either x86, meaning
# 32 bit Intel instruction set, or x64 for 64 bit Intel instruction set.

IF(${CMAKE_SIZEOF_VOID_P} EQUAL 8)
    SET(PLATFORM_ABI x64)
ELSE()
    SET(PLATFORM_ABI x86)
ENDIF()

SET(BUILD_PLATFORM "${PLATFORM_NAME}:${PLATFORM_ABI}" CACHE STRING
    "This is the platform and ABI we're building for. Not changeable here; use a different CMake generator instead."
    FORCE)

# 
# Make sure "sandbox" input & output directories are set. During
# SimTK Core build, the build system will set them. Otherwise, we'll
# set them to sensible local values.
# If CMAKE_INSTALL_PREFIX is set then
# it is a sandbox installation area, otherwise we want to install
# in /usr/local/simbody or %ProgramFiles%\Simbody. 
#

# As of Ubuntu 12, 64 bit libraries are in lib, not lib64, so we are 
# not setting this on any platform now. If you are building on an
# older system you might have to set it manually.
SET(LIB64)

# For backward compatibility, use SimTK_INSTALL_PREFIX if
# CMAKE_INSTALL_PREFIX hasn't been set.
IF(SimTK_INSTALL_PREFIX AND CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX ${SimTK_INSTALL_PREFIX} CACHE PATH
      "Install path prefix." FORCE)
ENDIF()

# C compiler is gcc on Linux, gcc or cc on Mac where command line tools
# are installed from XCode. May be clang in Mac 10.8 or later.

IF(${CMAKE_C_COMPILER} MATCHES "cc" OR ${CMAKE_C_COMPILER} MATCHES "clang")
    IF(NOT CMAKE_INSTALL_PREFIX)
        SET(CMAKE_INSTALL_PREFIX "/usr/local/simbody" 
            CACHE PATH "Install directory")
    ENDIF(NOT CMAKE_INSTALL_PREFIX)

    # As of Ubuntu 12, 64 bit libraries are in lib, not lib64
    #IF(${CMAKE_SIZEOF_VOID_P} EQUAL 8)
    #    SET( LIB64 64 )
    #ENDIF()

ELSE() # Windows
    # On Win64, there are two environment variables, ProgramFiles(x86) for
    # 32 bit binaries, and ProgramW6432 for 64 bit binaries. The variable
    # ProgramFiles returns one or the other depending on whether the requesting
    # executable is 32 or 64 bits. That's CMake in this case, but we want
    # the right directory for the target. CMAKE_INSTALL_DIR after 2.8.3 gets
    # this right but still isn't what we want.
    IF( ${CMAKE_SIZEOF_VOID_P} EQUAL 8 )
	# 64 bit target on Win64
	set(PROGFILEDIR "$ENV{ProgramW6432}")
    ELSE() # Target is 32 bit
	set(PROGFILEDIR "$ENV{ProgramFiles(x86)}") # present if 64bit Windows
	if (NOT PROGFILEDIR)
	    set(PROGFILEDIR "$ENV{ProgramFiles}") # on 32bit Windows
	endif()
    ENDIF()

    IF(NOT CMAKE_INSTALL_PREFIX)
        SET(CMAKE_INSTALL_PREFIX "${PROGFILEDIR}/Simbody" 
            CACHE PATH "Install directory")
    ENDIF(NOT CMAKE_INSTALL_PREFIX)

ENDIF()

# On Mac, we'll use the defaults for the architecture and deployment target.
# This works OK on OSX 10.8. For older systems, or if you want to build for
# an older deployment target, or 32 bit or univeral binaries, you may have 
# to set CMAKE_OSX_ARCHITECTURES and/or CMAKE_OSX_DEPLOYMENT_TARGET. 
# Architectures are "i386", "x86_64", or "i386;x86_64" for universal.
# Deployment targets are "10.6" or whatever.
IF(APPLE)
    SET(LIB64 ) # never use lib64 on Mac
ENDIF(APPLE)

IF(UNIX AND NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release CACHE STRING "Debug or Release build" FORCE)
ENDIF (UNIX AND NOT CMAKE_BUILD_TYPE)

## Choose the maximum level of x86 instruction set that the compiler is 
## allowed to use. SSE2 is ubiquitous enough now that we don't mind
## abandoning machines that can't handle those instructions. SSE3 might
## also be reasonable by now (April 2009) so this default should be
## revisited soon. This can be set to a different value by the person
## running CMake.
SET(BUILD_INST_SET "sse2"  # use SSE2 instruction set by default
    CACHE STRING "CPU instruction level compiler is permitted to use.")
MARK_AS_ADVANCED( BUILD_INST_SET )

## When building in any of the Release modes, tell gcc to use full optimization and
## to generate SSE2 floating point instructions. Here we are specifying *all* of the
## Release flags, overriding CMake's defaults.
## Watch out for optimizer bugs in particular gcc versions!

IF(${CMAKE_C_COMPILER} MATCHES "cc" OR ${CMAKE_C_COMPILER} MATCHES "clang")
    string(TOLOWER ${BUILD_INST_SET} GCC_INST_SET)

    # Get the gcc version number in major.minor.build format
    execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion
                    OUTPUT_VARIABLE GCC_VERSION)

    # Unrolling fixed-count loops is a useful optimization for Simmatrix.
    SET(GCC_OPT_ENABLE "-funroll-loops")

    # If you know of optimization bugs that affect SimTK in particular
    # gcc versions, this is the place to turn off those optimizations.
    SET(GCC_OPT_DISABLE)
    # We know Gcc 4.4.3 on Ubuntu 10 is buggy and that Snow Leopard's
    # 4.2.1 is not. To be safe for now we'll assume anything over 4.3
    # should have these disabled.
    if (GCC_VERSION VERSION_GREATER 4.3 OR GCC_VERSION VERSION_EQUAL 4.3)
        SET(GCC_OPT_DISABLE "-fno-strict-aliasing")
        # Clang does not recognize these two flags.
        if (NOT ${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
            SET(GCC_OPT_DISABLE "-fno-tree-vrp -fno-guess-branch-probability")
        endif()
    endif()

    # C++
    SET(BUILD_CXX_FLAGS_DEBUG          "-g -m${GCC_INST_SET}" 
      CACHE STRING "g++ Debug build compile flags")
    SET(BUILD_CXX_FLAGS_RELEASE        
      "-DNDEBUG -O3 ${GCC_OPT_ENABLE} ${GCC_OPT_DISABLE} -m${GCC_INST_SET}" 
      CACHE STRING "g++ Release build compile flags")
    SET(BUILD_CXX_FLAGS_RELWITHDEBINFO 
      "-DNDEBUG -O3 -g ${GCC_OPT_ENABLE} ${GCC_OPT_DISABLE} -m${GCC_INST_SET}" 
      CACHE STRING "g++ RelWithDebInfo build compile flags")
    SET(BUILD_CXX_FLAGS_MINSIZEREL     "-DNDEBUG -Os -m${GCC_INST_SET}" 
      CACHE STRING "g++ MinSizeRel build compile flags")

    # C
    SET(BUILD_C_FLAGS_DEBUG            "-g -m${GCC_INST_SET}" 
      CACHE STRING "gcc Debug build compile flags")
    SET(BUILD_C_FLAGS_RELEASE          
      "-DNDEBUG -O3 ${GCC_OPT_ENABLE} ${GCC_OPT_DISABLE} -m${GCC_INST_SET}" 
      CACHE STRING "gcc Release build compile flags")
    SET(BUILD_C_FLAGS_RELWITHDEBINFO   
      "-DNDEBUG -O3 -g ${GCC_OPT_ENABLE} ${GCC_OPT_DISABLE} -m${GCC_INST_SET}" 
      CACHE STRING "gcc RelWithDebInfo build compile flags")
    SET(BUILD_C_FLAGS_MINSIZEREL       "-DNDEBUG -Os -m${GCC_INST_SET}" 
      CACHE STRING "gcc MinSizeRel build compile flags")

    # C++
    SET(CMAKE_CXX_FLAGS_DEBUG ${BUILD_CXX_FLAGS_DEBUG}
        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)
    SET(CMAKE_CXX_FLAGS_RELEASE ${BUILD_CXX_FLAGS_RELEASE}
        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)
    SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO ${BUILD_CXX_FLAGS_RELWITHDEBINFO}
        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)
    SET(CMAKE_CXX_FLAGS_MINSIZEREL ${BUILD_CXX_FLAGS_MINSIZEREL}
        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)

    # C
    SET(CMAKE_C_FLAGS_DEBUG ${BUILD_C_FLAGS_DEBUG}
        CACHE STRING "Can't change here -- see BUILD_C..." FORCE)
    SET(CMAKE_C_FLAGS_RELEASE ${BUILD_C_FLAGS_RELEASE}         
        CACHE STRING "Can't change here -- see BUILD_C..." FORCE)
    SET(CMAKE_C_FLAGS_RELWITHDEBINFO ${BUILD_C_FLAGS_RELWITHDEBINFO}
        CACHE STRING "Can't change here -- see BUILD_C..." FORCE)
    SET(CMAKE_C_FLAGS_MINSIZEREL ${BUILD_C_FLAGS_MINSIZEREL}
        CACHE STRING "Can't change here -- see BUILD_C..." FORCE)

ENDIF()

## When building in any of the Release modes, tell VC++ cl compiler to use 
## intrinsics (i.e. sqrt instruction rather than sqrt subroutine) with 
## flag /Oi.

## Careful -- "cl" would match "clang" also. Only check on Windows.
IF(WIN32 AND ${CMAKE_C_COMPILER} MATCHES "cl")
    STRING(TOUPPER ${BUILD_INST_SET} CL_INST_SET)
    set(BUILD_LIMIT_PARALLEL_COMPILES "" CACHE STRING
        "Set a maximum number of simultaneous compilations.")
    mark_as_advanced(BUILD_LIMIT_PARALLEL_COMPILES)
    set(mxcpu ${BUILD_LIMIT_PARALLEL_COMPILES}) # abbreviation

    ## C++
    SET(BUILD_CXX_FLAGS_DEBUG        
	"/D _DEBUG /MDd /Od /Ob0 /RTC1 /Zi /GS- /arch:${CL_INST_SET}" 
        CACHE STRING "VC++ Debug build compile flags")
    SET(BUILD_CXX_FLAGS_RELEASE        
	"/D NDEBUG /MD  /O2 /Ob2 /Oi /GS- /arch:${CL_INST_SET}" 
        CACHE STRING "VC++ Release build compile flags")
    SET(BUILD_CXX_FLAGS_RELWITHDEBINFO 
	"/D NDEBUG /MD  /O2 /Ob2 /Oi /Zi /GS- /arch:${CL_INST_SET}" 
        CACHE STRING "VC++ RelWithDebInfo build compile flags")
    SET(BUILD_CXX_FLAGS_MINSIZEREL 
	"/D NDEBUG /MD  /O1 /Ob1 /Oi /GS- /arch:${CL_INST_SET}" 
        CACHE STRING "VC++ MinSizeRel build compile flags")

    ## C
    SET(BUILD_C_FLAGS_DEBUG        
	"/D _DEBUG /MDd /Od /Ob0 /RTC1 /Zi /GS- /arch:${CL_INST_SET}" 
        CACHE STRING "VC Debug build compile flags")
    SET(BUILD_C_FLAGS_RELEASE        
	"/D NDEBUG /MD  /O2 /Ob2 /Oi /GS- /arch:${CL_INST_SET}" 
        CACHE STRING "VC Release build compile flags")
    SET(BUILD_C_FLAGS_RELWITHDEBINFO 
	"/D NDEBUG /MD  /O2 /Ob2 /Oi /Zi /GS- /arch:${CL_INST_SET}" 
        CACHE STRING "VC RelWithDebInfo build compile flags")
    SET(BUILD_C_FLAGS_MINSIZEREL 
	"/D NDEBUG /MD  /O1 /Ob1 /Oi /GS- /arch:${CL_INST_SET}" 
        CACHE STRING "VC MinSizeRel build compile flags")

    ## C++
    SET(CMAKE_CXX_FLAGS_DEBUG "/MP${mxcpu} ${BUILD_CXX_FLAGS_DEBUG}"
        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)
    SET(CMAKE_CXX_FLAGS_RELEASE "/MP${mxcpu} ${BUILD_CXX_FLAGS_RELEASE}"
        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)
    SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "/MP${mxcpu} ${BUILD_CXX_FLAGS_RELWITHDEBINFO}"
        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)
    SET(CMAKE_CXX_FLAGS_MINSIZEREL "/MP${mxcpu} ${BUILD_CXX_FLAGS_MINSIZEREL}"
        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)

    ## C
    SET(CMAKE_C_FLAGS_DEBUG "/MP${mxcpu} ${BUILD_C_FLAGS_DEBUG}"
        CACHE STRING "Can't change here -- see BUILD_C_..." FORCE)
    SET(CMAKE_C_FLAGS_RELEASE "/MP${mxcpu} ${BUILD_C_FLAGS_RELEASE}"
        CACHE STRING "Can't change here -- see BUILD_C_..." FORCE)
    SET(CMAKE_C_FLAGS_RELWITHDEBINFO "/MP${mxcpu} ${BUILD_C_FLAGS_RELWITHDEBINFO}"
        CACHE STRING "Can't change here -- see BUILD_C_..." FORCE)
    SET(CMAKE_C_FLAGS_MINSIZEREL "/MP${mxcpu} ${BUILD_C_FLAGS_MINSIZEREL}"
        CACHE STRING "Can't change here -- see BUILD_C_..." FORCE)

ENDIF()

MARK_AS_ADVANCED(
    BUILD_CXX_FLAGS_DEBUG
    BUILD_CXX_FLAGS_RELEASE
    BUILD_CXX_FLAGS_RELWITHDEBINFO
    BUILD_CXX_FLAGS_MINSIZEREL
    BUILD_C_FLAGS_DEBUG
    BUILD_C_FLAGS_RELEASE
    BUILD_C_FLAGS_RELWITHDEBINFO
    BUILD_C_FLAGS_MINSIZEREL
)

# Collect up information about the version of the simbody library we're building
# and make it available to the code so it can be built into the binaries.

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

# Ensure that debug libraries have "_d" appended to their names.
# CMake gets this right on Windows automatically with this definition.
# But on Unix or Cygwin we have to add the suffix manually
IF (${CMAKE_GENERATOR} MATCHES "Visual Studio")
    SET(CMAKE_DEBUG_POSTFIX "_d" CACHE INTERNAL "" FORCE)
ENDIF (${CMAKE_GENERATOR} MATCHES "Visual Studio")

# Determine which math libraries to use for this platform.
# Intel MKL: mkl_intel_c_dll;mkl_sequential_dll;mkl_core_dll
SET(BUILD_USING_OTHER_LAPACK "" CACHE STRING
 "If you have your own Lapack and Blas, put library basenames here, separated by semicolons. See LAPACK_BEING_USED to see what's actually being used.")

if(WIN32)
    SET(LAPACK_PLATFORM_DEFAULT liblapack;libblas)
else()
    SET(LAPACK_PLATFORM_DEFAULT lapack;blas)
endif()

SET(LAPACK_BEING_USED ${LAPACK_PLATFORM_DEFAULT} CACHE STRING
    "Basename of the actual Lapack library we're depending on; can't change here; see variable BUILD_USING_OTHER_LAPACK." FORCE)

IF(BUILD_USING_OTHER_LAPACK)
    SET(LAPACK_BEING_USED ${BUILD_USING_OTHER_LAPACK} CACHE STRING
"Basename of the actual Lapack library we're depending on; can't change here; see variable BUILD_USING_OTHER_LAPACK." FORCE)
ENDIF(BUILD_USING_OTHER_LAPACK)

IF(${CMAKE_C_COMPILER} MATCHES "cc" OR ${CMAKE_C_COMPILER} MATCHES "clang")
    IF(APPLE)
	    SET(REALTIME_LIB)
    ELSE()
	    SET(REALTIME_LIB rt)
    ENDIF()
    SET(MATH_LIBS_TO_USE    ${LAPACK_BEING_USED} pthread ${REALTIME_LIB} dl m)
    SET(MATH_LIBS_TO_USE_VN ${LAPACK_BEING_USED} pthread ${REALTIME_LIB} dl m)
ELSE()
    ## Assume Microsoft Visual Studio
    if (${PLATFORM_ABI} MATCHES "x64")
        SET(MATH_LIBS_TO_USE    ${LAPACK_BEING_USED} pthreadVC2_x64)
        SET(MATH_LIBS_TO_USE_VN ${LAPACK_BEING_USED} pthreadVC2_x64)
    else()
        SET(MATH_LIBS_TO_USE    ${LAPACK_BEING_USED} pthreadVC2)
        SET(MATH_LIBS_TO_USE_VN ${LAPACK_BEING_USED} pthreadVC2)
    endif()
ENDIF()

#
# Allow automated build and dashboard.
#
INCLUDE (Dart)
## When in Debug mode and running valgrind, some of the test
## cases take longer than the default 1500 seconds.
SET(DART_TESTING_TIMEOUT 7200)

IF (BUILD_TESTING)

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

    # Make a RUN_TESTS_PARALLEL target (thanks, Kevin!)
    # Specify number of cores to run for testing
    SET(TESTING_PROCESSOR_COUNT 4 CACHE STRING 
        "Number of CPUs to be used by the RUN_TESTS_PARALLEL target.")
    MARK_AS_ADVANCED(TESTING_PROCESSOR_COUNT)
    SET (cmd ${CMAKE_CTEST_COMMAND} -j${TESTING_PROCESSOR_COUNT})
    IF (MSVC)
        SET (cmd ${cmd} -C ${CMAKE_CFG_INTDIR})
    ELSE (MSVC)
        SET (cmd ${cmd} -C ${CMAKE_BUILD_TYPE})
    ENDIF (MSVC)
    ADD_CUSTOM_TARGET (RUN_TESTS_PARALLEL
        COMMAND ${cmd}
    )

ENDIF (BUILD_TESTING)

INCLUDE(ApiDoxygen.cmake)

# Each build should look in our local binary directory to find the
# earlier-built libraries that it depends on.
LINK_DIRECTORIES(${BUILD_BINARY_DIR})

# Each of these returns a list of API include directories for
# use by the later builds.
ADD_SUBDIRECTORY( Platform )
# PLATFORM_INCLUDE_DIRECTORIES now set
ADD_SUBDIRECTORY( SimTKcommon )
# SimTKCOMMON_INCLUDE_DIRECTORIES now set
ADD_SUBDIRECTORY( SimTKmath )
# SimTKMATH_INCLUDE_DIRECTORIES now set
ADD_SUBDIRECTORY( Simbody )
# SimTKSIMBODY_INCLUDE_DIRECTORIES now set (but not used)

FILE(GLOB TOPLEVEL_DOCS doc/*.pdf doc/*.txt)
INSTALL(FILES ${TOPLEVEL_DOCS} DESTINATION doc)

FILE(GLOB INSTALLED_CMAKE_STUFF cmake/*.cmake cmake/*.txt)
INSTALL(FILES ${INSTALLED_CMAKE_STUFF} DESTINATION share/cmake)
