#-------------------------------------------------------------------------------
# Simbody
#
# This is the master CMake file that coordinates
# the build of Simbody. There are four steps:
#    (1) Get files needed for particular platform
#    (2) Build SimTKcommon library
#    (3) Build SimTKmath library
#    (4) Build SimTKsimbody library
#    (5) Build examples
#
#-------------------------------------------------------------------------------
if(MSVC)
    # To properly support Visual Studio 2015.
    set(SIMBODY_REQUIRED_CMAKE_VERSION 3.1.3)
else()
    set(SIMBODY_REQUIRED_CMAKE_VERSION 2.8.10)
endif()
cmake_minimum_required(VERSION ${SIMBODY_REQUIRED_CMAKE_VERSION})

if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
    cmake_policy(SET CMP0005 NEW)
endif()

project(Simbody)

# At this point CMake will have set CMAKE_INSTALL_PREFIX to /usr/local on Linux
# or appropriate ProgramFiles directory on Windows, for example
# C:/Program Files/Simbody, C:/Program Files (x86)/Simbody, or the local
# language equivalent.

if(WIN32)
    set(CMAKE_INSTALL_DOCDIR doc CACHE PATH "documentation root (DATAROOTDIR/doc/PROJECT_NAME)")
else()
    # Redefine DOCDIR to use the project name in lowercase to avoid
    # problems with some platforms: NTFS on Win, XFS or JFS variants
    set(CMAKE_INSTALL_DOCDIR share/doc/simbody CACHE PATH "documentation root (DATAROOTDIR/doc/PROJECT_NAME)")
endif()

# Include GNUInstallDirs to get canonical paths
include(GNUInstallDirs)

# Declare the option for code coverage
option(SIMBODY_COVERAGE
       "Adding ability to assess test coverage (requires gcc or clang)."
       OFF)

# Check compiler version
if(MSVC)
    if(MSVC_VERSION LESS 1800 OR MSVC_VERSION EQUAL 1800)
        message(FATAL_ERROR "\nMSVC does not support C++ 2011 features, for "
                            "example 'constexpr'. Update to at least MSVC 2015 "
                            "or use a MinGW version on Windows.\nLook at the "
                            "README.md for more information.\nIf you have the"
                            " 'Visual C++ Compiler Nov 2013 CTP (CTP_Nov2013)'"
                            " comment this test and configure normally.")
    endif()
    if(SIMBODY_COVERAGE)
        message(FATAL_ERROR "Code coverage is not possible with MSVC.")
    endif()
elseif(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
    set(SIMBODY_REQUIRED_GCC_VERSION 4.9.0)
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS ${SIMBODY_REQUIRED_GCC_VERSION})
        message(FATAL_ERROR "GNU GCC/G++ version is too old to compile Simbody.\n"
                            "Simbody requires at least version : "
                            "${SIMBODY_REQUIRED_GCC_VERSION}")
    endif()
elseif(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
    set(SIMBODY_REQUIRED_CLANG_VERSION 3.4)
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS ${SIMBODY_REQUIRED_CLANG_VERSION})
        message(FATAL_ERROR "Clang version is too old to compile Simbody.\n"
                            "Simbody requires at least version : "
                            "${SIMBODY_REQUIRED_CLANG_VERSION}")
    endif()
endif()

set(SIMBODY_MAJOR_VERSION 3)
set(SIMBODY_MINOR_VERSION 7)
set(SIMBODY_PATCH_VERSION 0)

set(SIMBODY_COPYRIGHT_YEARS "2005-18")

# underbar separated list of dotted authors, no spaces or commas
set(SIMBODY_AUTHORS "Michael.Sherman_Peter.Eastman")

# Report the version number to the CMake UI. Don't include the
# build version if it is zero.
set(PATCH_VERSION_STRING)
if(SIMBODY_PATCH_VERSION)
    set(PATCH_VERSION_STRING ".${SIMBODY_PATCH_VERSION}")
endif()

set(SIMBODY_VERSION
    "${SIMBODY_MAJOR_VERSION}.${SIMBODY_MINOR_VERSION}${PATCH_VERSION_STRING}"
    CACHE STRING
    "This is the version that will be built (can't be changed in GUI)."
    FORCE)

set(SIMBODY_SONAME_VERSION
    "${SIMBODY_MAJOR_VERSION}.${SIMBODY_MINOR_VERSION}"
    CACHE STRING
    "Soname version; appended to names of shared libs
    (can't be changed in GUI)."
    FORCE)

# This is the suffix if we're building and depending on versioned libraries.
set(VN "_${SIMBODY_VERSION}")

set(BUILD_BINARY_DIR ${CMAKE_BINARY_DIR}
    CACHE PATH
    "The Simbody build (not the install) puts all the libraries and executables together here (with /Release, etc. appended on some platforms).")

# Make everything go in the same binary directory. (These are CMake-defined
# variables.)
set(EXECUTABLE_OUTPUT_PATH ${BUILD_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${BUILD_BINARY_DIR})

# Dynamic libraries, tests, and examples
set(BUILD_DYNAMIC_LIBRARIES TRUE CACHE BOOL
"Build dynamic versions of libraries")

# Static libraries, tests, and examples won't be built unless this
# is set.
set(BUILD_STATIC_LIBRARIES FALSE CACHE BOOL
"Build '_static' versions of libraries")

if(NOT BUILD_DYNAMIC_LIBRARIES AND NOT BUILD_STATIC_LIBRARIES)
    message(SEND_ERROR "BUILD_DYNAMIC_LIBRARIES and BUILD_STATIC_LIBRARIES "
        "are both set to FALSE. No compilation is possible!")
endif()

# Use this to generate a private set of libraries whose names
# won't conflict with installed versions.
set(BUILD_USING_NAMESPACE "" CACHE STRING
    "All library names will be prefixed with 'xxx_' if this is set to xxx.")

set(BUILD_UNVERSIONED_LIBRARIES TRUE CACHE BOOL
 "Build library names, and assume dependency names, with no version numbers?")

set(BUILD_VERSIONED_LIBRARIES FALSE CACHE BOOL
 "Build library names, and assume dependency names, with version numbers?")

set(NS)
if(BUILD_USING_NAMESPACE)
    set(NS "${BUILD_USING_NAMESPACE}_")
endif()

# Ensure that debug libraries have "_d" appended to their names.
set(CMAKE_DEBUG_POSTFIX "_d")

#
# These are the names of all the libraries we may generate. These are
# target names so can be used to specify dependencies of one library
# on another. (In Debug mode the actual targets will have "_d" appended.)
#

set(SimTKSIMBODY_LIBRARY_NAME ${NS}SimTKsimbody CACHE STRING
"Base name of the library being built; can't be changed here; see BUILD_USING_NAMESPACE variable."
FORCE)
set(SimTKMATH_LIBRARY_NAME ${NS}SimTKmath CACHE STRING
"Base name of the library being built; can't be changed here; see BUILD_USING_NAMESPACE variable."
FORCE)
set(SimTKCOMMON_LIBRARY_NAME ${NS}SimTKcommon CACHE STRING
"Base name of the library being built; can't be changed here; see BUILD_USING_NAMESPACE variable."
FORCE)


set(SimTKCOMMON_SHARED_LIBRARY ${SimTKCOMMON_LIBRARY_NAME})
set(SimTKCOMMON_STATIC_LIBRARY ${SimTKCOMMON_LIBRARY_NAME}_static)

set(SimTKCOMMON_LIBRARY_NAME_VN ${NS}SimTKcommon${VN})
set(SimTKCOMMON_SHARED_LIBRARY_VN ${SimTKCOMMON_LIBRARY_NAME_VN})
set(SimTKCOMMON_STATIC_LIBRARY_VN ${SimTKCOMMON_LIBRARY_NAME_VN}_static)

set(SimTKMATH_SHARED_LIBRARY ${SimTKMATH_LIBRARY_NAME})
set(SimTKMATH_STATIC_LIBRARY ${SimTKMATH_LIBRARY_NAME}_static)

set(SimTKMATH_LIBRARY_NAME_VN ${NS}SimTKmath${VN})
set(SimTKMATH_SHARED_LIBRARY_VN ${SimTKMATH_LIBRARY_NAME_VN})
set(SimTKMATH_STATIC_LIBRARY_VN ${SimTKMATH_LIBRARY_NAME_VN}_static)

set(SimTKSIMBODY_SHARED_LIBRARY ${SimTKSIMBODY_LIBRARY_NAME})
set(SimTKSIMBODY_STATIC_LIBRARY ${SimTKSIMBODY_LIBRARY_NAME}_static)

set(SimTKSIMBODY_LIBRARY_NAME_VN ${NS}SimTKsimbody${VN})
set(SimTKSIMBODY_SHARED_LIBRARY_VN ${SimTKSIMBODY_LIBRARY_NAME_VN})
set(SimTKSIMBODY_STATIC_LIBRARY_VN ${SimTKSIMBODY_LIBRARY_NAME_VN}_static)


# Caution: this variable is automatically created by the CMake
# enable_testing() command, but we'll take it over here for
# our own purposes too.
set(BUILD_TESTING ON CACHE BOOL
    "Control building of Simbody test programs.
    To actually build tests, one
    or both of BUILD_TESTS_AND_EXAMPLES_STATIC and
    BUILD_TESTS_AND_EXAMPLES_SHARED must be ON.")

set(BUILD_EXAMPLES ON CACHE BOOL
    "Control building of Simbody example programs.
    To actually build examples, one
    or both of BUILD_TESTS_AND_EXAMPLES_STATIC and
    BUILD_TESTS_AND_EXAMPLES_SHARED must be ON.")

# Set whether to build the Visualizer code.
set(BUILD_VISUALIZER ON CACHE BOOL
    "Control building of the visualizer component.")

# Turning this off reduces the build time (and space) substantially,
# but you may miss the occasional odd bug. Also currently on Windows it
# is easier to debug the static tests than the DLL-linked ones.
set(BUILD_TESTS_AND_EXAMPLES_STATIC ON CACHE BOOL
    "If BUILDING_STATIC_LIBRARIES and BUILD_TESTING or BUILD_EXAMPLES, build
    statically-linked test and example programs too? On Windows,
    statically-linked tests may be easier to debug than DLL-linked tests.
    Statically-linked examples are never installed.")
mark_as_advanced(BUILD_TESTS_AND_EXAMPLES_STATIC)

if(BUILD_DYNAMIC_LIBRARIES)
    set(BUILD_TESTS_AND_EXAMPLES_SHARED ON CACHE BOOL
        "If BUILD_TESTING or BUILD_EXAMPLES, build dynamically-linked ones?")
else()
    set(BUILD_TESTS_AND_EXAMPLES_SHARED OFF)
endif()
mark_as_advanced(BUILD_TESTS_AND_EXAMPLES_SHARED)

if(BUILD_TESTING AND NOT (BUILD_TESTS_AND_EXAMPLES_STATIC OR
        BUILD_TESTS_AND_EXAMPLES_SHARED))
    message(SEND_ERROR "No tests would be built, despite BUILD_EXAMPLES"
        "being on, because BUILD_TESTS_AND_EXAMPLES_STATIC and "
        "BUILD_TESTS_AND_EXAMPLES_SHARED are both off.")
endif()

if(BUILD_EXAMPLES AND NOT (BUILD_TESTS_AND_EXAMPLES_STATIC OR
        BUILD_TESTS_AND_EXAMPLES_SHARED))
    message(SEND_ERROR
        "No examples would be built, despite BUILD_EXAMPLES being on, "
        "because BUILD_TESTS_AND_EXAMPLES_STATIC and "
        "BUILD_TESTS_AND_EXAMPLES_SHARED are both off.")
endif()

# In addition to the platform name we need to know the Application Binary
# Interface (ABI) we're building for. Currently that is either x86, meaning
# 32 bit Intel instruction set, or x64 for 64 bit Intel instruction set.

if(${CMAKE_SIZEOF_VOID_P} EQUAL 8)
    set(PLATFORM_ABI x64)
else()
    set(PLATFORM_ABI x86)
endif()

set(BUILD_PLATFORM "${CMAKE_HOST_SYSTEM_NAME}:${PLATFORM_ABI}" CACHE STRING
    "This is the platform and ABI we're building for. Not changeable here; use a different CMake generator instead."
    FORCE)

# If CMAKE_INSTALL_PREFIX is /usr/local, then the LIBDIR should necessarily be
# lib/. Sometimes (on Linux), LIBDIR is something like x86_64-linux-gnu. The
# linker will search /usr/lib/x86_64-linux-gnu (this path is in
# /etc/ld.so.conf.d), but it will NOT search /usr/local/lib/x86-64-linux-gnu.
# HOWEVER, it WILL search /usr/local/lib. So that Linux users needn't modify
# their LD_LIBRARY_PATH if installing to /usr/local, we force the LIBDIR to be
# lib/.
# Note: CMake 3.0 fixes this issue. When we move to CMake 3.0, we can
# remove this if-statement. See issue #151.
if("${CMAKE_INSTALL_PREFIX}" STREQUAL "/usr/local" OR
        "${CMAKE_INSTALL_PREFIX}" STREQUAL "/usr/local/")
    # Overwrite both of these variables; we use both of them.
    set(CMAKE_INSTALL_LIBDIR "lib")
    set(CMAKE_INSTALL_FULL_LIBDIR
        "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
endif()


if(NOT MSVC AND NOT XCODE AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
        "Debug, RelWithDebInfo (recommended), or Release build"
        FORCE)
endif()


## Choose the maximum level of x86 instruction set that the compiler is
## allowed to use.
## Was using sse2 but changed to let the compilers choose. Most will
## probably use sse2 or later by default.
## On 64 bit MSVC 2013, the default is sse2 and the argument
## isn't recognized so don't specify it.
if(CMAKE_CL_64)
    set(default_build_inst_set)
else()
    set(default_build_inst_set)
endif()

## This can be set to a different value by the person running CMake.
set(BUILD_INST_SET ""
    CACHE STRING
    "CPU instruction level compiler is permitted to use (default: let compiler decide).")
mark_as_advanced( BUILD_INST_SET )

if(BUILD_INST_SET)
    set(inst_set_to_use ${BUILD_INST_SET})
else()
    set(inst_set_to_use ${default_build_inst_set})
endif()


# RPATH
# -----
set(SIMBODY_USE_INSTALL_RPATH FALSE)
if(APPLE AND NOT (${CMAKE_VERSION} VERSION_LESS 2.8.12))
    # CMake 2.8.12 introduced the ability to set RPATH for shared libraries on
    # OSX. This helps executables find the libraries they depend on without
    # having to set the DYLD_LIBRARY_PATH environment variable.

    # Consider a library libfoo.dyld and an executable bar.
    # On OSX, libraries have an "install name" that, when linking, is copied
    # into the target (e.g., bar). The "install name" can be the full path to
    # libfoo.dylib, in which case bar will have no trouble finding libfoo.dylib
    # at runtime (since it has the full path to it). This doesn't work if you
    # want to be able to relocate your project. Therefore, it's possible to use
    # special tokens in the "install name" that are evaluated separately. The
    # token "@exectuable_path" is evaluated (at run time) to the full path of
    # the exectuable (e.g., bar) that is trying to load libfoo.dylib. An even
    # more flexible token is "@rpath", which is evaluated to a path (called
    # RPATH) that can be baked into the executable just after compiling or any
    # time before running the executable (using the executable
    # "install_name_tool"). The RPATH stored in executables can also contain
    # "@executable_path", etc.

    # Set the install name of libraries to contain "@rpath".
    # This allows clients of our libraries to point to them however they wish.
    set(CMAKE_MACOSX_RPATH ON)

    # We only need to set RPATH in executables if the libraries are installed
    # into directories that are not already searched by the linker.
    list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
        "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" isSystemDir)
    # CMake sets isSystemDir to -1 if the lib dir is NOT a system lib dir.
    if("${isSystemDir}" STREQUAL "-1")
        # This variable is used later on to toggle if RPATH should be set for
        # specific targets.
        set(SIMBODY_USE_INSTALL_RPATH TRUE)
    endif()
endif()


## When building in any of the Release modes, tell gcc/clang to use
## not-quite most agressive optimization.  Here we
## are specifying *all* of the Release flags, overriding CMake's defaults.
## Watch out for optimizer bugs in particular gcc versions!

if(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU" OR
   ${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")

    # C++11
    # -----
    # Using C++11 on OSX requires using libc++ instead of libstd++.
    # libc++ is an implementation of the C++ standard library for OSX.
    if(APPLE)
        if(XCODE)
            set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD "c++11")
            set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
        else()
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
            if(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
                set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
            endif()
        endif()
    else() # not APPLE
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    endif()

    if(inst_set_to_use)
        string(TOLOWER ${inst_set_to_use} GCC_INST_SET)
        set(GCC_INST_SET "-m${GCC_INST_SET}")
    else()
        set(GCC_INST_SET)
    endif()

    # Get the gcc or clang version number in major.minor.build format
    execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion
                    OUTPUT_VARIABLE GCC_VERSION)

    # Unrolling fixed-count loops was a useful optimization for Simmatrix
    # in earlier gcc versions.
    # Doesn't have a big effect for current compiler crop and may be
    # pushing our luck with optimizer bugs. So let the compilers decide
    # how to handle loops instead.
    ##set(GCC_OPT_ENABLE "-funroll-loops")

    # If you know of optimization bugs that affect Simbody in particular
    # gcc versions, this is the place to turn off those optimizations.
    set(GCC_OPT_DISABLE)

    # C++
    set(BUILD_CXX_FLAGS_DEBUG          "-g ${GCC_INST_SET}")
    set(BUILD_CXX_FLAGS_RELEASE
      "-DNDEBUG -O2 ${GCC_OPT_ENABLE} ${GCC_OPT_DISABLE} ${GCC_INST_SET}")
    set(BUILD_CXX_FLAGS_RELWITHDEBINFO
    "-DNDEBUG -O2 -g ${GCC_OPT_ENABLE} ${GCC_OPT_DISABLE} ${GCC_INST_SET}")
    set(BUILD_CXX_FLAGS_MINSIZEREL     "-DNDEBUG -Os ${GCC_INST_SET}")

    # C
    set(BUILD_C_FLAGS_DEBUG            "-g ${GCC_INST_SET}")
    set(BUILD_C_FLAGS_RELEASE
      "-DNDEBUG -O2 ${GCC_OPT_ENABLE} ${GCC_OPT_DISABLE} ${GCC_INST_SET}")
    set(BUILD_C_FLAGS_RELWITHDEBINFO
    "-DNDEBUG -O2 -g ${GCC_OPT_ENABLE} ${GCC_OPT_DISABLE} ${GCC_INST_SET}")
    set(BUILD_C_FLAGS_MINSIZEREL       "-DNDEBUG -Os ${GCC_INST_SET}")

    # C++
    set(CMAKE_CXX_FLAGS_DEBUG ${BUILD_CXX_FLAGS_DEBUG}
        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)
    set(CMAKE_CXX_FLAGS_RELEASE ${BUILD_CXX_FLAGS_RELEASE}
        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO ${BUILD_CXX_FLAGS_RELWITHDEBINFO}
        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)
    set(CMAKE_CXX_FLAGS_MINSIZEREL ${BUILD_CXX_FLAGS_MINSIZEREL}
        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)

    # C
    set(CMAKE_C_FLAGS_DEBUG ${BUILD_C_FLAGS_DEBUG}
        CACHE STRING "Can't change here -- see BUILD_C..." FORCE)
    set(CMAKE_C_FLAGS_RELEASE ${BUILD_C_FLAGS_RELEASE}
        CACHE STRING "Can't change here -- see BUILD_C..." FORCE)
    set(CMAKE_C_FLAGS_RELWITHDEBINFO ${BUILD_C_FLAGS_RELWITHDEBINFO}
        CACHE STRING "Can't change here -- see BUILD_C..." FORCE)
    set(CMAKE_C_FLAGS_MINSIZEREL ${BUILD_C_FLAGS_MINSIZEREL}
        CACHE STRING "Can't change here -- see BUILD_C..." FORCE)

    set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
    if(SIMBODY_COVERAGE)
        find_package(codecov)
    endif()

endif()

## When building in any of the Release modes, tell VC++ cl compiler to use
## intrinsics (i.e. sqrt instruction rather than sqrt subroutine) with
## flag /Oi.
## Caution: can't use CMAKE_CXX_COMPILER_ID MATCHES MSVC here because
## "MSVC" is a predefined CMAKE variable and will get expanded to 1 or 0
if(MSVC)
    if(inst_set_to_use)
        string(TOUPPER ${inst_set_to_use} CL_INST_SET)
        set(CL_INST_SET "/arch:${CL_INST_SET}")
    else()
        set(CL_INST_SET)
    endif()

    set(BUILD_LIMIT_PARALLEL_COMPILES "" CACHE STRING
        "Set a maximum number of simultaneous compilations.")
    mark_as_advanced(BUILD_LIMIT_PARALLEL_COMPILES)
    set(mxcpu ${BUILD_LIMIT_PARALLEL_COMPILES}) # abbreviation

    ## C++
    set(BUILD_CXX_FLAGS_DEBUG
    "/D _DEBUG /MDd /Od /Ob0 /RTC1 /Zi /GS- ${CL_INST_SET}") 
    set(BUILD_CXX_FLAGS_RELEASE
    "/D NDEBUG /MD  /O2 /Ob2 /Oi /GS- ${CL_INST_SET}")
    set(BUILD_CXX_FLAGS_RELWITHDEBINFO
    "/D NDEBUG /MD  /O2 /Ob2 /Oi /Zi /GS- ${CL_INST_SET}")
    set(BUILD_CXX_FLAGS_MINSIZEREL
    "/D NDEBUG /MD  /O1 /Ob1 /Oi /GS- ${CL_INST_SET}")

    ## C
    set(BUILD_C_FLAGS_DEBUG
    "/D _DEBUG /MDd /Od /Ob0 /RTC1 /Zi /GS- ${CL_INST_SET}")
    set(BUILD_C_FLAGS_RELEASE
    "/D NDEBUG /MD  /O2 /Ob2 /Oi /GS- ${CL_INST_SET}")
    set(BUILD_C_FLAGS_RELWITHDEBINFO
    "/D NDEBUG /MD  /O2 /Ob2 /Oi /Zi /GS- ${CL_INST_SET}")
    set(BUILD_C_FLAGS_MINSIZEREL
    "/D NDEBUG /MD  /O1 /Ob1 /Oi /GS- ${CL_INST_SET}")

    ## C++
    set(CMAKE_CXX_FLAGS_DEBUG "/MP${mxcpu} ${BUILD_CXX_FLAGS_DEBUG}"
        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)
    set(CMAKE_CXX_FLAGS_RELEASE "/MP${mxcpu} ${BUILD_CXX_FLAGS_RELEASE}"
        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "/MP${mxcpu} ${BUILD_CXX_FLAGS_RELWITHDEBINFO}"
        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)
    set(CMAKE_CXX_FLAGS_MINSIZEREL "/MP${mxcpu} ${BUILD_CXX_FLAGS_MINSIZEREL}"
        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)

    ## C
    set(CMAKE_C_FLAGS_DEBUG "/MP${mxcpu} ${BUILD_C_FLAGS_DEBUG}"
        CACHE STRING "Can't change here -- see BUILD_C_..." FORCE)
    set(CMAKE_C_FLAGS_RELEASE "/MP${mxcpu} ${BUILD_C_FLAGS_RELEASE}"
        CACHE STRING "Can't change here -- see BUILD_C_..." FORCE)
    set(CMAKE_C_FLAGS_RELWITHDEBINFO "/MP${mxcpu} ${BUILD_C_FLAGS_RELWITHDEBINFO}"
        CACHE STRING "Can't change here -- see BUILD_C_..." FORCE)
    set(CMAKE_C_FLAGS_MINSIZEREL "/MP${mxcpu} ${BUILD_C_FLAGS_MINSIZEREL}"
        CACHE STRING "Can't change here -- see BUILD_C_..." FORCE)

endif()

# Collect up information about the version of the simbody library we're building
# and make it available to the code so it can be built into the binaries.
# TODO removed SVN_REVSION; replace with GIT_SHA1
# http://stackoverflow.com/questions/1435953/how-can-i-pass-git-sha1-to-compiler-as-definition-using-cmake

# CMake quotes automatically when building Visual Studio projects but we need
# to add them ourselves for Linux or Cygwin. Two cases to avoid duplicate quotes
# in Visual Studio which end up in the binary.

if(${CMAKE_GENERATOR} MATCHES "Visual Studio")
    set(NEED_QUOTES FALSE)
else()
    set(NEED_QUOTES TRUE)
endif()

##TODO: doesn't work without quotes in nightly build
set(NEED_QUOTES TRUE)

if(NEED_QUOTES)
    add_definitions(-DSimTK_SIMBODY_COPYRIGHT_YEARS="${SIMBODY_COPYRIGHT_YEARS}"
                    -DSimTK_SIMBODY_AUTHORS="${SIMBODY_AUTHORS}")
else()
    add_definitions(-DSimTK_SIMBODY_COPYRIGHT_YEARS=${SIMBODY_COPYRIGHT_YEARS}
                    -DSimTK_SIMBODY_AUTHORS=${SIMBODY_AUTHORS})
endif()

# Determine which math libraries to use for this platform.
# Intel MKL: mkl_intel_c_dll;mkl_sequential_dll;mkl_core_dll
set(BUILD_USING_OTHER_LAPACK "" CACHE STRING
 "If you have your own Lapack and Blas, put libraries here, separated by semicolons (full paths or paths that are on the (DY)LD_LIBRARY_PATH (UNIX) or PATH (Windows)). See LAPACK_BEING_USED to see what's actually being used.")
if(WIN32)
    set(LAPACK_PLATFORM_DEFAULT lapack;blas)
else()
    find_package(BLAS)
    find_package(LAPACK)
    if(BLAS_FOUND AND LAPACK_FOUND)
        set(LAPACK_PLATFORM_DEFAULT ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    else()
        message(WARNING "Could not find blas/lapack")
    endif()
endif()
set(LAPACK_BEING_USED ${LAPACK_PLATFORM_DEFAULT} CACHE STRING
    "Basename of the actual Lapack library we're depending on; can't change here; see variable BUILD_USING_OTHER_LAPACK." FORCE)

if(BUILD_USING_OTHER_LAPACK)
    set(LAPACK_BEING_USED ${BUILD_USING_OTHER_LAPACK} CACHE STRING
"Basename of the actual Lapack library we're depending on; can't change here; see variable BUILD_USING_OTHER_LAPACK." FORCE)
endif()

if(UNIX)
    if(NOT APPLE)
        set(REALTIME_LIB rt)
        # We do not use pthreads directly anymore, but on Linux, C++11 threads
        # are implemented with pthreads, and we must link to it.
        set(PTHREAD_LIB pthread)
    endif()
    set(MATH_LIBRARY m)
    set(DL_LIBRARY dl)
endif()

set(MATH_LIBS_TO_USE    ${LAPACK_BEING_USED} ${PTHREAD_LIB}
                        ${REALTIME_LIB} ${DL_LIBRARY} ${MATH_LIBRARY})
set(MATH_LIBS_TO_USE_VN ${MATH_LIBS_TO_USE})

#
# Allow automated build and dashboard.
#
include(Dart)
## When in Debug mode and running valgrind, some of the test
## cases take longer than the default 1500 seconds.
set(DART_TESTING_TIMEOUT 7200)

if(BUILD_TESTING)
    #
    # Testing
    #
    enable_testing()

    # Make a RUN_TESTS_PARALLEL target (thanks, Kevin!)
    # Specify number of cores to run for testing
    set(TESTING_PROCESSOR_COUNT 4 CACHE STRING
        "Number of CPUs to be used by the RUN_TESTS_PARALLEL target.")
    mark_as_advanced(TESTING_PROCESSOR_COUNT)
    set(cmd ${CMAKE_CTEST_COMMAND} -j${TESTING_PROCESSOR_COUNT})
    if(MSVC OR XCODE)
        set(cmd ${cmd} -C ${CMAKE_CFG_INTDIR})
    else()
        set(cmd ${cmd} -C ${CMAKE_BUILD_TYPE})
    endif()
    add_custom_target(RUN_TESTS_PARALLEL COMMAND ${cmd})
endif()

# These are used in Doxyfile.in and SimbodyConfig.cmake.in.
set(SIMBODY_INSTALL_DOXYGENDIR   "${CMAKE_INSTALL_DOCDIR}/api")
set(SIMBODY_DOXYGEN_TAGFILE_NAME "SimbodyDoxygenTagfile")
add_subdirectory(doc)

# Specify where visualizer should be installed. This needs to be in the
# root CMakeLists.txt so the cmake config file can see this value.
#
# Also specify where include files are installed.
if(WIN32)
    # Install visualizer to bin, since it needs to be co-located with dll's
    set(SIMBODY_VISUALIZER_REL_INSTALL_DIR ${CMAKE_INSTALL_BINDIR})
    # Install include files into base include folder since it's a sandbox
    set(SIMBODY_INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_INCLUDEDIR})
else()
    # Visualizer is not intended to be a user executable. Proper place is
    # inside the lib directory
    set(SIMBODY_VISUALIZER_REL_INSTALL_DIR ${CMAKE_INSTALL_LIBEXECDIR}/simbody)
    # Install include files in simbody subfolder to avoid polluting the
    # global build folder
    set(SIMBODY_INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_INCLUDEDIR}/simbody)
endif()
set(SIMBODY_VISUALIZER_INSTALL_DIR
    ${CMAKE_INSTALL_PREFIX}/${SIMBODY_VISUALIZER_REL_INSTALL_DIR})

# Each of these returns a list of API include directories for
# use by the later builds.
add_custom_target(PlatformFiles ALL)
add_subdirectory( Platform )
# PLATFORM_INCLUDE_DIRECTORIES now set
add_subdirectory( SimTKcommon )
# SimTKCOMMON_INCLUDE_DIRECTORIES now set
add_subdirectory( SimTKmath )
# SimTKMATH_INCLUDE_DIRECTORIES now set
add_subdirectory( Simbody )
# SimTKSIMBODY_INCLUDE_DIRECTORIES now set(but not used)
# GUI_NAME now set

if( BUILD_EXAMPLES )
    add_subdirectory( examples )
endif()

file(GLOB TOPLEVEL_DOCS LICENSE.txt *.md doc/*.pdf doc/*.md)
install(FILES ${TOPLEVEL_DOCS} DESTINATION ${CMAKE_INSTALL_DOCDIR})

# Add uninstall target
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)
add_custom_target(uninstall
    "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake")

# Make the cmake config files
set(PKG_NAME ${PROJECT_NAME})
set(PKG_LIBRARIES
    ${SimTKSIMBODY_LIBRARY_NAME}
    ${SimTKMATH_LIBRARY_NAME}
    ${SimTKCOMMON_LIBRARY_NAME}
    )

if(WIN32)
    set(SIMBODY_CMAKE_DIR cmake)
elseif(UNIX)
    set(SIMBODY_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/simbody/)
endif()

# Configure SimbodyConfig.cmake in a way that allows the installation to be
# relocatable.
include(CMakePackageConfigHelpers)
configure_package_config_file(
    ${CMAKE_SOURCE_DIR}/cmake/SimbodyConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/SimbodyConfigForInstall.cmake
    INSTALL_DESTINATION "${SIMBODY_CMAKE_DIR}"
    PATH_VARS # Variables to edit in the SimbodyConfig.cmake.in.
        CMAKE_INSTALL_PREFIX
        SIMBODY_INCLUDE_INSTALL_DIR
        CMAKE_INSTALL_LIBDIR
        CMAKE_INSTALL_BINDIR
        SIMBODY_VISUALIZER_INSTALL_DIR
        SIMBODY_INSTALL_DOXYGENDIR
    )
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cmake/SimbodyConfigForInstall.cmake
        DESTINATION ${SIMBODY_CMAKE_DIR}
        RENAME SimbodyConfig.cmake)

# Create a file that allows clients to Simbody to ensure they have the version
# of Simbody they want.
# Requires CMake 2.8.8.
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/SimbodyConfigVersion.cmake
    VERSION "${SIMBODY_VERSION}"
    COMPATIBILITY SameMajorVersion)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cmake/SimbodyConfigVersion.cmake
        DESTINATION ${SIMBODY_CMAKE_DIR})

install(EXPORT SimbodyTargets DESTINATION
        "${SIMBODY_CMAKE_DIR}") #optionally: COMPONENT dev

# Install a sample CMakeLists.txt that uses SimbodyConfig.cmake.
install(FILES ${CMAKE_SOURCE_DIR}/cmake/SampleCMakeLists.txt
        DESTINATION ${SIMBODY_CMAKE_DIR})

# Make the pkgconfig file: select the proper flags depending on the platform
if(WIN32)
    if(${PLATFORM_ABI} MATCHES "x64")
        # win 64 bits
        set(PKGCONFIG_PLATFORM_LIBS "-lliblapack -llibblas")
    else()
        set(PKGCONFIG_PLATFORM_LIBS "-lliblapack -llibblas")
    endif()
elseif(APPLE)
    set(PKGCONFIG_PLATFORM_LIBS     "-llapack -lblas -ldl")
else()
    set(PKGCONFIG_PLATFORM_LIBS     "-llapack -lblas -lpthread -lrt -ldl -lm")
endif()

configure_file(${CMAKE_SOURCE_DIR}/cmake/pkgconfig/simbody.pc.in
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/pkgconfig/simbody.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cmake/pkgconfig/simbody.pc
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig/)
