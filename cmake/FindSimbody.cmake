# FindSimbody.cmake
#
# Simbios National Center for Physics Based Simulation of Biological Structures
# Stanford University
# This cmake file created 2012 by Michael Sherman and is in the public
# domain; Simbody itself is open source under the Apache 2.0 license.
#
# This is a CMake "find" module that will try to find the Simbody multibody
# dynamics (physics) package installed somewhere on your computer. Simbody
# is part of the SimTK biosimulation toolkit. For more information, see
# https://simtk.org/home/simbody.
#
# To use this in a find_package() command from your own CMakeLists.txt file,
# make sure this file is in a directory that is in the CMAKE_MODULE_PATH.
# You can add a directory to that path with a line like this:
#     list(APPEND CMAKE_MODULE_PATH "myModuleDir")
#
# Then, to use Simbody include these lines:
#
#     find_package(Simbody REQUIRED)
#     include_directories(${Simbody_INCLUDE_DIR})
#     link_directories(${Simbody_LIB_DIR})
#     add_executable(myexe ${my_source_files} ${my_header_files})
#     target_link_libraries(myexe ${Simbody_LIBRARIES})
#
# If you don't want to make it REQUIRED, you can check Simbody_FOUND after
# find_package() returns.
# TODO: no version selection is implemented here yet; if you provide it
# to find_package() it will be ignored.
#
# This includes several libraries:
#     SimTKsimbody
#     SimTKmath
#     SimTKcommon
#     Windows only: liblapack libblas pthreadVC2[_x64]
# The above will be in Simbody_ROOT_DIR/lib.
#
# On Mac and Linux we don't provide our own lapack but expect it to be
# available.
#     Mac/Linux only: lapack blas
#
# Once done this will define:
#
#   Simbody_FOUND - Whether search for Simbody libraries and headers succeeded.
#   Simbody_ROOT_DIR - the installation directory; all the pieces must be
#                      found together
#   Simbody_INCLUDE_DIR - location of Simbody.h
#   Simbody_LIB_DIR     - location of libSimTKsimbody.{a,so,dylib} or SimTKsimbody.lib
#   Simbody_BIN_DIR     - location of VisualizerGUI and .dll's on Windows
#   Simbody_LIBRARIES   - suitable for target_link_libraries(); includes
#                           both optimized and debug libraries if both are
#                           available
#   Simbody_STATIC_LIBRARIES - suitable for target_link_libraries(); includes
#                              both optimized and debug static libraries if
#                              both are available
#
# The following environment variables are used if set, in order of decreasing
# preference:
#   SIMBODY_HOME
#   SimTK_INSTALL_DIR
#
# Otherwise we look in platform-specific standard places, in this order:
#   <standardPlaces>/Simbody, /simbody
#   <standardPlaces>/SimTK, /simtk
#
# This module looks for certain CMake variables on input and behaves
# accordingly if they are present:
#
#   SimTK_INSTALL_DIR
#       This is commonly set by other SimTK software and overrides the
#       environment variables if present. Note that this has the same name
#       as one of the environment variables but is distinct.
#   SimTK_SDK
#       If this is set (probably by the nightly build system) then that
#       will be the only place we'll look for Simbody, overriding environment
#       variables if they are set, and overriding the above CMake variable.

cmake_minimum_required(VERSION 2.8)

# Get values of relevant environment variables for convenient testing.
set(ENV_SIMBODY_HOME_VALUE $ENV{SIMBODY_HOME})
set(ENV_SimTK_INSTALL_DIR_VALUE $ENV{SimTK_INSTALL_DIR})

if(SimTK_SDK) # a CMake variable, not an environment variable
    # This is the SimTK.org nightly builds location. If set it overrides
    # everything else.
    set(Simbody_SEARCH_PATHS "${SimTK_SDK}")
elseif(SimTK_INSTALL_DIR) # a CMake variable, not an environment variable
    set(Simbody_SEARCH_PATHS "${SimTK_INSTALL_DIR}")
elseif(ENV_SIMBODY_HOME_VALUE)
    set(Simbody_SEARCH_PATHS "${ENV_SIMBODY_HOME_VALUE}")
elseif(ENV_SimTK_INSTALL_DIR_VALUE) # != the eponymous CMake variable.
    set(Simbody_SEARCH_PATHS "${ENV_SimTK_INSTALL_DIR_VALUE}")
else()
    # Nothing is set so we'll have to hunt.
    set(Simbody_SEARCH_PATHS)

    # UNIX includes Mac, Linux, and Cygwin
    if(UNIX)
        list(APPEND Simbody_SEARCH_PATHS /usr/local)
    endif()

    if(APPLE) # Mac only
        list(APPEND Simbody_SEARCH_PATHS /Developer)
    endif()

    # WIN32 includes Windows 32 & 64 bit, and Cygwin
    if(WIN32)
        if( ${CMAKE_SIZEOF_VOID_P} EQUAL 8 )
	      # 64 bit target on Win64
	      set(PROGFILE_DIR "$ENV{ProgramW6432}")
        else() # Target is 32 bit
	      # present if 64bit Windows
	      set(PROGFILE_DIR "$ENV{ProgramFiles(x86)}")
	      if(NOT PROGFILE_DIR)
	        set(PROGFILE_DIR "$ENV{ProgramFiles}") # on 32bit Windows
	      endif()
        endif()
        list(APPEND Simbody_SEARCH_PATHS ${PROGFILE_DIR})
    endif()
endif(SimTK_SDK)

# We'll use the main Simbody header as the key to finding the installation
# root directory. We're assuming the header is in
#    Simbody_INCLUDE_DIR == ${Simbody_ROOT_DIR}/include
# So stripping off the "/include" should give the root directory.

# Force this to be recalculated every time.
set(Simbody_INCLUDE_DIR "Simbody_INCLUDE_DIR-NOTFOUND" CACHE PATH
    "The Simbody and SimTK include directory." FORCE)

message("${Simbody_SEARCH_PATHS}")


foreach(pth IN LISTS Simbody_SEARCH_PATHS)
  find_path(Simbody_INCLUDE_DIR
    NAMES "SimTKsimbody.h" "Simbody.h"
    PATHS "${pth}"
    PATH_SUFFIXES "include"
                  "Simbody/include" "simbody/include"
                  "SimTK/include" "simtk/include"
    NO_DEFAULT_PATH
    DOC "Location of top-level installed Simbody header files"
  )
endforeach()

# This will only be executed if the first loop fails. We're getting
# desperate!
find_path(Simbody_INCLUDE_DIR
    NAMES "SimTKsimbody.h" "Simbody.h"
    PATH_SUFFIXES "include"
                  "Simbody/include" "simbody/include"
                  "SimTK/include" "simtk/include"
    DOC "Location of top-level installed Simbody header files"
)

get_filename_component(Simbody_ROOT_DIR_TEMP "${Simbody_INCLUDE_DIR}" PATH)
set(Simbody_ROOT_DIR "${Simbody_ROOT_DIR_TEMP}" CACHE PATH
    "Where we found Simbody; use SimTK_INSTALL_DIR to change." FORCE)

set(Simbody_LIB_DIR ${Simbody_ROOT_DIR}/lib CACHE PATH
    "Location of Simbody and SimTK libraries." FORCE)
set(Simbody_BIN_DIR ${Simbody_ROOT_DIR}/bin CACHE PATH
    "Location of Simbody-related executables and (on Windows) dlls." FORCE)

set(Simbody_LIBRARY_LIST SimTKsimbody;SimTKmath;SimTKcommon)

if(WIN32)
    set(Simbody_LAPACK_LIBRARY_LIST liblapack;libblas)
else()
    set(Simbody_LAPACK_LIBRARY_LIST lapack;blas)
endif()


if(WIN32)
    if( ${CMAKE_SIZEOF_VOID_P} EQUAL 8 )
        set(Simbody_EXTRA_LIBRARY_LIST pthreadVC2_x64)
    else()
        set(Simbody_EXTRA_LIBRARY_LIST pthreadVC2)
    endif()
elseif(APPLE)
    set(Simbody_EXTRA_LIBRARY_LIST pthread;dl)
else()
    set(Simbody_EXTRA_LIBRARY_LIST pthread;rt;dl;m)
endif()


# Find out which of the libraries are available.
find_library(Simbody_LIBRARY NAMES SimTKsimbody
    PATHS ${Simbody_LIB_DIR}
    DOC "This is the main Simbody library."
    NO_DEFAULT_PATH)
find_library(Simbody_STATIC_LIBRARY NAMES SimTKsimbody_static
    PATHS ${Simbody_LIB_DIR}
    DOC "This is the main Simbody static library."
    NO_DEFAULT_PATH)
find_library(Simbody_DEBUG_LIBRARY NAMES SimTKsimbody_d
    PATHS ${Simbody_LIB_DIR}
    DOC "This is the main Simbody debug library."
    NO_DEFAULT_PATH)
find_library(Simbody_STATIC_DEBUG_LIBRARY NAMES SimTKsimbody_static_d
    PATHS ${Simbody_LIB_DIR}
    DOC "This is the main Simbody static debug library."
    NO_DEFAULT_PATH)


# Set composite Simbody_LIBRARIES variable
set(LIBS)
if(Simbody_LIBRARY AND Simbody_DEBUG_LIBRARY)
    foreach(libname IN LISTS Simbody_LIBRARY_LIST)
        set(LIBS ${LIBS} optimized "${libname}" debug "${libname}_d")
    endforeach()
elseif(Simbody_LIBRARY)
    foreach(libname IN LISTS Simbody_LIBRARY_LIST)
        set(LIBS ${LIBS} "${libname}")
    endforeach()
elseif(Simbody_DEBUG_LIBRARY)
    foreach(libname IN LISTS Simbody_LIBRARY_LIST)
        set(LIBS ${LIBS} "${libname}_d")
    endforeach()
endif()

if(LIBS)
    foreach(lapack_lib IN LISTS Simbody_LAPACK_LIBRARY_LIST)
        set(LIBS ${LIBS} "${lapack_lib}")
    endforeach()
    foreach(extra_lib IN LISTS Simbody_EXTRA_LIBRARY_LIST)
        set(LIBS ${LIBS} "${extra_lib}")
    endforeach()
    set(Simbody_LIBRARIES ${LIBS} CACHE STRING
        "Simbody dynamic libraries" FORCE)
else()
    set(Simbody_LIBRARIES Simbody_LIBRARIES-NOTFOUND CACHE STRING
        "Simbody dynamic libraries" FORCE)
endif()

# Static library
set(LIBS)
if(Simbody_STATIC_LIBRARY AND Simbody_STATIC_DEBUG_LIBRARY)
    foreach(libname IN LISTS Simbody_LIBRARY_LIST)
        set(LIBS ${LIBS} optimized "${libname}_static"
		         debug     "${libname}_static_d")
    endforeach()
elseif(Simbody_STATIC_LIBRARY)
    foreach(libname IN LISTS Simbody_LIBRARY_LIST)
        set(LIBS ${LIBS} "${libname}_static")
    endforeach()
elseif(Simbody_STATIC_DEBUG_LIBRARY)
    foreach(libname IN LISTS Simbody_LIBRARY_LIST)
        set(LIBS ${LIBS} "${libname}_static_d")
    endforeach()
endif()

if(LIBS)
    # these aren't available in static
    foreach(lapack_lib IN LISTS Simbody_LAPACK_LIBRARY_LIST)
        set(LIBS ${LIBS} "${lapack_lib}")
    endforeach()
    foreach(extra_lib IN LISTS Simbody_EXTRA_LIBRARY_LIST)
        set(LIBS ${LIBS} "${extra_lib}")
    endforeach()
    set(Simbody_STATIC_LIBRARIES "${LIBS}" CACHE STRING
        "Simbody static libraries" FORCE)
else()
    set(Simbody_STATIC_LIBRARIES Simbody_STATIC_LIBRARIES-NOTFOUND CACHE STRING
        "Simbody static libraries" FORCE)
endif()

# This CMake-supplied script provides standard error handling.
include(FindPackageHandleStandardArgs OPTIONAL)
find_package_handle_standard_args(Simbody DEFAULT_MSG
	Simbody_INCLUDE_DIR)

# Not all the variables we produced need be returned.

unset(Simbody_LIBRARY CACHE)
unset(Simbody_DEBUG_LIBRARY CACHE)
unset(Simbody_STATIC_LIBRARY CACHE)
unset(Simbody_STATIC_DEBUG_LIBRARY CACHE)
mark_as_advanced(Simbody_LIBRARIES Simbody_STATIC_LIBRARIES)

mark_as_advanced(
    Simbody_ROOT_DIR
    Simbody_INCLUDE_DIR
    Simbody_BIN_DIR
    Simbody_LIB_DIR
)
