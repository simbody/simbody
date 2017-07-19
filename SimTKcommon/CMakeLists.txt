#---------------------------------------------------
# SimTKcommon
#
# Creates SimTK Core library, base name=SimTKcommon.
# Default libraries are shared & optimized. Variants
# are created for static (_static) and debug (_d) and
# provision is made for an optional "namespace" (ns)
# and version number (vn).
#
# Windows:
#   [ns_]SimTKcommon[_vn][_d].dll
#   [ns_]SimTKcommon[_vn][_d].lib
#   [ns_]SimTKcommon[_vn]_static[_d].lib
# Linux/Cygwin:
#   lib[ns_]SimTKcommon[_vn][_d].so
#   lib[ns_]SimTKcommon[_vn]_static[_d].a
# Mac:
#   lib[ns_]SimTKcommon[_vn][_d].dylib
#   lib[ns_]SimTKcommon[_vn]_static[_d].a
#
# Targets are installed in
#   %ProgramFiles%\SimTK\lib,bin       (Win32)
#   %ProgramFiles(x86)%\SimTK\lib,bin  (32 bit target on Win64)
#   %ProgramW6432%\SimTK\lib,bin       (64 bit target on Win64)
#   /usr/local/SimTK/lib[64]           (Linux, Mac, Cygwin)
#
#----------------------------------------------------

project(SimTKcommon)

# SimTKcommon depends on the Platform files but nothing else.
include_directories(${PLATFORM_INCLUDE_DIRECTORIES})

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
set(SimTKCOMMON_SOURCE_SUBDIRS
    . Scalar SmallMatrix Mechanics BigMatrix
    Geometry Simulation Random Polynomial)

# Collect up information about the version of the SimTKcommon library
# we're building and make it available to the code so it can be built
# into the binaries. This also determines the versioned library names
# in which case all dependencies must use the same version.

set(SimTKCOMMON_MAJOR_VERSION ${SIMBODY_MAJOR_VERSION})
set(SimTKCOMMON_MINOR_VERSION ${SIMBODY_MINOR_VERSION})
set(SimTKCOMMON_PATCH_VERSION ${SIMBODY_PATCH_VERSION})

# Report the version number to the CMake UI. Don't include the
# build version if it is zero.
set(PATCH_VERSION_STRING)
if(SimTKCOMMON_PATCH_VERSION)
    set(PATCH_VERSION_STRING ".${SimTKCOMMON_PATCH_VERSION}")
endif()

set(SimTKCOMMON_VERSION
    "${SimTKCOMMON_MAJOR_VERSION}.${SimTKCOMMON_MINOR_VERSION}${PATCH_VERSION_STRING}")

# This is the suffix if we're building and depending on versioned libraries.
set(VN "_${SimTKCOMMON_VERSION}")


set(SimTKCOMMON_COPYRIGHT_YEARS "2005-10")

# underbar separated list of dotted authors, no spaces or commas
set(SimTKCOMMON_AUTHORS "Michael.Sherman_Peter.Eastman")

set(SimTKCOMMON_SVN_REVISION ${SIMBODY_SVN_REVISION})

add_definitions(-DSimTK_SimTKCOMMON_LIBRARY_NAME=${SimTKCOMMON_LIBRARY_NAME}
                -DSimTK_SimTKCOMMON_MAJOR_VERSION=${SimTKCOMMON_MAJOR_VERSION}
                -DSimTK_SimTKCOMMON_MINOR_VERSION=${SimTKCOMMON_MINOR_VERSION}
                -DSimTK_SimTKCOMMON_PATCH_VERSION=${SimTKCOMMON_PATCH_VERSION})

if(NEED_QUOTES)
   add_definitions(-DSimTK_SimTKCOMMON_SVN_REVISION="${SimTKCOMMON_SVN_REVISION}"
                   -DSimTK_SimTKCOMMON_COPYRIGHT_YEARS="${SimTKCOMMON_COPYRIGHT_YEARS}"
                   -DSimTK_SimTKCOMMON_AUTHORS="${SimTKCOMMON_AUTHORS}")
else(NEED_QUOTES)
   add_definitions(-DSimTK_SimTKCOMMON_SVN_REVISION=${SimTKCOMMON_SVN_REVISION}
                   -DSimTK_SimTKCOMMON_COPYRIGHT_YEARS=${SimTKCOMMON_COPYRIGHT_YEARS}
                   -DSimTK_SimTKCOMMON_AUTHORS=${SimTKCOMMON_AUTHORS})
endif(NEED_QUOTES)

# -DSimTK_SimTKCOMMON_LIBRARY_TYPE has to be defined in the target subdirectories.
# -DSimTKcommon_EXPORTS defined automatically when Windows DLL build is being done.

set(SHARED_TARGET ${SimTKCOMMON_LIBRARY_NAME})
set(STATIC_TARGET ${SimTKCOMMON_LIBRARY_NAME}_static)
set(SHARED_TARGET_VN ${SimTKCOMMON_LIBRARY_NAME}${VN})
set(STATIC_TARGET_VN ${SimTKCOMMON_LIBRARY_NAME}${VN}_static)

## Test against the unversioned libraries if they are being build;
## otherwise against the versioned libraries.
if(BUILD_UNVERSIONED_LIBRARIES)
    set(TEST_SHARED_TARGET ${SHARED_TARGET})
    set(TEST_STATIC_TARGET ${STATIC_TARGET})
else(BUILD_UNVERSIONED_LIBRARIES)
    set(TEST_SHARED_TARGET ${SHARED_TARGET_VN})
    set(TEST_STATIC_TARGET ${STATIC_TARGET_VN})
endif(BUILD_UNVERSIONED_LIBRARIES)


# These are all the places to search for header files which are
# to be part of the API.
set(API_INCLUDE_DIRS) # start empty
set(SimTKCOMMON_INCLUDE_DIRS) # start empty
foreach(subdir ${SimTKCOMMON_SOURCE_SUBDIRS})
    # append
    set(API_INCLUDE_DIRS ${API_INCLUDE_DIRS}
     ${PROJECT_SOURCE_DIR}/${subdir}/include
     ${PROJECT_SOURCE_DIR}/${subdir}/include/SimTKcommon
     ${PROJECT_SOURCE_DIR}/${subdir}/include/SimTKcommon/internal)
    # Referencing headers must always be done relative to this level.
    set(SimTKCOMMON_INCLUDE_DIRS ${SimTKCOMMON_INCLUDE_DIRS}
     ${PROJECT_SOURCE_DIR}/${subdir}/include)
endforeach(subdir)


# Include the SimTKcommon API include directories now so that SimTKcommon code
# can use them.
include_directories(${SimTKCOMMON_INCLUDE_DIRS})

# Pass up the include directories list to the parent so
# subsequent libraries can use them.
set(SimTKCOMMON_INCLUDE_DIRECTORIES ${SimTKCOMMON_INCLUDE_DIRS}
    PARENT_SCOPE)

# We'll need both *relative* path names, starting with their API_INCLUDE_DIRS,
# and absolute pathnames.
set(API_REL_INCLUDE_FILES)   # start these out empty
set(API_ABS_INCLUDE_FILES)

foreach(dir ${API_INCLUDE_DIRS})
    file(GLOB fullpaths ${dir}/*.h)    # returns full pathnames
    set(API_ABS_INCLUDE_FILES ${API_ABS_INCLUDE_FILES} ${fullpaths})

    foreach(pathname ${fullpaths})
        get_filename_component(filename ${pathname} NAME)
        set(API_REL_INCLUDE_FILES ${API_REL_INCLUDE_FILES} ${dir}/${filename})
    endforeach(pathname)
endforeach(dir)

# collect up source files
set(SOURCE_FILES) # empty
set(SOURCE_INCLUDE_FILES)

foreach(subdir ${SimTKCOMMON_SOURCE_SUBDIRS})
    file(GLOB src_files  ${subdir}/src/*.cpp)
    file(GLOB incl_files ${subdir}/src/*.h)
    set(SOURCE_FILES         ${SOURCE_FILES}         ${src_files})   #append
    set(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})
endforeach(subdir)

#
# Installation
#
# libraries are installed from their subdirectories; headers here

# install headers
file(GLOB CORE_HEADERS     include/*.h                      */include/*.h)
file(GLOB TOP_HEADERS      include/SimTKcommon/*.h          */include/SimTKcommon/*.h)
file(GLOB INTERNAL_HEADERS include/SimTKcommon/internal/*.h */include/SimTKcommon/internal/*.h)
install(FILES ${CORE_HEADERS} DESTINATION ${SIMBODY_INCLUDE_INSTALL_DIR})
install(FILES ${TOP_HEADERS} DESTINATION ${SIMBODY_INCLUDE_INSTALL_DIR}/SimTKcommon)
install(FILES ${INTERNAL_HEADERS} DESTINATION ${SIMBODY_INCLUDE_INSTALL_DIR}/SimTKcommon/internal)

file(GLOB SIMTKCOMMON_DOCS doc/*.pdf doc/*.txt doc/*.md)
install(FILES ${SIMTKCOMMON_DOCS} DESTINATION ${CMAKE_INSTALL_DOCDIR})

# These are at the end because we want them processed after
# all the various variables have been set above.

if(BUILD_STATIC_LIBRARIES)
    add_subdirectory( staticTarget )
endif()
if(BUILD_DYNAMIC_LIBRARIES)
    add_subdirectory( sharedTarget )
endif()

if( BUILD_TESTING )
    add_subdirectory( tests )
endif( BUILD_TESTING )

