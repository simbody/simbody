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

PROJECT (SimTKcommon)

# SimTKcommon depends on the Platform files but nothing else.
INCLUDE_DIRECTORIES(${PLATFORM_INCLUDE_DIRECTORIES})

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
SET(SimTKCOMMON_SOURCE_SUBDIRS
    . Scalar SmallMatrix Mechanics BigMatrix
    Geometry Simulation Random Polynomial)

# Collect up information about the version of the SimTKcommon library 
# we're building and make it available to the code so it can be built 
# into the binaries. This also determines the versioned library names
# in which case all dependencies must use the same version.

SET(SimTKCOMMON_MAJOR_VERSION ${SIMBODY_MAJOR_VERSION})
SET(SimTKCOMMON_MINOR_VERSION ${SIMBODY_MINOR_VERSION})
SET(SimTKCOMMON_PATCH_VERSION ${SIMBODY_PATCH_VERSION})

# Report the version number to the CMake UI. Don't include the 
# build version if it is zero.
SET(PATCH_VERSION_STRING)
IF(SimTKCOMMON_PATCH_VERSION)
    SET(PATCH_VERSION_STRING ".${SimTKCOMMON_PATCH_VERSION}")
ENDIF()

SET(SimTKCOMMON_VERSION 
    "${SimTKCOMMON_MAJOR_VERSION}.${SimTKCOMMON_MINOR_VERSION}${PATCH_VERSION_STRING}")

# This is the suffix if we're building and depending on versioned libraries.
SET(VN "_${SimTKCOMMON_VERSION}")


SET(SimTKCOMMON_COPYRIGHT_YEARS "2005-10")

# underbar separated list of dotted authors, no spaces or commas
SET(SimTKCOMMON_AUTHORS "Michael.Sherman_Peter.Eastman")

SET (SimTKCOMMON_SVN_REVISION ${SIMBODY_SVN_REVISION}) 

ADD_DEFINITIONS(-DSimTK_SimTKCOMMON_LIBRARY_NAME=${SimTKCOMMON_LIBRARY_NAME}
                -DSimTK_SimTKCOMMON_MAJOR_VERSION=${SimTKCOMMON_MAJOR_VERSION}
                -DSimTK_SimTKCOMMON_MINOR_VERSION=${SimTKCOMMON_MINOR_VERSION}
                -DSimTK_SimTKCOMMON_PATCH_VERSION=${SimTKCOMMON_PATCH_VERSION})

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

SET(SHARED_TARGET ${SimTKCOMMON_LIBRARY_NAME})
SET(STATIC_TARGET ${SimTKCOMMON_LIBRARY_NAME}_static)
SET(SHARED_TARGET_VN ${SimTKCOMMON_LIBRARY_NAME}${VN})
SET(STATIC_TARGET_VN ${SimTKCOMMON_LIBRARY_NAME}${VN}_static)

## Test against the unversioned libraries if they are being build;
## otherwise against the versioned libraries.
IF(BUILD_UNVERSIONED_LIBRARIES)
    SET(TEST_SHARED_TARGET ${SHARED_TARGET})
    SET(TEST_STATIC_TARGET ${STATIC_TARGET})
ELSE(BUILD_UNVERSIONED_LIBRARIES)
    SET(TEST_SHARED_TARGET ${SHARED_TARGET_VN})
    SET(TEST_STATIC_TARGET ${STATIC_TARGET_VN})
ENDIF(BUILD_UNVERSIONED_LIBRARIES)


# These are all the places to search for header files which are
# to be part of the API.
SET(API_INCLUDE_DIRS) # start empty
SET(SimTKCOMMON_INCLUDE_DIRS) # start empty
FOREACH(subdir ${SimTKCOMMON_SOURCE_SUBDIRS})
    # append
    SET(API_INCLUDE_DIRS ${API_INCLUDE_DIRS}
     ${PROJECT_SOURCE_DIR}/${subdir}/include 
     ${PROJECT_SOURCE_DIR}/${subdir}/include/SimTKcommon 
     ${PROJECT_SOURCE_DIR}/${subdir}/include/SimTKcommon/internal)
    # Referencing headers must always be done relative to this level.
    SET(SimTKCOMMON_INCLUDE_DIRS ${SimTKCOMMON_INCLUDE_DIRS}
     ${PROJECT_SOURCE_DIR}/${subdir}/include)
ENDFOREACH(subdir)


# Include the SimTKcommon API include directories now so that SimTKcommon code 
# can use them.
INCLUDE_DIRECTORIES(${SimTKCOMMON_INCLUDE_DIRS})

# Pass up the include directories list to the parent so
# subsequent libraries can use them.
SET(SimTKCOMMON_INCLUDE_DIRECTORIES ${SimTKCOMMON_INCLUDE_DIRS}
    PARENT_SCOPE)

# We'll need both *relative* path names, starting with their API_INCLUDE_DIRS,
# and absolute pathnames.
SET(API_REL_INCLUDE_FILES)   # start these out empty
SET(API_ABS_INCLUDE_FILES)

FOREACH(dir ${API_INCLUDE_DIRS})
    FILE(GLOB fullpaths ${dir}/*.h)    # returns full pathnames
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
ENDFOREACH(subdir)

#
# Installation
#
# libraries are installed from their subdirectories; headers here

# install headers
FILE(GLOB CORE_HEADERS     include/*.h                      */include/*.h)
FILE(GLOB TOP_HEADERS      include/SimTKcommon/*.h          */include/SimTKcommon/*.h)
FILE(GLOB INTERNAL_HEADERS include/SimTKcommon/internal/*.h */include/SimTKcommon/internal/*.h)
INSTALL_FILES(/${SIMBODY_INCLUDE_INSTALL_DIR}/                     FILES ${CORE_HEADERS})
INSTALL_FILES(/${SIMBODY_INCLUDE_INSTALL_DIR}/SimTKcommon/         FILES ${TOP_HEADERS})
INSTALL_FILES(/${SIMBODY_INCLUDE_INSTALL_DIR}/SimTKcommon/internal FILES ${INTERNAL_HEADERS})

FILE(GLOB SIMTKCOMMON_DOCS doc/*.pdf doc/*.txt doc/*.md)
INSTALL(FILES ${SIMTKCOMMON_DOCS} DESTINATION ${CMAKE_INSTALL_DOCDIR})

# These are at the end because we want them processed after
# all the various variables have been set above.

IF(BUILD_STATIC_LIBRARIES)
    ADD_SUBDIRECTORY( staticTarget )
ENDIF()
ADD_SUBDIRECTORY( sharedTarget )

IF( BUILD_TESTING )
  ADD_SUBDIRECTORY( tests )
ENDIF( BUILD_TESTING )

