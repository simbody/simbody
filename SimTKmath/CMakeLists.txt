#---------------------------------------------------
# Simmath 
#
# Creates SimTK Core library, base name=SimTKmath.
# Default libraries are shared & optimized. Variants
# are created for static (_static) and debug (_d) and
# provision is made for an optional "namespace" (ns)
# and version number (vn).
#
# Windows:
#   [ns_]SimTKmath[_vn][_d].dll
#   [ns_]SimTKmath[_vn][_d].lib
#   [ns_]SimTKmath[_vn]_static[_d].lib
# Unix:
#   lib[ns_]SimTKmath[_vn][_d].so
#   lib[ns_]SimTKmath[_vn]_static[_d].a
#
# All libraries are installed in 
#   %ProgramFiles%\SimTK\lib  (Windows)
#   /usr/local/SimTK/lib        (UNIX)
#
#----------------------------------------------------

PROJECT (SimTKmath)

# SimTKmath depends on PlatformFiles and SimTKcommon only.
INCLUDE_DIRECTORIES(${PLATFORM_INCLUDE_DIRECTORIES}
	            ${SimTKCOMMON_INCLUDE_DIRECTORIES})

# code is in "./src" and in each of these subdirectories/src
SET(SIMMATH_SOURCE_SUBDIRS LinearAlgebra Integrators Optimizers Geometry)

# Collect up information about the version of the SimTKmath library 
# we're building and make it available to the code so it can be built 
# into the binaries. This also determines the versioned library names
# in which case all dependencies must use the same version.

SET(SIMMATH_MAJOR_VERSION ${SIMBODY_MAJOR_VERSION})
SET(SIMMATH_MINOR_VERSION ${SIMBODY_MINOR_VERSION})
SET(SIMMATH_PATCH_VERSION ${SIMBODY_PATCH_VERSION})


# Report the version number to the CMake UI. Don't include the 
# patch version if it is zero.
SET(PATCH_VERSION_STRING)
IF(SIMMATH_PATCH_VERSION)
    SET(PATCH_VERSION_STRING ".${SIMMATH_PATCH_VERSION}")
ENDIF()

SET(SIMMATH_VERSION 
    "${SIMMATH_MAJOR_VERSION}.${SIMMATH_MINOR_VERSION}${PATCH_VERSION_STRING}")

# This is the suffix if we're building and depending on versioned libraries.
SET(VN "_${SIMMATH_VERSION}")

SET(SIMMATH_COPYRIGHT_YEARS "2005-10")

# underbar separated list of dotted authors, no spaces or commas
SET(SIMMATH_AUTHORS         "Jack.Middleton_Michael.Sherman")

SET (SIMMATH_SVN_REVISION ${SIMBODY_SVN_REVISION}) 

ADD_DEFINITIONS(-DSimTK_SIMMATH_LIBRARY_NAME=${SimTKMATH_LIBRARY_NAME}
                -DSimTK_SIMMATH_MAJOR_VERSION=${SIMMATH_MAJOR_VERSION}
                -DSimTK_SIMMATH_MINOR_VERSION=${SIMMATH_MINOR_VERSION}
		-DSimTK_SIMMATH_PATCH_VERSION=${SIMMATH_PATCH_VERSION})

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

SET(SHARED_TARGET ${SimTKMATH_LIBRARY_NAME})
SET(STATIC_TARGET ${SimTKMATH_LIBRARY_NAME}_static)
SET(SHARED_TARGET_VN ${SimTKMATH_LIBRARY_NAME}${VN})
SET(STATIC_TARGET_VN ${SimTKMATH_LIBRARY_NAME}${VN}_static)

# On Unix or Cygwin we have to add the debug suffix manually
IF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)
    SET(SHARED_TARGET ${SHARED_TARGET}_d)
    SET(STATIC_TARGET ${STATIC_TARGET}_d)
    SET(SHARED_TARGET_VN ${SHARED_TARGET_VN}_d)
    SET(STATIC_TARGET_VN ${STATIC_TARGET_VN}_d)
ENDIF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)

## Test against the unversioned libraries if they are being built;
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
SET(SimTKMATH_INCLUDE_DIRS) # start empty
FOREACH(subdir . ${SIMMATH_SOURCE_SUBDIRS})
    # append
    SET(API_INCLUDE_DIRS ${API_INCLUDE_DIRS}
	 ${PROJECT_SOURCE_DIR}/${subdir}/include 
	 ${PROJECT_SOURCE_DIR}/${subdir}/include/simmath 
	 ${PROJECT_SOURCE_DIR}/${subdir}/include/simmath/internal)

    # Referencing headers must always be done relative to this level.
    SET(SimTKMATH_INCLUDE_DIRS ${SimTKMATH_INCLUDE_DIRS}
	 ${PROJECT_SOURCE_DIR}/${subdir}/include)
ENDFOREACH(subdir)

# Include the Simmath API include directories now so that Simmath code 
# can use them.
INCLUDE_DIRECTORIES(${SimTKMATH_INCLUDE_DIRS})

# And pass API include directories up to the parent so subsequent libraries
# can find the headers too.
SET(SimTKMATH_INCLUDE_DIRECTORIES ${SimTKMATH_INCLUDE_DIRS}
    PARENT_SCOPE)


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
SET(SOURCE_INCLUDE_DIRS) # in case some low-level source needs a -I include dir
SET(SOURCE_GROUPS)
SET(SOURCE_GROUP_FILES)

# first process all the source subdirectories
FOREACH(subdir ${SIMMATH_SOURCE_SUBDIRS})
    add_subdirectory(${subdir})
ENDFOREACH(subdir)

if (SOURCE_GROUPS)
    list(LENGTH SOURCE_GROUPS NGROUPS)
    math(EXPR lastgrpnum ${NGROUPS}-1)
    foreach( grpnum RANGE 0 ${lastgrpnum} )
        list(GET SOURCE_GROUPS ${grpnum} grp)
        list(GET SOURCE_GROUP_FILES ${grpnum} grpfile)
        source_group("${grp}" FILES "${grpfile}")
    endforeach()
endif()

# then process ./src
FILE(GLOB src_files  ./src/*.cpp)
FILE(GLOB incl_files ./src/*.h)
SET(SOURCE_FILES         ${SOURCE_FILES}         ${src_files})   #append
SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})

# Add low-level source include directories if any.
INCLUDE_DIRECTORIES(${SOURCE_INCLUDE_DIRS})

#
# Installation
#

# libraries and examples are installed from their subdirectories; headers here

# install headers
FILE(GLOB CORE_HEADERS     include/*.h                  */include/*.h)
FILE(GLOB TOP_HEADERS      include/simmath/*.h          */include/simmath/*.h)
FILE(GLOB INTERNAL_HEADERS include/simmath/internal/*.h */include/simmath/internal/*.h)

INSTALL_FILES(/${SIMBODY_INCLUDE_INSTALL_DIR}/                   FILES ${CORE_HEADERS})
INSTALL_FILES(/${SIMBODY_INCLUDE_INSTALL_DIR}/SimTKmath/         FILES ${TOP_HEADERS})
INSTALL_FILES(/${SIMBODY_INCLUDE_INSTALL_DIR}/SimTKmath/internal FILES ${INTERNAL_HEADERS})

FILE(GLOB SIMMATH_DOCS doc/*.pdf doc/*.txt)
INSTALL(FILES ${SIMMATH_DOCS} DESTINATION ${CMAKE_INSTALL_DOCDIR})

# These are at the end because we want them processed after
# all the various variables have been set above.

IF(BUILD_STATIC_LIBRARIES)
    ADD_SUBDIRECTORY( staticTarget )
ENDIF()
ADD_SUBDIRECTORY( sharedTarget )

IF( BUILD_EXAMPLES )
  ADD_SUBDIRECTORY( examples )
ENDIF( BUILD_EXAMPLES )

IF( BUILD_TESTING )
  ADD_SUBDIRECTORY( tests )
ENDIF( BUILD_TESTING )
