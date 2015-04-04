#---------------------------------------------------
# Simbody 
#
# Creates SimTK Core library, base name=SimTKsimbody.
# Default libraries are shared & optimized. Variants
# are created for static (_static) and debug (_d) and
# provision is made for an optional "namespace" (ns)
# and version number (vn).
#
# Windows:
#   [ns_]SimTKsimbody[_vn][_d].dll
#   [ns_]SimTKsimbody[_vn][_d].lib
#   [ns_]SimTKsimbody[_vn]_static[_d].lib
# Unix:
#   lib[ns_]SimTKsimbody[_vn][_d].so
#   lib[ns_]SimTKsimbody[_vn]_static[_d].a
#
# All libraries are installed in 
#   %ProgramFiles%\SimTK\lib[64]  (Windows)
#   /usr/local/SimTK/lib[64]        (UNIX)
#
#----------------------------------------------------

PROJECT (SimTKsimbody)

# SimTKsimbody depends on PlatformFiles, SimTKcommon, and SimTKmath only.
INCLUDE_DIRECTORIES(${PLATFORM_INCLUDE_DIRECTORIES}
                    ${SimTKCOMMON_INCLUDE_DIRECTORIES}
                    ${SimTKMATH_INCLUDE_DIRECTORIES})

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
# The BUILD_VISUALIZER variable determines whether we actually build the
# simbody-visualizer, but we still have to build the library-side support for
# the Visualizer so that Simbody tests won't fail to build.
SET(SIMBODY_SOURCE_SUBDIRS . Visualizer)

ADD_DEFINITIONS(-DSimTK_SIMBODY_LIBRARY_NAME=${SimTKSIMBODY_LIBRARY_NAME}
                -DSimTK_SIMBODY_MAJOR_VERSION=${SIMBODY_MAJOR_VERSION}
                -DSimTK_SIMBODY_MINOR_VERSION=${SIMBODY_MINOR_VERSION}
                -DSimTK_SIMBODY_PATCH_VERSION=${SIMBODY_PATCH_VERSION})


SET(SHARED_TARGET ${SimTKSIMBODY_LIBRARY_NAME})
SET(STATIC_TARGET ${SimTKSIMBODY_LIBRARY_NAME}_static)
SET(SHARED_TARGET_VN ${SimTKSIMBODY_LIBRARY_NAME}${VN})
SET(STATIC_TARGET_VN ${SimTKSIMBODY_LIBRARY_NAME}${VN}_static)

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
SET(SimTKSIMBODY_INCLUDE_DIRS) # start empty
FOREACH(subdir ${SIMBODY_SOURCE_SUBDIRS})
    # append
    SET(API_INCLUDE_DIRS ${API_INCLUDE_DIRS}
     ${PROJECT_SOURCE_DIR}/${subdir}/include 
     ${PROJECT_SOURCE_DIR}/${subdir}/include/simbody 
     ${PROJECT_SOURCE_DIR}/${subdir}/include/simbody/internal)

    # Referencing headers must always be done relative to this level.
    SET(SimTKSIMBODY_INCLUDE_DIRS ${SimTKSIMBODY_INCLUDE_DIRS}
     ${PROJECT_SOURCE_DIR}/${subdir}/include)
ENDFOREACH(subdir)

# Include the Simbody API include directories now so that Simbody code 
# can use them.
INCLUDE_DIRECTORIES(${SimTKSIMBODY_INCLUDE_DIRS})

# And pass API include directories up to the parent so subsequent libraries
# can find the headers too.
SET(SimTKSIMBODY_INCLUDE_DIRECTORIES ${SimTKSIMBODY_INCLUDE_DIRS}
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

FOREACH(subdir ${SIMBODY_SOURCE_SUBDIRS})
    FILE(GLOB src_files  ${subdir}/src/*.cpp ${subdir}/src/*/*.cpp)
    FILE(GLOB incl_files ${subdir}/src/*.h)
    SET(SOURCE_FILES         ${SOURCE_FILES}         ${src_files})   #append
    SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})

    #INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/${subdir}/include)
ENDFOREACH(subdir)

#INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src)

# libraries are installed from their subdirectories; headers here

# install headers
FILE(GLOB CORE_HEADERS     include/*.h                  */include/*.h)
FILE(GLOB TOP_HEADERS      include/simbody/*.h          */include/simbody/*.h)
FILE(GLOB INTERNAL_HEADERS include/simbody/internal/*.h */include/simbody/internal/*.h)
INSTALL_FILES(/${SIMBODY_INCLUDE_INSTALL_DIR}/                 FILES ${CORE_HEADERS})
INSTALL_FILES(/${SIMBODY_INCLUDE_INSTALL_DIR}/simbody/         FILES ${TOP_HEADERS})
INSTALL_FILES(/${SIMBODY_INCLUDE_INSTALL_DIR}/simbody/internal FILES ${INTERNAL_HEADERS})

# Notice that the (deprecated) INSTALL_FILES(dir ...) command works differently than the
# newer (recommended) INSTALL(dir FILES ...) command
# when the target directory starts with a slash
FILE(GLOB SIMBODY_DOCS doc/*.pdf doc/*.txt doc/*.md)
INSTALL(FILES ${SIMBODY_DOCS} DESTINATION ${CMAKE_INSTALL_DOCDIR})

# SIMBODY_VISUALIZER_INSTALL_DIR defined in root CMakeLists.txt
ADD_DEFINITIONS(-DSIMBODY_VISUALIZER_INSTALL_DIR="${SIMBODY_VISUALIZER_INSTALL_DIR}/")
ADD_DEFINITIONS(-DSIMBODY_VISUALIZER_REL_INSTALL_DIR="${SIMBODY_VISUALIZER_REL_INSTALL_DIR}/")

# These are at the end because we want them processed after
# all the various variables have been set above.

IF(BUILD_STATIC_LIBRARIES)
    ADD_SUBDIRECTORY( staticTarget )
ENDIF()
ADD_SUBDIRECTORY( sharedTarget )

IF (BUILD_VISUALIZER)
  # NOTE: If you change GUI_NAME, you must also change the value of GuiAppName
  # in VisualizerProtocol.cpp
  # TODO even that is not sufficient.
  SET(GUI_NAME "simbody-visualizer")
  ADD_SUBDIRECTORY(Visualizer/simbody-visualizer)
ENDIF (BUILD_VISUALIZER)

IF( BUILD_TESTING )
  ADD_SUBDIRECTORY( tests )
ENDIF( BUILD_TESTING )

