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

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
SET(SIMBODY_SOURCE_SUBDIRS .)
IF (BUILD_VISUALIZER)
    SET(SIMBODY_SOURCE_SUBDIRS ${SIMBODY_SOURCE_SUBDIRS} Visualizer)
ENDIF (BUILD_VISUALIZER)

SET(SIMBODY_LIBRARY_NAME ${NS}SimTKsimbody CACHE STRING
"Base name of the library being built; can't be changed here; see BUILD_USING_NAMESPACE variable."
FORCE)

ADD_DEFINITIONS(-DSimTK_SIMBODY_LIBRARY_NAME=${SIMBODY_LIBRARY_NAME})

SET(SHARED_TARGET ${SIMBODY_LIBRARY_NAME})
SET(STATIC_TARGET ${SIMBODY_LIBRARY_NAME}_static)
SET(SHARED_TARGET_VN ${SIMBODY_LIBRARY_NAME}${VN})
SET(STATIC_TARGET_VN ${SIMBODY_LIBRARY_NAME}${VN}_static)

# On Unix or Cygwin we have to add the suffix manually
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

SET(SimTKCOMMON_LIBRARY_NAME ${NS}SimTKcommon)
SET(SimTKCOMMON_SHARED_LIBRARY ${SimTKCOMMON_LIBRARY_NAME})
SET(SimTKCOMMON_STATIC_LIBRARY ${SimTKCOMMON_LIBRARY_NAME}_static)

SET(SimTKCOMMON_LIBRARY_NAME_VN ${NS}SimTKcommon${VN})
SET(SimTKCOMMON_SHARED_LIBRARY_VN ${SimTKCOMMON_LIBRARY_NAME_VN})
SET(SimTKCOMMON_STATIC_LIBRARY_VN ${SimTKCOMMON_LIBRARY_NAME_VN}_static)

SET(SimTKMATH_LIBRARY_NAME ${NS}SimTKmath)
SET(SimTKMATH_SHARED_LIBRARY ${SimTKMATH_LIBRARY_NAME})
SET(SimTKMATH_STATIC_LIBRARY ${SimTKMATH_LIBRARY_NAME}_static)

SET(SimTKMATH_LIBRARY_NAME_VN ${NS}SimTKmath${VN})
SET(SimTKMATH_SHARED_LIBRARY_VN ${SimTKMATH_LIBRARY_NAME_VN})
SET(SimTKMATH_STATIC_LIBRARY_VN ${SimTKMATH_LIBRARY_NAME_VN}_static)

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

    ## Make sure we find these locally before looking in SimTK/include if
    ## Simbody was previously installed there.
    INCLUDE_DIRECTORIES(BEFORE ${PROJECT_SOURCE_DIR}/${subdir}/include)
ENDFOREACH(subdir)

INCLUDE_DIRECTORIES(BEFORE ${PROJECT_SOURCE_DIR}/src)

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_SOURCE_DIR}/SimTKmath/include)
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_SOURCE_DIR}/SimTKcommon/include)
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_SOURCE_DIR}/Prerequisites/include)

INCLUDE(ApiDoxygen.cmake)

# libraries and examples are installed from their subdirectories; headers here

# install headers
FILE(GLOB CORE_HEADERS     include/*.h                  */include/*.h)
FILE(GLOB TOP_HEADERS      include/simbody/*.h          */include/simbody/*.h)
FILE(GLOB INTERNAL_HEADERS include/simbody/internal/*.h */include/simbody/internal/*.h)
INSTALL_FILES(/include/                 FILES ${CORE_HEADERS})
INSTALL_FILES(/include/simbody/         FILES ${TOP_HEADERS})
INSTALL_FILES(/include/simbody/internal FILES ${INTERNAL_HEADERS})

# Notice that the (deprecated) INSTALL_FILES(dir ...) command works differently than the
# newer (recommended) INSTALL(dir FILES ...) command
# when the target directory starts with a slash
INSTALL( FILES ${PROJECT_SOURCE_DIR}/doc/SimbodyTheoryManual.pdf DESTINATION doc/simbody )

# These are at the end because we want them processed after
# all the various variables have been set above.

ADD_SUBDIRECTORY( staticTarget )
ADD_SUBDIRECTORY( sharedTarget )

IF (BUILD_VISUALIZER)
  SET(GUI_NAME "VisualizerGUI")
  ADD_SUBDIRECTORY(Visualizer/VisualizerGUI)
ENDIF (BUILD_VISUALIZER)

IF( BUILD_EXAMPLES )
  ADD_SUBDIRECTORY( examples )
ENDIF( BUILD_EXAMPLES )

IF( BUILD_TESTING )
  ADD_SUBDIRECTORY( tests )
ENDIF( BUILD_TESTING )

