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

project(SimTKsimbody)

# SimTKsimbody depends on PlatformFiles, SimTKcommon, and SimTKmath only.
include_directories(${PLATFORM_INCLUDE_DIRECTORIES}
                    ${SimTKCOMMON_INCLUDE_DIRECTORIES}
                    ${SimTKMATH_INCLUDE_DIRECTORIES})

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
# The BUILD_VISUALIZER variable determines whether we actually build the
# simbody-visualizer, but we still have to build the library-side support for
# the Visualizer so that Simbody tests won't fail to build.
set(SIMBODY_SOURCE_SUBDIRS . Visualizer)

add_definitions(-DSimTK_SIMBODY_LIBRARY_NAME=${SimTKSIMBODY_LIBRARY_NAME}
                -DSimTK_SIMBODY_MAJOR_VERSION=${SIMBODY_MAJOR_VERSION}
                -DSimTK_SIMBODY_MINOR_VERSION=${SIMBODY_MINOR_VERSION}
                -DSimTK_SIMBODY_PATCH_VERSION=${SIMBODY_PATCH_VERSION})


set(SHARED_TARGET ${SimTKSIMBODY_LIBRARY_NAME})
set(STATIC_TARGET ${SimTKSIMBODY_LIBRARY_NAME}_static)
set(SHARED_TARGET_VN ${SimTKSIMBODY_LIBRARY_NAME}${VN})
set(STATIC_TARGET_VN ${SimTKSIMBODY_LIBRARY_NAME}${VN}_static)

## Test against the unversioned libraries if they are being built;
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
set(SimTKSIMBODY_INCLUDE_DIRS) # start empty
foreach(subdir ${SIMBODY_SOURCE_SUBDIRS})
    # append
    set(API_INCLUDE_DIRS ${API_INCLUDE_DIRS}
     ${PROJECT_SOURCE_DIR}/${subdir}/include
     ${PROJECT_SOURCE_DIR}/${subdir}/include/simbody
     ${PROJECT_SOURCE_DIR}/${subdir}/include/simbody/internal)

    # Referencing headers must always be done relative to this level.
    set(SimTKSIMBODY_INCLUDE_DIRS ${SimTKSIMBODY_INCLUDE_DIRS}
     ${PROJECT_SOURCE_DIR}/${subdir}/include)
endforeach(subdir)

# Include the Simbody API include directories now so that Simbody code
# can use them.
include_directories(${SimTKSIMBODY_INCLUDE_DIRS})

# And pass API include directories up to the parent so subsequent libraries
# can find the headers too.
set(SimTKSIMBODY_INCLUDE_DIRECTORIES ${SimTKSIMBODY_INCLUDE_DIRS}
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

foreach(subdir ${SIMBODY_SOURCE_SUBDIRS})
    file(GLOB src_files  ${subdir}/src/*.cpp ${subdir}/src/*/*.cpp)
    file(GLOB incl_files ${subdir}/src/*.h)
    set(SOURCE_FILES         ${SOURCE_FILES}         ${src_files})   #append
    set(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})

    #include_directories(${PROJECT_SOURCE_DIR}/${subdir}/include)
endforeach(subdir)

#include_directories(${PROJECT_SOURCE_DIR}/src)

# libraries are installed from their subdirectories; headers here

# install headers
file(GLOB CORE_HEADERS     include/*.h                  */include/*.h)
file(GLOB TOP_HEADERS      include/simbody/*.h          */include/simbody/*.h)
file(GLOB INTERNAL_HEADERS include/simbody/internal/*.h */include/simbody/internal/*.h)
install(FILES ${CORE_HEADERS} DESTINATION ${SIMBODY_INCLUDE_INSTALL_DIR})
install(FILES ${TOP_HEADERS} DESTINATION ${SIMBODY_INCLUDE_INSTALL_DIR}/simbody)
install(FILES ${INTERNAL_HEADERS} DESTINATION ${SIMBODY_INCLUDE_INSTALL_DIR}/simbody/internal)

file(GLOB SIMBODY_DOCS doc/*.pdf doc/*.txt doc/*.md)
install(FILES ${SIMBODY_DOCS} DESTINATION ${CMAKE_INSTALL_DOCDIR})

# SIMBODY_VISUALIZER_INSTALL_DIR defined in root CMakeLists.txt
add_definitions(-DSIMBODY_VISUALIZER_INSTALL_DIR="${SIMBODY_VISUALIZER_INSTALL_DIR}/")
add_definitions(-DSIMBODY_VISUALIZER_REL_INSTALL_DIR="${SIMBODY_VISUALIZER_REL_INSTALL_DIR}/")
# Create a relative path from lib dir to dir containing simbody-visualizer;
# this is used in VisualizerProtocol.cpp (though, only useful on non-Windows).
file(RELATIVE_PATH lib_dir_to_install_dir
    "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" "${CMAKE_INSTALL_PREFIX}")
set(lib_dir_to_viz_dir
    "${lib_dir_to_install_dir}${SIMBODY_VISUALIZER_REL_INSTALL_DIR}")
add_definitions(
    -DSIMBODY_PATH_FROM_LIBDIR_TO_VIZ_DIR="${lib_dir_to_viz_dir}")

# These are at the end because we want them processed after
# all the various variables have been set above.

if(BUILD_STATIC_LIBRARIES)
    add_subdirectory( staticTarget )
endif()
if(BUILD_DYNAMIC_LIBRARIES)
    add_subdirectory( sharedTarget )
endif()

if(MINGW)
    if(${PLATFORM_ABI} MATCHES "x86")
        message(STATUS "Visualization is disabled because the shipped version of glut32.dll uses a")
        message(STATUS "different convention call than MinGW")
        message(STATUS "Use a 64 bit version of MinGW to build visualizer, or provide your own implementation of glut")
        set(BUILD_VISUALIZER FALSE)
    endif()
endif()
if(BUILD_VISUALIZER)
    # NOTE: If you change GUI_NAME, you must also change the value of GuiAppName
    # in VisualizerProtocol.cpp
    # TODO even that is not sufficient.
    set(GUI_NAME "simbody-visualizer")
    set(GUI_NAME ${GUI_NAME} PARENT_SCOPE)
    add_subdirectory(Visualizer/simbody-visualizer)
endif()

if( BUILD_TESTING )
    add_subdirectory( tests )
endif()

