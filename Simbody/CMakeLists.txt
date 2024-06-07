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

project(SimTKsimbody VERSION ${Simbody_VERSION})

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
# The BUILD_VISUALIZER variable determines whether we actually build the
# simbody-visualizer, but we still have to build the library-side support for
# the Visualizer so that Simbody tests won't fail to build.
set(SIMBODY_SOURCE_SUBDIRS . Visualizer)

set(SIMBODY_COPYRIGHT_YEARS "2005-20")

# underbar separated list of dotted authors, no spaces or commas
set(SIMBODY_AUTHORS "Michael.Sherman_Peter.Eastman")

# These are all the places to search for header files which are
# to be part of the API.
set(API_INCLUDE_DIRS) # start empty
set(SimTKSIMBODY_BUILD_INCLUDE_DIRS) # start empty
foreach(subdir ${SIMBODY_SOURCE_SUBDIRS})
    list(APPEND API_INCLUDE_DIRS
     ${PROJECT_SOURCE_DIR}/${subdir}/include
     ${PROJECT_SOURCE_DIR}/${subdir}/include/simbody
     ${PROJECT_SOURCE_DIR}/${subdir}/include/simbody/internal)

    # Referencing headers must always be done relative to this level.
    list(APPEND SimTKSIMBODY_BUILD_INCLUDE_DIRS
        "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/${subdir}/include>")
endforeach()

# We'll need both *relative* path names, starting with their API_INCLUDE_DIRS,
# and absolute pathnames.
set(API_REL_INCLUDE_FILES)   # start these out empty
set(API_ABS_INCLUDE_FILES)

foreach(dir ${API_INCLUDE_DIRS})
    file(GLOB fullpaths ${dir}/*.h)    # returns full pathnames
    list(APPEND API_ABS_INCLUDE_FILES ${fullpaths})

    foreach(pathname ${fullpaths})
        get_filename_component(filename ${pathname} NAME)
        list(APPEND API_REL_INCLUDE_FILES ${dir}/${filename})
    endforeach()
endforeach()

# collect up source files
set(SOURCE_FILES) # empty
set(SOURCE_INCLUDE_FILES)

foreach(subdir ${SIMBODY_SOURCE_SUBDIRS})
    file(GLOB src_files  ${subdir}/src/*.cpp ${subdir}/src/*/*.cpp)
    file(GLOB incl_files ${subdir}/src/*.h)
    list(APPEND SOURCE_FILES         ${src_files})
    list(APPEND SOURCE_INCLUDE_FILES ${incl_files})

    #include_directories(${PROJECT_SOURCE_DIR}/${subdir}/include)
endforeach()

#include_directories(${PROJECT_SOURCE_DIR}/src)

# libraries are installed from their subdirectories; headers here

# install headers
file(GLOB CORE_HEADERS     include/*.h                  */include/*.h)
file(GLOB TOP_HEADERS      include/simbody/*.h          */include/simbody/*.h)
file(GLOB INTERNAL_HEADERS include/simbody/internal/*.h */include/simbody/internal/*.h)
install(FILES ${CORE_HEADERS} DESTINATION ${SIMBODY_INCLUDE_INSTALL_DIR})
install(FILES ${TOP_HEADERS} DESTINATION ${SIMBODY_INCLUDE_INSTALL_DIR}/simbody)
install(FILES ${INTERNAL_HEADERS} DESTINATION ${SIMBODY_INCLUDE_INSTALL_DIR}/simbody/internal)

if( INSTALL_DOCS )
    file(GLOB SIMBODY_DOCS doc/*.pdf doc/*.txt doc/*.md)
    install(FILES ${SIMBODY_DOCS} DESTINATION ${CMAKE_INSTALL_DOCDIR})
endif()

# set locally-scoped BUILD_SHARED_LIBS to globally cache'd SIMBODY_BUILD_SHARED_LIBS which defaults to
# BUILD_SHARED_LIBS if defined, or ON otherwise
set(BUILD_SHARED_LIBS ${SIMBODY_BUILD_SHARED_LIBS})

add_library(SimTKsimbody
    ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES}
    ${API_ABS_INCLUDE_FILES})

target_link_libraries(SimTKsimbody PUBLIC
    SimTKmath SimTKcommon ${MATH_LIBS_TO_USE})

target_include_directories(SimTKsimbody
    PUBLIC
        "$<INSTALL_INTERFACE:${SIMBODY_INCLUDE_INSTALL_DIR}>"
        ${SimTKSIMBODY_BUILD_INCLUDE_DIRS}
    )

set_target_properties(SimTKsimbody PROPERTIES
    PROJECT_LABEL "Code - SimTKsimbody"
    SOVERSION ${SIMBODY_SONAME_VERSION}
    OUTPUT_NAME ${SimTKSIMBODY_LIBRARY_NAME}
    ARCHIVE_OUTPUT_NAME ${SimTKSIMBODY_LIBRARY_NAME}${SIMBODY_STATIC_LIBRARIES_POSTFIX})


# Create a relative path from lib dir to dir containing simbody-visualizer;
# this is used in VisualizerProtocol.cpp (though, only useful on non-Windows).
file(RELATIVE_PATH lib_dir_to_install_dir
    "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" "${CMAKE_INSTALL_PREFIX}")
set(lib_dir_to_viz_dir
    "${lib_dir_to_install_dir}${SIMBODY_VISUALIZER_REL_INSTALL_DIR}")

set(shared_lib "$<BOOL:${BUILD_SHARED_LIBS}>")
set(static_lib "$<NOT:${shared_lib}>")
target_compile_definitions(SimTKsimbody
    PUBLIC
        "$<${static_lib}:SimTK_USE_STATIC_LIBRARIES>"
    PRIVATE
        "$<${shared_lib}:SimTK_SIMBODY_BUILDING_SHARED_LIBRARY>"
        "$<${static_lib}:SimTK_SIMBODY_BUILDING_STATIC_LIBRARY>"
        SimTK_SIMBODY_LIBRARY_NAME=${SimTKSIMBODY_LIBRARY_NAME}
        SimTK_SIMBODY_MAJOR_VERSION=${SimTKsimbody_VERSION_MAJOR}
        SimTK_SIMBODY_MINOR_VERSION=${SimTKsimbody_VERSION_MINOR}
        SimTK_SIMBODY_PATCH_VERSION=${SimTKsimbody_VERSION_PATCH}
        SIMBODY_VISUALIZER_INSTALL_DIR="${SIMBODY_VISUALIZER_INSTALL_DIR}/"
        SIMBODY_VISUALIZER_REL_INSTALL_DIR="${SIMBODY_VISUALIZER_REL_INSTALL_DIR}/"
        SIMBODY_PATH_FROM_LIBDIR_TO_VIZ_DIR="${lib_dir_to_viz_dir}"
        SimTK_SIMBODY_COPYRIGHT_YEARS="${SIMBODY_COPYRIGHT_YEARS}"
        SimTK_SIMBODY_AUTHORS="${SIMBODY_AUTHORS}"
)

# install library; on Windows the .dll goes in the bin directory.
install(TARGETS SimTKsimbody EXPORT SimbodyTargets
    PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_WRITE GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    ARCHIVE
        PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ GROUP_WRITE WORLD_READ
        DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )

if(NOT SIMBODY_BUILD_SHARED_LIBS AND SIMBODY_COVERAGE)
    add_coverage(SimTKsimbody)
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

