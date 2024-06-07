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

project(SimTKcommon VERSION ${Simbody_VERSION})

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
set(SimTKCOMMON_SOURCE_SUBDIRS
    . Scalar SmallMatrix Mechanics BigMatrix
    Geometry Simulation Random Polynomial)

# Collect up information about the version of the SimTKcommon library
# we're building and make it available to the code so it can be built
# into the binaries. This also determines the versioned library names
# in which case all dependencies must use the same version.

set(SimTKCOMMON_COPYRIGHT_YEARS "2005-10")

# underbar separated list of dotted authors, no spaces or commas
set(SimTKCOMMON_AUTHORS "Michael.Sherman_Peter.Eastman")

# These are all the places to search for header files which are
# to be part of the API.
set(API_INCLUDE_DIRS) # start empty
set(SimTKCOMMON_BUILD_INCLUDE_DIRS) # start empty
foreach(subdir ${SimTKCOMMON_SOURCE_SUBDIRS})
    list(APPEND API_INCLUDE_DIRS
     ${PROJECT_SOURCE_DIR}/${subdir}/include
     ${PROJECT_SOURCE_DIR}/${subdir}/include/SimTKcommon
     ${PROJECT_SOURCE_DIR}/${subdir}/include/SimTKcommon/internal)
    # Referencing headers must always be done relative to this level.
    list(APPEND SimTKCOMMON_BUILD_INCLUDE_DIRS
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

foreach(subdir ${SimTKCOMMON_SOURCE_SUBDIRS})
    file(GLOB src_files  ${subdir}/src/*.cpp)
    file(GLOB incl_files ${subdir}/src/*.h)
    list(APPEND SOURCE_FILES         ${src_files})
    list(APPEND SOURCE_INCLUDE_FILES ${incl_files})
endforeach()

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

if( INSTALL_DOCS )
    file(GLOB SIMTKCOMMON_DOCS doc/*.pdf doc/*.txt doc/*.md)
    install(FILES ${SIMTKCOMMON_DOCS} DESTINATION ${CMAKE_INSTALL_DOCDIR})
endif()

# set locally-scoped BUILD_SHARED_LIBS to globally cache'd SIMBODY_BUILD_SHARED_LIBS which defaults to
# BUILD_SHARED_LIBS if defined, or ON otherwise
set(BUILD_SHARED_LIBS ${SIMBODY_BUILD_SHARED_LIBS})

add_library(SimTKcommon
    ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES}
    ${API_ABS_INCLUDE_FILES})

add_dependencies(SimTKcommon PlatformFiles)

target_link_libraries(SimTKcommon PUBLIC ${MATH_LIBS_TO_USE})

target_include_directories(SimTKcommon
    PUBLIC
        "$<INSTALL_INTERFACE:${SIMBODY_INCLUDE_INSTALL_DIR}>"
        ${SimTKCOMMON_BUILD_INCLUDE_DIRS}
    )

set_target_properties(SimTKcommon PROPERTIES
    PROJECT_LABEL "Code - SimTKcommon"
    SOVERSION "${SimTKcommon_VERSION_MAJOR}.${SimTKcommon_VERSION_MINOR}"
    OUTPUT_NAME ${SimTKCOMMON_LIBRARY_NAME}
    ARCHIVE_OUTPUT_NAME ${SimTKCOMMON_LIBRARY_NAME}${SIMBODY_STATIC_LIBRARIES_POSTFIX})

set(shared_lib "$<BOOL:${BUILD_SHARED_LIBS}>")
set(static_lib "$<NOT:${shared_lib}>")
target_compile_definitions(SimTKcommon
    PUBLIC
        "$<${static_lib}:SimTK_USE_STATIC_LIBRARIES>"
    PRIVATE
        "$<${shared_lib}:SimTK_SimTKCOMMON_BUILDING_SHARED_LIBRARY>"
        "$<${static_lib}:SimTK_SimTKCOMMON_BUILDING_STATIC_LIBRARY>"
        SimTK_SimTKCOMMON_LIBRARY_NAME=${SimTKCOMMON_LIBRARY_NAME}
        SimTK_SimTKCOMMON_MAJOR_VERSION=${SimTKcommon_VERSION_MAJOR}
        SimTK_SimTKCOMMON_MINOR_VERSION=${SimTKcommon_VERSION_MINOR}
        SimTK_SimTKCOMMON_PATCH_VERSION=${SimTKcommon_VERSION_PATCH}
        SimTK_SimTKCOMMON_COPYRIGHT_YEARS="${SimTKCOMMON_COPYRIGHT_YEARS}"
        SimTK_SimTKCOMMON_AUTHORS="${SimTKCOMMON_AUTHORS}")
# -DSimTKcommon_EXPORTS defined automatically when Windows DLL build is being done.

# Setting as PUBLIC means that the cxx_std_11 requirement will be propagated to dependents
# (maintaining intent of previous behavior, but more universally applied and effective)
# See https://cmake.org/cmake/help/latest/manual/cmake-compile-features.7.html#language-standard-flags
# for info re: ordering of the compiler flag(s) for language standards and CMAKE_<LANG>_FLAGS
target_compile_features(SimTKcommon PUBLIC cxx_std_11)

# Warning 4996 is not propagated to downstream targets
target_compile_options(SimTKcommon INTERFACE
    "$<$<AND:$<CXX_COMPILER_ID:MSVC>,$<CONFIG:Debug>>:/wd4996>")

# install library; on Windows .dll's go in bin directory
install(TARGETS SimTKcommon EXPORT SimbodyTargets
    PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_WRITE GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    ARCHIVE
        PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ GROUP_WRITE WORLD_READ
        DESTINATION ${CMAKE_INSTALL_LIBDIR})

if(NOT SIMBODY_BUILD_SHARED_LIBS AND SIMBODY_COVERAGE)
    add_coverage(SimTKcommon)
endif()

if( BUILD_TESTING )
    add_subdirectory( tests )
endif()

