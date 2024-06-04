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

project(SimTKmath VERSION ${Simbody_VERSION})

# code is in "./src" and in each of these subdirectories/src
set(SIMMATH_SOURCE_SUBDIRS LinearAlgebra Integrators Optimizers Geometry)

# Collect up information about the version of the SimTKmath library
# we're building and make it available to the code so it can be built
# into the binaries. This also determines the versioned library names
# in which case all dependencies must use the same version.

set(SIMMATH_COPYRIGHT_YEARS "2005-10")

# underbar separated list of dotted authors, no spaces or commas
set(SIMMATH_AUTHORS         "Jack.Middleton_Michael.Sherman")

# These are all the places to search for header files which are
# to be part of the API.
set(API_INCLUDE_DIRS) # start empty
set(SimTKMATH_INCLUDE_DIRS) # start empty
foreach(subdir . ${SIMMATH_SOURCE_SUBDIRS})
    list(APPEND API_INCLUDE_DIRS
     ${PROJECT_SOURCE_DIR}/${subdir}/include
     ${PROJECT_SOURCE_DIR}/${subdir}/include/simmath
     ${PROJECT_SOURCE_DIR}/${subdir}/include/simmath/internal)

    # Referencing headers must always be done relative to this level.
    list(APPEND SimTKMATH_INCLUDE_DIRS
        "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/${subdir}/include>")
endforeach()

# We'll need both *relative* path names, starting with
# their API_INCLUDE_DIRS, and absolute pathnames.
set(API_REL_INCLUDE_FILES)   # start these out empty
set(API_ABS_INCLUDE_FILES)

foreach(dir ${API_INCLUDE_DIRS})
    file(GLOB fullpaths ${dir}/*.h)    # returns full pathnames
    list(APPEND API_ABS_INCLUDE_FILES ${fullpaths})

    foreach(pathname ${fullpaths})
        get_filename_component(filename ${pathname} NAME)
        list(APPEND API_REL_INCLUDE_FILES
        ${dir}/${filename})
    endforeach()
endforeach()

# collect up source files
set(SOURCE_FILES) # empty
set(SOURCE_INCLUDE_FILES) # for dependency tracking only, I believe (sherm)
set(SOURCE_INCLUDE_DIRS) # in case some low-level source needs a -I include dir
set(SOURCE_GROUPS)
set(SOURCE_GROUP_FILES)

# first process all the source subdirectories
foreach(subdir ${SIMMATH_SOURCE_SUBDIRS})
    add_subdirectory(${subdir})
endforeach()

if(SOURCE_GROUPS)
    list(LENGTH SOURCE_GROUPS NGROUPS)
    math(EXPR lastgrpnum ${NGROUPS}-1)
    foreach( grpnum RANGE 0 ${lastgrpnum} )
        list(GET SOURCE_GROUPS ${grpnum} grp)
        list(GET SOURCE_GROUP_FILES ${grpnum} grpfile)
        source_group("${grp}" FILES "${grpfile}")
    endforeach()
endif()

# then process ./src
file(GLOB src_files  ./src/*.cpp)
file(GLOB incl_files ./src/*.h)
list(APPEND SOURCE_FILES         ${src_files})
list(APPEND SOURCE_INCLUDE_FILES ${incl_files})

#
# Installation
#

# libraries are installed from their subdirectories; headers here

# install headers
file(GLOB CORE_HEADERS     include/*.h                  */include/*.h)
file(GLOB TOP_HEADERS      include/simmath/*.h          */include/simmath/*.h)
file(GLOB INTERNAL_HEADERS include/simmath/internal/*.h */include/simmath/internal/*.h)

install(FILES ${CORE_HEADERS} DESTINATION ${SIMBODY_INCLUDE_INSTALL_DIR})
install(FILES ${TOP_HEADERS} DESTINATION ${SIMBODY_INCLUDE_INSTALL_DIR}/simmath)
install(FILES ${INTERNAL_HEADERS} DESTINATION ${SIMBODY_INCLUDE_INSTALL_DIR}/simmath/internal)

file(GLOB SIMMATH_DOCS doc/*.pdf doc/*.txt doc/*.md)
install(FILES ${SIMMATH_DOCS} DESTINATION ${CMAKE_INSTALL_DOCDIR})

#
# Set up file groups for better browsing in Visual Studio.
#

# Override Api headers to put them in the Api group.
source_group("Header Files\\Api" FILES ${API_ABS_INCLUDE_FILES})

# Now put things in their specified groups.
if(SOURCE_GROUPS)
    list(LENGTH SOURCE_GROUPS NGROUPS)
    math(EXPR lastgrpnum ${NGROUPS}-1)
    foreach( grpnum RANGE 0 ${lastgrpnum} )
        list(GET SOURCE_GROUPS ${grpnum} grp)
        list(GET SOURCE_GROUP_FILES ${grpnum} grpfile)
        source_group("${grp}" FILES "${grpfile}")
    endforeach()
endif()

# These are defaults for files that don't otherwise get mentioned.
source_group("Source Files\\Others" FILES ${SOURCE_FILES})
source_group("Header Files\\Others" FILES ${SOURCE_INCLUDE_FILES})

# set locally-scoped BUILD_SHARED_LIBS to globally cache'd SIMBODY_BUILD_SHARED_LIBS which defaults to
# BUILD_SHARED_LIBS if defined, or ON otherwise
set(BUILD_SHARED_LIBS ${SIMBODY_BUILD_SHARED_LIBS})

add_library(SimTKmath
    ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES}
    ${API_ABS_INCLUDE_FILES})

# SimTKmath depends on PlatformFiles and SimTKcommon only.
target_link_libraries(SimTKmath PUBLIC
    SimTKcommon ${MATH_LIBS_TO_USE}) # TODO: Do we really need to link to MATH_LIBS_TO_USE?

target_include_directories(SimTKmath
    PUBLIC
        "$<INSTALL_INTERFACE:${SIMBODY_INCLUDE_INSTALL_DIR}>"
        ${SimTKMATH_INCLUDE_DIRS}
    PRIVATE
        ${PLATFORM_INCLUDE_DIRECTORIES}
        ${SOURCE_INCLUDE_DIRS})

set_target_properties(SimTKmath PROPERTIES
    PROJECT_LABEL "Code - SimTKmath"
    SOVERSION "${SimTKmath_VERSION_MAJOR}.${SimTKmath_VERSION_MINOR}"
    OUTPUT_NAME ${SimTKMATH_LIBRARY_NAME}
    ARCHIVE_OUTPUT_NAME ${SimTKMATH_LIBRARY_NAME}${SIMBODY_STATIC_LIBRARIES_POSTFIX})

set(shared_lib "$<BOOL:${BUILD_SHARED_LIBS}>")
set(static_lib "$<NOT:${shared_lib}>")
target_compile_definitions(SimTKmath 
    PUBLIC
        "$<${static_lib}:SimTK_USE_STATIC_LIBRARIES>"
    PRIVATE
        "$<${shared_lib}:SimTK_SIMMATH_BUILDING_SHARED_LIBRARY>"
        "$<${static_lib}:SimTK_SIMMATH_BUILDING_STATIC_LIBRARY>"
        SimTK_SIMMATH_LIBRARY_NAME=${SimTKMATH_LIBRARY_NAME}
        SimTK_SIMMATH_MAJOR_VERSION=${SimTKmath_VERSION_MAJOR}
        SimTK_SIMMATH_MINOR_VERSION=${SimTKmath_VERSION_MINOR}
        SimTK_SIMMATH_PATCH_VERSION=${SimTKmath_VERSION_PATCH}
        SimTK_SIMMATH_COPYRIGHT_YEARS="${SIMMATH_COPYRIGHT_YEARS}"
        SimTK_SIMMATH_AUTHORS="${SIMMATH_AUTHORS}")

# install library; on Windows both .lib and .dll go in the lib directory.
install(TARGETS SimTKmath EXPORT SimbodyTargets
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
    add_coverage(SimTKmath)
endif()

if( BUILD_TESTING )
    add_subdirectory( tests )
endif()
