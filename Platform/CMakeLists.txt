#-----------------------------------------------------------
# Simbody/Platform
#
# Simbody has dependencies on several low-level packages
# for which the details differ on different platforms.
# Dependencies are:
#    - pthreads
#    - Lapack
# Pthreads is standard on Linux platforms, and close enough
# on Macs. On Windows we have to supply our own.
# We provide our own high-performance Lapack library
# for Windows; Linux and Macs are expected to have one.
#-----------------------------------------------------------

set(PLATFORM_ROOT ${CMAKE_SOURCE_DIR}/Platform/${CMAKE_HOST_SYSTEM_NAME})

if(WIN32 AND NOT CMAKE_HOST_SYSTEM_NAME MATCHES Windows)
 # only provided for native builds
 return()
endif()

set(LIB_ABI_DIR "${PLATFORM_ROOT}/lib_${PLATFORM_ABI}")
file(GLOB LIB_FILES
    "${LIB_ABI_DIR}/*.a"
    "${LIB_ABI_DIR}/*.so"
    "${LIB_ABI_DIR}/*.so.*"
    "${LIB_ABI_DIR}/*.dylib"
    "${LIB_ABI_DIR}/*.lib"
    "${LIB_ABI_DIR}/*.dll")

# This is where we're going to put these binaries.
set(COPIED_LIB_FILES)
foreach(LIBF ${LIB_FILES})
    get_filename_component(LIBF_ROOT ${LIBF} NAME)
    set(COPIED_LIB_FILES ${COPIED_LIB_FILES}
    "${BUILD_BINARY_DIR}/${CMAKE_CFG_INTDIR}/${LIBF_ROOT}")
endforeach()


add_custom_target(PlatformFiles_subdir DEPENDS ${COPIED_LIB_FILES})
add_dependencies(PlatformFiles PlatformFiles_subdir)
set_target_properties(PlatformFiles
    PROPERTIES PROJECT_LABEL "Code - Platform Files")

# During build, we want to copy all of the platform binaries
# into the same binary directory that we are using for all of
# the ones we build. That ensures that our build will always
# use the current versions of these.
#
# At installation, we copy both includes and binaries into
# the appropriate locations of the SimTK_INSTALL_DIR.

foreach(LIBF ${LIB_FILES})
    get_filename_component(LIBF_ROOT ${LIBF} NAME)
    get_filename_component(LIBF_SUFFIX ${LIBF} EXT)
    set(COPIED_LIBF "${BUILD_BINARY_DIR}/${CMAKE_CFG_INTDIR}/${LIBF_ROOT}")
    file(TO_NATIVE_PATH "${LIBF}" LIBF_SRC)
    file(TO_NATIVE_PATH "${COPIED_LIBF}" LIBF_DEST)

    add_custom_command(OUTPUT "${COPIED_LIBF}"
        COMMAND ${CMAKE_COMMAND} -E copy "${LIBF_SRC}" "${LIBF_DEST}"
        DEPENDS "${LIBF}"
        COMMENT "Copy ${LIBF_SRC} -> ${LIBF_DEST}"
        VERBATIM)

    if(LIBF_SUFFIX STREQUAL ".dll")
        install(FILES ${LIBF}
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                GROUP_READ GROUP_WRITE GROUP_EXECUTE
                    WORLD_READ WORLD_EXECUTE
        DESTINATION ${CMAKE_INSTALL_BINDIR})
    else()
        install(FILES ${LIBF}
        PERMISSIONS OWNER_READ OWNER_WRITE
                GROUP_READ GROUP_WRITE
                    WORLD_READ
        DESTINATION ${CMAKE_INSTALL_LIBDIR})
    endif()
endforeach()

# These are directories that all SimTK libraries should link with
# to find the appropriate Platform libraries they depend on.
set(PLATFORM_LINK_LIBRARIES "${PLATFORM_ROOT}/lib_${PLATFORM_ABI}"
    PARENT_SCOPE)

# There is just one item in PLATFORM_INCLUDE_DIRECTORIES.
set(PLATFORM_INCLUDE_DIRECTORIES "${PLATFORM_ROOT}/include")
# Copy to parent.
set(PLATFORM_INCLUDE_DIRECTORIES "${PLATFORM_INCLUDE_DIRECTORIES}"
    PARENT_SCOPE)

# This needs an outer loop if you add more include directories.
file(GLOB INCL_FILES "${PLATFORM_INCLUDE_DIRECTORIES}/*.h")
foreach(INCLF ${INCL_FILES})
    get_filename_component(INCLF_ROOT ${INCLF} NAME)
    install(FILES ${INCLF}
    PERMISSIONS OWNER_READ OWNER_WRITE
            GROUP_READ GROUP_WRITE
                WORLD_READ
        DESTINATION include )
endforeach()

