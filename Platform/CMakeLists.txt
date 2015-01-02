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

SET(PLATFORM_ROOT ${CMAKE_SOURCE_DIR}/Platform/${CMAKE_HOST_SYSTEM_NAME})

if(WIN32 AND NOT CMAKE_HOST_SYSTEM_NAME MATCHES Windows)
 # only provided for native builds
 return()
endif ()

SET(LIB_ABI_DIR "${PLATFORM_ROOT}/lib_${PLATFORM_ABI}")
FILE(GLOB LIB_FILES 
    "${LIB_ABI_DIR}/*.a"
    "${LIB_ABI_DIR}/*.so"
    "${LIB_ABI_DIR}/*.so.*"
    "${LIB_ABI_DIR}/*.dylib"
    "${LIB_ABI_DIR}/*.lib"
    "${LIB_ABI_DIR}/*.dll")

# This is where we're going to put these binaries.
SET(COPIED_LIB_FILES)
FOREACH(LIBF ${LIB_FILES})
    GET_FILENAME_COMPONENT(LIBF_ROOT ${LIBF} NAME)
    SET(COPIED_LIB_FILES ${COPIED_LIB_FILES}
    "${BUILD_BINARY_DIR}/${CMAKE_CFG_INTDIR}/${LIBF_ROOT}")
ENDFOREACH()


ADD_CUSTOM_TARGET(PlatformFiles ALL DEPENDS ${COPIED_LIB_FILES})
SET_TARGET_PROPERTIES(PlatformFiles
    PROPERTIES PROJECT_LABEL "Code - Platform Files")

# During build, we want to copy all of the platform binaries
# into the same binary directory that we are using for all of
# the ones we build. That ensures that our build will always
# use the current versions of these.
#
# At installation, we copy both includes and binaries into
# the appropriate locations of the SimTK_INSTALL_DIR.

FOREACH(LIBF ${LIB_FILES})
    GET_FILENAME_COMPONENT(LIBF_ROOT ${LIBF} NAME)
    GET_FILENAME_COMPONENT(LIBF_SUFFIX ${LIBF} EXT)
    SET(COPIED_LIBF "${BUILD_BINARY_DIR}/${CMAKE_CFG_INTDIR}/${LIBF_ROOT}")
    FILE(TO_NATIVE_PATH "${LIBF}" LIBF_SRC)
    FILE(TO_NATIVE_PATH "${COPIED_LIBF}" LIBF_DEST)

    ADD_CUSTOM_COMMAND(OUTPUT "${COPIED_LIBF}"
        COMMAND ${CMAKE_COMMAND} -E copy "${LIBF_SRC}" "${LIBF_DEST}"
        DEPENDS "${LIBF}"
        COMMENT "Copy ${LIBF_SRC} -> ${LIBF_DEST}"
        VERBATIM)

    IF (LIBF_SUFFIX STREQUAL ".dll")
        INSTALL(FILES ${LIBF} 
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE 
                GROUP_READ GROUP_WRITE GROUP_EXECUTE 
                    WORLD_READ WORLD_EXECUTE
        DESTINATION ${CMAKE_INSTALL_BINDIR})
    ELSE()
        INSTALL(FILES ${LIBF} 
        PERMISSIONS OWNER_READ OWNER_WRITE
                GROUP_READ GROUP_WRITE
                    WORLD_READ
        DESTINATION ${CMAKE_INSTALL_LIBDIR})
    ENDIF()
ENDFOREACH()

# These are directories that all SimTK libraries should link with
# to find the appropriate Platform libraries they depend on.
SET(PLATFORM_LINK_LIBRARIES "${PLATFORM_ROOT}/lib_${PLATFORM_ABI}"
    PARENT_SCOPE)

# There is just one item in PLATFORM_INCLUDE_DIRECTORIES.
SET(PLATFORM_INCLUDE_DIRECTORIES "${PLATFORM_ROOT}/include")
# Copy to parent.
SET(PLATFORM_INCLUDE_DIRECTORIES "${PLATFORM_INCLUDE_DIRECTORIES}" 
    PARENT_SCOPE)

# This needs an outer loop if you add more include directories.
FILE(GLOB INCL_FILES "${PLATFORM_INCLUDE_DIRECTORIES}/*.h")
FOREACH(INCLF ${INCL_FILES})
    GET_FILENAME_COMPONENT(INCLF_ROOT ${INCLF} NAME)
    INSTALL(FILES ${INCLF} 
    PERMISSIONS OWNER_READ OWNER_WRITE
            GROUP_READ GROUP_WRITE
                WORLD_READ
        DESTINATION ${CMAKE_INSTALL_PREFIX}/include )
ENDFOREACH()

