#-----------------------------------------------------------
# Simbody/Platform
#
# Simbody has dependencies on several low-level packages
# for which the details differ on different platforms.
# Dependencies are:
#    - pthreads
#    - Blas, Lapack
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

if(MINGW)
    # Caution: No cross-compiling is taken into account:
    # If we use a 32 bit compiler, we will get a 32 executable.
    # The same for 64 bit platform
    #
    # 1 Check the thread model
    include(CheckCXXSourceCompiles)
    CHECK_CXX_SOURCE_COMPILES("#include <mutex>\nint main(){std::mutex m;return 0;}" GPP_SUPPORTS_STD_MUTEX)
    if(NOT GPP_SUPPORTS_STD_MUTEX)
        message(STATUS "Your MinGW version does not support std::mutex.")
        message(STATUS "This is probably because of the win32 thread model used.")
        message(STATUS "Change to a posix thread model.")
        message(STATUS "Also note that exception mechanism is important before downloading another MinGW version")
        message(STATUS "    For a 32 bit instatallation, use a Posix thread with DWARF (DW2) exception mechanism")
        message(STATUS "    For a 64 bit instatallation, use a Posix thread with SJLJ exception mechanism")
        message(FATAL_ERROR "Your MinGW version does not support std::mutex. -> Change to a posix thread model")
    endif()
    # 2 Check the exception mechanism
    get_filename_component(MINGW_BIN_PATH ${CMAKE_C_COMPILER} PATH)
    message(STATUS "MINGW BIN PATH = ${MINGW_BIN_PATH}")
    set(MINGW_LIB_DIRECTORY "")
    if(${PLATFORM_ABI} MATCHES "x64")
        set(MINGW_LIB_DIRECTORY "x86_64-w64-mingw32")
    else()
        set(MINGW_LIB_DIRECTORY "i686-w64-mingw32")
    endif()
    if(EXISTS "${MINGW_BIN_PATH}/../${MINGW_LIB_DIRECTORY}/")
        file(GLOB MINGW_DLLS ${MINGW_BIN_PATH}/../${MINGW_LIB_DIRECTORY}/lib/*.dll)
    else()
        # Maybe Dlls are present in the mingw bin directory
        file(GLOB MINGW_DLLS ${MINGW_BIN_PATH}/*.dll)
    endif()
    if(NOT MINGW_DLLS)
        message(FATAL_ERROR "No MinGW DLL was found. CMake could not find the MinGW dll of your installation")
    endif()
    set(MINGW_GCC_EXCEPTION "")
    foreach(f ${MINGW_DLLS})
        if(${f} MATCHES "libgcc")
            message(STATUS "MINGW libgcc library = ${f}")
            if(${f} MATCHES "dw2")
                set(MINGW_GCC_EXCEPTION "dw2")
            elseif(${f} MATCHES "sjlj")
                set(MINGW_GCC_EXCEPTION "sjlj")
            elseif(${f} MATCHES "seh")
                set(MINGW_GCC_EXCEPTION "seh")
            else()
                message(FATAL_ERROR "MinGW exception was not identified amongst dw2, seh, sjlj")
            endif()
        endif()
    endforeach()
    if(${MINGW_GCC_EXCEPTION} STREQUAL "" )
        message(STATUS "MinGW exception mechanism could not be found")
        message(FATAL_ERROR "Please update your version of MinGW")
    endif()
    if(NOT BUILD_USING_OTHER_LAPACK)
        # We have to ensure that the exception mechanism of MinGW
        # is the same as the one used with the blas and lapack libraries provided
        if(${PLATFORM_ABI} MATCHES "x64")
            if(NOT ${MINGW_GCC_EXCEPTION} STREQUAL "sjlj" )
                message(STATUS "MinGW exception mechanism does not match the one used to generate Blas and Lapack")
                message(FATAL_ERROR "Please provide your own Blas and Lapack libraries or use a SJLJ version of MinGW")
            endif()
        else()
            if(NOT ${MINGW_GCC_EXCEPTION} STREQUAL "dw2" )
                message(STATUS "MinGW exception mechanism does not match the one used to generate Blas and Lapack")
                message(FATAL_ERROR "Please provide your own Blas and Lapack libraries or use a DWARF (DW2) version of MinGW")
            endif()
        endif()
    endif()
endif()

set(LIB_ABI_DIR "${PLATFORM_ROOT}/lib_${PLATFORM_ABI}")
file(GLOB LIB_FILES
    "${LIB_ABI_DIR}/*.a"
    "${LIB_ABI_DIR}/*.so"
    "${LIB_ABI_DIR}/*.so.*"
    "${LIB_ABI_DIR}/*.dylib"
    "${LIB_ABI_DIR}/*.lib"
    "${LIB_ABI_DIR}/*.dll")

if(BUILD_USING_OTHER_LAPACK)
    if(MINGW)
        # If we use MinGW, we bypass the classical *dll and use
        # the MinGW ones.
        # TODO: We should select only the ones needed...
        set(LIB_FILES ${MINGW_DLLS} ${LAPACK_BEING_USED})
    endif()
endif()

# message(STATUS "List of dependencies : ${LIB_FILES}" )

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
    IF(MINGW)
        set(LIBF_SRC ${LIBF})
        set(LIBF_DEST ${COPIED_LIBF})
    ELSE()
        file(TO_NATIVE_PATH "${LIBF}" LIBF_SRC)
        file(TO_NATIVE_PATH "${COPIED_LIBF}" LIBF_DEST)
    ENDIF()

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
if(MSVC)
    set(PLATFORM_INCLUDE_DIRECTORIES "${PLATFORM_ROOT}/include")
endif()
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

