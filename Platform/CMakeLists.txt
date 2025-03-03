#-----------------------------------------------------------
# Simbody/Platform
#
# Simbody has dependencies on several low-level packages
# for which the details differ on different platforms.
# Dependencies are:
#    - Blas, Lapack
# We provide our own high-performance Lapack library
# for Windows; Linux and Macs are expected to have one.
# We used to supply our own Pthreads on Windows, but this
# was replaced with C++11 threading in 2018.
#-----------------------------------------------------------

# Nothing here is needed done if Unix, WINDOWS_USE_EXTERNAL_LIBS or (non-native) cross-compiling for Windows
if(NOT WIN32 OR WINDOWS_USE_EXTERNAL_LIBS OR (WIN32 AND NOT CMAKE_HOST_SYSTEM_NAME MATCHES Windows))
    # Users/system must provide own BLAS/LAPACK
    return()
endif()

set(PLATFORM_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/${CMAKE_HOST_SYSTEM_NAME})

# Native cross-compiling (e.g. compiling on Windows using non-MSVC compiler, aka MinGW)
if(MINGW)
    # Caution: No cross-compiling is taken into account:
    # If we use a 32 bit compiler, we will get a 32 executable.
    # The same for 64 bit platform

    # 1 Check the thread model
    include(CheckCXXSourceCompiles)
    CHECK_CXX_SOURCE_COMPILES("#include <mutex>\nint main(){std::mutex m;return 0;}" GPP_SUPPORTS_STD_MUTEX)
    if(NOT GPP_SUPPORTS_STD_MUTEX)
        message(STATUS
            "Your MinGW version does not support std::mutex."
            "This is probably because of the win32 thread model used."
            "Change to a posix thread model."
            "Also note that exception mechanism is important before downloading another MinGW version"
            "    For a 32 bit instatallation, use a Posix thread with DWARF (DW2) exception mechanism"
            "    For a 64 bit instatallation, use a Posix thread with SJLJ exception mechanism")
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

# Create "fake" IMPORTED targets to represent the pre-built platform libraries.
# When CMake sees that a target links to, e.g., the blas target, CMake will use
# the appropriate library paths below.
# Note: we don't need to do this for the other platform libraries
# (gcc_s_sjlj-1, gfortran, quadmath, ...) since Simbody's libraries do not link
# to them when compiling; we only need to do this for import libraries (.lib).
# All the platform libraries (including quadmath, etc.) still must be on the
# PATH when running executables using Simbody.
if(WIN32 AND NOT BUILD_USING_OTHER_LAPACK)
    # declare transitive dependencies for blas/lapack (see simbody/simbody#771)
    #
    # the transitive dependencies don't, themselves, have IMPLIBs, because
    # they have already been linked by blas/lapack. However, cmake can become
    # confused when it sees a library with no IMPLIB, so we include blas's
    # here to satisfy cmake
    #
    # (if you link to blas or lapack, you're going to link to the blas IMPLIB
    # anyway, so it makes no difference)
    add_library(libgcc_s_sjlj-1 SHARED IMPORTED)
    set_target_properties(libgcc_s_sjlj-1 PROPERTIES
        IMPORTED_IMPLIB "${LIB_ABI_DIR}/libblas.lib"
        IMPORTED_LOCATION "${LIB_ABI_DIR}/libgcc_s_sjlj-1.dll"
        )
    add_library(libgfortran-3 SHARED IMPORTED)
    set_target_properties(libgfortran-3 PROPERTIES
        IMPORTED_IMPLIB "${LIB_ABI_DIR}/libblas.lib"
        IMPORTED_LOCATION "${LIB_ABI_DIR}/libgfortran-3.dll"
        )
    add_library(libquadmath-0 SHARED IMPORTED)
    set_target_properties(libquadmath-0 PROPERTIES
        IMPORTED_IMPLIB "${LIB_ABI_DIR}/libblas.lib"
        IMPORTED_LOCATION "${LIB_ABI_DIR}/libquadmath-0.dll"
        )

    # Without GLOBAL, this target is only available in this dir. and below,
    # but we want to use these targets everywhere in this project.
    add_library(blas SHARED IMPORTED GLOBAL)
    set_target_properties(blas PROPERTIES
        IMPORTED_IMPLIB "${LIB_ABI_DIR}/libblas.lib"
        IMPORTED_LOCATION "${LIB_ABI_DIR}/libblas.dll"
        )
    # link to transitive dependencies
    target_link_libraries(blas INTERFACE libgcc_s_sjlj-1 libgfortran-3 libquadmath-0)

    # CMake won't install import libraries, but downstream libraries that depend on our
    # vendored blas/lapack need the import libs to link
    install(FILES $<TARGET_PROPERTY:blas,IMPORTED_IMPLIB>
        PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ GROUP_WRITE WORLD_READ
        DESTINATION ${CMAKE_INSTALL_LIBDIR})

    add_library(lapack SHARED IMPORTED GLOBAL)
    set_target_properties(lapack PROPERTIES
        IMPORTED_IMPLIB "${LIB_ABI_DIR}/liblapack.lib"
        IMPORTED_LOCATION "${LIB_ABI_DIR}/liblapack.dll"
        # lapack depends on blas:
        INTERFACE_LINK_LIBRARIES blas
        )
    install(FILES $<TARGET_PROPERTY:lapack,IMPORTED_IMPLIB>
        PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ GROUP_WRITE WORLD_READ
        DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
