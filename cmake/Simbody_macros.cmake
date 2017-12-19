# These CMake functions serve to reduce duplication across CMakeLists.txt files.

include(CMakeParseArguments)

# Copy DLL files from a dependency's installation into Simbody's
# build and install directories. This is a Windows-specific function enabled 
# only on Windows. The intention is to allow the runtime loader to find all 
# the required DLLs without editing the PATH environment variable.
# Copied from OpenSimMacros.cmake.
# In the future, we could use the VS_USER_PROPS_CXX property instead
# https://gitlab.kitware.com/cmake/cmake/commit/ef121ca0c33fb4931007c38b22c046998694b052
function(simbody_copy_dlls DEP_NAME DEP_INSTALL_DIR)
    # On Windows, copy dlls into the Simbody binary directory.
    if(WIN32)
        file(GLOB_RECURSE DLLS ${DEP_INSTALL_DIR}/*.dll)
        if(NOT DLLS)
            message(FATAL_ERROR "Zero DLLs found in directory "
                                "${DEP_INSTALL_DIR}.")
        endif()
        set(DEST_DIR "${CMAKE_BINARY_DIR}/${CMAKE_CFG_INTDIR}")
        set(DLL_NAMES "")
        foreach(DLL IN LISTS DLLS)
            get_filename_component(DLL_NAME ${DLL} NAME)
            list(APPEND DLLS_DEST "${DEST_DIR}/${DLL_NAME}")
            set(DLL_NAMES "${DLL_NAMES} ${DLL_NAME}")
        endforeach()
        add_custom_command(OUTPUT ${DLLS_DEST}
            COMMAND ${CMAKE_COMMAND} -E make_directory ${DEST_DIR}
            COMMAND ${CMAKE_COMMAND} -E copy ${DLLS} ${DEST_DIR}
            DEPENDS ${DLLS}
            COMMENT "Copying DLLs from ${DEP_INSTALL_DIR} to ${DEST_DIR}:${DLL_NAMES}")
        add_custom_target(Copy_${DEP_NAME}_DLLs ALL DEPENDS ${DLLS_DEST})
        set_target_properties(Copy_${DEP_NAME}_DLLs PROPERTIES
            PROJECT_LABEL "Copy ${DEP_NAME} DLLs" FOLDER "simbody")
        install(FILES ${DLLS} DESTINATION ${CMAKE_INSTALL_BINDIR})
    endif()
endfunction()
