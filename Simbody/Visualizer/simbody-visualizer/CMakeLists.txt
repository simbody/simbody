# Generate the visualizer application and deal with its dependencies
# on OpenGL and glut. On Windows we're going to install our own glut; for
# other platforms we depend on their already being a glut (or freeglut)
# available.

find_package(OpenGL)
if(WIN32)
    set(glut32dir "${CMAKE_CURRENT_SOURCE_DIR}/glut32")
    set(GLUT32_HEADERS "${glut32dir}/glut.h" "${glut32dir}/glext.h")
    if(${PLATFORM_ABI} MATCHES "x64")
        set(glut32libdir "${glut32dir}/lib64")
    else()
        set(glut32libdir "${glut32dir}/lib")
    endif()
    set(GLUT_LIBRARIES "${glut32libdir}/glut32.lib")
    # Update the local copy of the dll; destination is Release or Debug.
    # Must fix slashes so that Windows "copy" command can work.
    file(TO_NATIVE_PATH  "${glut32libdir}/glut32.dll" gl32src)
    file(TO_NATIVE_PATH  "${BUILD_BINARY_DIR}/${CMAKE_CFG_INTDIR}/glut32.dll"
         gl32dest)
    set(SIMBODY_HAS_GLUT TRUE)
else()
    find_package(GLUT) # sets GLUT_LIBRARIES
    if(APPLE)
        if(GLUT_FOUND)
            set(SIMBODY_HAS_GLUT TRUE)
        else()
            set(SIMBODY_HAS_GLUT FALSE)
        endif()
    else()
        if(GLUT_FOUND AND GLUT_Xmu_LIBRARY AND GLUT_Xi_LIBRARY)
            set(SIMBODY_HAS_GLUT TRUE)
        else()
            set(SIMBODY_HAS_GLUT FALSE)
        endif()
    endif()
endif()

if(OPENGL_FOUND AND SIMBODY_HAS_GLUT)

add_executable(${GUI_NAME} 
    simbody-visualizer.cpp lodepng.cpp lodepng.h
    ${GLUT32_HEADERS}) # only on Windows

if(NOT WIN32)
    target_include_directories(${GUI_NAME} PRIVATE ${GLUT_INCLUDE_DIR})
endif()

# If building as debug, append the debug postfix to the name of the executable.
# CMAKE_DEBUG_POSTFIX only affects non-executable targets, but we use its value
# to set the postfix for this executable.
set_target_properties(${GUI_NAME} PROPERTIES
        PROJECT_LABEL "Code - ${GUI_NAME}"
        DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})

# On OSX, bake the relative path to the Simbody libraries into the visualizer
# executable. Then there's no need to set `DYLD_LIBRARY_PATH` to find the
# libraries when using the visualizer.
if(${SIMBODY_USE_INSTALL_RPATH})
    # @executable_path only makes sense on OSX, so if we use RPATH on
    # Linux we'll have to revisit.

    # vis_dir_to_install_dir is most likely "../"
    file(RELATIVE_PATH vis_dir_to_install_dir
        "${SIMBODY_VISUALIZER_INSTALL_DIR}"
        "${CMAKE_INSTALL_PREFIX}")
    set(vis_dir_to_lib_dir "${vis_dir_to_install_dir}${CMAKE_INSTALL_LIBDIR}")
    set_target_properties(${GUI_NAME} PROPERTIES
        INSTALL_RPATH "\@executable_path/${vis_dir_to_lib_dir}"
        )
endif()

if(BUILD_DYNAMIC_LIBRARIES)
    target_link_libraries(${GUI_NAME} ${TEST_SHARED_TARGET})
else()
    target_link_libraries(${GUI_NAME} ${TEST_STATIC_TARGET})
endif()
    
target_link_libraries(${GUI_NAME} ${GLUT_LIBRARIES} ${OPENGL_LIBRARIES})

# SIMBODY_VISUALIZER_INSTALL_DIR is set in the root CMakeLists.txt
install(TARGETS ${GUI_NAME}
         PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                             GROUP_READ GROUP_WRITE GROUP_EXECUTE
                             WORLD_READ WORLD_EXECUTE
         DESTINATION ${SIMBODY_VISUALIZER_INSTALL_DIR})

# on Windows we also have to copy and later install the glut32.dll
if(WIN32)
    # This gets executed whenever the named target gets built.
    add_custom_command(
        TARGET ${GUI_NAME}
        COMMAND ${CMAKE_COMMAND} -E copy ${gl32src} ${gl32dest}
        COMMENT "Copying glut32.dll" VERBATIM)

    install(FILES ${glut32libdir}/glut32.dll
        DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

else()
    message(WARNING "Visualizer will not be built because some of its dependencies (OpenGL or GLUT) are missing")
endif()
