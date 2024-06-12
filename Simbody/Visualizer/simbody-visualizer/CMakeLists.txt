# Generate the visualizer application and deal with its dependencies
# on OpenGL and glut. On Windows we're going to install our own glut; for
# other platforms we depend on their already being a glut (or freeglut)
# available.

set(MISSING_DEPENDENCIES "")

set(OpenGL_GL_PREFERENCE LEGACY)
find_package(OpenGL)
if (NOT OPENGL_FOUND)
    list(APPEND MISSING_DEPENDENCIES OpenGL)
endif()

if(WIN32)
    set(glut32dir "${CMAKE_CURRENT_SOURCE_DIR}/glut32")
    set(GLUT32_HEADERS "${glut32dir}/glut.h" "${glut32dir}/glext.h")
    if(${PLATFORM_ABI} MATCHES "x64")
        set(glut32libdir "${glut32dir}/lib64")
    else()
        set(glut32libdir "${glut32dir}/lib")
    endif()

    add_library(glut32 SHARED IMPORTED)
    set_target_properties(glut32 PROPERTIES
        IMPORTED_IMPLIB "${glut32libdir}/glut32.lib"
        IMPORTED_LOCATION "${glut32libdir}/glut32.dll"
        )

    set(GLUT_LIBRARIES glut32)
    set(SIMBODY_HAS_GLUT TRUE)
else()
    find_package(GLUT) # sets GLUT_LIBRARIES
    if(APPLE)
        if(GLUT_FOUND)
            set(SIMBODY_HAS_GLUT TRUE)
        else()
            set(SIMBODY_HAS_GLUT FALSE)
            list(APPEND MISSING_DEPENDENCIES GLUT)
        endif()
    else()
        if(GLUT_FOUND AND GLUT_Xmu_LIBRARY AND GLUT_Xi_LIBRARY)
            set(SIMBODY_HAS_GLUT TRUE)
        else()
            set(SIMBODY_HAS_GLUT FALSE)

            if (NOT GLUT_FOUND)
                list(APPEND MISSING_DEPENDENCIES GLUT)
            endif()
            if (NOT GLUT_Xmu_LIBRARY)
                list(APPEND MISSING_DEPENDENCIES libxmu)
            endif()
            if (NOT GLUT_Xi_LIBRARY)
                list(APPEND MISSING_DEPENDENCIES libxi)
            endif()
        endif()
    endif()
endif()

if(MISSING_DEPENDENCIES)
    message(WARNING "Visualizer will not be built because some of its dependencies (${MISSING_DEPENDENCIES}) are missing")
else()
    add_executable(simbody-visualizer MACOSX_BUNDLE
        simbody-visualizer.cpp lodepng.cpp lodepng.h
        ${GLUT32_HEADERS}) # only on Windows

    target_link_libraries(simbody-visualizer PUBLIC SimTKsimbody ${GLUT_LIBRARIES} ${OPENGL_LIBRARIES})

    target_include_directories(simbody-visualizer
        PRIVATE
            $<$<BOOL:${WIN32}>:${GLUT_INCLUDE_DIR}>)

    # If building as debug, append the debug postfix to the name of the executable.
    # CMAKE_DEBUG_POSTFIX only affects non-executable targets, but we use its value
    # to set the postfix for this executable.
    # Setting the target property DEBUG_POSTFIX does not work for MACOSX bundles,
    # so we use DEBUG_OUTPUT_NAME instead.
    set_target_properties(simbody-visualizer PROPERTIES
            PROJECT_LABEL "Code - simbody-visualizer"
            DEBUG_OUTPUT_NAME simbody-visualizer${CMAKE_DEBUG_POSTFIX})

    # Copy third-party libs (e.g. vendored GLUT) to the build directory
    add_custom_command(TARGET simbody-visualizer PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy -t $<TARGET_FILE_DIR:simbody-visualizer> $<TARGET_RUNTIME_DLLS:simbody-visualizer>
        COMMAND_EXPAND_LISTS
    )

    # On OSX, bake the relative path to the Simbody libraries into the visualizer
    # executable. Then there's no need to set `DYLD_LIBRARY_PATH` to find the
    # libraries when using the visualizer.
    if(${SIMBODY_USE_INSTALL_RPATH})
        # @executable_path only makes sense on OSX, so if we use RPATH on
        # Linux we'll have to revisit.

        # vis_dir_to_install_dir is most likely "../"
        if (APPLE)
            set(vis_install_dir
                    "${SIMBODY_VISUALIZER_INSTALL_DIR}/simbody-visualizer.app/Contents/MacOS")
        else()
            set(vis_install_dir "${SIMBODY_VISUALIZER_INSTALL_DIR}")
        endif()
        file(RELATIVE_PATH vis_dir_to_install_dir
                "${vis_install_dir}"
                "${CMAKE_INSTALL_PREFIX}")
        set(vis_dir_to_lib_dir "${vis_dir_to_install_dir}${CMAKE_INSTALL_LIBDIR}")
        set_target_properties(simbody-visualizer PROPERTIES
            INSTALL_RPATH "\@executable_path/${vis_dir_to_lib_dir}"
            )
    endif()

    # SIMBODY_VISUALIZER_INSTALL_DIR is set in the root CMakeLists.txt
    install(TARGETS simbody-visualizer
            PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                                GROUP_READ GROUP_WRITE GROUP_EXECUTE
                                WORLD_READ WORLD_EXECUTE
            DESTINATION ${SIMBODY_VISUALIZER_INSTALL_DIR})

    if(WIN32)
        # on Windows we also have to copy and later install the glut32.dll
        # no need to specify permissions on Windows
        install(IMPORTED_RUNTIME_ARTIFACTS simbody-visualizer
            RUNTIME_DEPENDENCY_SET glut-vendored-dlls
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
        install(RUNTIME_DEPENDENCY_SET glut-vendored-dlls
            PRE_EXCLUDE_REGEXES "^(api|ext)-ms"
            POST_EXCLUDE_REGEXES ".*system32/.*\\.dll"
            )
    endif()

endif()
