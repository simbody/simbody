# Build and install TaskSpaceControl-UR10 example, from within a
# Simbody source build. This must be invoked as a CMake subdirectory from
# the main examples CMakeLists.txt file.

set(header_files UR10.h)
set(shared_header_files
    ../shared/TaskSpace.h ../shared/SimbodyExampleHelper.h)

set(source_files TaskSpaceControl-UR10.cpp UR10.cpp)
set(shared_source_files ../shared/TaskSpace.cpp)

set(all_header_files ${header_files} ${shared_header_files})
set(all_source_files ${source_files} ${shared_source_files})

add_definitions("-DSIMBODY_EXAMPLE_NAME=\"${EX_NAME}\"")

if(BUILD_TESTS_AND_EXAMPLES_SHARED)
    # Link with shared library
    add_executable(${EX_NAME} ${all_source_files} ${all_header_files})
    if(GUI_NAME)
        add_dependencies(${EX_NAME} ${GUI_NAME})
    endif()
    set_target_properties(${EX_NAME}
                          PROPERTIES
                          PROJECT_LABEL "Example - ${EX_NAME}")
    target_link_libraries(${EX_NAME} ${EXAMPLES_SHARED_TARGET})
    # Don't install Debug examples
    install(TARGETS ${EX_NAME}
            DESTINATION ${EXAMPLES_INSTALL_BIN}
            CONFIGURATIONS Release RelWithDebInfo MinSizeRel)
endif()

if(BUILD_STATIC_LIBRARIES AND BUILD_TESTS_AND_EXAMPLES_STATIC)
    # Link with static library
    set(ex_static ${EX_NAME}Static)
    add_executable(${ex_static} ${all_source_files} ${all_header_files})
    if(GUI_NAME)
        add_dependencies(${ex_static} ${GUI_NAME})
    endif()
    set_target_properties(${ex_static}
                          PROPERTIES
                          COMPILE_FLAGS "-DSimTK_USE_STATIC_LIBRARIES"
                          PROJECT_LABEL "Example - ${ex_static}")
    target_link_libraries(${ex_static} ${EXAMPLES_STATIC_TARGET})
    # Don't install static examples
endif()

# Copy geometry to the binary directory that will
# be the working directory when the example is run from a source build;
# and install it in the examples installation.
foreach(extradir geometry models)
  file(GLOB extradir_file  "${extradir}/*")
  foreach(xfile ${extradir_file})
    get_filename_component(xfile_name ${xfile} NAME)
    configure_file(${xfile}
       ${CMAKE_CURRENT_BINARY_DIR}/${extradir}/${xfile_name} COPYONLY)
    install(FILES ${xfile} DESTINATION
            ${EXAMPLES_INSTALL_SRC}/${EX_NAME}/${extradir})
  endforeach()
endforeach()

# install source for example
install(FILES ${source_files} ${header_files}
        DESTINATION ${EXAMPLES_INSTALL_SRC}/${EX_NAME})

# install the installed version of CMakeLists.txt
install(FILES InstalledCMakeLists.txt
        DESTINATION ${EXAMPLES_INSTALL_SRC}/${EX_NAME}
        RENAME CMakeLists.txt)
