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

add_executable(${EX_NAME} ${all_source_files} ${all_header_files})
if(TARGET simbody-visualizer)
    add_dependencies(${EX_NAME} simbody-visualizer)
endif()
set_target_properties(${EX_NAME}
    PROPERTIES
    FOLDER "Examples" # Organize target list in IDE
    PROJECT_LABEL "Example - ${EX_NAME}")
target_link_libraries(${EX_NAME} SimTKsimbody)

if(SIMBODY_BUILD_SHARED_LIBS) # Don't install static examples
    install(TARGETS ${EX_NAME}
        DESTINATION ${EXAMPLES_INSTALL_BIN}
        # Don't install Debug examples
        CONFIGURATIONS Release RelWithDebInfo MinSizeRel)
endif()

# Where this example's source and data files are installed.
cmake_path(APPEND EXAMPLES_INSTALL_SRC ${EX_NAME} OUTPUT_VARIABLE example_src_dest)

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
            ${example_src_dest}/${extradir})
  endforeach()
endforeach()

# install source for example
install(FILES ${source_files} ${header_files}
        DESTINATION ${example_src_dest})

# install the installed version of CMakeLists.txt
install(FILES InstalledCMakeLists.txt 
        DESTINATION ${example_src_dest}
        RENAME CMakeLists.txt)
