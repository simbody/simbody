# Build and install TaskSpaceControl-Atlas example, from within a
# Simbody source build. This must be invoked as a CMake subdirectory from
# the main examples CMakeLists.txt file.

set(ex_name TaskSpaceControl-Atlas)
set(header_files Atlas.h URDFReader.h)
set(source_files TaskSpaceControl-Atlas.cpp Atlas.cpp URDFReader.cpp)

ADD_DEFINITIONS(
    -DSIMBODY_EXAMPLE_INSTALL_DIR="${EXAMPLES_INSTALL_SRC}/${ex_name}/")

IF (BUILD_TESTS_AND_EXAMPLES_SHARED)
    # Link with shared library
    ADD_EXECUTABLE(${ex_name} ${source_files} ${header_files})
    IF(GUI_NAME)
        ADD_DEPENDENCIES(${ex_name} ${GUI_NAME})
    ENDIF()
    SET_TARGET_PROPERTIES(${ex_name}
                          PROPERTIES
                          PROJECT_LABEL "Example - ${ex_name}")
    TARGET_LINK_LIBRARIES(${ex_name} ${EXAMPLES_SHARED_TARGET})
    INSTALL(TARGETS ${ex_name} DESTINATION ${EXAMPLES_INSTALL_BIN})
ENDIF()

IF (BUILD_STATIC_LIBRARIES AND BUILD_TESTS_AND_EXAMPLES_STATIC)
    # Link with static library
    SET(ex_static ${ex_name}Static)
    ADD_EXECUTABLE(${ex_static} ${source_files} ${header_files})
    IF(GUI_NAME)
        ADD_DEPENDENCIES(${ex_static} ${GUI_NAME})
    ENDIF()
    SET_TARGET_PROPERTIES(${ex_static}
                          PROPERTIES
                          COMPILE_FLAGS "-DSimTK_USE_STATIC_LIBRARIES"
                          PROJECT_LABEL "Example - ${ex_static}")
    TARGET_LINK_LIBRARIES(${ex_static} ${EXAMPLES_STATIC_TARGET})
    # Don't install static examples
ENDIF()

IF(UNIX)
  EXECUTE_PROCESS(COMMAND "${CMAKE_COMMAND}" -E create_symlink
                  ${CMAKE_BINARY_DIR}/${ex_name}
                  ${CMAKE_CURRENT_BINARY_DIR}/${ex_name}
  )
ENDIF()

# Copy geometry and models to the binary directory that will 
# be the working directory when the example is run from a source build;
# and install it in the examples installation.
FOREACH(extradir geometry models)
  FILE(GLOB extradir_file  "${extradir}/*")
  FOREACH(xfile ${extradir_file})
    GET_FILENAME_COMPONENT(xfile_name ${xfile} NAME)
    CONFIGURE_FILE(${xfile}
        ${CMAKE_CURRENT_BINARY_DIR}/${extradir}/${xfile_name} COPYONLY)
    INSTALL(FILES ${xfile} DESTINATION 
            ${EXAMPLES_INSTALL_SRC}/${ex_name}/${extradir})
    # This allows to run this example in the build tree, assuming
    # no files get overwritten.
    IF(UNIX)
      EXECUTE_PROCESS(COMMAND "${CMAKE_COMMAND}" -E create_symlink
          ${CMAKE_CURRENT_BINARY_DIR}/${extradir}/${xfile_name}
          ${CMAKE_BINARY_DIR}/${extradir}/${xfile_name}
          )
    ENDIF()
  ENDFOREACH()
ENDFOREACH()

# install source for example
INSTALL(FILES ${source_files} ${header_files}
        DESTINATION ${EXAMPLES_INSTALL_SRC}/${ex_name})

# install the installed version of CMakeLists.txt
#INSTALL(FILES InstalledCMakeLists.txt DESTINATION ${EXAMPLES_INSTALL_SRC}
#    RENAME CMakeLists.txt)
