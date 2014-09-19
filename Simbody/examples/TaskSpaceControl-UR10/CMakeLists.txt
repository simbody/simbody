# Build and install TaskSpaceControl-UR10 example, from within a
# Simbody source build. This must be invoked as a CMake subdirectory from
# the main examples CMakeLists.txt file.

set(ex_name TaskSpaceControl-UR10)
set(header_files UR10.h)
set(source_files TaskSpaceControl-UR10.cpp UR10.cpp)

IF (BUILD_TESTS_AND_EXAMPLES_SHARED)
    # Link with shared library
    ADD_EXECUTABLE(${ex_name} ${source_files} ${header_files})
    IF(GUI_NAME)
        ADD_DEPENDENCIES(${ex_name} ${GUI_NAME})
    ENDIF()
    SET_TARGET_PROPERTIES(${ex_name}
                          PROPERTIES
                          PROJECT_LABEL "Example - ${ex_name}")
    TARGET_LINK_LIBRARIES(${ex_name} ${TEST_SHARED_TARGET})
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
    TARGET_LINK_LIBRARIES(${ex_static} ${TEST_STATIC_TARGET})
    # Don't install static examples
ENDIF()

IF(NOT WIN32)
  EXECUTE_PROCESS(COMMAND "${CMAKE_COMMAND}" -E create_symlink
                  ${EXAMPLES_INSTALL_BIN}
                  ${CMAKE_CURRENT_BINARY_DIR}/bin
  )
  INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin
          DESTINATION ${EXAMPLES_SYMLINK_BIN}
  )
ENDIF()

# Copy geometry to the binary directory that will 
# be the working directory when the example is run from a source build;
# and install it in the examples installation (both for binaries and
# for source).
FOREACH(geodir geometry)
  FILE(GLOB geodir_geom  "${geodir}/*")
  FOREACH(geom ${geodir_geom})
    GET_FILENAME_COMPONENT(geom_name ${geom} NAME)
    set(relname ${geodir}/${geom_name})
    CONFIGURE_FILE(${geom}
                   ${CMAKE_CURRENT_BINARY_DIR}/${relname} COPYONLY)
    INSTALL(FILES ${geom} DESTINATION ${EXAMPLES_INSTALL_SRC}/${geodir})
    INSTALL(FILES ${geom} DESTINATION ${EXAMPLES_INSTALL_BIN}/${geodir})
  ENDFOREACH()
ENDFOREACH()

# install source for example
INSTALL(FILES ${source_files} ${header_files}
        DESTINATION ${EXAMPLES_INSTALL_SRC}/${ex_name})

# install the installed version of CMakeLists.txt
#INSTALL(FILES InstalledCMakeLists.txt DESTINATION ${EXAMPLES_INSTALL_SRC}
#    RENAME CMakeLists.txt)
