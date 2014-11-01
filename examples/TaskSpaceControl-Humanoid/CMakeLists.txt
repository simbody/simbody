# Build and install example, from within a
# Simbody source build. This must be invoked as a CMake subdirectory from
# the main examples CMakeLists.txt file.

set(header_files Humanoid.h)
set(source_files TaskSpaceHumanoid.cpp)

ADD_DEFINITIONS("-DSIMBODY_EXAMPLE_NAME=\"${EX_NAME}\"")

IF (BUILD_TESTS_AND_EXAMPLES_SHARED)
    # Link with shared library
    ADD_EXECUTABLE(${EX_NAME} ${source_files} ${header_files})
    IF(GUI_NAME)
        ADD_DEPENDENCIES(${EX_NAME} ${GUI_NAME})
    ENDIF()
    SET_TARGET_PROPERTIES(${EX_NAME}
                          PROPERTIES
                          PROJECT_LABEL "Example - ${EX_NAME}")
    TARGET_LINK_LIBRARIES(${EX_NAME} ${EXAMPLES_SHARED_TARGET})
    INSTALL(TARGETS ${EX_NAME} DESTINATION ${EXAMPLES_INSTALL_BIN})
ENDIF()

IF (BUILD_STATIC_LIBRARIES AND BUILD_TESTS_AND_EXAMPLES_STATIC)
    # Link with static library
    SET(ex_static ${EX_NAME}Static)
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

IF(NOT WIN32)
  EXECUTE_PROCESS(COMMAND "${CMAKE_COMMAND}" -E create_symlink
                  ${EXAMPLES_INSTALL_BIN}
                  ${CMAKE_CURRENT_BINARY_DIR}/bin
  )
  INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin
          DESTINATION ${EXAMPLES_SYMLINK_BIN}
  )
ENDIF()

# install source for example
INSTALL(FILES ${source_files} ${header_files}
        DESTINATION ${EXAMPLES_INSTALL_SRC}/${EX_NAME})

# install the installed version of CMakeLists.txt
#INSTALL(FILES InstalledCMakeLists.txt DESTINATION ${EXAMPLES_INSTALL_SRC}
#    RENAME CMakeLists.txt)
