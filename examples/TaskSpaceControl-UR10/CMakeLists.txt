# Build and install TaskSpaceControl-UR10 example, from within a
# Simbody source build. This must be invoked as a CMake subdirectory from
# the main examples CMakeLists.txt file.

set(header_files UR10.h)
set(source_files TaskSpaceControl-UR10.cpp UR10.cpp)

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

# Copy geometry to the binary directory that will 
# be the working directory when the example is run from a source build;
# and install it in the examples installation.
FOREACH(extradir geometry)
  FILE(GLOB extradir_file  "${extradir}/*")
  FOREACH(xfile ${extradir_file})
    GET_FILENAME_COMPONENT(xfile_name ${xfile} NAME)
    CONFIGURE_FILE(${xfile}
       ${CMAKE_CURRENT_BINARY_DIR}/${extradir}/${xfile_name} COPYONLY)
    INSTALL(FILES ${xfile} DESTINATION 
            ${EXAMPLES_INSTALL_SRC}/${EX_NAME}/${extradir})
  ENDFOREACH()
ENDFOREACH()

# install source for example
INSTALL(FILES ${source_files} ${header_files}
        DESTINATION ${EXAMPLES_INSTALL_SRC}/${EX_NAME})

# install the installed version of CMakeLists.txt
#INSTALL(FILES InstalledCMakeLists.txt DESTINATION ${EXAMPLES_INSTALL_SRC}
#    RENAME CMakeLists.txt)
