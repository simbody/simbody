# Build and install example, from within a
# Simbody source build. This must be invoked as a CMake subdirectory from
# the main examples CMakeLists.txt file.

set(ex_name TaskSpaceControl-Humanoid)
set(header_files Humanoid.h)
set(source_files TaskSpaceHumanoid.cpp)

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

# install source for example
INSTALL(FILES ${source_files} ${header_files}
        DESTINATION ${EXAMPLES_INSTALL_SRC}/${ex_name})

# install the installed version of CMakeLists.txt
#INSTALL(FILES InstalledCMakeLists.txt DESTINATION ${EXAMPLES_INSTALL_SRC}
#    RENAME CMakeLists.txt)
