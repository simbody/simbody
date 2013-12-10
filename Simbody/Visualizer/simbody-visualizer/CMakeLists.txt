# Generate the visualizer application and deal with its dependencies
# on OpenGL and glut. On Windows we're going to install our own glut; for
# other platforms we depend on their already being a glut (or freeglut)
# available.

INCLUDE (FindOpenGL)
IF(WIN32)
    SET(glut32dir "${CMAKE_CURRENT_SOURCE_DIR}/glut32")
    SET(GLUT32_HEADERS "${glut32dir}/glut.h" "${glut32dir}/glext.h")
    IF(${PLATFORM_ABI} MATCHES "x64")
        SET(glut32libdir "${glut32dir}/lib64")
    ELSE()
        SET(glut32libdir "${glut32dir}/lib")
    ENDIF()
    SET(GLUT_LIBRARIES "${glut32libdir}/glut32.lib")
    # Update the local copy of the dll; destination is Release or Debug.
    # Must fix slashes so that Windows "copy" command can work.
    FILE(TO_NATIVE_PATH  "${glut32libdir}/glut32.dll" gl32src)
    FILE(TO_NATIVE_PATH  "${BUILD_BINARY_DIR}/${CMAKE_CFG_INTDIR}/glut32.dll"
         gl32dest)
ELSE()
    INCLUDE(FindGLUT) # sets GLUT_LIBRARIES
ENDIF()

ADD_EXECUTABLE(${GUI_NAME} 
	simbody-visualizer.cpp lodepng.cpp lodepng.h
	${GLUT32_HEADERS}) # only on Windows

SET_TARGET_PROPERTIES(${GUI_NAME} PROPERTIES
        PROJECT_LABEL "Code - ${GUI_NAME}")

TARGET_LINK_LIBRARIES(${GUI_NAME} 
	${TEST_SHARED_TARGET} ${GLUT_LIBRARIES} ${OPENGL_LIBRARIES})

IF(WIN32)
  # Install visualizer to bin, since it needs to be co-located with dll's
  SET(SIMBODY_VISUALIZER_INSTALL_DIR ${CMAKE_INSTALL_FULL_BINDIR})
ELSE()
  # Visualizer is not intended to be a user executable. Proper place is
  # inside the lib directory
  SET(SIMBODY_VISUALIZER_INSTALL_DIR ${CMAKE_INSTALL_LIBEXECDIR}/simbody)
ENDIF()
INSTALL(TARGETS ${GUI_NAME}
         PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                             GROUP_READ GROUP_WRITE GROUP_EXECUTE
                             WORLD_READ WORLD_EXECUTE
         DESTINATION ${SIMBODY_VISUALIZER_INSTALL_DIR})

# on Windows we also have to copy and later install the glut32.dll
IF(WIN32)
    # This gets executed whenever the named target gets built.
    add_custom_command(
        TARGET ${GUI_NAME}
        COMMAND copy ${gl32src} ${gl32dest}
	COMMENT "Copying glut32.dll" VERBATIM)

    INSTALL(FILES ${glut32libdir}/glut32.dll
	    DESTINATION ${CMAKE_INSTALL_FULL_BINDIR}) #same as GUI
ENDIF()
