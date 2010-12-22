# Generate the VisualizerGUI application and deal with its dependencies
# on OpenGL and glut. On Windows we're going to install our own glut; for
# other platforms we depend on their already being a glut (or freeglut)
# available.

INCLUDE (FindOpenGL)
IF(WIN32)
    SET(glut32dir "${CMAKE_CURRENT_SOURCE_DIR}/glut32")
    SET(GLUT32_HEADERS "${glut32dir}/glut.h" "${glut32dir}/glext.h")
    SET(glut32libdir "${glut32dir}/lib${LIB64}")
    SET(GLUT_LIBRARIES "${glut32libdir}/glut32.lib")
    # Update the local copy of the dll; destination is Release or Debug.
    # Must fix slashes so that Windows "copy" command can work.
    FILE(TO_NATIVE_PATH  "${glut32libdir}/glut32.dll" gl32src)
    FILE(TO_NATIVE_PATH  "${PROJECT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/glut32.dll"
         gl32dest)
    # This gets executed whenever the named target gets built.
    add_custom_command(
        TARGET ${GUI_NAME}
        COMMAND copy ${gl32src} ${gl32dest}
	COMMENT "Copying glut32.dll" VERBATIM)
ELSE()
    INCLUDE(FindGLUT) # sets GLUT_LIBRARIES
ENDIF()

ADD_EXECUTABLE(${GUI_NAME} 
	VisualizerGUI.cpp lodepng.cpp lodepng.h
	${GLUT32_HEADERS}) # only on Windows

TARGET_LINK_LIBRARIES(${GUI_NAME} 
	${TEST_SHARED_TARGET} ${GLUT_LIBRARIES} ${OPENGL_LIBRARIES})

IF(WIN32)
    # Update the local copy of the dll; destination is Release or Debug.
    # Must fix slashes so that Windows "copy" command can work.
    FILE(TO_NATIVE_PATH  "${glut32libdir}/glut32.dll" gl32src)
    FILE(TO_NATIVE_PATH  "${PROJECT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/glut32.dll"
         gl32dest)
    # This gets executed whenever the named target gets built.
    add_custom_command(
        TARGET ${GUI_NAME}
        COMMAND copy ${gl32src} ${gl32dest}
	COMMENT "Copying glut32.dll" VERBATIM)
ENDIF()

INSTALL(TARGETS ${GUI_NAME}
         PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                             GROUP_READ GROUP_WRITE GROUP_EXECUTE
                             WORLD_READ WORLD_EXECUTE
         DESTINATION ${SimTK_INSTALL_PREFIX}/bin)

# on Windows we also have to install the glut32.dll
IF(WIN32)
    INSTALL(FILES ${glut32libdir}/glut32.dll
	    DESTINATION ${SimTK_INSTALL_PREFIX}/bin) #same as GUI
ENDIF()
