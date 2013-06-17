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
    FILE(TO_NATIVE_PATH  "${BUILD_BINARY_DIR}/${CMAKE_CFG_INTDIR}/glut32.dll"
         gl32dest)
ELSE()
    INCLUDE(FindGLUT) # sets GLUT_LIBRARIES
ENDIF()

SET(VISUALIZER_BASE_SOURCES 
	#	Scene.h 
	#	VisualizerBase.h 
	#VisualizerBase.cpp
	)

#Ogre

if(UNIX)
	if(EXISTS "/usr/local/lib/OGRE/cmake")

	  set(CMAKE_MODULE_PATH "/usr/local/lib/OGRE/cmake/;${CMAKE_MODULE_PATH}")
	  set(OGRE_SAMPLES_INCLUDEPATH "/usr/local/share/OGRE/samples/Common/include/") # We could just *assume* that developers uses this basepath : /usr/local

    elseif(EXISTS "/usr/share/OGRE/cmake")

	  set(CMAKE_MODULE_PATH "/usr/share/OGRE/cmake/;${CMAKE_MODULE_PATH}")
	  set(OGRE_SAMPLES_INCLUDEPATH "/usr/share/OGRE-1.7.4/Samples/Common/include/")

	elseif(EXISTS "/usr/lib/OGRE/cmake")

	  set(CMAKE_MODULE_PATH "/usr/lib/OGRE/cmake/;${CMAKE_MODULE_PATH}")
	  set(OGRE_SAMPLES_INCLUDEPATH "/usr/share/OGRE/samples/Common/include/") # Otherwise, this one

	else ()
	  message(SEND_ERROR "Failed to find module path.")
	endif()
endif(UNIX)

ADD_EXECUTABLE(${GUI_NAME} ${VISUALIZER_BASE_SOURCES}
	VisualizerGUI.cpp lodepng.cpp lodepng.h
	${GLUT32_HEADERS}) # only on Windows

SET_TARGET_PROPERTIES(${GUI_NAME} PROPERTIES
        PROJECT_LABEL "Code - ${GUI_NAME}")

TARGET_LINK_LIBRARIES(${GUI_NAME} 
	${TEST_SHARED_TARGET} ${GLUT_LIBRARIES} ${OPENGL_LIBRARIES} ${Boost_LIBRARIES})

INSTALL(TARGETS ${GUI_NAME}
         PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                             GROUP_READ GROUP_WRITE GROUP_EXECUTE
                             WORLD_READ WORLD_EXECUTE
         DESTINATION ${SimTK_INSTALL_PREFIX}/bin)

# on Windows we also have to copy and later install the glut32.dll
IF(WIN32)
    # This gets executed whenever the named target gets built.
    add_custom_command(
        TARGET ${GUI_NAME}
        COMMAND copy ${gl32src} ${gl32dest}
	COMMENT "Copying glut32.dll" VERBATIM)

    INSTALL(FILES ${glut32libdir}/glut32.dll
	    DESTINATION ${SimTK_INSTALL_PREFIX}/bin) #same as GUI
ENDIF()
