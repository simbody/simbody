# Generate the Visualization GUI application.

INCLUDE (FindOpenGL)
SET(GLUT_LIBRARIES)
SET(WINMM)
SET(freeglutsrc)
SET(freeglutincl)
IF(WIN32)
    #use local freeglut
    INCLUDE_DIRECTORIES(freeglut)
    FILE(GLOB freeglutsrc freeglut/*.c)
    FILE(GLOB freeglutincl freeglut/*.h freeglut/GL/*.h)
    ADD_DEFINITIONS(-DFREEGLUT_STATIC -DFREEGLUT_LIB_PRAGMAS=0)
    SET(WINMM winmm)
ELSE()
  INCLUDE (FindGLUT) # must be freeglut
ENDIF(WIN32)

SET(GUI_NAME "VisualizationGUI")
INCLUDE_DIRECTORIES(${SimTK_SDK}/include)
ADD_EXECUTABLE(${GUI_NAME} 
	VisualizationGUI.cpp lodepng.cpp lodepng.h
	${freeglutsrc} ${freeglutincl})

TARGET_LINK_LIBRARIES(${GUI_NAME} 
	${TEST_SHARED_TARGET} ${GLUT_LIBRARIES} ${OPENGL_LIBRARIES} ${WINMM})

INSTALL(TARGETS ${GUI_NAME}
         PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                             GROUP_READ GROUP_WRITE GROUP_EXECUTE
                             WORLD_READ WORLD_EXECUTE
         RUNTIME DESTINATION ${SimTK_INSTALL_PREFIX}/bin)
