# Generate the Visualization GUI application.

#IF (WIN32)
#  SET(GLUT_LIBRARIES ${CMAKE_CURRENT_SOURCE_DIR}/glut32.lib)
#  INSTALL(FILES glut32.dll glut32.lib DESTINATION ${SimTK_INSTALL_PREFIX}/lib)
#ELSE (WIN32)
#  INCLUDE (FindGLUT)
#ENDIF (WIN32)
#use local freeglut instead
SET(GLUT_LIBRARIES)
INCLUDE (FindOpenGL)
SET(WINMM)
IF(WIN32)
    SET(WINMM winmm)
ENDIF(WIN32)
SET(GUI_NAME "VisualizationGUI")
INCLUDE_DIRECTORIES(${SimTK_SDK}/include freeglut)
FILE(GLOB freeglutsrc freeglut/*.c)
FILE(GLOB freeglutincl freeglut/*.h freeglut/GL/*.h)
ADD_DEFINITIONS(-DFREEGLUT_STATIC -DFREEGLUT_LIB_PRAGMAS=0)
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
