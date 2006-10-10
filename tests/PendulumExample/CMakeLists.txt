
SET(VTK_DIR "C:/Program Files/VTK/lib/vtk-5.0" CACHE STRING "VTK installation")

FIND_PACKAGE(VTK)
IF (VTK_FOUND)
  INCLUDE(${VTK_USE_FILE})
ELSE(VTK_FOUND)
  MESSAGE(FATAL_ERROR
     "Cannot build without VTK.  Please set VTK_DIR.")
ENDIF (VTK_FOUND)

LINK_LIBRARIES(
  vtkHybrid
  vtkFiltering
  vtkGenericFiltering
  vtkVolumeRendering 
  vtkRendering
  vtkIO
  vtkGraphics
  vtkCommon
)

ADD_EXECUTABLE(PendulumExample       PendulumExample.cpp)
ADD_EXECUTABLE(PendulumExampleStatic PendulumExampleStatic.cpp)

TARGET_LINK_LIBRARIES(PendulumExample       ${SHARED_TARGET})
TARGET_LINK_LIBRARIES(PendulumExampleStatic ${STATIC_TARGET})