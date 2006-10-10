
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

ADD_EXECUTABLE(CoarseRNA       CoarseRNA.cpp)
ADD_EXECUTABLE(CoarseRNAStatic CoarseRNAStatic.cpp)

TARGET_LINK_LIBRARIES(CoarseRNA       ${SHARED_TARGET})
TARGET_LINK_LIBRARIES(CoarseRNAStatic ${STATIC_TARGET})