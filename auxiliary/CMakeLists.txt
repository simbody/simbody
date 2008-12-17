SUBDIRS (staticTarget sharedTarget)

IF (BUILD_TESTING)
    SUBDIRS (tests)
ENDIF (BUILD_TESTING)

SET(AUX_SHARED_TARGET ${SIMBODY_LIBRARY_NAME}_aux)
SET(AUX_STATIC_TARGET ${SIMBODY_LIBRARY_NAME}_aux_static)

IF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)
    SET(AUX_SHARED_TARGET ${AUX_SHARED_TARGET}_d)
    SET(AUX_STATIC_TARGET ${AUX_STATIC_TARGET}_d)
ENDIF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/include)

FILE(GLOB_RECURSE AUX_SOURCE_FILES ${PROJECT_SOURCE_DIR}/auxiliary/src/*.cpp )
FILE(GLOB_RECURSE AUX_INCLUDE_FILES ${PROJECT_SOURCE_DIR}/auxiliary/src/*.h 
${PROJECT_SOURCE_DIR}/auxiliary/include/*.h )


   ## Find VTK
   SET(VTK_DIR ${SimTK_SDK}/core/lib/vtk-5.0 CACHE PATH "VTK installation")

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

