SUBDIRS (staticTarget sharedTarget tests)

SET(AUX_SHARED_TARGET ${SIMBODY_LIBRARY_NAME}_aux)
SET(AUX_STATIC_TARGET ${SIMBODY_LIBRARY_NAME}_aux_static)

IF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)
    SET(AUX_SHARED_TARGET ${AUX_SHARED_TARGET}_d)
    SET(AUX_STATIC_TARGET ${AUX_STATIC_TARGET}_d)
ENDIF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)

FILE(GLOB AUX_SOURCE_FILES ${PROJECT_SOURCE_DIR}/auxiliary/src/*.cpp )
FILE(GLOB AUX_INCLUDE_FILES ${PROJECT_SOURCE_DIR}/auxiliary/src/*.h )


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



ADD_TEST(ContactTest                  ${TEST_PATH}/ContactTest)
ADD_TEST(ContactTestStatic            ${TEST_PATH}/ContactTestStatic)
ADD_TEST(PendulumExample              ${TEST_PATH}/PendulumExample)
ADD_TEST(PendulumExampleStatic        ${TEST_PATH}/PendulumExampleStatic)
ADD_TEST(TwoPendulums                 ${TEST_PATH}/TwoPendulums)
ADD_TEST(TwoPendulumsStatic           ${TEST_PATH}/TwoPendulumsStatic)
ADD_TEST(SBPendulum1                  ${TEST_PATH}/SBPendulum1)
ADD_TEST(SBPendulum1Static            ${TEST_PATH}/SBPendulum1Static)
ADD_TEST(CoarseRNA                    ${TEST_PATH}/CoarseRNA)
ADD_TEST(CoarseRNAStatic              ${TEST_PATH}/CoarseRNAStatic)

