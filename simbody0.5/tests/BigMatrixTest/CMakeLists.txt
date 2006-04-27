ADD_EXECUTABLE(BigMatrixTest       BigMatrixTest.cpp)
ADD_EXECUTABLE(BigMatrixTestStatic BigMatrixTest.cpp)

IF( WIN32 )
   TARGET_LINK_LIBRARIES(BigMatrixTest       ${SHARED_TARGET})
   TARGET_LINK_LIBRARIES(BigMatrixTestStatic ${STATIC_TARGET})
ELSE( WIN32 )
   TARGET_LINK_LIBRARIES(BigMatrixTest       ${SHARED_TARGET} pthread )
   TARGET_LINK_LIBRARIES(BigMatrixTestStatic ${STATIC_TARGET} pthread )
ENDIF( WIN32 )

