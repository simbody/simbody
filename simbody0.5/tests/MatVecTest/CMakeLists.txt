ADD_EXECUTABLE(MatVecTest       MatVecTest.cpp)
ADD_EXECUTABLE(MatVecTestStatic MatVecTest.cpp)

IF( WIN32 )
   TARGET_LINK_LIBRARIES(MatVecTest       ${SHARED_TARGET})
   TARGET_LINK_LIBRARIES(MatVecTestStatic ${STATIC_TARGET})
ELSE( WIN32 )
   TARGET_LINK_LIBRARIES(MatVecTest       ${SHARED_TARGET} pthread )
   TARGET_LINK_LIBRARIES(MatVecTestStatic ${STATIC_TARGET} pthread )
ENDIF( WIN32 )

