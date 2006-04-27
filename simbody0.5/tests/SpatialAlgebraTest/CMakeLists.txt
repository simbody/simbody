ADD_EXECUTABLE(SpatialAlgebraTest       SpatialAlgebraTest.cpp)
ADD_EXECUTABLE(SpatialAlgebraTestStatic SpatialAlgebraTest.cpp)

IF( WIN32 )
   TARGET_LINK_LIBRARIES(SpatialAlgebraTest       ${SHARED_TARGET})
   TARGET_LINK_LIBRARIES(SpatialAlgebraTestStatic ${STATIC_TARGET})
ELSE( WIN32 )
   TARGET_LINK_LIBRARIES(SpatialAlgebraTest       ${SHARED_TARGET} pthread )
   TARGET_LINK_LIBRARIES(SpatialAlgebraTestStatic ${STATIC_TARGET} pthread )
ENDIF( WIN32 )

