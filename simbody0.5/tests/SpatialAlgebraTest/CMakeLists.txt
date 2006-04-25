ADD_EXECUTABLE(SpatialAlgebraTest       SpatialAlgebraTest.cpp)
ADD_EXECUTABLE(SpatialAlgebraTestStatic SpatialAlgebraTest.cpp)

TARGET_LINK_LIBRARIES(SpatialAlgebraTest       ${SHARED_TARGET})
TARGET_LINK_LIBRARIES(SpatialAlgebraTestStatic ${STATIC_TARGET})