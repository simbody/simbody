ADD_EXECUTABLE(SpatialAlgebraTest       SpatialAlgebraTest.cpp)
ADD_EXECUTABLE(SpatialAlgebraTestStatic SpatialAlgebraTestStatic.cpp)

TARGET_LINK_LIBRARIES(SpatialAlgebraTest       ${SHARED_TARGET})
TARGET_LINK_LIBRARIES(SpatialAlgebraTestStatic ${STATIC_TARGET})

SET_TARGET_PROPERTIES(SpatialAlgebraTest PROPERTIES
	PROJECT_LABEL "Test - SpatialAlgebraTest")

SET_TARGET_PROPERTIES(SpatialAlgebraTestStatic PROPERTIES
	PROJECT_LABEL "Test - SpatialAlgebraTestStatic")