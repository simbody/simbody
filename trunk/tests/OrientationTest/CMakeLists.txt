
ADD_EXECUTABLE(OrientationTest       OrientationTest.cpp)
ADD_EXECUTABLE(OrientationTestStatic OrientationTestStatic.cpp)

TARGET_LINK_LIBRARIES(OrientationTest       ${SHARED_TARGET})
TARGET_LINK_LIBRARIES(OrientationTestStatic ${STATIC_TARGET})