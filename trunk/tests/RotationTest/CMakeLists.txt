ADD_EXECUTABLE(RotationTest       RotationTest.cpp)
ADD_EXECUTABLE(RotationTestStatic RotationTestStatic.cpp)

TARGET_LINK_LIBRARIES(RotationTest       ${SHARED_TARGET})
TARGET_LINK_LIBRARIES(RotationTestStatic ${STATIC_TARGET})
