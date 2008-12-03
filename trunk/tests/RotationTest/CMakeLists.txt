ADD_EXECUTABLE(RotationTest       RotationTest.cpp)
ADD_EXECUTABLE(RotationTestStatic RotationTestStatic.cpp)

TARGET_LINK_LIBRARIES(RotationTest       ${SHARED_TARGET})
TARGET_LINK_LIBRARIES(RotationTestStatic ${STATIC_TARGET})

SET_TARGET_PROPERTIES(RotationTest PROPERTIES
	PROJECT_LABEL "Test - RotationTest")

SET_TARGET_PROPERTIES(RotationTestStatic PROPERTIES
	PROJECT_LABEL "Test - RotationTestStatic")
