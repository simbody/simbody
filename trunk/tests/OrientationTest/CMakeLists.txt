
ADD_EXECUTABLE(OrientationTest       OrientationTest.cpp)
ADD_EXECUTABLE(OrientationTestStatic OrientationTestStatic.cpp)

TARGET_LINK_LIBRARIES(OrientationTest       ${SHARED_TARGET})
TARGET_LINK_LIBRARIES(OrientationTestStatic ${STATIC_TARGET})

SET_TARGET_PROPERTIES(OrientationTest PROPERTIES
	PROJECT_LABEL "Test - OrientationTest")

SET_TARGET_PROPERTIES(OrientationTestStatic PROPERTIES
	PROJECT_LABEL "Test - OrientationTestStatic")