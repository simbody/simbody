ADD_EXECUTABLE(RandomTest       RandomTest.cpp)
ADD_EXECUTABLE(RandomTestStatic RandomTestStatic.cpp)

TARGET_LINK_LIBRARIES(RandomTest       ${SHARED_TARGET})
TARGET_LINK_LIBRARIES(RandomTestStatic ${STATIC_TARGET})

SET_TARGET_PROPERTIES(RandomTest PROPERTIES
	PROJECT_LABEL "Test - RandomTest")

SET_TARGET_PROPERTIES(RandomTestStatic PROPERTIES
	PROJECT_LABEL "Test - RandomTestStatic")
