ADD_EXECUTABLE(EnumerationTest       EnumerationTest.cpp)
ADD_EXECUTABLE(EnumerationTestStatic EnumerationTestStatic.cpp)

TARGET_LINK_LIBRARIES(EnumerationTest       ${SHARED_TARGET})
TARGET_LINK_LIBRARIES(EnumerationTestStatic ${STATIC_TARGET})

SET_TARGET_PROPERTIES(EnumerationTest PROPERTIES
	PROJECT_LABEL "Test - EnumerationTest")

SET_TARGET_PROPERTIES(EnumerationTestStatic PROPERTIES
	PROJECT_LABEL "Test - EnumerationTestStatic")
