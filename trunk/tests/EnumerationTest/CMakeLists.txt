ADD_EXECUTABLE(EnumerationTest       EnumerationTest.cpp)
ADD_EXECUTABLE(EnumerationTestStatic EnumerationTestStatic.cpp)

TARGET_LINK_LIBRARIES(EnumerationTest       ${SHARED_TARGET})
TARGET_LINK_LIBRARIES(EnumerationTestStatic ${STATIC_TARGET})
