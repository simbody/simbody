ADD_EXECUTABLE(RandomTest       RandomTest.cpp)
ADD_EXECUTABLE(RandomTestStatic RandomTestStatic.cpp)

TARGET_LINK_LIBRARIES(RandomTest       ${SHARED_TARGET})
TARGET_LINK_LIBRARIES(RandomTestStatic ${STATIC_TARGET})
