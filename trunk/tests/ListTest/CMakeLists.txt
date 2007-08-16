ADD_EXECUTABLE(ListTest       ListTest.cpp)
ADD_EXECUTABLE(ListTestStatic ListTestStatic.cpp)

TARGET_LINK_LIBRARIES(ListTest       ${SHARED_TARGET})
TARGET_LINK_LIBRARIES(ListTestStatic ${STATIC_TARGET})
