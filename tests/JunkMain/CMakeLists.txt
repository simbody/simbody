
ADD_EXECUTABLE(JunkMain       JunkMain.cpp)
ADD_EXECUTABLE(JunkMainStatic JunkMainStatic.cpp)

TARGET_LINK_LIBRARIES(JunkMain       ${SHARED_TARGET})
TARGET_LINK_LIBRARIES(JunkMainStatic ${STATIC_TARGET})