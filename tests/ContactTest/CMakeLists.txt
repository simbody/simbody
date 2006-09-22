
ADD_EXECUTABLE(ContactTest       ContactTest.cpp)
ADD_EXECUTABLE(ContactTestStatic ContactTestStatic.cpp)

TARGET_LINK_LIBRARIES(ContactTest       ${SHARED_TARGET})
TARGET_LINK_LIBRARIES(ContactTestStatic ${STATIC_TARGET})