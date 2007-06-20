ADD_EXECUTABLE(IntegratorTest       IntegratorTest.cpp)
ADD_EXECUTABLE(IntegratorTestStatic IntegratorTestStatic.cpp)

TARGET_LINK_LIBRARIES(IntegratorTest       ${SHARED_TARGET})
TARGET_LINK_LIBRARIES(IntegratorTestStatic ${STATIC_TARGET})