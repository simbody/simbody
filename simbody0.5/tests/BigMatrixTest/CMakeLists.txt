ADD_EXECUTABLE(BigMatrixTest       BigMatrixTest.cpp)
ADD_EXECUTABLE(BigMatrixTestStatic BigMatrixTest.cpp)

TARGET_LINK_LIBRARIES(BigMatrixTest       ${SHARED_TARGET})
TARGET_LINK_LIBRARIES(BigMatrixTestStatic ${STATIC_TARGET})