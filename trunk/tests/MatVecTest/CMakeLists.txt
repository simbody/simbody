ADD_EXECUTABLE(MatVecTest       MatVecTest.cpp)
ADD_EXECUTABLE(MatVecTestStatic MatVecTestStatic.cpp)

TARGET_LINK_LIBRARIES(MatVecTest       ${SHARED_TARGET})
TARGET_LINK_LIBRARIES(MatVecTestStatic ${STATIC_TARGET})