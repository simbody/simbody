ADD_EXECUTABLE(TypesafeEnumTest       TypesafeEnumTest.cpp)
ADD_EXECUTABLE(TypesafeEnumTestStatic TypesafeEnumTestStatic.cpp)

TARGET_LINK_LIBRARIES(TypesafeEnumTest       ${SHARED_TARGET})
TARGET_LINK_LIBRARIES(TypesafeEnumTestStatic ${STATIC_TARGET})
