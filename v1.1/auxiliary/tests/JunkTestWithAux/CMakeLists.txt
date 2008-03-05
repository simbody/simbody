
ADD_EXECUTABLE(JunkTestWithAux       JunkTestWithAux.cpp)
ADD_EXECUTABLE(JunkTestWithAuxStatic JunkTestWithAuxStatic.cpp)

TARGET_LINK_LIBRARIES(JunkTestWithAux       ${SHARED_TARGET} ${AUX_SHARED_TARGET})
TARGET_LINK_LIBRARIES(JunkTestWithAuxStatic ${STATIC_TARGET} ${AUX_STATIC_TARGET})