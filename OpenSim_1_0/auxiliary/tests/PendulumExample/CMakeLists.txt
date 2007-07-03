
ADD_EXECUTABLE(PendulumExample       PendulumExample.cpp)
ADD_EXECUTABLE(PendulumExampleStatic PendulumExampleStatic.cpp)

TARGET_LINK_LIBRARIES(PendulumExample       ${SHARED_TARGET} ${AUX_SHARED_TARGET} )
TARGET_LINK_LIBRARIES(PendulumExampleStatic ${STATIC_TARGET} ${AUX_STATIC_TARGET})
