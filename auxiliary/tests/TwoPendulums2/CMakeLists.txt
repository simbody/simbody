
ADD_EXECUTABLE(TwoPendulums2       TwoPendulums2.cpp)
ADD_EXECUTABLE(TwoPendulums2Static TwoPendulums2Static.cpp)

TARGET_LINK_LIBRARIES(TwoPendulums2       ${SHARED_TARGET} ${AUX_SHARED_TARGET})
TARGET_LINK_LIBRARIES(TwoPendulums2Static ${STATIC_TARGET} ${AUX_STATIC_TARGET})
