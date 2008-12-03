
ADD_EXECUTABLE(TwoPendulums       TwoPendulums.cpp)
ADD_EXECUTABLE(TwoPendulumsStatic TwoPendulumsStatic.cpp)

TARGET_LINK_LIBRARIES(TwoPendulums       ${SHARED_TARGET} ${AUX_SHARED_TARGET})
TARGET_LINK_LIBRARIES(TwoPendulumsStatic ${STATIC_TARGET} ${AUX_STATIC_TARGET})


SET_TARGET_PROPERTIES(TwoPendulums PROPERTIES
	PROJECT_LABEL "Adhoc - TwoPendulums")
SET_TARGET_PROPERTIES(TwoPendulumsStatic PROPERTIES
	PROJECT_LABEL "Adhoc - TwoPendulumsStatic")
