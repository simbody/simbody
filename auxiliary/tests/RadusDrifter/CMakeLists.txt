
ADD_EXECUTABLE(RadusDrifter       RadusDrifter.cpp)
ADD_EXECUTABLE(RadusDrifterStatic RadusDrifterStatic.cpp)

TARGET_LINK_LIBRARIES(RadusDrifter       ${SHARED_TARGET} ${AUX_SHARED_TARGET})
TARGET_LINK_LIBRARIES(RadusDrifterStatic ${STATIC_TARGET} ${AUX_STATIC_TARGET})

SET_TARGET_PROPERTIES(RadusDrifter PROPERTIES
	PROJECT_LABEL "Adhoc - RadusDrifter")
SET_TARGET_PROPERTIES(RadusDrifterStatic PROPERTIES
	PROJECT_LABEL "Adhoc - RadusDrifterStatic")
