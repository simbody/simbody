
ADD_EXECUTABLE(RadusDrifter       RadusDrifter.cpp)
ADD_EXECUTABLE(RadusDrifterStatic RadusDrifterStatic.cpp)

TARGET_LINK_LIBRARIES(RadusDrifter       ${SHARED_TARGET} ${AUX_SHARED_TARGET})
TARGET_LINK_LIBRARIES(RadusDrifterStatic ${STATIC_TARGET} ${AUX_STATIC_TARGET})
