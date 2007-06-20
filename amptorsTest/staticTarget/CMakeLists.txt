## This whole directory exists just so I could define this extra preprocessor value.
ADD_DEFINITIONS(-DSimTK_SIMBODY_BUILDING_STATIC_LIBRARY -DSimTK_USE_STATIC_LIBRARIES)

ADD_LIBRARY(${STATIC_TARGET} STATIC ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_ABS_INCLUDE_FILES})
TARGET_LINK_LIBRARIES(${STATIC_TARGET} 
    debug ${SimTKCPODES_STATIC_LIBRARY}_d optimized ${SimTKCPODES_STATIC_LIBRARY}
    debug ${SimTKCOMMON_STATIC_LIBRARY}_d optimized ${SimTKCOMMON_STATIC_LIBRARY}
    ${MATH_LIBS_TO_USE})

# install the library; on Windows both .lib and .dll go in the lib directory.
INSTALL_TARGETS(/core/lib${LIB64} RUNTIME_DIRECTORY /core/lib${LIB64} ${STATIC_TARGET})