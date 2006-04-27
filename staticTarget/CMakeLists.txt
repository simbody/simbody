## This whole directory exists just so I could define this extra preprocessor value.

ADD_DEFINITIONS(-DSimTK_SIMBODY_TYPE=static)
ADD_LIBRARY(${STATIC_TARGET} STATIC ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_ABS_INCLUDE_FILES})
TARGET_LINK_LIBRARIES(${STATIC_TARGET} 
    ${SIMMATRIX_STATIC_LIBRARY} 
    ${SimTKCOMMON_STATIC_LIBRARY} 
    ${MATH_LIBS_TO_USE})

# install the library; on Windows both .lib and .dll go in the lib directory.
INSTALL_TARGETS(/core/lib${LIB64} RUNTIME_DIRECTORY /core/lib${LIB64} ${STATIC_TARGET})