## This whole directory exists just so I could define this extra preprocessor value.

ADD_DEFINITIONS(-DSimTK_SimTKCOMMON_BUILDING_SHARED_LIBRARY)
ADD_LIBRARY(${SHARED_TARGET} SHARED ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_ABS_INCLUDE_FILES})

# install the library; on Windows both .lib and .dll go in the lib directory.
INSTALL_TARGETS(/core/lib${LIB64} RUNTIME_DIRECTORY /core/lib${LIB64} ${SHARED_TARGET})



