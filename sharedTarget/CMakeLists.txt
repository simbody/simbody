## This whole directory exists just so I could define this extra preprocessor value.
ADD_DEFINITIONS(-DSimTK_SIMBODY_BUILDING_SHARED_LIBRARY)

ADD_LIBRARY(${SHARED_TARGET} SHARED ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_ABS_INCLUDE_FILES})
TARGET_LINK_LIBRARIES(${SHARED_TARGET} 
		debug ${SimTKCPODES_SHARED_LIBRARY}_d optimized ${SimTKCPODES_SHARED_LIBRARY}
		debug ${SIMMATRIX_SHARED_LIBRARY}_d   optimized ${SIMMATRIX_SHARED_LIBRARY}
		debug ${SimTKCOMMON_SHARED_LIBRARY}_d optimized ${SimTKCOMMON_SHARED_LIBRARY}
		${MATH_LIBS_TO_USE})

# install the library; on Windows both .lib and .dll go in the lib directory.
INSTALL_TARGETS(/core/lib${LIB64} RUNTIME_DIRECTORY /core/lib${LIB64} ${SHARED_TARGET})



