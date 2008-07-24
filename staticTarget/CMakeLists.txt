## This whole directory exists just so I could define this extra preprocessor value.
ADD_DEFINITIONS(-DSimTK_SIMBODY_BUILDING_STATIC_LIBRARY -DSimTK_USE_STATIC_LIBRARIES)

ADD_LIBRARY(${STATIC_TARGET} STATIC ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_ABS_INCLUDE_FILES})

IF( CMAKE_BUILD_TYPE MATCHES Release )
    TARGET_LINK_LIBRARIES(${STATIC_TARGET} ${SimTKMATH_STATIC_LIBRARY} 
                                           ${SimTKCPODES_STATIC_LIBRARY}  
                                           ${SimTKCOMMON_STATIC_LIBRARY} ${MATH_LIBS_TO_USE})
ELSE( CMAKE_BUILD_TYPE MATCHES Release )
    TARGET_LINK_LIBRARIES(${STATIC_TARGET} ${SimTKMATH_STATIC_LIBRARY}_d 
                                           ${SimTKCPODES_STATIC_LIBRARY}_d  
                                           ${SimTKCOMMON_STATIC_LIBRARY}_d ${MATH_LIBS_TO_USE})
ENDIF( CMAKE_BUILD_TYPE MATCHES Release )


# install the library; on Windows both .lib and .dll go in the lib directory.
INSTALL_TARGETS(/core/lib${LIB64} RUNTIME_DIRECTORY /core/lib${LIB64} ${STATIC_TARGET})
