## This whole directory exists just so I could define this extra preprocessor value.
ADD_DEFINITIONS(-DSimTK_SIMBODY_BUILDING_SHARED_LIBRARY)

ADD_LIBRARY(${SHARED_TARGET} SHARED ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_ABS_INCLUDE_FILES})

IF( CMAKE_BUILD_TYPE MATCHES Release )
    TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${SimTKMATH_SHARED_LIBRARY} 
                                           ${SimTKCPODES_SHARED_LIBRARY}  
                                           ${SimTKCOMMON_SHARED_LIBRARY}  ${MATH_LIBS_TO_USE})
ELSE( CMAKE_BUILD_TYPE MATCHES Release )
    TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${SimTKMATH_SHARED_LIBRARY}_d 
                                           ${SimTKCPODES_SHARED_LIBRARY}_d  
                                           ${SimTKCOMMON_SHARED_LIBRARY}_d  ${MATH_LIBS_TO_USE})
ENDIF( CMAKE_BUILD_TYPE MATCHES Release )


# install the library; on Windows both .lib and .dll go in the lib directory.
INSTALL_TARGETS(/core/lib${LIB64} RUNTIME_DIRECTORY /core/lib${LIB64} ${SHARED_TARGET})



