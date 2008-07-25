ADD_DEFINITIONS(-DSimTK_SIMBODY_BUILDING_STATIC_LIBRARY -DSimTK_USE_STATIC_LIBRARIES)

ADD_LIBRARY(${AUX_STATIC_TARGET} STATIC ${AUX_SOURCE_FILES} ${AUX_INCLUDE_FILES} ${API_ABS_INCLUDE_FILES})


IF( CMAKE_BUILD_TYPE MATCHES Release )
    TARGET_LINK_LIBRARIES(${AUX_STATIC_TARGET} ${STATIC_TARGET}                        
                                           ${SimTKCPODES_STATIC_LIBRARY}
                                           ${SimTKCOMMON_STATIC_LIBRARY}  
                                           ${MATH_LIBS_TO_USE})
ELSE( CMAKE_BUILD_TYPE MATCHES Release )
    TARGET_LINK_LIBRARIES(${AUX_STATIC_TARGET} ${STATIC_TARGET}                          
                                           ${SimTKCPODES_STATIC_LIBRARY}_d
                                           ${SimTKCOMMON_STATIC_LIBRARY}_d  
                                           ${MATH_LIBS_TO_USE})
ENDIF( CMAKE_BUILD_TYPE MATCHES Release )


# install the library; on Windows both .lib and .dll go in the lib directory.
INSTALL(TARGETS ${AUX_STATIC_TARGET}
                 ARCHIVE DESTINATION ${SimTK_INSTALL_PREFIX}/core/lib${LIB64}
                 LIBRARY DESTINATION ${SimTK_INSTALL_PREFIX}/core/lib${LIB64}
                 RUNTIME DESTINATION ${SimTK_INSTALL_PREFIX}/core/lib${LIB64}
)


