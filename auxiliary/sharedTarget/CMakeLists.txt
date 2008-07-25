## This whole directory exists just so I could define this extra preprocessor value.

ADD_DEFINITIONS(-DSimTK_SIMBODY_BUILDING_SHARED_LIBRARY)

ADD_LIBRARY(${AUX_SHARED_TARGET} SHARED ${AUX_SOURCE_FILES} ${AUX_INCLUDE_FILES} ${API_ABS_INCLUDE_FILES})


IF( CMAKE_BUILD_TYPE MATCHES Release )
    TARGET_LINK_LIBRARIES(${AUX_SHARED_TARGET} ${SHARED_TARGET}  
                                           ${SimTKCPODES_SHARED_LIBRARY}  
                                           ${SimTKCOMMON_SHARED_LIBRARY}  
                                           ${MATH_LIBS_TO_USE})
ELSE( CMAKE_BUILD_TYPE MATCHES Release )
    TARGET_LINK_LIBRARIES(${AUX_SHARED_TARGET} ${SHARED_TARGET} 
                                           ${SimTKCPODES_SHARED_LIBRARY}_d  
                                           ${SimTKCOMMON_SHARED_LIBRARY}_d  
                                           ${MATH_LIBS_TO_USE})
ENDIF( CMAKE_BUILD_TYPE MATCHES Release )


# install the library; on Windows both .lib and .dll go in the lib directory.
INSTALL(TARGETS ${AUX_SHARED_TARGET}
                 PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
                 ARCHIVE DESTINATION ${SimTK_INSTALL_PREFIX}/core/lib${LIB64}
                 LIBRARY DESTINATION ${SimTK_INSTALL_PREFIX}/core/lib${LIB64}
                 RUNTIME DESTINATION ${SimTK_INSTALL_PREFIX}/core/lib${LIB64}
)

