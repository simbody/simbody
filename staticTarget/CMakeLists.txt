## This whole directory exists just so I could define these extra 
## preprocessor values.

ADD_DEFINITIONS(-DSimTK_SIMBODY_BUILDING_STATIC_LIBRARY 
		        -DSimTK_USE_STATIC_LIBRARIES)

IF(BUILD_UNVERSIONED_LIBRARIES)

    ADD_LIBRARY(${STATIC_TARGET} STATIC 
			    ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} 
				${API_ABS_INCLUDE_FILES})
    
    TARGET_LINK_LIBRARIES(${STATIC_TARGET} 
                 debug ${SimTKMATH_STATIC_LIBRARY}_d 
				 optimized ${SimTKMATH_STATIC_LIBRARY}
                 debug ${SimTKCOMMON_STATIC_LIBRARY}_d 
				 optimized ${SimTKCOMMON_STATIC_LIBRARY}
                 ${MATH_LIBS_TO_USE})
    
    SET_TARGET_PROPERTIES(${STATIC_TARGET} PROPERTIES
    	PROJECT_LABEL "Library - ${STATIC_TARGET}")
    
    # install library; on Windows both .lib and .dll go in the lib directory.
    
    INSTALL(TARGETS ${STATIC_TARGET}  
                     PERMISSIONS OWNER_READ OWNER_WRITE 
					 GROUP_READ GROUP_WRITE 
					 WORLD_READ 
                     ARCHIVE DESTINATION ${SimTK_INSTALL_PREFIX}/lib${LIB64} 
                     LIBRARY DESTINATION ${SimTK_INSTALL_PREFIX}/lib${LIB64} 
                     RUNTIME DESTINATION ${SimTK_INSTALL_PREFIX}/bin)

ENDIF(BUILD_UNVERSIONED_LIBRARIES)



IF(BUILD_VERSIONED_LIBRARIES)

    ADD_LIBRARY(${STATIC_TARGET_VN} STATIC 
			    ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} 
				${API_ABS_INCLUDE_FILES})
    
    TARGET_LINK_LIBRARIES(${STATIC_TARGET_VN} 
                 debug ${SimTKMATH_STATIC_LIBRARY_VN}_d 
				 optimized ${SimTKMATH_STATIC_LIBRARY_VN}
                 debug ${SimTKCOMMON_STATIC_LIBRARY_VN}_d 
				 optimized ${SimTKCOMMON_STATIC_LIBRARY_VN}
                 ${MATH_LIBS_TO_USE_VN})
    
    SET_TARGET_PROPERTIES(${STATIC_TARGET_VN} PROPERTIES
    	PROJECT_LABEL "Library - ${STATIC_TARGET_VN}")
    
    # install library; on Windows both .lib and .dll go in the lib directory.
    
    INSTALL(TARGETS ${STATIC_TARGET_VN}  
                     PERMISSIONS OWNER_READ OWNER_WRITE 
					 GROUP_READ GROUP_WRITE 
					 WORLD_READ 
                     ARCHIVE DESTINATION ${SimTK_INSTALL_PREFIX}/lib${LIB64} 
                     LIBRARY DESTINATION ${SimTK_INSTALL_PREFIX}/lib${LIB64} 
                     RUNTIME DESTINATION ${SimTK_INSTALL_PREFIX}/bin)

ENDIF(BUILD_VERSIONED_LIBRARIES)
