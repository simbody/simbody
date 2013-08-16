## This whole directory exists just so I could define these extra 
## preprocessor values.

ADD_DEFINITIONS(-DSimTK_SIMMATH_BUILDING_STATIC_LIBRARY 
		        -DSimTK_USE_STATIC_LIBRARIES)

#
# Set up file groups for better browsing in Visual Studio.
#

# Override Api headers to put them in the Api group.
source_group("Header Files\\Api" FILES ${API_ABS_INCLUDE_FILES})

# Now put things in their specified groups.
if (SOURCE_GROUPS)
    list(LENGTH SOURCE_GROUPS NGROUPS)
    math(EXPR lastgrpnum ${NGROUPS}-1)
    foreach( grpnum RANGE 0 ${lastgrpnum} )
        list(GET SOURCE_GROUPS ${grpnum} grp)
        list(GET SOURCE_GROUP_FILES ${grpnum} grpfile)
        source_group("${grp}" FILES "${grpfile}")
    endforeach()
endif()

# These are defaults for files that don't otherwise get mentioned.
source_group("Source Files\\Others" FILES ${SOURCE_FILES})
source_group("Header Files\\Others" FILES ${SOURCE_INCLUDE_FILES})

IF(BUILD_UNVERSIONED_LIBRARIES)

    ADD_LIBRARY(${STATIC_TARGET} STATIC 
			    ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} 
				${API_ABS_INCLUDE_FILES})
    
    TARGET_LINK_LIBRARIES(${STATIC_TARGET} 
                 debug ${SimTKCOMMON_STATIC_LIBRARY}_d 
				 optimized ${SimTKCOMMON_STATIC_LIBRARY}
                 ${MATH_LIBS_TO_USE})
    
    SET_TARGET_PROPERTIES(${STATIC_TARGET} PROPERTIES
    	PROJECT_LABEL "Code - ${STATIC_TARGET}")
    
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
                 debug ${SimTKCOMMON_STATIC_LIBRARY_VN}_d 
				 optimized ${SimTKCOMMON_STATIC_LIBRARY_VN}
                 ${MATH_LIBS_TO_USE_VN})
    
    SET_TARGET_PROPERTIES(${STATIC_TARGET_VN} PROPERTIES
    	PROJECT_LABEL "Code - ${STATIC_TARGET_VN}")
    
    # install library; on Windows both .lib and .dll go in the lib directory.
    INSTALL(TARGETS ${STATIC_TARGET_VN}
                     PERMISSIONS OWNER_READ OWNER_WRITE 
					 GROUP_READ GROUP_WRITE 
					 WORLD_READ
                     ARCHIVE DESTINATION ${SimTK_INSTALL_PREFIX}/lib${LIB64}
                     LIBRARY DESTINATION ${SimTK_INSTALL_PREFIX}/lib${LIB64}
                     RUNTIME DESTINATION ${SimTK_INSTALL_PREFIX}/bin)

ENDIF(BUILD_VERSIONED_LIBRARIES)
