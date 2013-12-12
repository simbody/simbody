## This whole directory exists just so I could define this extra preprocessor value.

ADD_DEFINITIONS(-DSimTK_SIMMATH_BUILDING_SHARED_LIBRARY)

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

    ADD_LIBRARY(${SHARED_TARGET} SHARED 
			    ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} 
				${API_ABS_INCLUDE_FILES})
    
    TARGET_LINK_LIBRARIES(${SHARED_TARGET} 
                 debug ${SimTKCOMMON_SHARED_LIBRARY}_d 
				 optimized ${SimTKCOMMON_SHARED_LIBRARY}
                 ${MATH_LIBS_TO_USE})
    
    SET_TARGET_PROPERTIES(${SHARED_TARGET} PROPERTIES
    	PROJECT_LABEL "Code - ${SHARED_TARGET}"
    	SOVERSION ${SIMBODY_SONAME_VERSION})
    
    # install library; on Windows both .lib and .dll go in the lib directory.
    INSTALL(TARGETS ${SHARED_TARGET}
                     PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE 
					 GROUP_READ GROUP_WRITE GROUP_EXECUTE 
					 WORLD_READ WORLD_EXECUTE
                     LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
                     ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
                     RUNTIME DESTINATION ${CMAKE_INSTALL_FULL_BINDIR})

ENDIF(BUILD_UNVERSIONED_LIBRARIES)


IF(BUILD_VERSIONED_LIBRARIES)

    ADD_LIBRARY(${SHARED_TARGET_VN} SHARED 
			    ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} 
				${API_ABS_INCLUDE_FILES})
    
    TARGET_LINK_LIBRARIES(${SHARED_TARGET_VN} 
                 debug ${SimTKCOMMON_SHARED_LIBRARY_VN}_d 
				 optimized ${SimTKCOMMON_SHARED_LIBRARY_VN}
                 ${MATH_LIBS_TO_USE_VN})
    
    SET_TARGET_PROPERTIES(${SHARED_TARGET_VN} PROPERTIES
    	PROJECT_LABEL "Code - ${SHARED_TARGET_VN}" 
        SOVERSION ${SIMBODY_SONAME_VERSION})
    
    # install library; on Windows both .lib and .dll go in the lib directory.
    INSTALL(TARGETS ${SHARED_TARGET_VN}
                     PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE 
					 GROUP_READ GROUP_WRITE GROUP_EXECUTE 
					 WORLD_READ WORLD_EXECUTE
                     LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
                     ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
                     RUNTIME DESTINATION ${CMAKE_INSTALL_FULL_BINDIR})

ENDIF(BUILD_VERSIONED_LIBRARIES)
