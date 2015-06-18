## This whole directory exists just so I could define this extra preprocessor value.

add_definitions(-DSimTK_SIMMATH_BUILDING_SHARED_LIBRARY)

#
# Set up file groups for better browsing in Visual Studio.
#

# Override Api headers to put them in the Api group.
source_group("Header Files\\Api" FILES ${API_ABS_INCLUDE_FILES})

# Now put things in their specified groups.
if(SOURCE_GROUPS)
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

if(BUILD_UNVERSIONED_LIBRARIES)

    add_library(${SHARED_TARGET} SHARED 
                ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} 
                ${API_ABS_INCLUDE_FILES})
    
    target_link_libraries(${SHARED_TARGET} 
                 ${SimTKCOMMON_SHARED_LIBRARY}
                 ${MATH_LIBS_TO_USE})
    
    set_target_properties(${SHARED_TARGET} PROPERTIES
        PROJECT_LABEL "Code - ${SHARED_TARGET}"
        SOVERSION ${SIMBODY_SONAME_VERSION})
    
    # install library; on Windows both .lib and .dll go in the lib directory.
    install(TARGETS ${SHARED_TARGET} EXPORT SimbodyTargets
                     PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE 
                     GROUP_READ GROUP_WRITE GROUP_EXECUTE 
                     WORLD_READ WORLD_EXECUTE
                     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
                     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

endif(BUILD_UNVERSIONED_LIBRARIES)


if(BUILD_VERSIONED_LIBRARIES)

    add_library(${SHARED_TARGET_VN} SHARED 
                ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} 
                ${API_ABS_INCLUDE_FILES})
    
    target_link_libraries(${SHARED_TARGET_VN} 
                 ${SimTKCOMMON_SHARED_LIBRARY_VN}
                 ${MATH_LIBS_TO_USE_VN})
    
    set_target_properties(${SHARED_TARGET_VN} PROPERTIES
        PROJECT_LABEL "Code - ${SHARED_TARGET_VN}" 
        SOVERSION ${SIMBODY_SONAME_VERSION})
    
    # install library; on Windows both .lib and .dll go in the lib directory.
    install(TARGETS ${SHARED_TARGET_VN} EXPORT SimbodyTargets
                     PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE 
                     GROUP_READ GROUP_WRITE GROUP_EXECUTE 
                     WORLD_READ WORLD_EXECUTE
                     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
                     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

endif(BUILD_VERSIONED_LIBRARIES)
