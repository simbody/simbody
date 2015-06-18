## This whole directory exists just so I could define these extra 
## preprocessor values.

add_definitions(-DSimTK_SIMMATH_BUILDING_STATIC_LIBRARY 
                -DSimTK_USE_STATIC_LIBRARIES)

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

    add_library(${STATIC_TARGET} STATIC 
                ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} 
                ${API_ABS_INCLUDE_FILES})
    
    target_link_libraries(${STATIC_TARGET} 
                 ${SimTKCOMMON_STATIC_LIBRARY}
                 ${MATH_LIBS_TO_USE})
    
    set_target_properties(${STATIC_TARGET} PROPERTIES
        PROJECT_LABEL "Code - ${STATIC_TARGET}")
    
    # install library; on Windows both .lib and .dll go in the lib directory.
    install(TARGETS ${STATIC_TARGET} EXPORT SimbodyTargets
                     PERMISSIONS OWNER_READ OWNER_WRITE 
                     GROUP_READ GROUP_WRITE 
                     WORLD_READ
                     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
                     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

endif(BUILD_UNVERSIONED_LIBRARIES)


if(BUILD_VERSIONED_LIBRARIES)

    add_library(${STATIC_TARGET_VN} STATIC 
                ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} 
                ${API_ABS_INCLUDE_FILES})
    
    target_link_libraries(${STATIC_TARGET_VN} 
                 ${SimTKCOMMON_STATIC_LIBRARY_VN}
                 ${MATH_LIBS_TO_USE_VN})
    
    set_target_properties(${STATIC_TARGET_VN} PROPERTIES
        PROJECT_LABEL "Code - ${STATIC_TARGET_VN}")
    
    # install library; on Windows both .lib and .dll go in the lib directory.
    install(TARGETS ${STATIC_TARGET_VN} EXPORT SimbodyTargets
                     PERMISSIONS OWNER_READ OWNER_WRITE 
                     GROUP_READ GROUP_WRITE 
                     WORLD_READ
                     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
                     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

endif(BUILD_VERSIONED_LIBRARIES)
