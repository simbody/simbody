## This whole directory exists just so I could define this extra 
## preprocessor value.

add_definitions(-DSimTK_SIMBODY_BUILDING_SHARED_LIBRARY)

if(BUILD_UNVERSIONED_LIBRARIES)

    add_library(${SHARED_TARGET} SHARED 
                ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} 
                ${API_ABS_INCLUDE_FILES})
    
    target_link_libraries(${SHARED_TARGET} 
                 ${SimTKMATH_SHARED_LIBRARY}
                 ${SimTKCOMMON_SHARED_LIBRARY}
                 ${MATH_LIBS_TO_USE})
    
    set_target_properties(${SHARED_TARGET} PROPERTIES
        PROJECT_LABEL "Code - ${SHARED_TARGET}"
        SOVERSION ${SIMBODY_SONAME_VERSION})
    
    # install library; on Windows the .dll goes in the bin directory.
    
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
                 ${SimTKMATH_SHARED_LIBRARY_VN}
                 ${SimTKCOMMON_SHARED_LIBRARY_VN}
                 ${MATH_LIBS_TO_USE_VN})
    
    set_target_properties(${SHARED_TARGET_VN} PROPERTIES
        PROJECT_LABEL "Code - ${SHARED_TARGET_VN}"
        SOVERSION ${SIMBODY_SONAME_VERSION})
    
    # install library; on Windows the .dll goes in the bin directory.
    
    install(TARGETS ${SHARED_TARGET_VN} EXPORT SimbodyTargets
                     PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE 
                     GROUP_READ GROUP_WRITE GROUP_EXECUTE 
                     WORLD_READ WORLD_EXECUTE
                     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
                     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
        
endif(BUILD_VERSIONED_LIBRARIES)
