## This whole directory exists just so I could define this extra 
## preprocessor value.

add_definitions(-DSimTK_SimTKCOMMON_BUILDING_SHARED_LIBRARY)

if(BUILD_UNVERSIONED_LIBRARIES)
    add_library(${SHARED_TARGET} SHARED 
        ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} 
        ${API_ABS_INCLUDE_FILES})

    add_dependencies(${SHARED_TARGET} PlatformFiles)

    target_link_libraries(${SHARED_TARGET} ${MATH_LIBS_TO_USE})

    set_target_properties(${SHARED_TARGET} PROPERTIES
        PROJECT_LABEL "Code - ${SHARED_TARGET}"
        SOVERSION ${SIMBODY_SONAME_VERSION})
    
    # install library; on Windows .dll's go in bin directory

    install(TARGETS ${SHARED_TARGET} EXPORT SimbodyTargets
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE 
                    GROUP_READ GROUP_WRITE GROUP_EXECUTE 
                    WORLD_READ WORLD_EXECUTE
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} )

    if(NOT ${CMAKE_VERSION} VERSION_LESS 3.0.2)
        # This removes the need for client projects to use
        # `include_directories(${Simbody_INCLUDE_DIR})`. However, this ability
        # was only introduced in CMake 3.0.
        target_include_directories(${SHARED_TARGET} INTERFACE
            $<INSTALL_INTERFACE:${SIMBODY_INCLUDE_INSTALL_DIR}>
            )
    endif()

    if(SIMBODY_COVERAGE)
        add_coverage(${SHARED_TARGET})
    endif()

endif(BUILD_UNVERSIONED_LIBRARIES)

if(BUILD_VERSIONED_LIBRARIES)
    add_library(${SHARED_TARGET_VN} SHARED 
        ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} 
        ${API_ABS_INCLUDE_FILES})

    add_dependencies(${SHARED_TARGET_VN} PlatformFiles)

    target_link_libraries(${SHARED_TARGET_VN} ${MATH_LIBS_TO_USE_VN})

    set_target_properties(${SHARED_TARGET_VN} PROPERTIES
        PROJECT_LABEL "Code - ${SHARED_TARGET_VN}"
        SOVERSION ${SIMBODY_SONAME_VERSION})
           
    # install library; on Windows .dll's go in bin directory

    install(TARGETS ${SHARED_TARGET_VN} EXPORT SimbodyTargets
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE 
                    GROUP_READ GROUP_WRITE GROUP_EXECUTE 
                    WORLD_READ WORLD_EXECUTE
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} )

    if(NOT ${CMAKE_VERSION} VERSION_LESS 3.0.2)
        target_include_directories(${SHARED_TARGET_VN} INTERFACE
            $<INSTALL_INTERFACE:${SIMBODY_INCLUDE_INSTALL_DIR}>
            )
    endif()

    if(SIMBODY_COVERAGE)
        add_coverage(${SHARED_TARGET_VN})
    endif()

endif(BUILD_VERSIONED_LIBRARIES)
