## This whole directory exists just so I could define this extra 
## preprocessor values.

add_definitions(-DSimTK_SimTKCOMMON_BUILDING_STATIC_LIBRARY 
                -DSimTK_USE_STATIC_LIBRARIES)

if(BUILD_UNVERSIONED_LIBRARIES)
    add_library(${STATIC_TARGET} STATIC 
        ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} 
        ${API_ABS_INCLUDE_FILES})

    add_dependencies(${STATIC_TARGET} PlatformFiles)

    target_link_libraries(${STATIC_TARGET} ${MATH_LIBS_TO_USE})

    set_target_properties(${STATIC_TARGET} PROPERTIES
        PROJECT_LABEL "Code - ${STATIC_TARGET}")
        
    # Warning 4996 is not propagated to downstream targets
    if(WIN32 AND NOT ${CMAKE_VERSION} VERSION_LESS 3.0.2) 
        target_compile_options(${STATIC_TARGET} INTERFACE 
            "$<$<CONFIG:Debug>:/wd4996>")
    endif()

    install(TARGETS ${STATIC_TARGET} EXPORT SimbodyTargets
        PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ GROUP_WRITE WORLD_READ 
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} )

    if(NOT ${CMAKE_VERSION} VERSION_LESS 3.0.2)
        # This removes the need for client projects to use
        # `include_directories(${Simbody_INCLUDE_DIR})`. However, this ability
        # was only introduced in CMake 3.0.
        target_include_directories(${STATIC_TARGET} INTERFACE
            $<INSTALL_INTERFACE:${SIMBODY_INCLUDE_INSTALL_DIR}>
            )
    endif()
endif(BUILD_UNVERSIONED_LIBRARIES)

if(BUILD_VERSIONED_LIBRARIES)
    add_library(${STATIC_TARGET_VN} STATIC 
        ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} 
        ${API_ABS_INCLUDE_FILES})

    add_dependencies(${STATIC_TARGET_VN} PlatformFiles)

    target_link_libraries(${STATIC_TARGET_VN} ${MATH_LIBS_TO_USE_VN})

    set_target_properties(${STATIC_TARGET_VN} PROPERTIES
        PROJECT_LABEL "Code - ${STATIC_TARGET_VN}")
    
    # Warning 4996 is not propagated to downstream targets  
    if(WIN32 AND NOT ${CMAKE_VERSION} VERSION_LESS 3.0.2) 
        target_compile_options(${STATIC_TARGET_VN} INTERFACE 
            "$<$<CONFIG:Debug>:/wd4996>")
    endif()

    install(TARGETS ${STATIC_TARGET_VN} EXPORT SimbodyTargets
        PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ GROUP_WRITE WORLD_READ 
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} )

    if(NOT ${CMAKE_VERSION} VERSION_LESS 3.0.2)
        target_include_directories(${STATIC_TARGET_VN} INTERFACE
            $<INSTALL_INTERFACE:${SIMBODY_INCLUDE_INSTALL_DIR}>
            )
    endif()
endif(BUILD_VERSIONED_LIBRARIES)
