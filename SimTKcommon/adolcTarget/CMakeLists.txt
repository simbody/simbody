## This whole directory exists just so I could define this extra
## preprocessor value.

add_definitions(-DSimTK_SimTKCOMMON_BUILDING_SHARED_LIBRARY)

if(BUILD_UNVERSIONED_LIBRARIES)
    add_library(${ADOLC_TARGET} SHARED
            ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES}
            ${API_ABS_INCLUDE_FILES})

    add_dependencies(${ADOLC_TARGET} PlatformFiles)

    # TODO
    target_include_directories(${ADOLC_TARGET} PRIVATE "${ADOLC_INCLUDES}")

    target_link_libraries(${ADOLC_TARGET} ${MATH_LIBS_TO_USE}
            ${ADOLC_LIBRARIES})

    set_target_properties(${ADOLC_TARGET} PROPERTIES
            PROJECT_LABEL "Code - ${ADOLC_TARGET}"
            SOVERSION ${SIMBODY_SONAME_VERSION})

    # "PUBLIC" means this definition will be passed to the compiler when
    # compiling any target that links to ADOLC_TARGET.
    target_compile_definitions(${ADOLC_TARGET} PUBLIC SimTK_REAL_IS_ADOUBLE)
    
    # Warning 4996 is not propagated to downstream targets  
    if(WIN32 AND NOT ${CMAKE_VERSION} VERSION_LESS 3.0.2) 
        target_compile_options(${ADOLC_TARGET} INTERFACE 
            "$<$<CONFIG:Debug>:/wd4996>")
    endif()

    # install library; on Windows .dll's go in bin directory

    install(TARGETS ${ADOLC_TARGET} EXPORT SimbodyTargets
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
            GROUP_READ GROUP_WRITE GROUP_EXECUTE
            WORLD_READ WORLD_EXECUTE
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} )
endif(BUILD_UNVERSIONED_LIBRARIES)

if(BUILD_VERSIONED_LIBRARIES)
    add_library(${ADOLC_TARGET_VN} SHARED
            ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES}
            ${API_ABS_INCLUDE_FILES})

    add_dependencies(${ADOLC_TARGET_VN} PlatformFiles)

    target_include_directories(${ADOLC_TARGET_VN} "${ADOLC_INCLUDES}")

    target_link_libraries(${ADOLC_TARGET_VN} ${MATH_LIBS_TO_USE_VN}
            ${ADOLC_LIBRARIES})

    # "PUBLIC" means this definition will be passed to the compiler when
    # compiling any target that links to ADOLC_TARGET.
    target_compile_definitions(${ADOLC_TARGET_VN} PUBLIC SimTK_REAL_IS_ADOUBLE)

    set_target_properties(${ADOLC_TARGET_VN} PROPERTIES
            PROJECT_LABEL "Code - ${ADOLC_TARGET_VN}"
            SOVERSION ${SIMBODY_SONAME_VERSION})
            
    # Warning 4996 is not propagated to downstream targets
    if(WIN32 AND NOT ${CMAKE_VERSION} VERSION_LESS 3.0.2) 
        target_compile_options(${ADOLC_TARGET_VN} INTERFACE 
            "$<$<CONFIG:Debug>:/wd4996>")
    endif()

    # install library; on Windows .dll's go in bin directory

    install(TARGETS ${ADOLC_TARGET_VN} EXPORT SimbodyTargets
            PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
            GROUP_READ GROUP_WRITE GROUP_EXECUTE
            WORLD_READ WORLD_EXECUTE
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} )
endif(BUILD_VERSIONED_LIBRARIES)

