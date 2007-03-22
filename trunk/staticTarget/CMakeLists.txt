## This whole directory exists just so I could define this extra preprocessor value.

ADD_DEFINITIONS(-DSimTK_SIMMATH_BUILDING_STATIC_LIBRARY -DSimTK_USE_STATIC_LIBRARIES)
ADD_LIBRARY(${STATIC_TARGET} STATIC ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_ABS_INCLUDE_FILES})

TARGET_LINK_LIBRARIES(${STATIC_TARGET} 
    debug ${SimTKCOMMON_STATIC_LIBRARY}_d optimized ${SimTKCOMMON_STATIC_LIBRARY}
    ${MATH_LIBS_TO_USE})

# install the library; on Windows both .lib and .dll go in the lib directory.
INSTALL(TARGETS ${STATIC_TARGET}
                 PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ GROUP_WRITE WORLD_READ
                 ARCHIVE DESTINATION ${SimTK_INSTALL_PREFIX}/core/lib${LIB64}
                 LIBRARY DESTINATION ${SimTK_INSTALL_PREFIX}/core/lib${LIB64}
                 RUNTIME DESTINATION ${SimTK_INSTALL_PREFIX}/core/lib${LIB64}
)

