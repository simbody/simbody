## This whole directory exists just so I could define this extra preprocessor value.

ADD_DEFINITIONS(-DSIMTK_SIMBODY_TYPE=shared)
ADD_LIBRARY(${SHARED_TARGET} SHARED ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_ABS_INCLUDE_FILES})
TARGET_LINK_LIBRARIES(${SHARED_TARGET} simmatrix simtkcommon ${MATH_LIBS_TO_USE})

# install the library; on Windows both .lib and .dll go in the lib directory.
INSTALL_TARGETS(/core/lib RUNTIME_DIRECTORY /core/lib ${SHARED_TARGET})



