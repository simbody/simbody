PROJECT (Simbody)
SUBDIRS (Geometry Mechanics tests)

## Pick up any top-level stuff
FILE(GLOB SOURCE_FILES src/*.cpp )
FILE(GLOB INCLUDE_FILES include/simbody/*.h include/simbody/internal/*.h src/*.h)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include
                    ${PROJECT_SOURCE_DIR}/Geometry/include
                    ${PROJECT_SOURCE_DIR}/Mechanics/include)

ADD_LIBRARY(simbody ${SOURCE_FILES} ${INCLUDE_FILES})


#
# Find the right place to look for dependency resolution (which
# can be the default install directory or specified SIMTK_SDK directory).
#
# Also, etermine which math libraries to use for this platform.
# TODO: change ACML to simtklapack when it is ready.
#

IF(${CMAKE_C_COMPILER} MATCHES "gcc")

    ## A gcc build
    SET(MATH_LIBS_TO_USE acml c m g2c)
    LINK_DIRECTORIES(/usr/local/lib /lib)

    IF(NOT SIMTK_SDK)
        SET(SIMTK_SDK /usr/local/SimTK/core)
    ENDIF(NOT SIMTK_SDK)

ELSE(${CMAKE_C_COMPILER} MATCHES "gcc")

    ## Assume Microsoft Visual Studio
    SET(MATH_LIBS_TO_USE libacml_dll)
    LINK_DIRECTORIES("c:/Program Files/SimTK/core/lib")

    IF(NOT SIMTK_SDK)
        SET(SIMTK_SDK "c:/Program Files/SimTK/core")
    ENDIF(NOT SIMTK_SDK)

ENDIF(${CMAKE_C_COMPILER} MATCHES "gcc")

SET(CMAKE_INSTALL_PREFIX ${SIMTK_SDK})

LINK_DIRECTORIES(${SIMTK_SDK}/lib/simtk)
LINK_DIRECTORIES(${SIMTK_SDK}/lib/simmatrix)
LINK_DIRECTORIES(${SIMTK_SDK}/lib)
INCLUDE_DIRECTORIES(${SIMTK_SDK}/include)



#
# Testing
#

IF (UNIX)
  ADD_DEFINITIONS(-fprofile-arcs -ftest-coverage)
  LINK_LIBRARIES(gcov)
ENDIF (UNIX)

INCLUDE (Dart)

ENABLE_TESTING()

IF (EXECUTABLE_OUTPUT_PATH)
  SET (TEST_PATH ${EXECUTABLE_OUTPUT_PATH})
ELSE (EXECUTABLE_OUTPUT_PATH)
  SET (TEST_PATH .)
ENDIF (EXECUTABLE_OUTPUT_PATH)

ADD_TEST(OrientationTest       ${TEST_PATH}/Geometry/tests/OrientationTest/OrientationTest)
ADD_TEST(SBPendulum1           ${TEST_PATH}/tests/SBPendulum1/SBPendulum1)