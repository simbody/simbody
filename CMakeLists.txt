PROJECT (Simbody)
SUBDIRS (Geometry Mechanics MultibodyComputation tests)

## Pick up any top-level stuff
FILE(GLOB SOURCE_FILES src/*.cpp )
FILE(GLOB INCLUDE_FILES include/simbody/*.h include/simbody/internal/*.h src/*.h)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include)
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
    LINK_DIRECTORIES(${PROJECT_SOURCE_DIR}/bin)

    IF(NOT SIMTK_SDK)
        SET(SIMTK_SDK "c:/Program Files/SimTK/core")
    ENDIF(NOT SIMTK_SDK)

ENDIF(${CMAKE_C_COMPILER} MATCHES "gcc")

SET(CMAKE_INSTALL_PREFIX ${SIMTK_SDK})

LINK_DIRECTORIES(${SIMTK_SDK}/lib/simtk)
LINK_DIRECTORIES(${SIMTK_SDK}/lib/simmatrix)
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

ADD_TEST(testClass             ${TEST_PATH}/MultibodyComputation/CDSlib/testClass)
ADD_TEST(IVMMoleculeTest       ${TEST_PATH}/MultibodyComputation/IVMMoleculeModeling/IVMMoleculeTest)
ADD_TEST(OrientationTest       ${TEST_PATH}/Geometry/tests/OrientationTest/OrientationTest)
