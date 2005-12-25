PROJECT (Simbody)
SUBDIRS (Geometry Mechanics MultibodyComputation MultibodyModeling)

IF(${CMAKE_C_COMPILER} MATCHES "gcc")
    ## A gcc build
    SET(MATH_LIBS_TO_USE acml c m g2c)
    LINK_DIRECTORIES(/usr/local/lib /lib)
    SET(SIMTK_FRAMEWORK_INSTALL_DIR /usr/local/SimTK/framework)
ELSE(${CMAKE_C_COMPILER} MATCHES "gcc")
    ## Assume Microsoft Visual Studio
    SET(MATH_LIBS_TO_USE libacml_dll)
    LINK_DIRECTORIES(${PROJECT_SOURCE_DIR}/bin)
    SET(SIMTK_FRAMEWORK_INSTALL_DIR "c:/Program Files/SimTK/framework")
ENDIF(${CMAKE_C_COMPILER} MATCHES "gcc")

LINK_DIRECTORIES(${SIMTK_FRAMEWORK_INSTALL_DIR}/lib/simtk)
LINK_DIRECTORIES(${SIMTK_FRAMEWORK_INSTALL_DIR}/lib/simmatrix)

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
ADD_TEST(MultibodyModelingTest ${TEST_PATH}/MultibodyModeling/MultibodyModelingTest)
