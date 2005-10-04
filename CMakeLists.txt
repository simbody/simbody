PROJECT (Simbody)
SUBDIRS (CDSlib common RigidBodies Molecules)

IF(${CMAKE_C_COMPILER} MATCHES "gcc")
    ## A gcc build
    SET(MATH_LIBS_TO_USE acml c m g2c)
    LINK_DIRECTORIES(/usr/local/lib /lib)
ELSE(${CMAKE_C_COMPILER} MATCHES "gcc")
    ## Assume Microsoft Visual Studio
    SET(MATH_LIBS_TO_USE libacml_dll)
    LINK_DIRECTORIES(${PROJECT_SOURCE_DIR}/bin)
ENDIF(${CMAKE_C_COMPILER} MATCHES "gcc")



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

ADD_TEST(testClass ${TEST_PATH}/CDSlib/testClass)
ADD_TEST(MoleculeTest ${TEST_PATH}/Molecules/MoleculeTest)
