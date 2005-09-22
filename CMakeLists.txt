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


