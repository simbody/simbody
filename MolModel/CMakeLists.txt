################################################################################################
### The following simtk_sdk configuration will be unnecessary when this is called from above ###
################################################################################################
IF(${CMAKE_C_COMPILER} MATCHES "gcc")

    IF(NOT SimTK_INSTALL_PREFIX)
        SET(SimTK_INSTALL_PREFIX "/usr/local/SimTK")
    ENDIF(NOT SimTK_INSTALL_PREFIX)

    IF(NOT SimTK_SDK)
        SET(SimTK_SDK "/usr/local/SimTK")
    ENDIF(NOT SimTK_SDK)

ELSE(${CMAKE_C_COMPILER} MATCHES "gcc")

    IF(NOT SimTK_INSTALL_PREFIX)
        SET(SimTK_INSTALL_PREFIX "C:/Program Files/SimTK")
    ENDIF(NOT SimTK_INSTALL_PREFIX)

    IF(NOT SimTK_SDK)
        SET(SimTK_SDK "C:/Program Files/SimTK")
    ENDIF(NOT SimTK_SDK)

ENDIF(${CMAKE_C_COMPILER} MATCHES "gcc")

INCLUDE_DIRECTORIES(${SimTK_SDK}/core/include)
################################################################################################
################################################################################################
################################################################################################


INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/../Chemistry/include)
FILE(GLOB INC1 ../Chemistry/include/*.h)
FILE(GLOB INC2 include/*.h)
SET(INCLUDES ${INC1} ${INC2})

SET(AAT_SRCS ../Chemistry/src/AminoAcidType.cpp)
ADD_EXECUTABLE(TestAminoAcidType tests/TestAminoAcidType.cpp ${AAT_SRCS} ${INCLUDES})

SET(AA_SRCS ${AAT_SRCS} ../Chemistry/src/AminoAcid.cpp)
ADD_EXECUTABLE(TestAminoAcid tests/TestAminoAcid.cpp ${AA_SRCS})

SET(PROTEIN_SRCS ${AA_SRCS} ../Chemistry/src/Protein.cpp)
ADD_EXECUTABLE(TestProtein tests/TestProtein.cpp ${PROTEIN_SRCS})

# FILE(GLOB SIMBODY_SRCS ../src/*.cpp)
SET(MOL_SIM_SRCS ${PROTEIN_SRCS} ../Chemistry/src/Ethane.cpp)
ADD_EXECUTABLE(TestSimMol1 tests/TestSimMol1.cpp ${MOL_SIM_SRCS} ${INCLUDES})

TARGET_LINK_LIBRARIES(TestSimMol1 "C:/Program Files/SIMTK/core/lib/simtkcommon.lib")
TARGET_LINK_LIBRARIES(TestSimMol1 "C:/Program Files/SIMTK/core/lib/simmatrix.lib")
TARGET_LINK_LIBRARIES(TestSimMol1 "C:/Program Files/SIMTK/core/lib/simbody.lib")
