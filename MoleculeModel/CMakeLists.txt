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

## Watch out -- this next line will cause CMake to find installed versions of
## simmatrix libraries before the local one. Be sure to delete installed ones first.
LINK_DIRECTORIES(${SimTK_SDK}/core/lib${LIB64})

################################################################################################
################################################################################################
################################################################################################


INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/../Chemistry/include ${CMAKE_CURRENT_SOURCE_DIR}/../include)
# INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/../Chemistry/include)
FILE(GLOB INC1 ../Chemistry/include/*.h)
FILE(GLOB INC2 include/*.h)
SET(INCLUDES ${INC1} ${INC2})

SET(ATOM_SRCS ../Chemistry/src/Atom.cpp ../Chemistry/src/ChemicalElement.cpp)
ADD_EXECUTABLE(TestAtom tests/TestAtom.cpp ${ATOM_SRCS} ${INCLUDES})

SET(AAT_SRCS ../Chemistry/src/AminoAcidType.cpp ${ATOM_SRCS})
ADD_EXECUTABLE(TestAminoAcidType tests/TestAminoAcidType.cpp ${AAT_SRCS} ${INCLUDES})

SET(AA_SRCS ${AAT_SRCS} ../Chemistry/src/AminoAcid.cpp)
ADD_EXECUTABLE(TestAminoAcid tests/TestAminoAcid.cpp ${AA_SRCS})

SET(PROTEIN_SRCS ${AA_SRCS} ../Chemistry/src/Protein.cpp ../Chemistry/src/Molecule.cpp)
ADD_EXECUTABLE(TestProtein tests/TestProtein.cpp ${PROTEIN_SRCS})

FILE(GLOB MODEL_SRCS ../MoleculeModel/src/*.cpp)
SET(MOL_SIM_SRCS ${PROTEIN_SRCS} ${MODEL_SRCS} ../Chemistry/src/Ethane.cpp)
ADD_EXECUTABLE(TestSimMol0 tests/TestSimMol0.cpp ${MOL_SIM_SRCS} ${INCLUDES})
ADD_EXECUTABLE(TestSimMol1 tests/TestSimMol1.cpp ${MOL_SIM_SRCS} ${INCLUDES})

TARGET_LINK_LIBRARIES(TestSimMol0 simtkcommon)
TARGET_LINK_LIBRARIES(TestSimMol0 simmatrix)
TARGET_LINK_LIBRARIES(TestSimMol0 simbody)

TARGET_LINK_LIBRARIES(TestSimMol1 simtkcommon)
TARGET_LINK_LIBRARIES(TestSimMol1 simmatrix)
TARGET_LINK_LIBRARIES(TestSimMol1 simbody)
