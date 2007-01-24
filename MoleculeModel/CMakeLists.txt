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
FILE(GLOB INC1 ../Chemistry/include/chemistry/*.h)
FILE(GLOB INC2 include/*.h)
FILE(GLOB INC3 include/internal/*.h)
SET(INCLUDES ${INC1} ${INC2} ${INC3})

FILE(GLOB CHEM_SRCS ../Chemistry/src/*.cpp)

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

ADD_EXECUTABLE(TestSimMol0 tests/TestSimMol0.cpp ${MODEL_SRCS} ${CHEM_SRCS} ${INCLUDES})
ADD_EXECUTABLE(TestSimMol1 tests/TestSimMol1.cpp ${MODEL_SRCS} ${CHEM_SRCS} ${INCLUDES})
TARGET_LINK_LIBRARIES(TestSimMol0 SimTKcommon)
TARGET_LINK_LIBRARIES(TestSimMol0 SimTKsimbody)
TARGET_LINK_LIBRARIES(TestSimMol1 SimTKcommon)
TARGET_LINK_LIBRARIES(TestSimMol1 SimTKsimbody)

ADD_EXECUTABLE(TestArgonModel tests/TestArgonModel.cpp ${MODEL_SRCS} ${CHEM_SRCS} ${INCLUDES})
TARGET_LINK_LIBRARIES(TestArgonModel SimTKcommon_d SimTKsimbody_d)

ADD_EXECUTABLE(TestOxygenModel tests/TestOxygenModel.cpp ${MODEL_SRCS} ${CHEM_SRCS} ${INCLUDES})
TARGET_LINK_LIBRARIES(TestOxygenModel SimTKcommon_d SimTKsimbody_d)

ADD_EXECUTABLE(TestAtomCopyConstructor tests/TestAtomCopyConstructor.cpp ${MODEL_SRCS} ${CHEM_SRCS} ${INCLUDES})
TARGET_LINK_LIBRARIES(TestAtomCopyConstructor SimTKcommon_d SimTKsimbody_d)

ADD_EXECUTABLE(TestAtomModelCopyConstructor tests/TestAtomModelCopyConstructor.cpp ${MODEL_SRCS} ${CHEM_SRCS} ${INCLUDES})
TARGET_LINK_LIBRARIES(TestAtomModelCopyConstructor SimTKcommon_d SimTKsimbody_d)

ADD_EXECUTABLE(TestRealize tests/TestRealize.cpp  ${MODEL_SRCS} ${CHEM_SRCS} ${INCLUDES})
TARGET_LINK_LIBRARIES(TestRealize SimTKcommon_d SimTKsimbody_d)
