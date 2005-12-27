FILE(GLOB SOURCE_FILES *.cpp)
FILE(GLOB INCLUDE_FILES *.h)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/MultibodyComputation/IVMXplorInterface 
                    ${PROJECT_SOURCE_DIR}/MultibodyComputation/IVMMoleculeModeling 
                    ${PROJECT_SOURCE_DIR}/MultibodyComputation/IVMRigidBodies 
                    ${PROJECT_SOURCE_DIR}/MultibodyComputation/CDSlib)
ADD_LIBRARY(IVMXplorInterface ${SOURCE_FILES} ${INCLUDE_FILES})
