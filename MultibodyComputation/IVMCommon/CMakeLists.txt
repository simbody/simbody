FILE(GLOB SOURCE_FILES *.cpp)
FILE(GLOB INCLUDE_FILES *.h)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/MultibodyComputation/IVMCommon 
                    ${PROJECT_SOURCE_DIR}/MultibodyComputation/IVMMoleculeModeling 
                    ${PROJECT_SOURCE_DIR}/MultibodyComputation/IVMRigidBodies 
                    ${PROJECT_SOURCE_DIR}/MultibodyComputation/CDSlib)
ADD_LIBRARY(IVMCommon ${SOURCE_FILES} ${INCLUDE_FILES})
