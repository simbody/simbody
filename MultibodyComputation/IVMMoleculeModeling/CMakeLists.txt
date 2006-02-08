FILE(GLOB SOURCE_FILES *.cpp)
FILE(GLOB INCLUDE_FILES *.h)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/MultibodyComputation/IVMMoleculeModeling
                    ${PROJECT_SOURCE_DIR}/MultibodyComputation/IVMRigidBodies
                    ${PROJECT_SOURCE_DIR}/MultibodyComputation/CDSlib
                    ${PROJECT_SOURCE_DIR}/MultibodyComputation/SimbodyRigidBodies
                    ${PROJECT_SOURCE_DIR}/MultibodyModeling/include
                    ${PROJECT_SOURCE_DIR}/include 
                    ${PROJECT_SOURCE_DIR}/Geometry/include
                    ${PROJECT_SOURCE_DIR}/Mechanics/include)

ADD_LIBRARY(IVMMoleculeModeling ${SOURCE_FILES} ${INCLUDE_FILES})

ADD_EXECUTABLE(IVMMoleculeTest ${SOURCE_FILES} ${INCLUDE_FILES})
TARGET_LINK_LIBRARIES(IVMMoleculeTest IVMRigidBodies CDSlib
                      SimbodyRigidBodies simmatrix
                      ${MATH_LIBS_TO_USE})
