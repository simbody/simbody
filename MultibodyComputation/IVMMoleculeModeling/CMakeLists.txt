FILE(GLOB SOURCE_FILES *.cpp)
FILE(GLOB INCLUDE_FILES *.h)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/MultibodyComputation/IVMMoleculeModeling
                    ${PROJECT_SOURCE_DIR}/MultibodyComputation/IVMRigidBodies
                    ${PROJECT_SOURCE_DIR}/MultibodyComputation/CDSlib)
ADD_LIBRARY(IVMMoleculeModeling ${SOURCE_FILES} ${INCLUDE_FILES})

ADD_EXECUTABLE(IVMMoleculeTest ${SOURCE_FILES} ${INCLUDE_FILES})
TARGET_LINK_LIBRARIES(IVMMoleculeTest IVMRigidBodies CDSlib ${MATH_LIBS_TO_USE})
