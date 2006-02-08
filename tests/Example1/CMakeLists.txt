FILE(GLOB SOURCE_FILES *.cpp)
FILE(GLOB INCLUDE_FILES *.h)


INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include 
                    ${PROJECT_SOURCE_DIR}/Geometry/include
                    ${PROJECT_SOURCE_DIR}/Mechanics/include
                    ${PROJECT_SOURCE_DIR}/MultibodyComputation/IVMSimbodyInterface/include
                    ${PROJECT_SOURCE_DIR}/MultibodyModeling/include)


ADD_EXECUTABLE(Example1 ${SOURCE_FILES} ${INCLUDE_FILES})
TARGET_LINK_LIBRARIES(Example1 IVMSimbodyInterface SimbodyRigidBodies MultibodyModeling Mechanics Geometry IVMRigidBodies CDSlib
                      simmatrix simtk ${MATH_LIBS_TO_USE})