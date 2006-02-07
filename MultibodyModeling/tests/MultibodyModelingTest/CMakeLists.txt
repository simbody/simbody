FILE(GLOB SOURCE_FILES *.cpp)
FILE(GLOB INCLUDE_FILES *.h)


INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include 
                    ${PROJECT_SOURCE_DIR}/Geometry/include
                    ${PROJECT_SOURCE_DIR}/Mechanics/include
                    ${PROJECT_SOURCE_DIR}/MultibodyModeling/include
                    ${SIMTK_CORE_INSTALL_DIR}/include)


ADD_EXECUTABLE(MultibodyModelingTest ${SOURCE_FILES} ${INCLUDE_FILES})
TARGET_LINK_LIBRARIES(MultibodyModelingTest MultibodyModeling Mechanics Geometry
                      BigMatrix SmallMatrix Scalar simtk ${MATH_LIBS_TO_USE})