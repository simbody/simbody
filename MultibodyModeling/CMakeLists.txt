SUBDIRS (tests)

FILE(GLOB SOURCE_FILES src/*.cpp)
FILE(GLOB INCLUDE_FILES include/simbody/*.h include/simbody/internal/*.h src/*.h)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/MultibodyModeling
                    ${PROJECT_SOURCE_DIR}/Mechanics/include
                    ${PROJECT_SOURCE_DIR}/Geometry/include
		    ${PROJECT_SOURCE_DIR}/include
                    ${PROJECT_SOURCE_DIR}/MultibodyComputation/IVMRigidBodies
                    ${SIMTK_FRAMEWORK_INSTALL_DIR}/include)
ADD_LIBRARY(MultibodyModeling ${SOURCE_FILES} ${INCLUDE_FILES})

ADD_EXECUTABLE(MultibodyModelingTest ${SOURCE_FILES} ${INCLUDE_FILES})
TARGET_LINK_LIBRARIES(MultibodyModelingTest MultibodyModeling Mechanics Geometry IVMRigidBodies CDSlib
                      BigMatrix SmallMatrix Scalar simtk ${MATH_LIBS_TO_USE})