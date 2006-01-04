FILE(GLOB SOURCE_FILES *.cpp)
FILE(GLOB INCLUDE_FILES *.h)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/MultibodyComputation/SimbodyRigidBodies
                    ${PROJECT_SOURCE_DIR}/Mechanics/include
                    ${PROJECT_SOURCE_DIR}/Geometry/include
		    ${PROJECT_SOURCE_DIR}/include
                    ${SIMTK_FRAMEWORK_INSTALL_DIR}/include 
                    ${PROJECT_SOURCE_DIR}/MultibodyComputation/CDSlib) ### this is going away
ADD_LIBRARY(SimbodyRigidBodies ${SOURCE_FILES} ${INCLUDE_FILES})
