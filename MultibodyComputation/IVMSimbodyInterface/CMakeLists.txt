SUBDIRS(tests)

FILE(GLOB SOURCE_FILES src/*.cpp )
FILE(GLOB INCLUDE_FILES 
    include/*.h
    include/simbody/*.h
    include/simbody/internal/*.h
    src/*.h)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/MultibodyComputation/IVMSimbodyInterface/include
                    ${PROJECT_SOURCE_DIR}/MultibodyComputation/IVMRigidBodies
                    ${PROJECT_SOURCE_DIR}/MultibodyComputation/CDSLib
                    ${PROJECT_SOURCE_DIR}/MultibodyModeling/include
                    ${PROJECT_SOURCE_DIR}/Mechanics/include
                    ${PROJECT_SOURCE_DIR}/Geometry/include  
                    ${PROJECT_SOURCE_DIR}/include
                    ${SIMTK_FRAMEWORK_INSTALL_DIR}/include)


ADD_LIBRARY(IVMSimbodyInterface ${SOURCE_FILES} ${INCLUDE_FILES})
