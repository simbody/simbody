
ADD_EXECUTABLE(TestForOpenSim      ScrewJointExample.cpp)

TARGET_LINK_LIBRARIES(TestForOpenSim       ${SHARED_TARGET} ${AUX_SHARED_TARGET} )

