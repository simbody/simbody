
ADD_EXECUTABLE(Ethane       Ethane.cpp)
ADD_EXECUTABLE(EthaneStatic EthaneStatic.cpp)

TARGET_LINK_LIBRARIES(Ethane       ${SHARED_TARGET})
TARGET_LINK_LIBRARIES(EthaneStatic ${STATIC_TARGET})