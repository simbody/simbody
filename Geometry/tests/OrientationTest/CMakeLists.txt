FILE(GLOB SOURCE_FILES *.cpp)
FILE(GLOB INCLUDE_FILES *.h)


INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include 
                    ${PROJECT_SOURCE_DIR}/Geometry/include)


ADD_EXECUTABLE(OrientationTest ${SOURCE_FILES} ${INCLUDE_FILES})
TARGET_LINK_LIBRARIES(OrientationTest Geometry simmatrix simtk ${MATH_LIBS_TO_USE})