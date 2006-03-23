FILE(GLOB SOURCE_FILES *.cpp)
FILE(GLOB INCLUDE_FILES *.h)


INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include 
                    ${PROJECT_SOURCE_DIR}/Geometry/include
                    ${PROJECT_SOURCE_DIR}/Mechanics/include)


ADD_EXECUTABLE(SBPendulum1 ${SOURCE_FILES} ${INCLUDE_FILES})
TARGET_LINK_LIBRARIES(SBPendulum1 simbody Mechanics Geometry
                      simmatrix simtk ${MATH_LIBS_TO_USE})