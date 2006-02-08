SUBDIRS(tests)

FILE(GLOB SOURCE_FILES src/*.cpp )
FILE(GLOB INCLUDE_FILES 
    include/simbody/*.h
    include/simbody/internal/*.h
    src/*.h)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/Geometry/include 
                    ${PROJECT_SOURCE_DIR}/include)

ADD_LIBRARY(Geometry ${SOURCE_FILES} ${INCLUDE_FILES})
