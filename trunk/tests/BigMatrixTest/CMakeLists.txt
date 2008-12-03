ADD_EXECUTABLE(BigMatrixTest       BigMatrixTest.cpp)
ADD_EXECUTABLE(BigMatrixTestStatic BigMatrixTestStatic.cpp)

TARGET_LINK_LIBRARIES(BigMatrixTest       ${SHARED_TARGET})
TARGET_LINK_LIBRARIES(BigMatrixTestStatic ${STATIC_TARGET})

SET_TARGET_PROPERTIES(BigMatrixTest PROPERTIES
	PROJECT_LABEL "Adhoc - BigMatrixTest")

SET_TARGET_PROPERTIES(BigMatrixTestStatic PROPERTIES
	PROJECT_LABEL "Adhoc - BigMatrixTestStatic")