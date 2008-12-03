ADD_EXECUTABLE(PolynomialTest       PolynomialTest.cpp)
ADD_EXECUTABLE(PolynomialTestStatic PolynomialTestStatic.cpp)

TARGET_LINK_LIBRARIES(PolynomialTest       ${SHARED_TARGET})
TARGET_LINK_LIBRARIES(PolynomialTestStatic ${STATIC_TARGET})

SET_TARGET_PROPERTIES(PolynomialTest PROPERTIES
	PROJECT_LABEL "Test - PolynomialTest")

SET_TARGET_PROPERTIES(PolynomialTestStatic PROPERTIES
	PROJECT_LABEL "Test - PolynomialTestStatic")
