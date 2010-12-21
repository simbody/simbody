all: UserGuideLimits UserGuide ConstrainedNumericalDiffOptimization ConstrainedOptimization ParameterConstrainedOptimization SimpleDifferentiator UnconstrainedNumericalDiffOptimization UnconstrainedOptimization IPOPTParameterConstrainedOptimization

inc_path = -I../include -I/usr/local/SimTK/core/include -I../include/simmath -I../include/simmath/internal 
lib_path = -L/usr/local/SimTK/core/lib64
libs     = -lSimTKmath_d -lSimTKlapack 


IPOPTParameterConstrainedOptimization : IPOPTParameterConstrainedOptimization.cpp
	g++ $(inc_path) $(lib_path)  $(libs) IPOPTParameterConstrainedOptimization.cpp -o IPOPTParameterConstrainedOptimization
UserGuideLimits : UserGuideLimits.cpp
	g++ $(inc_path) $(lib_path)  $(libs) UserGuideLimits.cpp -o UserGuideLimits

UserGuide : UserGuide.cpp
	g++ $(inc_path) $(lib_path)  $(libs) UserGuide.cpp -o UserGuide

ConstrainedNumericalDiffOptimization : ConstrainedNumericalDiffOptimization.cpp
	g++ $(inc_path) $(lib_path)  $(libs) ConstrainedNumericalDiffOptimization.cpp -o ConstrainedNumericalDiffOptimization

ConstrainedOptimization : ConstrainedOptimization.cpp
	g++ $(inc_path) $(lib_path)  $(libs) ConstrainedOptimization.cpp -o ConstrainedOptimization

ParameterConstrainedOptimization : ParameterConstrainedOptimization.cpp
	g++ $(inc_path) $(lib_path)  $(libs) ParameterConstrainedOptimization.cpp -o ParameterConstrainedOptimization

SimpleDifferentiator : SimpleDifferentiator.cpp
	g++ $(inc_path) $(lib_path)  $(libs) SimpleDifferentiator.cpp -o SimpleDifferentiator

UnconstrainedNumericalDiffOptimization : UnconstrainedNumericalDiffOptimization.cpp
	g++ $(inc_path) $(lib_path)  $(libs) UnconstrainedNumericalDiffOptimization.cpp -o UnconstrainedNumericalDiffOptimization

UnconstrainedOptimization : UnconstrainedOptimization.cpp
	g++ $(inc_path) $(lib_path)  $(libs) UnconstrainedOptimization.cpp -o UnconstrainedOptimization

clean: 
	rm ConstrainedNumericalDiffOptimization 
	rm ConstrainedOptimization 
	rm ParameterConstrainedOptimization 
	rm SimpleDifferentiator 
	rm UnconstrainedNumericalDiffOptimization 
	rm UnconstrainedOptimization
