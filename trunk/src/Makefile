all:
	g++ -g -m32 -fPIC -shared -I/usr/local/SimTK/core/include/lapack  lbfgs.c  optimizerImplementation.cpp optimizer.cpp -o libSimmath.so 
	gcc -g -m32 ctest.c -L. -lSimmath -L/usr/local/SimTK/core/lib -lSimTKlapack -o ctest
	g++ -g -m32 cpptest.cpp -L. -lSimmath -L/usr/local/SimTK/core/lib -lSimTKlapack -o cpptest
#	g++ -g -m32 lbfgs.c ctest.c  -I/usr/local/SimTK/core/include/lapack -L/usr/local/SimTK/core/lib -lSimTKlapack -o ctest
	/usr/bin/g77 -g -m32 ftest.f -L. -lSimmath -L/usr/local/SimTK/core/lib -lSimTKlapack -o ftest
