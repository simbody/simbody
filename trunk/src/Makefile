all:
	g++ -g -m32 -fPIC -shared  -I /usr/local/SimTK/core/include -I/usr/local/SimTK/core/include/lapack optimizer_C.c lbfgs.c optimizerImplementation.cpp optimizer.cpp -o libSimmath.so 
	g++ -g -c -m32 cpptest.cpp -I /usr/local/SimTK/core/include 
	g++ -g -m32 cpptest.o -I /usr/local/SimTK/core/include -L. -lSimmath -L/usr/local/SimTK/core/lib -lsimmatrix_d -lSimTKlapack -o cpptest
	gcc -g -m32 ctest.c -L. -lSimmath -L/usr/local/SimTK/core/lib -lsimmatrix_d -lSimTKlapack -o ctest
#	g++ -g -m32 lbfgs.c ctest.c  -I/usr/local/SimTK/core/include/lapack -L/usr/local/SimTK/core/lib -lSimTKlapack -o ctest
	/usr/bin/g77 -g -m32 ftest.f -L. -lSimmath -L/usr/local/SimTK/core/lib -lsimmatrix_d -lSimTKlapack -o ftest
