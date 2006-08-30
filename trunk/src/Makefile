all: cpptest ctest ftest
libSimmath.so: optimizer_C.c lbfgs.c OptimizerImplementation.cpp Optimizer.cpp OptimizerImplementation.h Optimizer.h OptimizerInterface.h
	g++ -g -m32 -fPIC -shared  -I /usr/local/SimTK/core/include -I/usr/local/SimTK/core/include/lapack optimizer_C.c lbfgs.c OptimizerImplementation.cpp Optimizer.cpp -o libSimmath.so 

cpptest.o: cpptest.cpp
	g++ -g -c -m32 cpptest.cpp -I /usr/local/SimTK/core/include 

cpptest: cpptest.o  libSimmath.so
	g++ -m32 cpptest.o -I /usr/local/SimTK/core/include -L. -lSimmath -L/usr/local/SimTK/core/lib -lsimmatrix_d -lSimTKlapack -o cpptest

ctest: libSimmath.so ctest.c
	gcc -g -m32  -I /usr/local/SimTK/core/include ctest.c -L. -lSimmath -L/usr/local/SimTK/core/lib -lsimmatrix_d -lSimTKlapack -o ctest

ftest: ftest.f  libSimmath.so
	/usr/bin/g77 -g -m32 ftest.f -L. -lSimmath -L/usr/local/SimTK/core/lib -lsimmatrix_d -lSimTKlapack -o ftest
