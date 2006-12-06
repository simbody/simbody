iptest: ipopt_test.o Optimizer.o LBFGSOptimizer.o LBFGSBOptimizer.o lbfgs.o lbfgsb.o InteriorPointOptimizer.o
	g++ -g -m32 ipopt_test.o LBFGSOptimizer.o LBFGSBOptimizer.o lbfgsb.o lbfgs.o  InteriorPointOptimizer.o Optimizer.o -L/usr/local/SimTK/core/lib -I/usr/local/SimTK/core/include -L../sharedTarget -lSimTKlapack -lsimmatrix_d -lsimmath_d -o iptest
lbfgsb: lbfgsb_test.o Optimizer.o LBFGSOptimizer.o LBFGSBOptimizer.o lbfgs.o lbfgsb.o InteriorPointOptimizer.o
	g++ -g -m32 lbfgsb_test.o LBFGSOptimizer.o LBFGSBOptimizer.o lbfgsb.o lbfgs.o  InteriorPointOptimizer.o Optimizer.o -L/usr/local/SimTK/core/lib -I/usr/local/SimTK/core/include -L../sharedTarget -lSimTKlapack -lsimmatrix_d -lsimmath_d -o lbfgsb
lbfgs: lbfgs_test.o Optimizer.o LBFGSOptimizer.o LBFGSBOptimizer.o lbfgs.o lbfgsb.o InteriorPointOptimizer.o
	g++ -g -m32 lbfgs_test.o LBFGSOptimizer.o LBFGSBOptimizer.o lbfgsb.o lbfgs.o  InteriorPointOptimizer.o Optimizer.o -L/usr/local/SimTK/core/lib -I/usr/local/SimTK/core/include -L../sharedTarget -lSimTKlapack -lsimmatrix_d -lsimmath_d -o lbfgs

ctest: ctest.o optimizer_C.o Optimizer.o LBFGSOptimizer.o LBFGSBOptimizer.o InteriorPointOptimizer.o lbfgs.o lbfgsb.o
	gcc -g -m32 ctest.o optimizer_C.o Optimizer.o LBFGSOptimizer.o LBFGSBOptimizer.o InteriorPointOptimizer.o lbfgsb.o lbfgs.o -L/usr/local/SimTK/core/lib -L../sharedTarget -I/usr/local/SimTK/core/include -lSimTKlapack -lsimmath_d  -o ctest

lbfgsb_test: lbfgsb_test.o Optimizer.o LBFGSOptimizer.o LBFGSBOptimizer.o lbfgs.o lbfgsb.o InteriorPointOptimizer.o
	g++ -g -m32 lbfgsb_test.o Optimizer.o LBFGSOptimizer.o LBFGSBOptimizer.o InteriorPointOptimizer.o lbfgsb.o lbfgs.o -L/usr/local/SimTK/core/lib -L../sharedTarget -I/usr/local/SimTK/core/include -lSimTKlapack -lsimmath_d  -o lbfgsb_test

lbfgs_test.o: lbfgs_test.cpp Simmath.h
	g++ -g -m32  -I/usr/local/SimTK/core/include -c lbfgs_test.cpp

ctest.o: ctest.c Simmath.h
	gcc -g -m32  -I/usr/local/SimTK/core/include -c ctest.c
ipopt_test.o: ipopt_test.cpp
	g++ -g -m32  -I/usr/local/SimTK/core/include -c ipopt_test.cpp
lbfgsb_test.o: lbfgsb_test.cpp
	g++ -g -m32  -I/usr/local/SimTK/core/include -c lbfgsb_test.cpp
optimizer_C.o: optimizer_C.cpp
	g++ -g -m32  -I/usr/local/SimTK/core/include -c optimizer_C.cpp
InteriorPointOptimizer.o: InteriorPointOptimizer.cpp
	g++ -g -m32  -I/usr/local/SimTK/core/include -c InteriorPointOptimizer.cpp
Optimizer.o: Optimizer.cpp
	g++ -g -m32  -I/usr/local/SimTK/core/include -c Optimizer.cpp
LBFGSOptimizer.o: LBFGSOptimizer.cpp
	g++ -g -m32  -I/usr/local/SimTK/core/include -c LBFGSOptimizer.cpp
LBFGSBOptimizer.o:  LBFGSBOptimizer.cpp
	g++ -g -m32  -I/usr/local/SimTK/core/include -c LBFGSBOptimizer.cpp
lbfgs.o: lbfgs.cpp
	g++ -g -m32  -I/usr/local/SimTK/core/include -c lbfgs.cpp
lbfgsb.o: lbfgsb.cpp
	g++ -g -m32  -I/usr/local/SimTK/core/include -c lbfgsb.cpp
