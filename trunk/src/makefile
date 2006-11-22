iptest: iptest.o Optimizer.o LBFGSOptimizer.o LBFGSBOptimizer.o InteriorPointOptimizer.o lbfgs.o lbfgsb.o
	g++ -g -m32 iptest.o Optimizer.o LBFGSOptimizer.o LBFGSBOptimizer.o InteriorPointOptimizer.o lbfgsb.o lbfgs.o -L/usr/local/SimTK/core/lib -L../sharedTarget -I/usr/local/SimTK/core/include -lSimTKlapack -lsimmath_d  -o iptest

lbfgsb_test: lbfgsb_test.o Optimizer.o LBFGSOptimizer.o LBFGSBOptimizer.o lbfgs.o lbfgsb.o InteriorPointOptimizer.o
	g++ -g -m32 lbfgsb_test.o Optimizer.o LBFGSOptimizer.o LBFGSBOptimizer.o InteriorPointOptimizer.o lbfgsb.o lbfgs.o -L/usr/local/SimTK/core/lib -L../sharedTarget -I/usr/local/SimTK/core/include -lSimTKlapack -lsimmath_d  -o lbfgsb_test

cpptest: cpptest.o Optimizer.o LBFGSOptimizer.o LBFGSBOptimizer.o lbfgs.o lbfgsb.o
	g++ -g -m32 cpptest.o Optimizer.o LBFGSOptimizer.o LBFGSBOptimizer.o lbfgsb.o lbfgs.o -L/usr/local/SimTK/core/lib -I/usr/local/SimTK/core/include -lSimTKlapack -lsimmatrix_d -o driver1

iptest.o: iptest.cc
	g++ -g -m32  -I/usr/local/SimTK/core/include -c iptest.cc
lbfgsb_test.o: lbfgsb_test.cc
	g++ -g -m32  -I/usr/local/SimTK/core/include -c lbfgsb_test.cc
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
