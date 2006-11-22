iptest: iptest.o Optimizer.o LBFGSOptimizer.o LBFGSBOptimizer.o InteriorPointOptimizer.o lbfgs.o lbfgsb.o
	g++ -g -m32 iptest.o Optimizer.o LBFGSOptimizer.o LBFGSBOptimizer.o InteriorPointOptimizer.o lbfgsb.o lbfgs.o -L/usr/local/SimTK/core/lib -L../sharedTarget -I/usr/local/SimTK/core/include -lSimTKlapack -lsimmath_d  -o iptest

driver1: driver1.o Optimizer.o LBFGSOptimizer.o LBFGSBOptimizer.o lbfgs.o lbfgsb.o
	g++ -g -m32 driver1.o Optimizer.o LBFGSOptimizer.o LBFGSBOptimizer.o lbfgsb.o lbfgs.o -L/usr/local/SimTK/core/lib -I/usr/local/SimTK/core/include -lSimTKlapack -lsimmatrix_d -o driver1

cpptest: cpptest.o Optimizer.o LBFGSOptimizer.o LBFGSBOptimizer.o lbfgs.o lbfgsb.o
	g++ -g -m32 cpptest.o Optimizer.o LBFGSOptimizer.o LBFGSBOptimizer.o lbfgsb.o lbfgs.o -L/usr/local/SimTK/core/lib -I/usr/local/SimTK/core/include -lSimTKlapack -lsimmatrix_d -o driver1

iptest.o: iptest.cc
	g++ -g -m32  -I/usr/local/SimTK/core/include -c iptest.cc
driver1.o: driver1.cc
	g++ -g -m32  -I/usr/local/SimTK/core/include -c driver1.cc
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
