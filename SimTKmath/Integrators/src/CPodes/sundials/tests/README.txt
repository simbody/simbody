List of serial CPODES examples

  cpdenx       : dense example
  cvbanx       : banded example
  cpsadamsx    : demonstration program for non-stiff problems 
  pend_test    : demonstration program with and without projection
  cpsdenx_lap  : dense example using Lapack linear algebra
  cpsbanx_lap  : banded example using Lapack linear algebra
  pendLr       : pendulum example using Lapack linear algebra

Sample results:

  SUNDIALS was built with the following options:

  ./configure CC=gcc F77=g77 CFLAGS="-g3 -O0" FFLAGS="-g3 -O0" --with-blas-lapack-libs="-L/home/radu/apps/lib -lSimTKlapack" --enable-examples 

  System Architecture: IA-32
  Processor Type: Intel Pentium 4 Xeon DP (i686)
  Operating System: Red Hat Enterprise Linux WS 3 (Taroon Update 7)
  C/Fortran Compilers: gcc/gfortran v4.1.0

  The SimTKlapack library provides ATLAS-tunned Blas and Lapack functions
