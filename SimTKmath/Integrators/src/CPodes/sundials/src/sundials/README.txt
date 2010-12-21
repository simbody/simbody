                           SUNDIALS 
                        Shared Module
                  Release 2.4.0, *** 2007


The family of solvers referred to as SUNDIALS consists of solvers CVODE 
(for ODE systems), CVODES (ODE with sensitivity analysis capabilities),
IDA (for differential-algebraic systems), and KINSOL (for nonlinear 
algebraic systems), 

The various solvers of this family share many subordinate modules contained
in this module:
- generic NVECTOR module
- generic linear solver modules (band, dense, smalldense, spgmr, bcg, tfqmr)
- definitions of SUNDIALS types (realtype, booleantype)
- common math functions (RpowerI, RPowerR, RSqrt, RAbs,...)


A. Documentation
----------------
All shared submodules are fully described in the user documentation for any of 
the SUNDIALS solvers [1,2,3,4]. PostScript and PDF files for the uer guide for 
a particular solver are available in the solver's directory.


B. Installation
---------------

For basic installation instructions see the file /sundials/INSTALL. 
For complete installation instructions see any of the user guides.


C. References
-------------

[1] A. C. Hindmarsh and R. Serban, "User Documentation for CVODE v2.4.0," 
    LLLNL technical report UCRL-MA-208108, November 2004.

[2] A. C. Hindmarsh and R. Serban, "User Documentation for CVODES v2.4.0," 
    LLNL technical report UCRL-MA-208111, November 2004.

[3] A. C. Hindmarsh and R. Serban, "User Documentation for IDA v2.4.0," 
    LLNL technical report UCRL-MA-208112, November 2004.

[4] A. M. Collier, A. C. Hindmarsh, R. Serban,and C. S. Woodward, "User 
    Documentation for KINSOL v2.4.0," LLNL technical report UCRL-MA-208116, 
    November 2004.


D. Releases
-----------

v. 2.4.0 - ***. 2007
v. 2.3.0 - Nov. 2006
v. 2.2.0 - Mar. 2006
v. 2.1.1 - May. 2005
v. 2.1.0 - Apr. 2005
v. 2.0.2 - Mar. 2005
v. 2.0.1 - Jan. 2005
v. 2.0   - Dec. 2004
v. 1.0   - Jul. 2002 (first SUNDIALS release)
v. 0.0   - Mar. 2002


E. Revision History
-------------------

v. 2.3.0 (Nov. 2006) ---> v. 2.4.0 (***. 2007)
---------------------------------------------------------

- New features
   - added a new generic linear solver module based on Blas + Lapack
     for both dense and banded matrices.

- Changes to user interface
   - common functionality for all direct linear solvers (dense, band, and
     the new Lapack solver) has been collected into the DLS (Direct Linear
     Solver) module, implemented in the files sundials_direct.h and
     sundials_direct.c (similar to the SPILS module for the iterative linear 
     solvers).
   - in order to include the new Lapack-based linear solver, all dimensions
     for the above linear solvers (problem sizes, bandwidths,... including 
     the underlying matrix data types) are now of type 'int' (and not 'long int').


v. 2.2.0 (Mar. 2006) ---> v. 2.3.0 (Nov. 2006)
----------------------------------------------

- Changes to the user interface
   - modified sundials_dense and sundials_smalldense to work with 
     rectangular m by n matrices (m <= n). 

- Changes related to the build system
   - reorganized source tree
   - exported header files are installed in solver-specific subdirectories
     of ${includedir}
   - sundialsTB is distributed only as part of the SUNDIALS tarball

v. 2.1.1 (May 2005) ---> v. 2.2.0 (Mar. 2006)
----------------------------------------------

- New features
   - added SPBCG (scaled preconditioned Bi-CGStab) linear solver module
   - added SPTFQMR (scaled preconditioned TFQMR) linear solver module

- Changes related to the build system
   - updated configure script and Makefiles for Fortran examples to avoid C++
     compiler errors (now use CC and MPICC to link only if necessary)
   - SUNDIALS shared header files are installed under a 'sundials' subdirectory
     of the install include directory
   - the shared object files are now linked into each SUNDIALS library rather
     than into a separate libsundials_shared library

- Changes to the user interface
   - added prefix 'sundials_' to all shared header files

v. 2.1.0 (Apr. 2005) ---> v. 2.1.1 (May.2005)
----------------------------------------------

- Changes to data structures
   - added N_VCloneEmpty to global vector operations table

v. 2.0.2 (Mar. 2005) ---> v. 2.1.0 (Apr. 2005)
----------------------------------------------

- none

v. 2.0.1 (Jan. 2005) ---> v. 2.0.2 (Mar. 2005)
----------------------------------------------

- Changes related to the build system
   - fixed autoconf-related bug to allow configuration with the PGI Fortran compiler
   - modified to use customized detection of the Fortran name mangling scheme 
     (autoconf's AC_F77_WRAPPERS routine is problematic on some platforms)
   - added --with-mpi-flags as a configure option to allow user to specify
     MPI-specific flags
   - updated Makefiles for Fortran examples to avoid C++ compiler errors (now use
     CC and MPICC to link)

v. 2.0 (Dec. 2004) ---> v. 2.0.1 (Jan. 2005)
--------------------------------------------

- Changes related to the build system
   - changed order of compiler directives in header files to avoid compilation
     errors when using a C++ compiler.

v. 1.0 (Jul. 2002) ---> v. 2.0 (Dec. 2004)
------------------------------------------

- Changes to the generic NVECTOR module
   - removed machEnv, redefined table of vector operations (now contained
     in the N_Vector structure itself).
   - all SUNDIALS functions create new N_Vector variables through cloning, using
     an N_Vector passed by the user as a template.
   - a particular NVECTOR implementation is supposed to provide user-callable 
     constructor and destructor functions.
   - removed from structure of vector operations the following functions:
     N_VNew, N_VNew_S, N_VFree, N_VFree_S, N_VMake, N_VDispose, N_VGetData,
     N_VSetData, N_VConstrProdPos, and N_VOneMask.
   - added in structure of vector operations the following functions:
     N_VClone, N_VDestroy, N_VSpace, N_VGetArrayPointer, N_VSetArrayPointer,
     and N_VWrmsNormMask.
   - Note that nvec_ser and nvec_par are now separate modules outside the 
     shared SUNDIALS module.

- Changes to the generic linear solvers
   - in SPGMR, added a dummy N_Vector argument to be used as a template
     for cloning.
   - in SPGMR, removed N (problem dimension) from argument list of SpgmrMalloc.
   - iterative.{c,h} replace iterativ.{c,h}
   - modified constant names in iterative.h (preconditioner types are prefixed
     with 'PREC_').
   - changed numerical values for MODIFIED_GS (from 0 to 1) and CLASSICAL_GS
     (from 1 to 2).

- Changes to sundialsmath submodule
   - replaced internal routine for estimation of unit roundoff with definition
     of unit roundoff from float.h
   - modified functions to call appropriate math routines given the precision 
     level specified by the user.

- Changes to sundialstypes submodule
   - removed type 'integertype'.
   - added definitions for 'BIG_REAL', 'SMALL_REAL', and 'UNIT_ROUNDOFF' using
     values from float.h based on the precision.
   - changed definition of macro RCONST to depend on precision.

v 0.0 (Mar. 2002) ---> v. 1.0 (Jul. 2002)
-----------------------------------------

20020321  Defined and implemented generic NVECTOR module, and separate serial/
          parallel NVECTOR modules, including serial/parallel F/C interfaces.
          Modified dense and band backsolve routines to take real* type for
          RHS and solution vector.
20020329  Named the DenseMat, BandMat, and SpgmrMemRec structures.
20020626  Changed type names to realtype, integertype, booleantype.
          Renamed llnltypes and llnlmath files.

