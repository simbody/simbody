                     NVECTOR_SERIAL
                Release 2.3.0, November 2006

Serial implementation of the NVECTOR module for SUNDIALS. 

NVECTOR_SERIAL defines the content field of N_Vector to be a structure 
containing the length of the vector, a pointer to the beginning of a 
contiguous data array, and a boolean flag indicating ownership of the 
data array.

NVECTOR_SERIAL defines five macros to provide access to the content of 
a serial N_Vector, several constructors for variables of type N_Vector,
a constructor for an array of variables of type N_Vector, and destructors
for N_Vector and N_Vector array.

NVECTOR_SERIAL provides implementations for all vector operations defined
by the generic NVECTOR module in the table of operations.


A. Documentation
----------------

The serial NVECTOR implementation is fully described in the user documentation 
for any of the SUNDIALS solvers [1,2,3,4]. PostScript and PDF files for the user
guide for a particular solver are available in the solver's directory.


B. Installation
---------------

For basic installation instructions see /sundials/INSTALL. 
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

v. 2.3.0 - Nov. 2006
v. 2.2.0 - Mar. 2006
v. 2.1.1 - May. 2005
v. 2.1.0 - Apr. 2005
v. 2.0.2 - Mar. 2005
v. 2.0.1 - Jan. 2005
v. 2.0   - Dec. 2004
v. 1.0   - Jul. 2002 (first SUNDIALS release)


E. Revision History
-------------------

v. 2.2.0 (Mar. 2006) ---> v. 2.3.0 (Nov. 2006)
----------------------------------------------

- Changes related to the build system
   - reorganized source tree. Header files in ${srcdir}/include/nvector; 
     sources in ${srcdir}/src/nvec_ser
   - exported header files in ${includedir}/sundials


v. 2.1.1 (May. 2005) ---> v. 2.2.0 (Mar. 2006)
----------------------------------------------

- none

v. 2.1.0 (Apr. 2005) ---> v. 2.1.1 (May. 2005)
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

- Revised to correspond to new generic NVECTOR module
  (see sundials/shared/README).
- Extended the list of user-callable functions provided by NVECTOR_SERIAL
  outside the table of vector operations.
- Revised the F/C interface to use underscore flags for name mapping
  and to use precision flag from configure.
- Revised F/C routine NVECTOR names for uniformity.
