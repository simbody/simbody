                              CPODES
                    Release 0.0.1, December 2006
                            Radu Serban
              Center for Applied Scientific Computing, LLNL

CPODES is a solver for stiff and nonstiff ODE IVP with invariants,
given in either explicit or implict form. In other words, CPODES solves:

       Explicit-form   or   Implicit-form

DE:     y' = f(t,y)          F(t,y,y') = 0
CN:     c(t,y) = 0           c(t,y) = 0
IC:     y(t0) = y0           y(t0) = y0; y'(t0) = yp0; F(t0,t0,yp0) = 0

where t0, y0, yp0 in R^N, f: R x R^N -> R^N, F: R x R^N x R^N -> R^N
and c: R x R^N -> R^M. The initial conditions must be lie on the
invariant manifold (i.e. c(t0, y0) = 0).

CPODES is part of a software family called SUNDIALS: SUite of Nonlinear
and DIfferential/ALgebraic equation Solvers.  This suite consists of
CVODE, CVODES, KINSOL, IDA, and CPODES. The directory structure of the
package reflects this family relationship.

CPODES can be used both on serial and parallel (MPI) computers, the main
difference being in the NVECTOR module of vector kernels.  The desired
version is obtained when compiling the example files by linking the
appropriate library of NVECTOR kernels.  In the parallel version,
communication between processors is done with the MPI (Message Passage
Interface) system.

The integration method in CPODES is based on variable-order, variable-step
implicit linear multistep methods, implemented with predictor-corrector
polynomials represented through Nordsieck history arrays.
CPODES provides both a nonstiff option (Adams method) and a stiff method
(Backward Differentiation Formula). The method order for the Adams option
varies between 1 and 12, while for the BDF method it varies between 1 and 5.
The resulting nonlinear systems can be solved either with a functional
iteration (appropriate for nonstiff problems and available ONLY for ODEs
in explicit-form) or with a (modified) Newton method (appropriate for stiff
problems). Therefore, it is recommended to use BDF+Newton for stiff problems
and Adams+Functional for nonstiff problems (if possible).

When the nonlinear systems are solved using the Newton method option,
a linear solver module must be attached to CPODES.  CPODEs provides both
direct and preconditioned Krylov (iterative) linear solvers. The following
direct solvers are currently available: dense and banded direct linear solver
based on internal SUNDIALS implementations and Lapack-based linear solvers
(both dense and banded).  Three different iterative solvers are available:
scaled preconditioned GMRES (SPGMR), scaled preconditioned BiCGStab (SPBCG),
and scaled preconditioned TFQMR (SPTFQMR). When CPODES is used in conjunction
with the serial NVECTOR module, any of the above linear solvers can be employed.
However, when CPODES is used with the parallel NVECTOR module, only the Krylov
linear solvers are available. For the serial version, CPODES also provides
a banded preconditioner module called CPBANDPRE available for use with the
Krylov solvers, while for the parallel version there is a preconditioner
module called CPBBDPRE which provides a band-block-diagonal preconditioner.

For any of the above linear solvers, the user has the option of providing
Jacobian information, but all of them also provide a default difference
quotient Jacobian approximation. For explicit-form ODEs, the Jacobian
required by CPODES is J = df/dy, while for implicit-form ODEs, it is
J = J = dF/dy' + gamma*dF/dy, where gamma is a scalar proportional to
the current integration step size. The value gamma is set by CPODES and
it should not be modified by a function evaluationg the Jacobian. Note also
that the iterative linear solvers are matrix-free, in the sense that they
only require Jacobian-vector products.

CPODES ensures that the solution lies on the invariant manifold by performing
a so-called coordinate projection.  The user has the option of providing their
own projection function, but CPODES provides an internal projection algorithm.
The CPODES internal projection algorithm applies a pseudo-inverse (if the user
indicates that the constraints are linear) or a performs Moore-Penrose iterations
(for nonlinear constraints).  In either case, the pseudo inverses are never
computed explicitly, but rather as solutions of linear systems in so-called
optimization matrix form. For the solution of these special-form linear systems,
CPODES offers several different methods, based either on a Schur-complement
approach or on various decompositions of the transposed of the constraint
Jacobian: LU, QR, or QR with column pivoting, the latter being applicable to
the case in which some of the constraints are known (or suspected) to be
redundant. Currently, only the dense and Lapack dense linear solvers in
SUNDIALS provide support for the linear algebra required in the projection
algorithm.

Therefore, in an attempt to advance the solution to a new point in time (i.e.
taking a new integration step), CPODES performs the following operations:
   - predict solution
   - solve nonlinear system and correct solution
   - project solution
   - test error
   - select order and step size for next step
and includes several recovery attempts in case there are convergence failures
(or difficulties) in the nonlinear solver or in the projection step, or if the
solution fails to satisfy the error test.

Other features offered by CPODES include:
   - dense output: an approximate solution (or derivatives of the solution up to
     the current method order) can be obtained for any time within the last step.
   - integration of pure quadrature equations: CPODES can treat more efficiently
     the integration of equations of the type yq' = q(t,y) if these are specified
     as such. This can be used to evaluate integrals I = \int_t0^t q(t,y) dt
     as I = yq(t), where yq' = q(t,y), yq(t0) = 0.
   - rootfinding: additional user-specified functions can be monitored and the
     solver will return whenever a root of such a function is encountered.
   - stability limit detection: since BDF methods of order higher than 2 are not
     absolutely stable. CPODES provides a mechanism which attempts to detect,
     in a direct manner, the presence of a stability region boundary that is
     limiting the step sizes in the presence of a weakly damped oscillation
     and reduce the BDF order.


A. Documentation
----------------

/sundials/doc/cpodes contains PDF files for the CPODES User Guide (cps_guide.pdf)
and the CPODES Examples (cps_examples.pdf) documents.

B. Installation
---------------

For basic installation instructions see the file /sundials/INSTALL_NOTES.
For complete installation instructions see the "CPODES Installation Procedure"
chapter in the CPODES User Guide.

C. References
-------------

[1] R. Serban, "User Documentation for CPODES"
    LLLNL technical report UCRL-SM-xxxxxx, Month Year.


D. Releases
-----------

v. 0.0.1       - Dec. 2006


E. Revision History
-------------------

