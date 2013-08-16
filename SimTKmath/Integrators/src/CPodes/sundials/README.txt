                            SUNDIALS 
    SUite of Nonlinear and DIfferential/ALgebraic equation Solvers
                   Release 2.3.0, November 2006
          Aaron Collier, Alan Hindmarsh, and Radu Serban
           Center for Applied Scientific Computing, LLNL


The family of solvers referred to as SUNDIALS consists of solvers CVODE 
(for ODE systems), CVODES (ODE with sensitivity analysis capabilities),
IDA (for differential-algebraic systems), and KINSOL (for nonlinear 
algebraic systems).

The various solvers of this family share many subordinate modules.
For this reason, it is organized as a family, with a directory structure 
that exploits that sharing. Each individual solver includes documentation 
on installation, along with full usage documentation.

The following is a list of the solver packages presently available.

CVODE:  A solver for stiff and nonstiff ODE systems y' = f(t,y).

CVODES: A solver for stiff and non-stiff ODE systems with sensitivity
        analysis capabilities.

IDA:    A solver for differential-algebraic systems F(t,y,y') = 0.

KINSOL: A solver for nonlinear algebraic systems F(u) = 0.

Warning to users who receive more than one of these individual solvers
at different times: The mixing of old and new versions SUNDIALS may fail.  
To avoid such failures, obtain all desired solvers at the same time.

For installation directions see the file INSTALL_NOTES.

For additional information on a particular solver, see the README file
in the solver directory (e.g. src/cvode/README).

                        Release history

+----------+------------------------------------------------------+
|          | SUNDIALS |            Solver version                 |
|   Date   |          +----------+----------+----------+----------+
|          | release  |   CVODE  | CVODES   |   IDA    |  KINSOL  |
+----------+----------+----------+----------+----------+----------+
| Jul 2002 |   1.0    |    2.0   |    1.0   |    2.0   |    2.0   |
| Dec 2004 |   2.0    |  2.2.0   |  2.1.0   |  2.2.0   |  2.2.0   |
| Jan 2005 |   2.0.1  |  2.2.1   |  2.1.1   |  2.2.1   |  2.2.1   |
| Mar 2005 |   2.0.2  |  2.2.2   |  2.1.2   |  2.2.2   |  2.2.2   |
| Apr 2005 |   2.1.0  |  2.3.0   |  2.2.0   |  2.3.0   |  2.3.0   |
| May 2005 |   2.1.1  |  2.3.0   |  2.3.0   |  2.3.0   |  2.3.0   |
| Mar 2006 |   2.2.0  |  2.4.0   |  2.4.0   |  2.4.0   |  2.4.0   |
| Nov 2006 |   2.3.0  |  2.5.0   |  2.5.0   |  2.5.0   |  2.5.0   |
+----------+----------+----------+----------+----------+----------+
