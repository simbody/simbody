/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2006/11/29 00:05:06 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the common header file for the Scaled, Preconditioned
 * Iterative Linear Solvers in CPODES.
 * -----------------------------------------------------------------
 */

#ifndef _CPSPILS_H
#define _CPSPILS_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <sundials/sundials_iterative.h>
#include <sundials/sundials_nvector.h>

/*
 * =================================================================
 *   C P S P I L S    C O N S T A N T S
 * =================================================================
 */

/* CPSPILS return values */

#define CPSPILS_SUCCESS          0
#define CPSPILS_MEM_NULL        -1
#define CPSPILS_LMEM_NULL       -2
#define CPSPILS_ILL_INPUT       -3
#define CPSPILS_MEM_FAIL        -4

/*
 * =================================================================
 *              F U N C T I O N   T Y P E S
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Type : CPSpilsPrecSetupExplFn and CPSpilsPrecSetupImplFn
 * -----------------------------------------------------------------
 *
 * Explicit ODE case
 *
 * The user-supplied preconditioner setup function PrecSetup and
 * the user-supplied preconditioner solve function PrecSolve
 * together must define left and right preconditoner matrices
 * P1 and P2 (either of which may be trivial), such that the
 * product P1*P2 is an approximation to the Newton matrix
 * M = I - gamma*J.  Here J is the system Jacobian J = df/dy,
 * and gamma is a scalar proportional to the integration step
 * size h.  The solution of systems P z = r, with P = P1 or P2,
 * is to be carried out by the PrecSolve function, and PrecSetup
 * is to do any necessary setup operations.
 *
 * The user-supplied preconditioner setup function PrecSetup
 * is to evaluate and preprocess any Jacobian-related data
 * needed by the preconditioner solve function PrecSolve.
 * This might include forming a crude approximate Jacobian,
 * and performing an LU factorization on the resulting
 * approximation to M.  This function will not be called in
 * advance of every call to PrecSolve, but instead will be called
 * only as often as necessary to achieve convergence within the
 * Newton iteration.  If the PrecSolve function needs no
 * preparation, the PrecSetup function can be NULL.
 *
 * For greater efficiency, the PrecSetup function may save
 * Jacobian-related data and reuse it, rather than generating it
 * from scratch.  In this case, it should use the input flag jok
 * to decide whether to recompute the data, and set the output
 * flag *jcurPtr accordingly.
 *
 * Each call to the PrecSetup function is preceded by a call to
 * the RhsFn f with the same (t,y) arguments.  Thus the PrecSetup
 * function can use any auxiliary data that is computed and
 * saved by the f function and made accessible to PrecSetup.
 *
 * A function PrecSetup must have the prototype given below.
 * Its parameters are as follows:
 *
 * t       is the current value of the independent variable.
 *
 * y       is the current value of the dependent variable vector,
 *          namely the predicted value of y(t).
 *
 * fy      is the vector f(t,y).
 *
 * jok     is an input flag indicating whether Jacobian-related
 *         data needs to be recomputed, as follows:
 *           jok == FALSE means recompute Jacobian-related data
 *                  from scratch.
 *           jok == TRUE  means that Jacobian data, if saved from
 *                  the previous PrecSetup call, can be reused
 *                  (with the current value of gamma).
 *         A Precset call with jok == TRUE can only occur after
 *         a call with jok == FALSE.
 *
 * jcurPtr is a pointer to an output integer flag which is
 *         to be set by PrecSetup as follows:
 *         Set *jcurPtr = TRUE if Jacobian data was recomputed.
 *         Set *jcurPtr = FALSE if Jacobian data was not recomputed,
 *                        but saved data was reused.
 *
 * gamma   is the scalar appearing in the Newton matrix.
 *
 * P_data  is a pointer to user data - the same as the P_data
 *         parameter passed to the CP*SetPreconditioner function.
 *
 * tmp1, tmp2, and tmp3 are pointers to memory allocated
 *                      for N_Vectors which can be used by
 *                      CPSpilsPrecSetupFn as temporary storage or
 *                      work space.
 *
 * NOTE: If the user's preconditioner needs other quantities,
 *       they are accessible as follows: hcur (the current stepsize)
 *       and ewt (the error weight vector) are accessible through
 *       CPodeGetCurrentStep and CPodeGetErrWeights, respectively).
 *       The unit roundoff is available as UNIT_ROUNDOFF defined in
 *       sundials_types.h.
 *
 * Returned value:
 * The value to be returned by the PrecSetup function is a flag
 * indicating whether it was successful.  This value should be
 *   0   if successful,
 *   > 0 for a recoverable error (step will be retried),
 *   < 0 for an unrecoverable error (integration is halted).
 *
 * -----------------------------------------------------------------
 *
 * Implicit ODE case
 *
 *
 * -----------------------------------------------------------------
 */

typedef int (*CPSpilsPrecSetupExplFn)(realtype t, N_Vector y, N_Vector fy,
                      booleantype jok, booleantype *jcurPtr,
                      realtype gamma, void *P_data,
                      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

typedef int (*CPSpilsPrecSetupImplFn)(realtype t, N_Vector y, N_Vector yp, N_Vector r,
                      realtype gamma, void *P_data,
                      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/*
 * -----------------------------------------------------------------
 * Type : CPSpilsPrecSolveExplFn and CPSpilsPrecSolvImplFn
 * -----------------------------------------------------------------
 *
 * Explicit ODE case
 *
 * The user-supplied preconditioner solve function PrecSolve
 * is to solve a linear system P z = r in which the matrix P is
 * one of the preconditioner matrices P1 or P2, depending on the
 * type of preconditioning chosen.
 *
 * A function PrecSolve must have the prototype given below.
 * Its parameters are as follows:
 *
 * t      is the current value of the independent variable.
 *
 * y      is the current value of the dependent variable vector.
 *
 * fy     is the vector f(t,y).
 *
 * b      is the right-hand side vector of the linear system.
 *
 * x      is the output vector computed by PrecSolve.
 *
 * gamma  is the scalar appearing in the Newton matrix.
 *
 * delta  is an input tolerance for use by PSolve if it uses
 *        an iterative method in its solution.  In that case,
 *        the residual vector Res = r - P z of the system
 *        should be made less than delta in weighted L2 norm,
 *        i.e., sqrt [ Sum (Res[i]*ewt[i])^2 ] < delta.
 *        Note: the error weight vector ewt can be obtained
 *        through a call to the routine CPodeGetErrWeights.
 *
 * lr     is an input flag indicating whether PrecSolve is to use
 *        the left preconditioner P1 or right preconditioner
 *        P2: lr = 1 means use P1, and lr = 2 means use P2.
 *
 * P_data is a pointer to user data - the same as the P_data
 *        parameter passed to the CP*SetPreconditioner function.
 *
 * tmp    is a pointer to memory allocated for an N_Vector
 *        which can be used by PSolve for work space.
 *
 * Returned value:
 * The value to be returned by the PrecSolve function is a flag
 * indicating whether it was successful.  This value should be
 *   0 if successful,
 *   positive for a recoverable error (step will be retried),
 *   negative for an unrecoverable error (integration is halted).
 *
 * -----------------------------------------------------------------
 *
 * Implicit case ODE
 *
 * -----------------------------------------------------------------
 */

typedef int (*CPSpilsPrecSolveExplFn)(realtype t, N_Vector y, N_Vector fy,
                      N_Vector b, N_Vector x,
                      realtype gamma, realtype delta,
                      int lr, void *P_data, N_Vector tmp);

typedef int (*CPSpilsPrecSolveImplFn)(realtype t, N_Vector y, N_Vector yp, N_Vector r,
                      N_Vector b, N_Vector x,
                      realtype gamma, realtype delta,
                      void *P_data, N_Vector tmp);

/*
 * -----------------------------------------------------------------
 * Type : CPSpilsJacTimesVecExplFn and CPSpilsJacTimesVecImplFn
 * -----------------------------------------------------------------
 *
 * Explicit ODE case
 *
 * The user-supplied function jtimes is to generate the product
 * J*v for given v, where J is the Jacobian df/dy, or an
 * approximation to it, and v is a given vector. It should return
 * 0 if successful a positive value for a recoverable error or
 * a negative value for an unrecoverable failure.
 *
 * A function jtimes must have the prototype given below. Its
 * parameters are as follows:
 *
 *   v        is the N_Vector to be multiplied by J.
 *
 *   Jv       is the output N_Vector containing J*v.
 *
 *   t        is the current value of the independent variable.
 *
 *   y        is the current value of the dependent variable
 *            vector.
 *
 *   fy       is the vector f(t,y).
 *
 *   jac_data is a pointer to user Jacobian data, the same as the
 *            jac_data parameter passed to the CP*SetJacTimesVecFn
 *            function.
 *
 *   tmp      is a pointer to memory allocated for an N_Vector
 *            which can be used by Jtimes for work space.
 *
 * -----------------------------------------------------------------
 *
 * Implicit ODE case
 *
 * -----------------------------------------------------------------
 */

typedef int (*CPSpilsJacTimesVecExplFn)(realtype t, N_Vector y, N_Vector fy,
                    N_Vector v, N_Vector Jv, void *jac_data,
                    N_Vector tmp);

typedef int (*CPSpilsJacTimesVecImplFn)(realtype t, realtype gm,
                    N_Vector y, N_Vector yp, N_Vector r,
                    N_Vector v, N_Vector Jv, void *jac_data,
                    N_Vector tmp1, N_Vector tmp2);

/*
 * =================================================================
 *          U S E R - C A L L A B L E   F U N C T I O N S
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Optional inputs to the CPSPILS linear solver
 * -----------------------------------------------------------------
 *
 * CPSpilsSetPrecType resets the type of preconditioner, pretype,
 *                from the value previously set.
 *                This must be one of PREC_NONE, PREC_LEFT,
 *                PREC_RIGHT, or PREC_BOTH.
 *
 * CPSpilsSetGSType specifies the type of Gram-Schmidt
 *                orthogonalization to be used. This must be one of
 *                the two enumeration constants MODIFIED_GS or
 *                CLASSICAL_GS defined in iterative.h. These correspond
 *                to using modified Gram-Schmidt and classical
 *                Gram-Schmidt, respectively.
 *                Default value is MODIFIED_GS.
 *
 * CPSpilsSetMaxl resets the maximum Krylov subspace size, maxl,
 *                from the value previously set.
 *                An input value <= 0, gives the default value.
 *
 * CPSpilsSetDelt specifies the factor by which the tolerance on
 *                the nonlinear iteration is multiplied to get a
 *                tolerance on the linear iteration.
 *                Default value is 0.05.
 *
 * CPSpilsSetPreconditioner specifies the PrecSetup and PrecSolve functions.
 *                as well as a pointer to user preconditioner data.
 *                This pointer is passed to PrecSetup and PrecSolve
 *                every time these routines are called.
 *                Default is NULL for al three arguments.
 *
 * CPSpilsSetJacTimesVecFn specifies the jtimes function and a pointer to
 *                user Jacobian data. This pointer is passed to jtimes every
 *                time the jtimes routine is called.
 *                Default is to use an internal finite difference
 *                approximation routine.
 *
 * The return value of CPSpilsSet* is one of:
 *    CPSPILS_SUCCESS   if successful
 *    CPSPILS_MEM_NULL  if the CPODES memory was NULL
 *    CPSPILS_LMEM_NULL if the CPSPILS memory was NULL
 *    CPSPILS_ILL_INPUT if an input has an illegal value
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPSpilsSetPrecType(void *cpode_mem, int pretype);
SUNDIALS_EXPORT int CPSpilsSetGSType(void *cpode_mem, int gstype);
SUNDIALS_EXPORT int CPSpilsSetMaxl(void *cpode_mem, int maxl);
SUNDIALS_EXPORT int CPSpilsSetDelt(void *cpode_mem, realtype delt);
SUNDIALS_EXPORT int CPSpilsSetPreconditioner(void *cpode_mem, void *pset, void *psolve, void *P_data);
SUNDIALS_EXPORT int CPSpilsSetJacTimesVecFn(void *cpode_mem, void *jtimes, void *jac_data);

/*
 * -----------------------------------------------------------------
 * Optional outputs from the CPSPILS linear solver
 * -----------------------------------------------------------------
 * CPSpilsGetWorkSpace returns the real and integer workspace used
 *                by the SPILS module.
 *
 * CPSpilsGetNumPrecEvals returns the number of preconditioner
 *                 evaluations, i.e. the number of calls made
 *                 to PrecSetup with jok==FALSE.
 *
 * CPSpilsGetNumPrecSolves returns the number of calls made to
 *                 PrecSolve.
 *
 * CPSpilsGetNumLinIters returns the number of linear iterations.
 *
 * CPSpilsGetNumConvFails returns the number of linear
 *                 convergence failures.
 *
 * CPSpilsGetNumJtimesEvals returns the number of calls to jtimes.
 *
 * CPSpilsGetNumRhsEvals returns the number of calls to the user
 *                 f routine due to finite difference Jacobian
 *                 times vector evaluation.
 *
 * CPSpilsGetLastFlag returns the last error flag set by any of
 *                 the CPSPILS interface functions.
 *
 * The return value of CPSpilsGet* is one of:
 *    CPSPILS_SUCCESS   if successful
 *    CPSPILS_MEM_NULL  if the CPODES memory was NULL
 *    CPSPILS_LMEM_NULL if the CPSPILS memory was NULL
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPSpilsGetWorkSpace(void *cpode_mem, long int *lenrwLS, long int *leniwLS);
SUNDIALS_EXPORT int CPSpilsGetNumPrecEvals(void *cpode_mem, long int *npevals);
SUNDIALS_EXPORT int CPSpilsGetNumPrecSolves(void *cpode_mem, long int *npsolves);
SUNDIALS_EXPORT int CPSpilsGetNumLinIters(void *cpode_mem, long int *nliters);
SUNDIALS_EXPORT int CPSpilsGetNumConvFails(void *cpode_mem, long int *nlcfails);
SUNDIALS_EXPORT int CPSpilsGetNumJtimesEvals(void *cpode_mem, long int *njvevals);
SUNDIALS_EXPORT int CPSpilsGetNumFctEvals(void *cpode_mem, long int *nfevalsLS);
SUNDIALS_EXPORT int CPSpilsGetLastFlag(void *cpode_mem, int *flag);

/*
 * -----------------------------------------------------------------
 * The following function returns the name of the constant
 * associated with a CPSPILS return flag
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT char *CPSpilsGetReturnFlagName(int flag);

#ifdef __cplusplus
}
#endif

#endif
