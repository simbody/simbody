/*
 * -----------------------------------------------------------------
 * $Revision: 1.6 $
 * $Date: 2007/03/20 14:33:17 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the interface file for the main CPODES integrator.
 * -----------------------------------------------------------------
 *
 * CPODES is used to solve numerically the ordinary initial value
 * problem with invariants:
 *
 *    y' = f(t,y)   or     F(t,y,y') = 0
 *    c(t,y) = 0           c(t,y) = 0
 *    y(t0) = y0           y(t0) = y0; y'(t0) = yp0 
 *
 * where t0, y0, yp0 in R^N, f: R x R^N -> R^N, F: R x R^N x R^N -> R^N 
 * and c: R x R^N -> R^M.
 *
 * -----------------------------------------------------------------
 */

#ifndef _CPODES_H
#define _CPODES_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <stdio.h>
#include <sundials/sundials_nvector.h>

/*
 * =================================================================
 *              C P O D E     C O N S T A N T S
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Inputs to CPodeCreate, CPodeInit, CPodeReInit, and CPode.
 * -----------------------------------------------------------------
 * Symbolic constants for the lmm_type and nls_type inputs to 
 * CPodeCreate, the ode_type and tol_type inputs to CPodeInit and 
 * CPodeReInit, the cnstr_type input to CPodeProjInit, as well as 
 * the input mode to CPode, are given below.
 *
 * lmm_type: The user of the CPODES package specifies whether to use
 *    the CP_ADAMS (Adams-Moulton) or CP_BDF (Backward Differentiation
 *    Formula) linear multistep method. The BDF method is recommended 
 *    for stiff problems, and the CP_ADAMS method is recommended for 
 *    nonstiff problems.
 *
 * nls_type: At each internal time step, a nonlinear equation must
 *    be solved. The user can specify either CP_FUNCTIONAL iteration, 
 *    which does not require linear algebra, or a CP_NEWTON iteration, 
 *    which requires the solution of linear systems. In the CP_NEWTON 
 *    case, the user also specifies a CPODES linear solver. 
 *    CP_NEWTON is recommended in case of stiff problems.
 *
 * ode_type: The ODE system can be given either in explicit form,
 *    y' = f(t,y) (in which case ode = CP_EXPL) or in implicit form
 *    (in which case ode = CP_IMPL).
 *
 * tol_type: This parameter specifies the relative and absolute
 *    tolerance types to be used. The CP_SS tolerance type means a
 *    scalar relative and absolute tolerance. The CP_SV tolerance type 
 *    means a scalar relative tolerance and a vector absolute tolerance 
 *    (a potentially different absolute tolerance for each vector 
 *    component). The CP_WF tolerance type means that the user provides 
 *    a function (of type CPEwtFn) to set the error weight vector.
 *
 * proj_type: If performing projection on the invariant manifold,
 *    this parameter specifies whether to use the internal algorithm
 *    or a user-provided projection function. The valid values are 
 *    CP_PROJ_USER and CP_PROJ_INERNAL.
 *    A value proj_type = CP_PROJ_USER indicates that a function to 
 *    perform the projection is supplied by the user.
 *    A value proj_type = CP_PROJ_INTERNAL indicates that the internal
 *    CPODES projection algorithm is to be used. In this case, the 
 *    user must specify the constraint function and a linear solver
 *    to be used.
 *
 * proj_norm: Type of the norm in which projection is to be performed.
 *    The valid values are CP_PROJ_L2NORM (in which case projection
 *    is done in Euclidian norm) and CP_PROJ_ERRNORM (in which case
 *    the projection is done in WRMS norm).
 *
 * cnstr_type: If the internal projection algorithm is used, then
 *    cnstr_type specifies the type of constraints.
 *    If the constraints are linear (cnstr_type = CP_CNSTR_LIN), 
 *    then an improved error estimate for the projected method 
 *    (obtained by projecting the error estimate for the original 
 *    method) is available at no additional cost and will always be used.
 *    If the constraints are nonlinear (cnstr_type = CP_CNSTR_NONLIN), 
 *    obtaining the global projection operator requires an additional 
 *    evaluation of the constraint Jacobian and constructing its
 *    More-Penrose pseudo-inverse. In this case, the default is
 *    to use the error estimate of the unprojected method (unless
 *    indicated otherwise by the user).
 *
 * mode: The mode input parameter to CPode indicates the job of the
 *    solver for the next user step. The CP_NORMAL mode is to have 
 *    the solver take internal steps until it has reached or just 
 *    passed the user specified tout parameter. The solver then 
 *    interpolates in order to return an approximate value of y(tout). 
 *    The CP_ONE_STEP option tells the solver to just take one internal 
 *    step and return the solution at the point reached by that step. 
 *    The CP_NORMAL_TSTOP and CP_ONE_STEP_TSTOP modes are similar to 
 *    CP_NORMAL and CP_ONE_STEP, respectively, except that the 
 *    integration never proceeds past the value tstop (specified 
 *    through the routine CPodeSetStopTime).
 * -----------------------------------------------------------------
 */

/* lmm_type */
#define CP_ADAMS          1
#define CP_BDF            2

/* nls_type */
#define CP_FUNCTIONAL     1
#define CP_NEWTON         2

/* ode_type */
#define CP_EXPL           1
#define CP_IMPL           2

/* tol_type */
#define CP_SS             1
#define CP_SV             2
#define CP_WF             3

/* proj_type */
#define CP_PROJ_USER      1
#define CP_PROJ_INTERNAL  2

/* proj_norm */
#define CP_PROJ_L2NORM    1
#define CP_PROJ_ERRNORM   2

/* cnstr_type */
#define CP_CNSTR_LIN      1
#define CP_CNSTR_NONLIN   2

/* mode */
#define CP_NORMAL         1
#define CP_ONE_STEP       2
#define CP_NORMAL_TSTOP   3
#define CP_ONE_STEP_TSTOP 4

/*
 * ----------------------------------------
 * CPODES return flags
 * ----------------------------------------
 */

#define CP_SUCCESS               0
#define CP_TSTOP_RETURN          1
#define CP_ROOT_RETURN           2

#define CP_WARNING              99

#define CP_TOO_MUCH_WORK        -1
#define CP_TOO_MUCH_ACC         -2
#define CP_ERR_FAILURE          -3
#define CP_CONV_FAILURE         -4

#define CP_LINIT_FAIL           -5
#define CP_LSETUP_FAIL          -6
#define CP_LSOLVE_FAIL          -7
#define CP_ODEFUNC_FAIL         -8
#define CP_FIRST_ODEFUNC_ERR    -9
#define CP_REPTD_ODEFUNC_ERR    -10
#define CP_UNREC_ODEFUNC_ERR    -11
#define CP_RTFUNC_FAIL          -12

#define CP_MEM_FAIL             -20
#define CP_MEM_NULL             -21
#define CP_ILL_INPUT            -22
#define CP_NO_MALLOC            -23
#define CP_BAD_K                -24
#define CP_BAD_T                -25
#define CP_BAD_DKY              -26
#define CP_TOO_CLOSE            -27

#define CP_NO_QUAD              -30
#define CP_QUADFUNC_FAIL        -31
#define CP_FIRST_QUADFUNC_ERR   -32
#define CP_REPTD_QUADFUNC_ERR   -33
#define CP_UNREC_QUADFUNC_ERR   -34

#define CP_BAD_IS               -40
#define CP_NO_SENS              -41
#define CP_SENSFUNC_FAIL        -42
#define CP_FIRST_SENSFUNC_ERR   -43
#define CP_REPTD_SENSFUNC_ERR   -44
#define CP_UNREC_SENSFUNC_ERR   -45

#define CP_PLINIT_FAIL          -50
#define CP_PLSETUP_FAIL         -51
#define CP_PLSOLVE_FAIL         -52
#define CP_CNSTRFUNC_FAIL       -53
#define CP_PROJFUNC_FAIL        -54
#define CP_PROJ_FAILURE         -55
#define CP_REPTD_CNSTRFUNC_ERR  -56
#define CP_REPTD_PROJFUNC_ERR   -57

#define CP_FIRST_CNSTRFUNC_ERR  -60
#define CP_NO_RECOVERY          -61
#define CP_LINESEARCH_FAIL      -62

/*
 * =================================================================
 *              F U N C T I O N   T Y P E S
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Type : CPRhsFn
 * -----------------------------------------------------------------
 * If the ODE is given in explicit form, the function which defines 
 * the right-hand side of the ODE system y'=f(t,y) must have type 
 * CPRhsFn.
 * Such a function takes as input the independent variable value t, 
 * and the dependent variable vector y.  It stores the result of 
 * f(t,y) in the vector fout.  The y and fout arguments are of type
 * N_Vector. (Allocation of memory for ydot is handled within CODES)
 * The f_data parameter is the same as the f_data parameter set by 
 * the user through the CPodeInit or CPodeReInit functions. This 
 * user-supplied pointer is passed to the user's fun function every 
 * time it is called.
 *
 * A CPRhsFn should return 0 if successful, a negative value if
 * an unrecoverable error occured, and a positive value if a 
 * recoverable error (e.g. invalid y values) occured. 
 * If an unrecoverable occured, the integration is halted. 
 * If a recoverable error occured, then (in most cases) CPODES
 * will try to correct and retry.
 * -----------------------------------------------------------------
 */

typedef int (*CPRhsFn)(realtype t, N_Vector y,
		       N_Vector fout, void *f_data);

/*
 * ----------------------------------------------------------------
 * Type : CPResFn                                                   
 * ----------------------------------------------------------------
 * If the ODE is given in implicit form, the function which 
 * defines the ODE system F(t,y,y')=0 must have type CPResFn.
 *
 * A CPResFn takes as input the independent variable value t,    
 * the dependent variable vector y, and the derivative (with     
 * respect to t) of the y vector, yp.  It stores the result of   
 * F(t,y,y') in the vector fout. The y, yp, and fout arguments are 
 * of type N_Vector. The f_data parameter is the pointer f_data 
 * passed by the user to the CPodeInit or CPodeReInit functions. 
 * This user-supplied pointer is passed to the user's fun function 
 * every time it is called.
 *                                                                
 * A CPResFn function should return a value of 0 if successful, a 
 * positive value if a recoverable error occured (e.g. y has an 
 * illegal value), or a negative value if a nonrecoverable error 
 * occured. In the latter case, the program halts. If a recoverable 
 * error occured, then (in most cases) the integrator will attempt 
 * to correct and retry.
 * ----------------------------------------------------------------
 */

typedef int (*CPResFn)(realtype t, N_Vector y, N_Vector yp,
		       N_Vector fout, void *f_data);

/*
 * -----------------------------------------------------------------
 * Type : CPCnstrFn
 * -----------------------------------------------------------------
 * The function cfun defining the invariant constraints c(t,y) = 0
 * must have type CPCnstrFn.
 *
 * A CPCnstrFn takes as input the independent variable value t and
 * the dependent variable vector y.  It stores the result of c(t,y)
 * in the vector cout.  The y and cout arguments are of type
 * N_Vector. (Allocation of memory for cout is handled within CPODES)
 * The c_data parameter is the same as the c_data parameter set by 
 * the user through the CPodeProjDefineConstraints routine.
 * This user-supplied pointer is passed to the user's cfun function
 * every time it is called.
 *
 * A CPCnstrFn should return 0 if successful, a negative value if
 * an unrecoverable error occured, and a positive value if a 
 * recoverable error (e.g. invalid y values) occured. 
 * If an unrecoverable occured, the integration is halted. 
 * If a recoverable error occured, then (in most cases) CPODES
 * will try to correct and retry.
 * -----------------------------------------------------------------
 */

typedef int (*CPCnstrFn)(realtype t, N_Vector y,
			 N_Vector cout, void *c_data);

/*
 * -----------------------------------------------------------------
 * Type : CPProjFn
 * -----------------------------------------------------------------
 * A user-supplied function to performs the projection onto the
 * invariant manifold c(t,y)=0, must have type CPProjFn. 
 *
 * A CPProjFn takes as input the independent variable t and the
 * current (corrected) variable vector ycur.
 * It must compute a correction vector corr such that 
 * y = ycurr + corr lies on the manifold (i.e. c(t,y)=0).
 * The value epsProj is provided to be used in the stopping test
 * of a nonlinear solver iteration (the iterations should be 
 * terminated when the WRMS norm of the current iterate update 
 * is below epsProj).
 *
 * Note that, if the projection is orthogonal then it can be written as 
 *    y = P * ycur + alpha(t), with P^2 = P,
 * then the projected error estimate is
 *    errP = P * err
 * and ||errP|| <= ||err||.
 * The vector err contains a (scaled) version of the current error
 * estimate. CPProjFn may also compute the projected error estimate 
 * and overwrite the vector err with errP. Otherwise, it should leave
 * err unchanged.
 * 
 * If err is NULL, a CPProjFn function should attempt NO projection
 * of the error estimate (in this case, it was called from within
 * the computation of consistent initial conditions).
 *
 * A CPProjFn should return 0 if successful, a negative value if
 * an unrecoverable error occured, and a positive value if a 
 * recoverable error (e.g. invalid y values) occured. 
 * If an unrecoverable occured, the integration is halted. 
 * If a recoverable error occured, then (in most cases) CPODES
 * will try to correct and retry.
 * -----------------------------------------------------------------
 * NOTE: If the user's projection routine needs other quantities,   
 *       they are accessible as follows: the error weight vector
 *       can be obtained by calling CPodeGetErrWeights. The unit
 *       roundoff is available as UNIT_ROUNDOFF defined in 
 *       sundials_types.h.
 * -----------------------------------------------------------------
 */

typedef int (*CPProjFn)(realtype t, N_Vector ycur, N_Vector corr,
			realtype epsProj, N_Vector err, void *pdata);

/*
 * -----------------------------------------------------------------
 * Type : CPQuadFn
 * -----------------------------------------------------------------
 * The qfun function which defines the right hand side of the
 * quadrature equations q' = fQ(t,y) must have type CPQuadFn.
 * It takes as input the value of the independent variable t and
 * the vector of states y and must store the result of fQ in qout.
 * (Allocation of memory for qout is handled by CPODES).
 * The q_data parameter is the same as the q_data parameter
 * set by the user through the CPodeQuadInit function and is
 * passed to the qfun function every time it is called.
 *
 * A CPQuadFn should return 0 if successful, a negative value if
 * an unrecoverable error occured, and a positive value if a 
 * recoverable error (e.g. invalid y values) occured. 
 * If an unrecoverable occured, the integration is halted. 
 * If a recoverable error occured, then (in most cases) CPODES
 * will try to correct and retry.
 * -----------------------------------------------------------------
 */

typedef int (*CPQuadFn)(realtype t, N_Vector y, 
			N_Vector qout, void *q_data);

/*
 * -----------------------------------------------------------------
 * Type : CPRootFn
 * -----------------------------------------------------------------
 * A function g, which defines a set of functions g_i(t,y,y') whose
 * roots are sought during the integration, must have type CPRootFn.
 * The function g takes as input the independent variable value
 * t, the dependent variable vector y, and its derivative yp=y'.
 * It stores the nrtfn values g_i(t,y,y') in the realtype array gout.
 * (Allocation of memory for gout is handled within CPODES.)
 * The g_data parameter is the same as that passed by the user
 * to the CPodeRootInit routine.  This user-supplied pointer is
 * passed to the user's g function every time it is called.
 *
 * A CPRootFn should return 0 if successful or a non-zero value
 * if an error occured (in which case the integration will be halted).
 * -----------------------------------------------------------------
 */

typedef int (*CPRootFn)(realtype t, N_Vector y, N_Vector yp,
			realtype *gout, void *g_data);

/*
 * -----------------------------------------------------------------
 * Type : CPEwtFn
 * -----------------------------------------------------------------
 * A function e, which sets the error weight vector ewt, must have
 * type CPEwtFn.
 * The function e takes as input the current dependent variable y.
 * It must set the vector of error weights used in the WRMS norm:
 * 
 *   ||y||_WRMS = sqrt [ 1/N * sum ( ewt_i * y_i)^2 ]
 *
 * Typically, the vector ewt has components:
 * 
 *   ewt_i = 1 / (reltol * |y_i| + abstol_i)
 *
 * The e_data parameter is the same as that passed by the user
 * to the CPodeSetEwtFn routine.  This user-supplied pointer is
 * passed to the user's e function every time it is called.
 * A CPEwtFn e must return 0 if the error weight vector has been
 * successfuly set and a non-zero value otherwise.
 * -----------------------------------------------------------------
 */

typedef int (*CPEwtFn)(N_Vector y, N_Vector ewt, void *e_data);

/*
 * -----------------------------------------------------------------
 * Type : CPErrHandlerFn
 * -----------------------------------------------------------------
 * A function eh, which handles error messages, must have type
 * CPErrHandlerFn.
 * The function eh takes as input the error code, the name of the
 * module reporting the error, the error message, and a pointer to
 * user data, the same as that passed to CPodeSetErrHandlerFn.
 * 
 * All error codes are negative, except CP_WARNING which indicates 
 * a warning (the solver continues).
 *
 * A CPErrHandlerFn has no return value.
 * -----------------------------------------------------------------
 */

typedef void (*CPErrHandlerFn)(int error_code, 
			       const char *module, const char *function, 
			       char *msg, void *eh_data); 

/*
 * =================================================================
 *          U S E R - C A L L A B L E   F U N C T I O N S
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Function : CPodeCreate
 * -----------------------------------------------------------------
 * CPodeCreate creates an internal memory block for a problem to
 * be solved by CPODES.
 *
 * ode_type  - form in which the ODE system is provided.
 *             The legal values are CP_EXPL or CP_IMPL (see above).
 *
 * lmm_type  - type of linear multistep method to be used.
 *             The legal values are CP_ADAMS and CP_BDF (see above).
 *
 * nls_type  - type of iteration used to solve the nonlinear
 *             system that arises during each internal time step.
 *             The legal values are CP_FUNCTIONAL and CP_NEWTON
 *             for ode_type = CP_EXPL and only CP_NEWTON for
 *             ode_type = CP_IMPL.
 *
 * If successful, CPodeCreate returns a pointer to initialized
 * problem memory. This pointer should be passed to CPodeInit.
 * If an initialization error occurs, CPodeCreate prints an error
 * message to standard err and returns NULL.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void *CPodeCreate(int ode_type, int lmm_type, int nls_type);

/*
 * -----------------------------------------------------------------
 * Function : CPodeInit
 * -----------------------------------------------------------------
 * CPodeInit allocates and initializes memory for a problem to be
 * solved by CPODES.
 *
 * cpode_mem - pointer to CPODES memory returned by CPodeCreate.
 *
 * fun       - name of the C function defining the ODE system.
 *             Depending on the ODE type (see CPodeCreate), fun 
 *             should be of type CPRhsFn (ode_type=CP_EXPL) or of
 *             type CPResFn (ode_type=CP_IMPL).
 *
 * f_data    - pointer to user data that will be passed to the 
 *             fun function every time it is called.
 *
 * t0        - initial value of t.
 *
 * y0        - initial condition vector y(t0).
 *
 * yp0       - initial condition vector y'(t0).
 *
 * tol_type  - type of tolerances to be used. The legal values are:
 *             CP_SS (scalar relative and absolute tolerances),
 *             CP_SV (scalar relative tolerance and vector
 *                    absolute tolerance).
 *             CP_WF (indicates that the user will provide a
 *                    function to evaluate the error weights.
 *                    In this case, reltol and abstol are ignored.)
 *
 * reltol    - scalar relative tolerance scalar.
 *
 * abstol    - pointer to the absolute tolerance scalar or
 *             an N_Vector of absolute tolerances.
 *
 * The parameters tol_type, reltol, and abstol define a vector of
 * error weights, ewt, with components
 *   ewt[i] = 1/(reltol*abs(y[i]) + abstol)     if tol_type = CP_SS
 *   ewt[i] = 1/(reltol*abs(y[i]) + abstol[i])  if tol_type = CP_SV.
 * This vector is used in all error and convergence tests, which
 * use a weighted RMS norm on all error-like vectors v:
 *    WRMSnorm(v) = sqrt( (1/N) sum(i=1..N) (v[i]*ewt[i])^2 ),
 * where N is the problem dimension.
 *
 * Return flag:
 *  CP_SUCCESS if successful
 *  CP_MEM_NULL if the CPODES memory was NULL
 *  CP_MEM_FAIL if a memory allocation failed
 *  CP_ILL_INPUT f an argument has an illegal value.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPodeInit(void *cpode_mem, 
			      void *fun, void *f_data,
			      realtype t0, N_Vector y0, N_Vector yp0,
			      int tol_type, realtype reltol, void *abstol);

/*
 * -----------------------------------------------------------------
 * Function : CPodeReInit
 * -----------------------------------------------------------------
 * CPodeReInit re-initializes CPODES for the solution of a problem,
 * where a prior call to CPodeInit has been made with the same
 * problem size N. CPodeReInit performs the same input checking
 * and initializations that CPodeInit does.
 * But it does no memory allocation, assuming that the existing
 * internal memory is sufficient for the new problem.
 *
 * The use of CPodeReInit requires that the maximum method order,
 * maxord, is no larger for the new problem than for the problem
 * specified in the last call to CPodeInit.  This condition is
 * automatically fulfilled if the multistep method parameter lmm_type
 * is unchanged (or changed from CP_ADAMS to CP_BDF) and the default
 * value for maxord is specified.
 *
 * All of the arguments to CPodeReInit have names and meanings
 * identical to those of CPodeInit.
 *
 * The return value of CPodeReInit is equal to CP_SUCCESS = 0 if
 * there were no errors; otherwise it is a negative int equal to:
 *   CP_MEM_NULL      indicating cpode_mem was NULL (i.e.,
 *                    CPodeCreate has not been called).
 *   CP_NO_MALLOC     indicating that cpode_mem has not been
 *                    allocated (i.e., CPodeInit has not been
 *                    called).
 *   CP_ILL_INPUT     indicating an input argument was illegal
 *                    (including an attempt to increase maxord).
 * In case of an error return, an error message is also printed.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPodeReInit(void *cpode_mem, 
				void *fun, void *f_data,
				realtype t0, N_Vector y0, N_Vector yp0,
				int tol_type, realtype reltol, void *abstol);

/*
 * -----------------------------------------------------------------
 * Function : CPodeProjInit
 * -----------------------------------------------------------------
 * CPodeProjInit initializes the internal CPODES coordinate 
 * projection algorithm. It must be called after CPodeCreate and
 * CPodeInit.
 *
 * The arguments are as follows:
 *
 * cpode_mem  - pointer to CPODES memory returned by CPodeCreate.
 *
 * proj_norm  - the norm in which projection is to be done. Legal
 *              values are CP_PROJ_L2NORM and CP_PROJ_ERRNORM.
 *
 * cnstr_type - the type of constraints. 
 *                CP_CNSTR_LIN :   linear constraints.
 *                CP_CNSTR_NONLIN: nonlinear constraints.
 *
 * cfun       - name of the user-supplied function (type CPCnstrFn)
 *              defining the invariant manifold.
 * 
 * c_data     - pointer to user data that will be passed to the 
 *              cfun function every time it is called.
 *
 * ctol       - a vector of "absolute tolerances" for the constraints.
 *              In default operation, this vector is only used as
 *              a template for cloning other vectors. However, if 
 *              enabled through CPodeSetProjTestCnstr, these values,
 *              together with reltol, are used to compute the 
 *              constraint WL2 norm and a projection will be 
 *              perfomed only if ||c(t)|_WL2 >= 1.0
 *
 * The return value of CPodeProjInit is equal to CP_SUCCESS = 0 if
 * there were no errors; otherwise it is a negative int equal to:
 *   CP_MEM_NULL  - cpode_mem was NULL 
 *                  (i.e., CPodeCreate has not been called).
 *   CP_NO_MALLOC - cpode_mem has not been allocated 
 *                  (i.e., CPodeInit has not been called).
 *   CP_MEM_FAIL  - a memory allocation failed.
 *   CP_ILL_INPUT - an input argument was illegal.
 * In case of an error return, an error message is also printed.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPodeProjInit(void *cpode_mem, int proj_norm, 
				  int cnstr_type, CPCnstrFn cfun, void *c_data, 
				  N_Vector ctol);

/*
 * -----------------------------------------------------------------
 * Function : CPodeProjDefine
 * -----------------------------------------------------------------
 * CPodeProjDefine initializes coordinate projection using a funciton
 * provided by the user. It must be called after CPodeInit. 
 *
 * The arguments are as follows:
 *
 * cpode_mem  - pointer to CPODES memory returned by CPodeCreate.
 *
 * pfun       - name of the user-supplied function (type CPProjFn)
 *              which will perform the projection.
 * 
 * p_data     - pointer to user data that will be passed to the 
 *              pfun function every time it is called.
 *
 * The return value of CPodeProjDefine is CP_SUCCESS if there were 
 * no errors, or CP_MEM_NULL if the cpode_mem argument was NULL 
 * (i.e., CPodeCreate has not been called).
 * In case of an error return, an error message is also printed.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPodeProjDefine(void *cpode_mem, CPProjFn pfun, void *p_data);

/*
 * -----------------------------------------------------------------
 * Function : CPodeQuadInit and CPodeQuadReInit
 * -----------------------------------------------------------------
 * CVodeQuadInit allocates and initializes memory related to
 * quadrature integration.
 *
 * cpode_mem - pointer to CPODES memory returned by CPodeCreate
 *
 * qfun      - the user-provided integrand routine.
 *
 * q_data    - pointer to user data that will be passed to the 
 *             qfun function every time it is called.
 *
 * q0        - N_Vector with initial values for quadratures
 *             (typically q0 has all zero components).
 *
 * CPodeQuadReInit re-initializes CPODES's quadrature related memory
 * for a problem, assuming it has already been allocated in prior calls 
 * to CPodeInit and CPodeQuadInit.  All problem specification inputs 
 * are checked for errors. The number of quadratures Nq is assumed to 
 * be unchanged since the previous call to CPodeQuadInit.
 *
 * Return values:
 *  CP_SUCCESS if successful
 *  CP_MEM_NULL if the CPODES memory was NULL
 *  CP_MEM_FAIL if a memory allocation failed
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPodeQuadInit(void *cpode_mem, CPQuadFn qfun, void *q_data, N_Vector q0);
SUNDIALS_EXPORT int CPodeQuadReInit(void *cpode_mem, CPQuadFn qfun, void *q_data, N_Vector q0);

/*
 * -----------------------------------------------------------------
 * Function : CPodeRootInit
 * -----------------------------------------------------------------
 * CPodeRootInit initializes a rootfinding problem to be solved
 * during the integration of the ODE system.  It must be called
 * after CPodeCreate, and before CPode.  The arguments are:
 *
 * cpode_mem - pointer to CPODES memory returned by CPodeCreate.
 *
 * nrtfn     - number of functions g_i, an int >= 0.
 *
 * gfun      - name of user-supplied function, of type CPRootFn,
 *             defining the functions g_i whose roots are sought.
 *
 * g_data    - pointer to user data that will be passed to the 
 *             gfun function every time it is called.
 *
 * If a new problem is to be solved with a call to CPodeReInit,
 * where the new problem has no root functions but the prior one
 * did, then call CPodeRootInit with nrtfn = 0.
 *
 * The return value of CPodeRootInit is CP_SUCCESS = 0 if there were
 * no errors; otherwise it is a negative int equal to:
 *   CP_MEM_NULL  - indicating cpode_mem was NULL.
 *   CP_MEM_FAIL  - indicating a memory allocation failed.
 *   CP_ILL_INPUT - indicating nrtfn > 0 but gfun = NULL.
 * In case of an error return, an error message is also printed.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPodeRootInit(void *cpode_mem, int nrtfn, CPRootFn gfun, void *g_data);

/*
 * -----------------------------------------------------------------
 * Function : CPodeCalcIC
 * -----------------------------------------------------------------
 * CPodeCalcIC calculates corrected initial conditions that are 
 * consistent with the invariant constraints and (for implicit-form
 * ODEs) with the ODE system itself. It first projects the initial
 * guess for the state vector (given by the user through CPodeInit
 * or CPodeReInit) and then, if necessary, computes a state derivative
 * vector as solution of F(t0, y0, y0') = 0.
 *
 * Note: If the initial conditions satisfy both the constraints and
 * (for implicit-form ODEs) the ODE itself, a call to CPodeCalcIC
 * is NOT necessary.
 *
 * A call to CpodeCalcIC must be preceded by a successful call to   
 * CPodeInit or CPodeReInit for the given ODE problem, and by a     
 * successful call to the linear system solver specification      
 * routine. 
 *
 * A call to CPodeCalcIC should precede the call(s) to CPode for
 * the given problem.
 *
 * If successful, CPodeCalcIC stores internally the corrected 
 * initial conditions which will be used to start the integration
 * at the first call to CPode.
 *
 * The only argument to CPodeCalcIC is the pointer to the CPODE
 * memory block returned by CPodeCreate.
 *
 * The return value of CPodeCalcIC is one of the following:
 *
 * CP_SUCCESS
 * CP_MEM_NULL
 * CP_NO_MALLOC
 * CP_ILL_INPUT
 * CP_LINIT_FAIL
 * CP_PLINIT_FAIL
 * CP_FIRST_CNSTRFUNC_ERR
 * CP_PROJ_FAILURE
 * CP_CNSTRFUNC_FAIL
 * CP_PROJFUNC_FAIL
 * CP_PLSETUP_FAIL
 * CP_PLSOLVE_FAIL
 * CP_NO_RECOVERY
 *
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPodeCalcIC(void *cpode_mem);

/*
 * -----------------------------------------------------------------
 * Integrator optional input specification functions
 * -----------------------------------------------------------------
 * The following functions can be called to set optional inputs
 * to values other than the defaults given below:
 *
 * Function                |  Optional input / [ default value ]
 * -----------------------------------------------------------------
 *                         |
 * CPodeSetErrHandlerFn    | user-provided ErrHandler function.
 *                         | [internal]
 *                         |
 * CPodeSetErrFile         | the file pointer for an error file
 *                         | where all CPODES warning and error
 *                         | messages will be written if the default
 *                         | internal error handling function is used. 
 *                         | This parameter can be stdout (standard 
 *                         | output), stderr (standard error), or a 
 *                         | file pointer (corresponding to a user 
 *                         | error file opened for writing) returned 
 *                         | by fopen.
 *                         | If not called, then all messages will
 *                         | be written to the standard error stream.
 *                         | [stderr]
 *                         |
 * -----------------------------------------------------------------
 *                         |
 * CPodeSetEwtFn           | user-provided EwtSet function efun and 
 *                         | a pointer to user data that will be
 *                         | passed to the user's efun function
 *                         | every time efun is called.
 *                         | [internal]
 *                         | [NULL]
 *                         |
 * CPodeSetMaxOrd          | maximum LMM order to be used by the
 *                         | solver.
 *                         | [12 for Adams , 5 for BDF]
 *                         |
 * CPodeSetMaxNumSteps     | maximum number of internal steps to be
 *                         | taken by the solver in its attempt to
 *                         | reach tout.
 *                         | [500]
 *                         |
 * CPodeSetMaxHnilWarns    | maximum number of warning messages
 *                         | issued by the solver that t+h==t on the
 *                         | next internal step. A value of -1 means
 *                         | no such messages are issued.
 *                         | [10]
 *                         |
 * CPodeSetStabLimDet      | flag to turn on/off stability limit
 *                         | detection (TRUE = on, FALSE = off).
 *                         | When BDF is used and order is 3 or
 *                         | greater, CPsldet is called to detect
 *                         | stability limit.  If limit is detected,
 *                         | the order is reduced.
 *                         | [FALSE]
 *                         |
 * CPodeSetInitStep        | initial step size.
 *                         | [estimated by CPODES]
 *                         |
 * CPodeSetMinStep         | minimum absolute value of step size
 *                         | allowed.
 *                         | [0.0]
 *                         |
 * CPodeSetMaxStep         | maximum absolute value of step size
 *                         | allowed.
 *                         | [infinity]
 *                         |
 * CPodeSetStopTime        | the independent variable value past
 *                         | which the solution is not to proceed.
 *                         | [infinity]
 *                         |
 * CPodeSetMaxErrTestFails | Maximum number of error test failures
 *                         | in attempting one step.
 *                         | [7]
 *                         |
 * -----------------------------------------------------------------
 *                         |
 * CPodeSetMaxNonlinIters  | Maximum number of nonlinear solver
 *                         | iterations at one solution.
 *                         | [3]
 *                         |
 * CPodeSetMaxConvFails    | Maximum number of convergence failures
 *                         | allowed in attempting one step.
 *                         | [10]
 *                         |
 * CPodeSetNonlinConvCoef  | Coefficient in the nonlinear
 *                         | convergence test.
 *                         | [0.1]
 *                         |
 * -----------------------------------------------------------------
 *                         |
 * CPodeSetProjUpdateErrEst| toggles ON/OFF projection of the
 *                         | error estimate.
 *                         | [TRUE]
 *                         |
 * CPodeSetProjFrequency   | frequency with which the projection
 *                         | step is performed. A value of 1 
 *                         | indicates that the projection step
 *                         | will be performed at every step.
 *                         | A value of 0 will disable projection.
 *                         | [1]
 *                         |
 * CPodeSetProjTestCnstr   | if TRUE, the internal projection 
 *                         | function will be performed only if
 *                         | the constraint violation is larger
 *                         | than the prescribed tolerances. 
 *                         | Otherwise, the tolerances are ignored
 *                         | and projection is always performed.
 *                         | [FALSE]
 *                         | 
 * CPodeSetProjLsetupFreq  | frequency with which the linear
 *                         | solver setup function is called
 *                         | (i.e. frequency of constraint 
 *                         | Jacobian evaluations). The default
 *                         | value of 1 forces a Jacobian
 *                         | evaluation before every single 
 *                         | projection step.
 *                         | [1]
 *                         |
 * CPodeSetProjNonlinConvCoef | Coefficient in the nonlinear
 *                         | convergence test (for projection).
 *                         | [0.1]
 *                         |
 * -----------------------------------------------------------------
 *                         |
 * CPodeSetQuadErrCon      | are quadrature variables considered in
 *                         | the error control?
 *                         | If yes, set tolerances for quadrature
 *                         | integration. 
 *                         | [errconQ = FALSE]
 *                         | [no tolerances]
 *                         | 
 * -----------------------------------------------------------------
 *                         |
 * CPodeSetTolerances      | Changes the integration tolerances
 *                         | between calls to CPode.
 *                         | [set by CPodeInit/CPodeReInit]
 *                         |
 * ---------------------------------------------------------------- 
 *                         |
 * CPodeSetRootDirection   | Specifies the direction of zero
 *                         | crossings to be monitored
 *                         | [both directions]
 *                         |
 * -----------------------------------------------------------------
 * Return flag:
 *   CP_SUCCESS   if successful
 *   CP_MEM_NULL  if the CPODES memory is NULL
 *   CP_ILL_INPUT if an argument has an illegal value
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPodeSetErrHandlerFn(void *cpode_mem, CPErrHandlerFn ehfun, void *eh_data);
SUNDIALS_EXPORT int CPodeSetErrFile(void *cpode_mem, FILE *errfp);

SUNDIALS_EXPORT int CPodeSetEwtFn(void *cpode_mem, CPEwtFn efun, void *e_data);
SUNDIALS_EXPORT int CPodeSetMaxOrd(void *cpode_mem, int maxord);
SUNDIALS_EXPORT int CPodeSetMaxNumSteps(void *cpode_mem, long int mxsteps);
SUNDIALS_EXPORT int CPodeSetMaxHnilWarns(void *cpode_mem, int mxhnil);
SUNDIALS_EXPORT int CPodeSetStabLimDet(void *cpode_mem, booleantype stldet);
SUNDIALS_EXPORT int CPodeSetInitStep(void *cpode_mem, realtype hin);
SUNDIALS_EXPORT int CPodeSetMinStep(void *cpode_mem, realtype hmin);
SUNDIALS_EXPORT int CPodeSetMaxStep(void *cpode_mem, realtype hmax);
SUNDIALS_EXPORT int CPodeSetStopTime(void *cpode_mem, realtype tstop);
SUNDIALS_EXPORT int CPodeSetMaxErrTestFails(void *cpode_mem, int maxnef);

SUNDIALS_EXPORT int CPodeSetMaxNonlinIters(void *cpode_mem, int maxcor);
SUNDIALS_EXPORT int CPodeSetMaxConvFails(void *cpode_mem, int maxncf);
SUNDIALS_EXPORT int CPodeSetNonlinConvCoef(void *cpode_mem, realtype nlscoef);

SUNDIALS_EXPORT int CPodeSetProjUpdateErrEst(void *cpode_mem, booleantype proj_err);
SUNDIALS_EXPORT int CPodeSetProjFrequency(void *cpode_mem, long int proj_freq);
SUNDIALS_EXPORT int CPodeSetProjTestCnstr(void *cpode_mem, booleantype test_cnstr);
SUNDIALS_EXPORT int CPodeSetProjLsetupFreq(void *cpode_mem, long int proj_lset_freq);
SUNDIALS_EXPORT int CPodeSetProjNonlinConvCoef(void *cpode_mem, realtype prjcoef);

SUNDIALS_EXPORT int CPodeSetQuadErrCon(void *cpode_mem, booleantype errconQ, 
				       int tol_typeQ, realtype reltolQ, void *abstolQ);

SUNDIALS_EXPORT int CPodeSetTolerances(void *cpode_mem,
				       int tol_type, realtype reltol, void *abstol);

SUNDIALS_EXPORT int CPodeSetRootDirection(void *cpode_mem, int *rootdir);

/*
 * -----------------------------------------------------------------
 * Function : CPode
 * -----------------------------------------------------------------
 * CPode integrates the ODE over an interval in t.
 * If mode is CP_NORMAL, then the solver integrates from its
 * current internal t value to a point at or beyond tout, then
 * interpolates to t = tout and returns y(tout) in the user-
 * allocated vector yout. If mode is CP_ONE_STEP, then the solver
 * takes one internal time step and returns in yout the value of
 * y at the new internal time. In this case, tout is used only
 * during the first call to CPode to determine the direction of
 * integration and the rough scale of the t variable.  If mode is
 * CP_NORMAL_TSTOP or CP_ONE_STEP_TSTOP, then CPode returns the
 * solution at tstop if that comes sooner than tout or the end of
 * the next internal step, respectively.  In any case,
 * the time reached by the solver is placed in (*tret). The
 * user is responsible for allocating the memory for this value.
 *
 * cpode_mem is the pointer to CPODE memory returned by
 *           CPodeCreate.
 *
 * tout  is the next time at which a computed solution is desired.
 *
 * tret  is a pointer to a real location. CPode sets (*tret) to
 *       the time reached by the solver and returns yout=y(*tret)
 *       and ypout=y'(*tret).
 *
 * yout  is the computed solution vector. In CP_NORMAL mode with no
 *       errors and no roots found, yout=y(tout).
 *
 * ypout is the computed derivative vector.
 *
 * mode  is CP_NORMAL, CP_ONE_STEP, CP_NORMAL_TSTOP, or 
 *       CP_ONE_STEP_TSTOP. These four modes are described above.
 *
 * Here is a brief description of each return value:
 *
 * CP_SUCCESS:      CPode succeeded and no roots were found.
 *
 * CP_ROOT_RETURN:  CPode succeeded, and found one or more roots.
 *                  If nrtfn > 1, call CPodeGetRootInfo to see
 *                  which g_i were found to have a root at (*tret).
 *
 * CP_TSTOP_RETURN: CPode succeeded and returned at tstop.
 *
 * CP_MEM_NULL:     The cpode_mem argument was NULL.
 *
 * CP_NO_MALLOC:    cpode_mem was not allocated.
 *
 * CP_ILL_INPUT:    One of the inputs to CPode is illegal. This
 *                  includes the situation when a component of the
 *                  error weight vectors becomes < 0 during
 *                  internal time-stepping.  It also includes the
 *                  situation where a root of one of the root
 *                  functions was found both at t0 and very near t0.
 *                  The ILL_INPUT flag will also be returned if the
 *                  linear solver routine CP--- (called by the user
 *                  after calling CPodeCreate) failed to set one of
 *                  the linear solver-related fields in cpode_mem or
 *                  if the linear solver's init routine failed. In
 *                  any case, the user should see the printed
 *                  error message for more details.
 *
 * CP_TOO_MUCH_WORK: The solver took mxstep internal steps but
 *                  could not reach tout. The default value for
 *                  mxstep is MXSTEP_DEFAULT = 500.
 *
 * CP_TOO_MUCH_ACC: The solver could not satisfy the accuracy
 *                  demanded by the user for some internal step.
 *
 * CP_ERR_FAILURE:  Error test failures occurred too many times
 *                  (= MXNEF = 7) during one internal time step or
 *                  occurred with |h| = hmin.
 *
 * CP_CONV_FAILURE: Convergence test failures occurred too many
 *                  times (= MXNCF = 10) during one internal time
 *                  step or occurred with |h| = hmin.
 *
 * CP_LINIT_FAIL:   The linear solver's initialization function 
 *                  failed.
 *
 * CP_LSETUP_FAIL:  The linear solver's setup routine failed in an
 *                  unrecoverable manner.
 *
 * CP_LSOLVE_FAIL:  The linear solver's solve routine failed in an
 *                  unrecoverable manner.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPode(void *cpode_mem, realtype tout, realtype *tret, 
			  N_Vector yout, N_Vector ypout, int mode);

/*
 * -----------------------------------------------------------------
 * Function : CPodeGetDky
 * -----------------------------------------------------------------
 * CPodeGetDky computes the kth derivative of the y function at
 * time t, where tn-hu <= t <= tn, tn denotes the current
 * internal time reached, and hu is the last internal step size
 * successfully used by the solver. The user may request
 * k=0, 1, ..., qu, where qu is the order last used. The
 * derivative vector is returned in dky. This vector must be
 * allocated by the caller. It is only legal to call this
 * function after a successful return from CPode.
 *
 * cpode_mem is the pointer to CPODES memory returned by
 *           CPodeCreate.
 *
 * t   is the time at which the kth derivative of y is evaluated.
 *     The legal range for t is [tn-hu,tn] as described above.
 *
 * k   is the order of the derivative of y to be computed. The
 *     legal range for k is [0,qu] as described above.
 *
 * dky is the output derivative vector [((d/dy)^k)y](t).
 *
 * The return value for CPodeGetDky is one of:
 *
 *   CP_SUCCESS:  CPodeGetDky succeeded.
 *
 *   CP_BAD_K:    k is not in the range 0, 1, ..., qu.
 *
 *   CP_BAD_T:    t is not in the interval [tn-hu,tn].
 *
 *   CP_BAD_DKY:  The dky argument was NULL.
 *
 *   CP_MEM_NULL: The cpode_mem argument was NULL.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPodeGetDky(void *cpode_mem, realtype t, int k, N_Vector dky);

/*
 * -----------------------------------------------------------------
 * Quadrature integration solution extraction routines
 * -----------------------------------------------------------------
 * The following functions can be called to obtain the quadrature
 * variables (or derivatives of them) after a successful integration
 * step.  If quadratures were not computed, they return CP_NO_QUAD.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPodeGetQuad(void *cpode_mem, realtype t, N_Vector yQout);
SUNDIALS_EXPORT int CPodeGetQuadDky(void *cpode_mem, realtype t, int k, N_Vector dky);


/*
 * -----------------------------------------------------------------
 * IC calculation optional output extraction functions
 * -----------------------------------------------------------------
 * CPodeGetConsistentIC returns the consistent initial conditions
 *       computed by CPodeCalcIC
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPodeGetConsistentIC(void *cpode_mem, N_Vector yy0, N_Vector yp0);

/*
 * -----------------------------------------------------------------
 * Integrator optional output extraction functions
 * -----------------------------------------------------------------
 * The following functions can be called to get optional outputs
 * and statistics related to the main integrator.
 * -----------------------------------------------------------------
 * CPodeGetWorkSpace returns the CPODES real and integer workspaces
 * CPodeGetNumSteps returns the cumulative number of internal
 *    steps taken by the solver
 * CPodeGetNumFctEvals returns the number of calls to the user's
 *    fun function
 * CPodeGetNumLinSolvSetups returns the number of calls made to
 *    the linear solver's setup routine
 * CPodeGetNumErrTestFails returns the number of local error test
 *    failures that have occured
 * CPodeGetLastOrder returns the order used during the last
 *    internal step
 * CPodeGetCurrentOrder returns the order to be used on the next
 *    internal step
 * CPodeGetNumStabLimOrderReds returns the number of order 
 *    reductions due to stability limit detection
 * CPodeGetActualInitStep returns the actual initial step size
 *    used by CPODES
 * CPodeGetLastStep returns the step size for the last internal step
 * CPodeGetCurrentStep returns the step size to be attempted on
 *    the next internal step
 * CPodeGetCurrentTime returns the current internal time reached
 *    by the solver
 * CPodeGetTolScaleFactor returns a suggested factor by which the
 *    user's tolerances should be scaled when too much accuracy has 
 *    been requested for some internal step
 * CPodeGetErrWeights returns the current error weight vector.
 *    The user must allocate space for eweight.
 * CPodeGetEstLocalErrors returns the vector of estimated local
 *    errors. The user must allocate space for ele.
 * CPodeGetNumGEvals returns the number of calls to the user's
 *    gfun function (for rootfinding)
 * CPodeGetRootInfo returns the indices for which g_i was found to 
 *    have a root. The user must allocate space for rootsfound. 
 *    For i = 0 ... nrtfn-1, rootsfound[i] = 1 if g_i has a root, 
 *    and = 0 if not.
 * CPodeGetIntegratorStats retruns most of the optional outputs as
 *    a group.
 * CPodeGet* return values:
 *   CP_SUCCESS   if succesful
 *   CP_MEM_NULL  if the CPODES memory was NULL
 *   CP_NO_SLDET  if stability limit was not turned on
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPodeGetWorkSpace(void *cpode_mem, long int *lenrw, long int *leniw);
SUNDIALS_EXPORT int CPodeGetNumSteps(void *cpode_mem, long int *nsteps);
SUNDIALS_EXPORT int CPodeGetNumFctEvals(void *cpode_mem, long int *nfevals);
SUNDIALS_EXPORT int CPodeGetNumLinSolvSetups(void *cpode_mem, long int *nlinsetups);
SUNDIALS_EXPORT int CPodeGetNumErrTestFails(void *cpode_mem, long int *netfails);
SUNDIALS_EXPORT int CPodeGetLastOrder(void *cpode_mem, int *qlast);
SUNDIALS_EXPORT int CPodeGetCurrentOrder(void *cpode_mem, int *qcur);
SUNDIALS_EXPORT int CPodeGetNumStabLimOrderReds(void *cpode_mem, long int *nslred);
SUNDIALS_EXPORT int CPodeGetActualInitStep(void *cpode_mem, realtype *hinused);
SUNDIALS_EXPORT int CPodeGetLastStep(void *cpode_mem, realtype *hlast);
SUNDIALS_EXPORT int CPodeGetCurrentStep(void *cpode_mem, realtype *hcur);
SUNDIALS_EXPORT int CPodeGetCurrentTime(void *cpode_mem, realtype *tcur);
SUNDIALS_EXPORT int CPodeGetTolScaleFactor(void *cpode_mem, realtype *tolsfac);
SUNDIALS_EXPORT int CPodeGetErrWeights(void *cpode_mem, N_Vector eweight);
SUNDIALS_EXPORT int CPodeGetEstLocalErrors(void *cpode_mem, N_Vector ele);
SUNDIALS_EXPORT int CPodeGetNumGEvals(void *cpode_mem, long int *ngevals);
SUNDIALS_EXPORT int CPodeGetRootInfo(void *cpode_mem, int *rootsfound);
SUNDIALS_EXPORT int CPodeGetIntegratorStats(void *cpode_mem, long int *nsteps,
					    long int *nfevals, long int *nlinsetups,
					    long int *netfails, int *qlast,
					    int *qcur, realtype *hinused, realtype *hlast,
					    realtype *hcur, realtype *tcur);

/*
 * -----------------------------------------------------------------
 * Nonlinear solver optional output extraction functions
 * -----------------------------------------------------------------
 * The following functions can be called to get optional outputs
 * and statistics related to the nonlinear solver.
 * -----------------------------------------------------------------
 * CPodeGetNumNonlinSolvIters returns the number of nonlinear
 *    solver iterations performed.
 * CPodeGetNumNonlinSolvConvFails returns the number of nonlinear
 *    convergence failures.
 * CPodeGetNonlinSolvStats returns the nonlinear solver optional 
 *    outputs in a group.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPodeGetNumNonlinSolvIters(void *cpode_mem, long int *nniters);
SUNDIALS_EXPORT int CPodeGetNumNonlinSolvConvFails(void *cpode_mem, long int *nncfails);
SUNDIALS_EXPORT int CPodeGetNonlinSolvStats(void *cpode_mem, long int *nniters,
					    long int *nncfails);

  
/*
 * -----------------------------------------------------------------
 * Projection optional output extraction functions
 * -----------------------------------------------------------------
 * The following functions can be called to get optional outputs
 * and statistics related to the projection step.
 * -----------------------------------------------------------------
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPodeGetProjNumProj(void *cpode_mem, long int *nproj);
SUNDIALS_EXPORT int CPodeGetProjNumCnstrEvals(void *cpode_mem, long int *nce);
SUNDIALS_EXPORT int CPodeGetProjNumLinSolvSetups(void *cpode_mem, long int *nsetupsP);
SUNDIALS_EXPORT int CPodeGetProjNumFailures(void *cpode_mem, long int *nprf);
SUNDIALS_EXPORT int CPodeGetProjStats(void *cpode_mem, long int *nproj,
				      long int *nce, long int *nsetupsP,
				      long int *nprf);

/*
 * -----------------------------------------------------------------
 * Quadrature integration optional output extraction functions
 * -----------------------------------------------------------------
 * The following functions can be called to get optional outputs
 * and statistics related to the integration of quadratures.
 * -----------------------------------------------------------------
 * CPodeGetQuadNumFunEvals returns the number of calls to the
 *    user function qfun defining the integrand
 * CPodeGetQuadErrWeights returns the vector of error weights for
 *    the quadrature variables. The user must allocate space for ewtQ.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPodeGetQuadNumFunEvals(void *cpode_mem, long int *nqevals);
SUNDIALS_EXPORT int CPodeGetQuadErrWeights(void *cpode_mem, N_Vector eQweight);

/*
 * -----------------------------------------------------------------
 * The following function returns the name of the constant 
 * associated with a CPODES return flag
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT char *CPodeGetReturnFlagName(int flag);

/*
 * -----------------------------------------------------------------
 * Function : CPodeFree
 * -----------------------------------------------------------------
 * CPodeFree frees the problem memory cpode_mem allocated by
 * CPodeCreate and CPodeInit.  Its only argument is the pointer
 * cpode_mem returned by CPodeCreate.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void CPodeFree(void **cpode_mem);

/*
 * -----------------------------------------------------------------
 * Function : CPodeQuadFree
 * -----------------------------------------------------------------
 * CPodeQuadFree frees the problem memory in cpode_mem allocated
 * for quadrature integration. Its only argument is the pointer
 * cpode_mem returned by CPodeCreate.
 * Note that CPodeQuadFree is called by CPodeFree.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void CPodeQuadFree(void *cpode_mem);

#ifdef __cplusplus
}
#endif

#endif
