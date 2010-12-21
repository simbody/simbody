/*
 * -----------------------------------------------------------------
 * $Revision: 1.7 $
 * $Date: 2007/10/26 21:51:29 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Private header file for CPODES.
 * -----------------------------------------------------------------
 */

#ifndef _CPODES_PRIVATE_H
#define _CPODES_PRIVATE_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <stdarg.h>
#include "cpodes_impl.h"

/*
 * =================================================================
 *   C P O D E S    I N T E R N A L    C O N S T A N T S
 * =================================================================
 */

/*
 * Default values for algorithmic constants
 * ---------------------------------------------------------------
 *
 * Step size and number of steps
 *
 *    HMIN_DEFAULT      default value for minimum step size
 *    HMAX_INV_DEFAULT  default value for inverse of maximum step size
 *    MXHNIL_DEFAULT    default value for mxhnil
 *    MXSTEP_DEFAULT    default value for maximum number of steps
 *
 * Nonlinear solver
 *    
 *    MXNCF         max no. of convergence failures during one step try
 *    NLS_MAXCOR    maximum no. of corrector iterations for the nonlinear solver
 *    NLS_TEST_COEF constant in error test (used in the test quantity tq[4])
 *    NLS_CRDOWN    constant used in the estimation of the convergence rate (crate)
 *    NLS_RDIV      declare divergence if ratio del/delp > NLS_RDIV
 *    NLS_MSBLS     max no. of steps between lsetup calls
 *    DGMAX         when nls_type=NEWTON, |gamma/gammap-1| > DGMAX => call lsetup
 *    ETACF         maximum step size decrease on convergence failure
 *
 * Error test
 *
 *    MXNEF       max no. of error test failures during one step try
 *    MXNEF1      max no. of error test failures before forcing a reduction of order
 *    SMALL_NEF   if an error failure occurs and SMALL_NEF <= nef <= MXNEF1, then
 *                reset eta =  MIN(eta, ETAMXF)
 *    LONG_WAIT   number of steps to wait before considering an order change when
 *                q==1 and MXNEF1 error test failures have occurred   
 *    SMALL_NST   nst > SMALL_NST => use ETAMX3 
 *
 * Projection
 *
 *    MXNPF         max no. of projection failures during one step try
 *    PRJ_MAXCOR    maximum no. of iterations for the nonlinear projection
 *    PRJ_TEST_COEF constant used in projection tolerances
 *    PRJ_CRDOWN    constant used in the estimation of the convergence rate (crateP)
 *    PRJ_RDIV      declare divergence if ratio del/delp > PRJ_RDIV
 *    PRJ_MSBLS     max no. of steps between lsetupP calls
 *
 * 
 * Other
 *
 *    FUZZ_FACTOR   factor used in defining an infinitesimal time interval
 */

#define HMIN_DEFAULT     RCONST(0.0)
#define HMAX_INV_DEFAULT RCONST(0.0)
#define MXHNIL_DEFAULT   10
#define MXSTEP_DEFAULT   500

#define MXNCF         10
#define NLS_MAXCOR    3
#define NLS_TEST_COEF RCONST(0.1)
#define NLS_CRDOWN    RCONST(0.3)
#define NLS_RDIV      RCONST(2.0)
#define NLS_MSBLS     20
#define DGMAX         RCONST(0.3)
#define ETACF         RCONST(0.25)

#define MXNEF         7
#define MXNEF1        3
#define SMALL_NEF     2
#define LONG_WAIT     10
#define SMALL_NST     10

#define MXNPF         10
#define PRJ_MAXCOR    3
#define PRJ_TEST_COEF RCONST(0.1)
#define PRJ_CRDOWN    RCONST(0.3)
#define PRJ_RDIV      RCONST(2.0)
#define PRJ_MSBLS     1

#define FUZZ_FACTOR RCONST(100.0)


/* 
 * Control constants for communication between main integrator and
 * lower level functions in cpStep
 * ---------------------------------------------------------------
 *
 * cpNls input nflag values:
 *    FIRST_CALL
 *    PREV_CONV_FAIL
 *    PREV_PROJ_FAIL
 *    PREV_ERR_FAIL
 *    
 * cpNls return values: 
 *    CP_SUCCESS,
 *    CP_LSETUP_FAIL, CP_LSOLVE_FAIL, CP_ODEFUNC_FAIL,
 *    CP_CONV_FAILURE, CP_REPTD_ODEFUNC_ERR,
 *    PREDICT_AGAIN
 *
 * cpDoProjection return values: 
 *    CP_SUCCESS,
 *    CP_PLSETUP_FAIL, CP_PLSOLVE_FAIL, CP_CNSTRFUNC_FAIL, CP_PROJFUNC_FAIL
 *    CP_PROJ_FAILURE, CP_REPTD_CNSTRFUNC_ERR, CP_REPTD_PROJFUNC_ERR,
 *    PREDICT_AGAIN
 * 
 * cpDoErrorTest return values: 
 *    CP_SUCCESS,
 *    CP_ERR_FAILURE,
 *    PREDICT_AGAIN
 * 
 * cpRcheck* return values:
 *    CP_RTFUNC_FAIL,
 *    RTFOUND,
 *    CP_SUCCESS
 */

#define PREDICT_AGAIN    +3

#define FIRST_CALL       +101
#define PREV_CONV_FAIL   +102
#define PREV_PROJ_FAIL   +103
#define PREV_ERR_FAIL    +104

#define CONV_FAIL        +110 
#define ODEFUNC_RECVR    +111
#define CNSTRFUNC_RECVR  +112
#define PROJFUNC_RECVR   +113
#define QUADFUNC_RECVR   +114

#define RTFOUND          +1

/*
 * CPODES Private Constants
 * ---------------------------------------------------------------
 */
  
#define ZERO    RCONST(0.0)
#define TINY    RCONST(1.0e-10)
#define PT001   RCONST(0.001)
#define PT01    RCONST(0.01)
#define PT1     RCONST(0.1)
#define PT2     RCONST(0.2)
#define PT25    RCONST(0.25)
#define HALF    RCONST(0.5)
#define ONE     RCONST(1.0)
#define ONEPSM  RCONST(1.000001)
#define TWO     RCONST(2.0)
#define THREE   RCONST(3.0)
#define FOUR    RCONST(4.0)
#define FIVE    RCONST(5.0)
#define TWELVE  RCONST(12.0)
#define HUNDRED RCONST(100.0)

/*
 * =================================================================
 *   C P O D E S   I N T E R N A L   F U N C T I O N S
 *         S H A R E D   A M O N G   M O D U L E S
 * =================================================================
 */

/* Nonlinear solver function */
int cpNls(CPodeMem cp_mem, int nflag, realtype saved_t, int *ncfPtr);

/* Projection step function */
int cpDoProjection(CPodeMem cp_mem, realtype saved_t, int *npfPtr);

/* Internal error weight computation */
int cpEwtSet(N_Vector ycur, N_Vector weight, void *edata);

/* Error return handler */
void cpProcessError(CPodeMem cp_mem, 
		    int error_code, const char *module, const char *fname, 
		    const char *msgfmt, ...);

/* Functions acting on the Nordsieck history array */ 
void cpRestore(CPodeMem cp_mem, realtype saved_t);
void cpRescale(CPodeMem cp_mem);

/* Function evaluating solution at given time */
int cpGetSolution(void *cpode_mem, realtype t, N_Vector yret, N_Vector ypret);

/* Root finding functions */
int cpRcheck1(CPodeMem cp_mem);
int cpRcheck2(CPodeMem cp_mem);
int cpRcheck3(CPodeMem cp_mem);
void cpRootFree(CPodeMem cp_mem);

/*
 * =================================================================
 *   M A C R O
 * =================================================================
 */

#define loop for(;;)

/*
 * =================================================================
 *   C V O D E    E R R O R    M E S S A G E S
 * =================================================================
 */

#if defined(SUNDIALS_EXTENDED_PRECISION)

#define MSG_TIME      "t = %Lg"
#define MSG_TIME_H    "t = %Lg and h = %Lg"
#define MSG_TIME_INT  "t = %Lg is not between tcur - hu = %Lg and tcur = %Lg."
#define MSG_TIME_TOUT "tout = %Lg"

#elif defined(SUNDIALS_DOUBLE_PRECISION)

#define MSG_TIME      "t = %lg"
#define MSG_TIME_H    "t = %lg and h = %lg"
#define MSG_TIME_INT  "t = %lg is not between tcur - hu = %lg and tcur = %lg."
#define MSG_TIME_TOUT "tout = %lg"

#else

#define MSG_TIME      "t = %g"
#define MSG_TIME_H    "t = %g and h = %g"
#define MSG_TIME_INT  "t = %g is not between tcur - hu = %g and tcur = %g."
#define MSG_TIME_TOUT "tout = %g"

#endif

/* Initialization and I/O error messages */

#define MSGCP_NO_MEM "cpode_mem = NULL illegal."
#define MSGCP_CPMEM_FAIL "Allocation of cpode_mem failed."
#define MSGCP_MEM_FAIL "A memory request failed."
#define MSGCP_BAD_ODE  "Illegal value for ode. The legal values are CP_EXPL and CP_IMPL."
#define MSGCP_BAD_LMM  "Illegal value for lmm_type. The legal values are CP_ADAMS and CP_BDF."
#define MSGCP_BAD_NLS  "Illegal value for nls_type. The legal values are CP_FUNCTIONAL and CP_NEWTON."
#define MSGCP_BAD_ODE_NLS  "Illegal combination ode_type=CP_IMPL nls_type=CP_FUNCTIONAL."
#define MSGCP_BAD_ITOL "Illegal value for tol_type. The legal values are CP_SS, CP_SV, and CP_WF."
#define MSGCP_NO_MALLOC "Attempt to call before CPodeInit."
#define MSGCP_BAD_PROJ "Illegal value for proj_type."
#define MSGCP_BAD_NORM "Illegal value for proj_norm."
#define MSGCP_BAD_CNSTR "Illegal value for cnstr_type."
#define MSGCP_CTOL_NULL "ctol = NULL illegal."
#define MSGCP_NEG_MAXORD "maxord <= 0 illegal."
#define MSGCP_BAD_MAXORD  "Illegal attempt to increase maximum method order."
#define MSGCP_NEG_MXSTEPS "mxsteps < 0 illegal."
#define MSGCP_SET_SLDET  "Attempt to use stability limit detection with the CP_ADAMS method illegal."
#define MSGCP_NEG_HMIN "hmin < 0 illegal."
#define MSGCP_NEG_HMAX "hmax < 0 illegal."
#define MSGCP_BAD_HMIN_HMAX "Inconsistent step size limits: hmin > hmax."
#define MSGCP_BAD_RELTOL "reltol < 0 illegal."
#define MSGCP_BAD_ABSTOL "abstol has negative component(s) (illegal)."
#define MSGCP_NULL_ABSTOL "abstol = NULL illegal."
#define MSGCP_NULL_Y0 "y0 = NULL illegal."
#define MSGCP_NULL_YP0 "yp0 = NULL illegal for ode_type=CP_IMPL."
#define MSGCP_NULL_F "fun = NULL illegal."
#define MSGCP_NULL_G "gfun = NULL illegal."
#define MSGCP_NULL_C "cfun = NULL illegal."
#define MSGCP_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGCP_BAD_K "Illegal value for k."
#define MSGCP_NULL_DKY "dky = NULL illegal."
#define MSGCP_BAD_T "Illegal value for t." MSG_TIME_INT
#define MSGCP_BAD_FREQ "proj_freq < 0 illegal."
#define MSGCP_BAD_LSFREQ "lset_freq <= 0 illegal."
#define MSGCP_TOO_LATE "CPodeGetConsistentIC can only be called before the first call to CPode."


/* CPode error messages */

#define MSGCP_LSOLVE_NULL "The linear solver's solve routine is NULL."
#define MSGCP_YOUT_NULL "yout = NULL illegal."
#define MSGCP_YPOUT_NULL "ypout = NULL illegal."
#define MSGCP_TRET_NULL "tret = NULL illegal."
#define MSGCP_BAD_EWT "Initial ewt has component(s) equal to zero (illegal)."
#define MSGCP_EWT_NOW_BAD "At " MSG_TIME ", a component of ewt has become <= 0."
#define MSGCP_BAD_MODE "Illegal value for mode."
#define MSGCP_BAD_H0 "h0 and tout - t0 inconsistent."
#define MSGCP_BAD_TOUT "Trouble interpolating at " MSG_TIME_TOUT ". tout too far back in direction of integration"
#define MSGCP_NO_EFUN "tol_type = CP_WF but no error weight function was provided."
#define MSGCP_NO_PFUN "proj_type = CP_PROJ_USER but no projection function was provided."
#define MSGCP_NO_CFUN "proj_type = CP_PROJ_INTERNAL but no constraint function was provided."
#define MSGCP_NO_TSTOP "mode = CP_NORMAL_TSTOP or mode = CP_ONE_STEP_TSTOP but tstop was not set."
#define MSGCP_EWT_FAIL "The user-provided EwtSet function failed."
#define MSGCP_EWT_NOW_FAIL "At " MSG_TIME ", the user-provides EwtSet function failed."
#define MSGCP_LINIT_FAIL "The linear solver's init routine failed."
#define MSGCP_HNIL_DONE "The above warning has been issued mxhnil times and will not be issued again for this problem."
#define MSGCP_TOO_CLOSE "tout too close to t0 to start integration."
#define MSGCP_MAX_STEPS "At " MSG_TIME ", mxstep steps taken before reaching tout."
#define MSGCP_TOO_MUCH_ACC "At " MSG_TIME ", too much accuracy requested."
#define MSGCP_HNIL "Internal " MSG_TIME_H " are such that t + h = t on the next step. The solver will continue anyway."
#define MSGCP_ERR_FAILS "At " MSG_TIME_H ", the error test failed repeatedly or with |h| = hmin."
#define MSGCP_CONV_FAILS "At " MSG_TIME_H ", the corrector convergence test failed repeatedly or with |h| = hmin."
#define MSGCP_SETUP_FAILED "At " MSG_TIME ", the setup routine failed in an unrecoverable manner."
#define MSGCP_SOLVE_FAILED "At " MSG_TIME ", the solve routine failed in an unrecoverable manner."
#define MSGCP_ODEFUNC_FAILED "At " MSG_TIME ", the ODE funciton failed in an unrecoverable manner."
#define MSGCP_ODEFUNC_UNREC "At " MSG_TIME ", the ODE function failed in a recoverable manner, but no recovery is possible."
#define MSGCP_ODEFUNC_REPTD "At " MSG_TIME "repeated recoverable ODE function errors."
#define MSGCP_ODEFUNC_FIRST "The ODE function failed at the first call."
#define MSGCP_RTFUNC_FAILED "At " MSG_TIME ", the rootfinding routine failed in an unrecoverable manner."
#define MSGCP_BAD_TSTOP "tstop is behind current " MSG_TIME "in the direction of integration."
#define MSGCP_NO_ROOT "Rootfinding was not initialized."

/* Projection error messages */

#define MSGCP_PLSOLVE_NULL "The projection linear solver's solve function is NULL."
#define MSGCP_PLINIT_FAIL "The projection linear solver's init function failed."
#define MSGCP_PLSETUP_FAILED "At " MSG_TIME ", the projection linear solver's setup function failed in an unrecoverable manner."
#define MSGCP_PLSOLVE_FAILED "At " MSG_TIME ", the projection linear solver's solve function failed in an unrecoverable manner."
#define MSGCP_PROJ_FAILS "At " MSG_TIME_H ", the projection failed repeatedly or with |h| = hmin."
#define MSGCP_CNSTRFUNC_FAILED "At " MSG_TIME ", the constraint function failed in an unrecoverable manner."
#define MSGCP_CNSTRFUNC_REPTD "At " MSG_TIME "repeated recoverable constraint function errors."
#define MSGCP_PROJFUNC_FAILED "At " MSG_TIME ", the projection function failed in an unrecoverable manner."
#define MSGCP_PROJFUNC_REPTD "At " MSG_TIME "repeated recoverable projection function errors."

/* Quadrature integration error messages */

#define MSGCP_NO_QUAD  "Illegal attempt to call before calling CPodeQuadInit."
#define MSGCP_BAD_ITOLQ "Illegal value for tol_typeQ. The legal values are CP_SS and CP_SV."
#define MSGCP_NULL_ABSTOLQ "abstolQ = NULL illegal."
#define MSGCP_BAD_RELTOLQ "reltolQ < 0 illegal."
#define MSGCP_BAD_ABSTOLQ "abstolQ has negative component(s) (illegal)."  
#define MSGCP_BAD_EWTQ "Initial ewtQ has component(s) equal to zero (illegal)."
#define MSGCP_EWTQ_NOW_BAD "At " MSG_TIME ", a component of ewtQ has become <= 0."
#define MSGCP_QUADFUNC_FAILED "At " MSG_TIME ", the quadrature function failed in an unrecoverable manner."
#define MSGCP_QUADFUNC_UNREC "At " MSG_TIME ", the quadrature function failed in a recoverable manner, but no recovery is possible."
#define MSGCP_QUADFUNC_REPTD "At " MSG_TIME "repeated recoverable quadrature function errors."
#define MSGCP_QUADFUNC_FIRST "The quadrature function failed at the first call."

/* IC calculation error messages */
#define MSGCP_IC_NO_WORK "Nothing to do for an explicit-form ODE if no constraints are defined."
#define MSGCP_IC_CNSTRFUNC_FIRST "The constraint function failed at the first call."
#define MSGCP_IC_CNSTRFUNC_FAILED "The constraint function failed in an unrecoverable manner."
#define MSGCP_IC_CNSTRFUNC_REPTD "Repeated recoverable constraint function errors."
#define MSGCP_IC_PROJFUNC_FAILED "The projection function failed in an unrecoverable manner."
#define MSGCP_IC_PLSETUP_FAILED "The projection linear solver's setup function failed in an unrecoverable manner."
#define MSGCP_IC_PLSOLVE_FAILED "The projection linear solver's solve function failed in an unrecoverable manner."
#define MSGCP_IC_NO_RECOVERY "The linear solver's solve function failed recoverably, but the Jacobian data is already current."
#define MSGCP_IC_PROJ_FAILS "WTF?"
#define MSGCP_IC_EWT_BAD "A component of ewt has become <= 0."
#define MSGCP_IC_EWT_FAIL "The user-provided EwtSet function failed."
#define MSGCP_IC_ODEFUNC_FIRST "The ODE function failed at the first call."
#define MSGCP_IC_ODEFUNC_FAILED "The ODE funciton failed in an unrecoverable manner."

#ifdef __cplusplus
}
#endif

#endif
