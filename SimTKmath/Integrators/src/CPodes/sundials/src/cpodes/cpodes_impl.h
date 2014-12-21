/*
 * -----------------------------------------------------------------
 * $Revision: 1.5 $
 * $Date: 2007/10/26 21:51:29 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Implementation header file for the main CPODES integrator.
 * -----------------------------------------------------------------
 */

#ifndef _CPODES_IMPL_H
#define _CPODES_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <cpodes/cpodes.h>

/*
 * =================================================================
 *   M A I N    I N T E G R A T O R    M E M O R Y    B L O C K
 * =================================================================
 */

/* Basic LMM constants */

#define ADAMS_Q_MAX 12     /* max value of q for ADAMS        */
#define BDF_Q_MAX    5     /* max value of q for BDF          */
#define Q_MAX  ADAMS_Q_MAX /* max value of q for either LMM   */
#define L_MAX  (Q_MAX+1)   /* max value of L for either LMM   */
#define NUM_TESTS    5     /* number of error test quantities */

/*
 * -----------------------------------------------------------------
 * Types : struct CPodeMemRec, CPodeMem
 * -----------------------------------------------------------------
 * The type CPodeMem is type pointer to struct CPodeMemRec.
 * This structure contains fields to keep track of problem state.
 * -----------------------------------------------------------------
 */

typedef struct CPodeMemRec {

  realtype cp_uround;    /* machine unit roundoff */

  /*-------------------------- 
    Problem Specification Data 
    --------------------------*/

  int cp_ode_type;             /* ODE type: ode = CP_IMPL or CP_EXPL          */
  CPResFn cp_fi;               /* F(t,y'(t),y(t)) = 0                         */
  CPRhsFn cp_fe;               /* y' = f(t,y(t))                              */
  void *cp_f_data;             /* user pointer passed to f                    */

  int cp_lmm_type;             /* lmm_type = CP_ADAMS or CP_BDF               */
  int cp_nls_type;             /* nls_type = CP_FUNCTIONAL or CP_NEWTON       */
  int cp_tol_type;             /* tol_type = CP_SS or CP_SV or CP_WF          */

  realtype cp_reltol;          /* relative tolerance                          */
  realtype cp_Sabstol;         /* scalar absolute tolerance                   */
  N_Vector cp_Vabstol;         /* vector absolute tolerance                   */
  CPEwtFn cp_efun;             /* function to set ewt                         */
  void *cp_e_data;             /* user pointer passed to efun                 */

  /*-----------------------
    Nordsieck History Array 
    -----------------------*/

  /* 
   * Nordsieck array, of size N x (q+1).
   * zn[j] is a vector of length N (j=0,...,q)
   * zn[j] = [1/factorial(j)] * h^j * (jth derivative of the interpolating polynomial)
   */

  N_Vector cp_zn[L_MAX];

  /*------------------------------------------------
    Other vectors of same length as the state vector 
    ------------------------------------------------*/

  N_Vector cp_ewt;             /* error weight vector.                        */
  N_Vector cp_y;               /* work space for y vector (= yout)            */
  N_Vector cp_yp;              /* work space for yp vector (= ypout)          */
  N_Vector cp_acor;            /* accumulated correction acor = y_n(m)-y_n(0) */
  N_Vector cp_tempv;           /* temporary storage vector                    */
  N_Vector cp_ftemp;           /* temporary storage vector                    */

  /*---------
    Step Data 
    ---------*/  

  int cp_q;                    /* current order                               */
  int cp_qprime;               /* order to be used on the next step           */ 
  int cp_next_q;               /* order to be used on the next step           */
  int cp_qwait;                /* no. of steps to wait before a change in q   */
  int cp_L;                    /* L = q + 1                                   */

  realtype cp_hin;             /* initial step size                           */
  realtype cp_h;               /* current step size                           */
  realtype cp_hprime;          /* step size to be used on the next step       */ 
  realtype cp_next_h;          /* step size to be used on the next step       */ 
  realtype cp_eta;             /* eta = hprime / h                            */
  realtype cp_hscale;          /* value of h used in zn                       */
  realtype cp_tn;              /* current internal value of t                 */
  realtype cp_tretlast;        /* value of tret last returned by CPode        */

  realtype cp_tau[L_MAX+1];    /* array of previous q+1 successful step sizes */
  realtype cp_tq[NUM_TESTS+1]; /* array of test quantities                    */
  realtype cp_l[L_MAX];        /* coefficients of Lambda(x) (degree q poly)   */
  realtype cp_p[L_MAX];        /* coefficients of Phi(x) (degree q poly)      */

  realtype cp_rl1;             /* the scalar 1/l[1]                           */
  realtype cp_gamma;           /* gamma = h * rl1                             */
  realtype cp_gammap;          /* gamma at the last setup call                */
  realtype cp_gamrat;          /* gamma / gammap                              */

  realtype cp_crate;           /* estimated corrector convergence rate        */
  realtype cp_acnrm;           /* ||acor||_wrms                               */
  realtype cp_nlscoef;         /* coeficient in nonlinear convergence test    */
  int  cp_mnewt;               /* Newton iteration counter                    */

  /*---------------
    Projection Data
    ---------------*/

  int cp_proj_type;            /* user vs. internal projection algorithm      */
  booleantype cp_proj_enabled; /* is projection enabled?                      */
  booleantype cp_applyProj;    /* was projection performed at current step?   */

  long int cp_proj_freq;       /* projection frequency                        */
  long int cp_nstlprj;         /* step number of last projection              */

  int cp_proj_norm;            /* type of projection norm (L2 or WRMS)        */
  int cp_cnstr_type;           /* type of constraints (lin. or nonlin.)       */
  CPCnstrFn cp_cfun;           /* 0  = c(t,y(t))                              */
  void *cp_c_data;             /* user pointer passed to cfun                 */

  CPProjFn cp_pfun;            /* function to perform projection              */
  void *cp_p_data;             /* user pointer passed to pfun                 */

  realtype cp_prjcoef;         /* coefficient in projection convergence test  */

  N_Vector cp_acorP;           /* projection correction (length N)            */
  N_Vector cp_errP;            /* projected error estimate (length N)         */
  N_Vector cp_yC;              /* saved corrected state (length N)            */

  long int cp_lsetupP_freq;    /* frequency of cnstr. Jacobian evaluation     */
  long int cp_nstlsetP;        /* step number of last lsetupP call            */

  realtype cp_crateP;          /* estimated conv. rate in nonlin. projection  */
  int cp_maxcorP;              /* maximum number of nonlin. proj. iterations  */

  booleantype cp_project_err;  /* should we project the error estimate?       */

  booleantype cp_test_cnstr;   /* test constraint norm before projection?     */
  N_Vector cp_ctol;            /* vector of constraint "absolute tolerances"  */

  N_Vector cp_ctemp;           /* temporary vectors (length M)                */
  N_Vector cp_tempvP1;
  N_Vector cp_tempvP2;
    
  booleantype cp_first_proj;   /* is this the first time we project?          */

  long int cp_nproj;           /* number of projection steps performed        */
  long int cp_nprf;            /* number of projection failures               */
  long int cp_nce;             /* number of calls to cfun                     */
  long int cp_nsetupsP;        /* number of calls to lsetupP                  */       
  int cp_maxnpf;               /* maximum number of projection failures       */

  /*-----------------------
    Quadrature Related Data 
    -----------------------*/

  booleantype cp_quadr;        /* are we integrating quadratures?             */

  CPQuadFn cp_qfun;            /* function defining the integrand             */
  void *cp_q_data;             /* user pointer passed to fQ                   */
  int cp_tol_typeQ;            /* tol_typeQ = CP_SS or CP_SV                  */
  booleantype cp_errconQ;      /* are quadratures included in error test?     */

  realtype cp_reltolQ;         /* relative tolerance for quadratures          */
  realtype cp_SabstolQ;        /* scalar absolute tolerance for quadratures   */
  N_Vector cp_VabstolQ;        /* vector absolute tolerance for quadratures   */

  N_Vector cp_znQ[L_MAX];      /* Nordsieck arrays for sensitivities          */
  N_Vector cp_ewtQ;            /* error weight vector for quadratures         */
  N_Vector cp_yQ;              /* Unlike y, yQ is not allocated by the user   */
  N_Vector cp_acorQ;           /* acorQ = yQ_n(m) - yQ_n(0)                   */
  N_Vector cp_tempvQ;          /* temporary storage vector (~ tempv)          */

  realtype cp_acnrmQ;          /* acnrmQ = ||acorQ||_WRMS                     */

  long int cp_lrw1Q;           /* no. of realtype words in 1 N_Vector yQ      */ 
  long int cp_liw1Q;           /* no. of integer words in 1 N_Vector yQ       */ 

  int cp_qmax_allocQ;          /* qmax used when allocating quad. memory      */

  booleantype cp_VabstolQMallocDone;  /* did we allocate memory for abstolQ?  */
  booleantype cp_quadMallocDone;      /* was quadrature memory allocated?     */

  long int cp_nqe;             /* number of calls to qfun                     */
  long int cp_netfQ;           /* number of quadr. error test failures        */

  /*-----------------
    Tstop information
    -----------------*/

  booleantype cp_tstopset;
  booleantype cp_istop;
  realtype cp_tstop;

  /*------
    Limits 
    ------*/

  int cp_qmax;          /* q <= qmax                                          */
  long int cp_mxstep;   /* maximum number of internal steps for one user call */
  int cp_maxcor;        /* maximum number of corrector iterations for the     */
  /* solution of the nonlinear equation                 */
  int cp_mxhnil;        /* maximum number of warning messages issued to the   */
  /* user that t + h == t for the next internal step    */
  int cp_maxnef;        /* maximum number of error test failures              */
  int cp_maxncf;        /* maximum number of nonlinear convergence failures   */

  realtype cp_hmin;     /* |h| >= hmin                                        */
  realtype cp_hmax_inv; /* |h| <= 1/hmax_inv                                  */
  realtype cp_etamax;   /* eta <= etamax                                      */

  /*--------
    Counters 
    --------*/

  long int cp_nst;             /* number of internal steps taken              */
  long int cp_nfe;             /* number of f calls                           */
  long int cp_ncfn;            /* number of corrector convergence failures    */
  long int cp_netf;            /* number of error test failures               */
  long int cp_nni;             /* number of nonlinear iterations performed    */
  long int cp_nsetups;         /* number of calls to lsetup                   */
  int cp_nhnil;                /* number of messages saying that t+h==t       */

  realtype cp_etaqm1;          /* ratio of new to old h for order q-1         */
  realtype cp_etaq;            /* ratio of new to old h for order q           */
  realtype cp_etaqp1;          /* ratio of new to old h for order q+1         */

  /*----------------------------
    Space requirements for CPODES 
    ----------------------------*/

  long int cp_lrw1;            /* no. of realtype words in 1 state N_Vector   */ 
  long int cp_liw1;            /* no. of integer words in 1 state N_Vector    */ 
  long int cp_lrw2;            /* no. of realtype words in 1 cnstr. N_Vector  */ 
  long int cp_liw2;            /* no. of integer words in 1 cnstr. N_Vector   */ 
  long int cp_lrw;             /* no. of realtype words in CPODES work vectors*/
  long int cp_liw;             /* no. of integer words in CPODES work vectors */

  /*---------------------------------------
    Implicit Integration Linear Solver Data 
    ---------------------------------------*/

  /* Linear Solver functions to be called */
  int (*cp_linit)(struct CPodeMemRec *cp_mem);
  int (*cp_lsetup)(struct CPodeMemRec *cp_mem, int convfail,
                   N_Vector yP, N_Vector ypP, N_Vector fctP,
                   booleantype *jcurPtr,
                   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3); 
  int (*cp_lsolve)(struct CPodeMemRec *cp_mem, N_Vector b, N_Vector weight,
                   N_Vector yC, N_Vector ypC, N_Vector fctC);
  void (*cp_lfree)(struct CPodeMemRec *cp_mem);

  /* Linear Solver specific memory */
  void *cp_lmem;           

  /* Does lsetup do anything? */
  booleantype cp_lsetup_exists;

  /*-----------------------------
    Projection Linear Solver Data 
    -----------------------------*/

  /* Linear Solver functions to be called */
  int (*cp_linitP)(struct CPodeMemRec *cp_mem);
  int (*cp_lsetupP)(struct CPodeMemRec *cp_mem, 
                    N_Vector y, N_Vector cy,
                    N_Vector c_tmp1, N_Vector c_tmp2, N_Vector s_tmp1); 
  int (*cp_lsolveP)(struct CPodeMemRec *cp_mem, N_Vector b, N_Vector x,
                    N_Vector y, N_Vector cy,
                    N_Vector c_tmp1, N_Vector s_tmp1);
  void (*cp_lmultP)(struct CPodeMemRec *cp_mem, N_Vector x, N_Vector Gx);
  void (*cp_lfreeP)(struct CPodeMemRec *cp_mem);

  /* Linear Solver specific memory */
  void *cp_lmemP;

  /* Does lsetupP do anything? */
  booleantype cp_lsetupP_exists;

  /*------------
    Saved Values
    ------------*/

  int cp_qu;                    /* last successful q value used               */
  long int cp_nstlset;          /* step number of last lsetup call            */
  realtype cp_h0u;              /* actual initial stepsize                    */
  realtype cp_hu;               /* last successful h value used               */
  realtype cp_saved_tq5;        /* saved value of tq[5]                       */
  booleantype cp_jcur;          /* Is Jacobian info current?                  */
  realtype cp_tolsf;            /* tolerance scale factor                     */
  int cp_qmax_alloc;            /* value of qmax used when allocating memory  */
  int cp_indx_acor;             /* index of the zn vector where acor is saved */

  booleantype cp_MallocDone;  
  booleantype cp_VabstolMallocDone;
  booleantype cp_projMallocDone;
  booleantype cp_rootMallocDone;


  /*-------------------------------------------
    Error handler function and error ouput file 
    -------------------------------------------*/

  CPErrHandlerFn cp_ehfun;     /* Error messages are handled by ehfun         */
  void *cp_eh_data;            /* user pointer passed to ehfun                */
  FILE *cp_errfp;              /* CPODES error messages are sent to errfp     */

  /*-------------------------
    Stability Limit Detection
    -------------------------*/

  booleantype cp_sldeton;      /* Is Stability Limit Detection on?            */
  realtype cp_ssdat[6][4];     /* scaled data array for STALD                 */
  int cp_nscon;                /* counter for STALD method                    */
  long int cp_nor;             /* counter for number of order reductions      */

  /*----------------
    Rootfinding Data
    ----------------*/

  booleantype cp_doRootfinding;/* Is rootfinding enabled?                     */
  CPRootFn cp_gfun;            /* Function g for roots sought                 */
  int cp_nrtfn;                /* number of components of g                   */
  void *cp_g_data;             /* pointer to user data for g                  */
  booleantype *cp_gactive;     /* flags for active/inactive g functions       */
  int *cp_iroots;              /* array for root information                  */
  int *cp_rootdir;             /* array specifying direction of zero-crossing */

  realtype cp_tlo;             /* nearest endpoint of interval in root search */
  realtype cp_thi;             /* farthest endpoint of interval in root search*/
  realtype cp_trout;           /* t value returned by rootfinding routine     */
  realtype *cp_glo;            /* saved array of g values at t = tlo          */
  realtype *cp_ghi;            /* saved array of g values at t = thi          */
  realtype *cp_grout;          /* array of g values at t = trout              */
  realtype cp_toutc;           /* copy of tout (if NORMAL mode)               */
  int cp_taskc;                /* copy of parameter task                      */
  int cp_irfnd;                /* flag showing whether last step had a root   */
  long int cp_nge;             /* counter for g evaluations                   */

  /*------------------------------
    Consistent IC calculation data
    ------------------------------*/

  realtype cp_icprj_convcoef;
  realtype cp_icprj_normtol;
  int cp_icprj_maxrcvr;
  int cp_icprj_maxiter;
  N_Vector cp_yy0;
  N_Vector cp_yp0;

} *CPodeMem;

/*
 * =================================================================
 *     I N T E R F A C E   T O    L I N E A R   S O L V E R S
 *       F O R   I M P L I C I T    I N T E G R A T I O N
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * int (*cp_linit)(CPodeMem cp_mem);
 * -----------------------------------------------------------------
 * The purpose of cp_linit is to complete initializations for a
 * specific linear solver, such as counters and statistics.
 * An LInitFn should return 0 if it has successfully initialized the
 * CPODES linear solver and a negative value otherwise.
 * If an error does occur, an appropriate message should be sent to
 * the error handler function.
 * -----------------------------------------------------------------
 */
  
/*
 * -----------------------------------------------------------------
 * int (*cp_lsetup)(CPodeMem cp_mem, int convfail,
 *                  N_Vector yP, N_Vector ypP, N_Vector fctP, 
 *                  booleantype *jcurPtr,
 *                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
 * -----------------------------------------------------------------
 * The job of cp_lsetup is to prepare the linear solver for subsequent 
 * calls to cp_lsolve. It may recompute Jacobian-related data is it 
 * deems necessary. Its parameters are as follows:
 *
 * cp_mem   - problem memory pointer of type CPodeMem. See the
 *            typedef earlier in this file.
 *
 * convfail - a flag to indicate any problem that occurred during
 *            the solution of the nonlinear equation on the
 *            current time step for which the linear solver is
 *            being used. This flag can be used to help decide
 *            whether the Jacobian data kept by a CPODES linear
 *            solver needs to be updated or not.
 *            Its possible values are documented below.
 *
 *    NOTE: convfail is UNDEFINED for implicit-form ODE
 *
 * yP       - the predicted y vector for the current CPODES
 *            internal step (i.e. zn[0])
 *
 * ypP      - the predicted y' vector for the current CPODES
 *            internal step (i.e. h*zn[1] for CP_IMPL)
 *
 *    NOTE: ypP = NULL for explicit-form ODE
 *
 * fctP     - f(yP) or F(tn, yP, ypP).
 *
 * jcurPtr  - a pointer to a boolean to be filled in by cp_lsetup.
 *            The function should set *jcurPtr=TRUE if its Jacobian
 *            data is current after the call and should set
 *            *jcurPtr=FALSE if its Jacobian data is not current.
 *            Note: If cp_lsetup calls for re-evaluation of
 *            Jacobian data (based on convfail and CPODES state
 *            data), it should return *jcurPtr=TRUE always;
 *            otherwise an infinite loop can result.
 *
 *    NOTE: jcurPtr is IGNORED for implicit-form ODE
 *
 * tmp1, tmp2, tmp3 - temporary N_Vectors provided for use by cp_lsetup.
 *
 * The cp_lsetup routine should return 0 if successful, a positive
 * value for a recoverable error, and a negative value for an
 * unrecoverable error.
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * int (*cp_lsolve)(CPodeMem cp_mem, N_Vector b, N_Vector weight,
 *                  N_Vector yC, N_Vector ypC, N_Vector fctC);
 * -----------------------------------------------------------------
 * cp_lsolve must solve the linear equation P x = b, where
 *   Explicit-form ODE:
 *         P is some approximation to (I - gamma * df/dy)
 *         and the RHS vector b is input. 
 *   Implicit-form ODE:
 *         P is some approximation to (dF/dy' + gamma * dF/dy)
 *         and b is the current residual.
 * Its arguments are as follows:
 *
 * cp_mem   - problem memory pointer of type CPodeMem. See the big
 *            typedef earlier in this file.
 *
 * b        - on input the right-hand side of the linear system
 *            to be solved. On output, the solution.
 *
 * weight   - te current weights in the WRMS norm.
 *
 * yC       - the solver's current approximation to y(tn)
 *
 * ypC      - the solver's current approximation to y'(tn)
 *
 *    NOTE: ypC = NULL for explicit-form ODE
 *
 * fctC     - f(tn, yC) or F(tn, yC, ypc).
 * 
 * The solution is to be returned in the vector b. 
 * The cp_lsolve function should return a positive value for a 
 * recoverable error and a negative value for an unrecoverable error. 
 * Success is indicated by a 0 return value.
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * void (*cp_lfree)(CPodeMem cp_mem);
 * -----------------------------------------------------------------
 * cp_lfree should free up any memory allocated by the linear
 * solver. This routine is called once a problem has been
 * completed and the linear solver is no longer needed.
 * -----------------------------------------------------------------
 */
  
/*
 * -----------------------------------------------------------------
 * Communication between CPODES and a CPODES Linear Solver
 * -----------------------------------------------------------------
 * convfail is passed as an input to cp_lsetup. When dealing with
 * an ODE in explicit form, this can be used to test whether the Jacobian
 * must be reevakluated or whether it is enough to use the up-to-date gamma
 * value. 
 * 
 * convfail can be one of:
 *
 * CP_NO_FAILURES : Either this is the first cp_setup call for this
 *                  step, or the local error test failed on the
 *                  previous attempt at this step (but the Newton
 *                  iteration converged).
 *
 * CP_FAIL_BAD_J  : This value is passed to cp_lsetup if
 *
 *                  (a) The previous Newton corrector iteration
 *                      did not converge and the linear solver's
 *                      setup routine indicated that its Jacobian-
 *                      related data is not current
 *                                   or
 *                  (b) During the previous Newton corrector
 *                      iteration, the linear solver's solve routine
 *                      failed in a recoverable manner and the
 *                      linear solver's setup routine indicated that
 *                      its Jacobian-related data is not current.
 *
 * CP_FAIL_OTHER  : During the current internal step try, the
 *                  previous Newton iteration failed to converge
 *                  even though the linear solver was using current
 *                  Jacobian-related data.
 * -----------------------------------------------------------------
 */


/*
 * =================================================================
 *     I N T E R F A C E   T O    L I N E A R   S O L V E R S
 *                 F O R    P R O J E C T I O N
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * int (*cp_linitP)(CPodeMem cp_mem);
 * -----------------------------------------------------------------
 * The purpose of cp_linitP is to complete initializations for a
 * specific linear solver, such as counters and statistics.
 * An linit function should return 0 if it has successfully 
 * initialized the CPODES linear solver and a negative value otherwise.
 * If an error does occur, an appropriate message should be sent to
 * the error handler function.
 * -----------------------------------------------------------------
 */
  
/*
 * -----------------------------------------------------------------
 * int (*cp_lsetupP)(CPodeMem cp_mem,
 *                   N_Vector y, N_Vector cy,
 *                   N_Vector c_tmp1, N_Vector c_tmp2, N_Vector s_tmp1);
 * -----------------------------------------------------------------
 * The job of cp_lsetupP is to prepare the linear solver for subsequent 
 * calls to cp_lsolveP. It may recompute Jacobian-related data if it 
 * deems necessary. Its parameters are as follows:
 *
 * cp_mem   - problem memory pointer of type CPodeMem. See the 
 *            typedef earlier in this file.
 *
 * y        - the corrected y vector for the current CPODES
 *            internal step (i.e. zn[0])
 *
 * cy      - c(tn, y)
 *
 * c_tmp1, c_tmp2 - temporary N_Vectors (of same length as c) provided 
 *                  for use by cp_lsetupP.
 * s_tmp1         - temporary N_Vectors (of same length as y) provided 
 *                  for use by cp_lsetupP.
 *
 * The cp_lsetupP routine should return 0 if successful, a positive
 * value for a recoverable error, and a negative value for an
 * unrecoverable error.
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * int (*cp_lsolveP)(CPodeMem cp_mem, N_Vector b, N_Vector x,
 *                   N_Vector y, N_Vector cy, 
 *                   N_Vector c_tmp1, N_Vector s_tmp1);
 * -----------------------------------------------------------------
 * cp_lsolveP must solve the linear equation:
 *
 *            
 *   Explicit-form ODE:
 *         P is some approximation to (I - gamma * df/dy)
 *         and the RHS vector b is input. 
 *   Implicit-form ODE:
 *         P is some approximation to (dF/dy' + gamma * dF/dy)
 *         and b is the current residual.
 * Its arguments are as follows:
 *
 * cp_mem   - problem memory pointer of type CPodeMem. See the big
 *            typedef earlier in this file.
 *
 * b        - on input the right-hand side of the linear system
 *            to be solved. On output, the solution.
 *
 * weight   - te current weights in the WRMS norm.
 *
 * yC       - the solver's current approximation to y(tn)
 *
 * ypC      - the solver's current approximation to y'(tn)
 *
 *    NOTE: ypC = NULL for explicit-form ODE
 *
 * fctC     - f(tn, yC) or F(tn, yC, ypc).
 * 
 * The solution is to be returned in the vector b. 
 * The cp_lsolve function should return a positive value for a 
 * recoverable error and a negative value for an unrecoverable error. 
 * Success is indicated by a 0 return value.
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * void (*cp_lmultP)(CPodeMem cp_mem, N_Vector x, N_Vector Gx);
 * -----------------------------------------------------------------
 * cp_lmultP should compute the Jacobian - vector product, for a
 * given vector x and return the result in Gx.
 * -----------------------------------------------------------------
 */
  
/*
 * -----------------------------------------------------------------
 * void (*cp_lfreeP)(CPodeMem cp_mem);
 * -----------------------------------------------------------------
 * cp_lfreeP should free up any memory allocated by the linear
 * solver. This routine is called once a problem has been
 * completed and the linear solver is no longer needed.
 * -----------------------------------------------------------------
 */
  
/*
 * -----------------------------------------------------------------
 * Communication between CPODES and a CPODES Linear Solver
 * -----------------------------------------------------------------
 * convfail is passed as an input to cp_lsetup. When dealing with
 * an ODE in explicit form, this can be used to test whether the Jacobian
 * must be reevakluated or whether it is enough to use the up-to-date gamma
 * value. 
 * 
 * convfail can be one of:
 *
 * CP_NO_FAILURES : Either this is the first cp_setup call for this
 *                  step, or the local error test failed on the
 *                  previous attempt at this step (but the Newton
 *                  iteration converged).
 *
 * CP_FAIL_BAD_J  : This value is passed to cp_lsetup if
 *
 *                  (a) The previous Newton corrector iteration
 *                      did not converge and the linear solver's
 *                      setup routine indicated that its Jacobian-
 *                      related data is not current
 *                                   or
 *                  (b) During the previous Newton corrector
 *                      iteration, the linear solver's solve routine
 *                      failed in a recoverable manner and the
 *                      linear solver's setup routine indicated that
 *                      its Jacobian-related data is not current.
 *
 * CP_FAIL_OTHER  : During the current internal step try, the
 *                  previous Newton iteration failed to converge
 *                  even though the linear solver was using current
 *                  Jacobian-related data.
 * -----------------------------------------------------------------
 */

  
#define CP_NO_FAILURES 0
#define CP_FAIL_BAD_J  1
#define CP_FAIL_OTHER  2
    

#ifdef __cplusplus
}
#endif

#endif
