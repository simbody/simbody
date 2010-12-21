/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2006/12/01 22:48:57 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban  @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Implementation of cmputation of consistent initial conditions
 * -----------------------------------------------------------------
 */

/*
 * TODO:
 *
 * cpicImplComputeYp
 */

/* 
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include "cpodes_private.h"
#include <sundials/sundials_math.h>

/* 
 * =================================================================
 * ALGORITHMIC CONSTANTS
 * =================================================================
 */

/* Maximum number of attempts to correct a recoverable
 * constraint function error */
#define MAX_RECVR     5
#define MAX_ITERS    10

/* 
 * =================================================================
 * PROTOTYPES FOR PRIVATE FUNCTIONS
 * =================================================================
 */

/* Initial consistency checks */
static int cpicInitialSetup(CPodeMem cp_mem);

/* Projection onto invariant manifold */
static int cpicDoProjection(CPodeMem cp_mem);

/* Fucntions called by cpicDoProjection */
static int cpicProjLinear(CPodeMem cp_mem);
static int cpicProjNonlinear(CPodeMem cp_mem);

/* Calculation of consistent y' */
static int cpicImplComputeYp(CPodeMem cp_mem);

static void cpicFailFlag(CPodeMem, int flag);

/* 
 * =================================================================
 * READIBILITY REPLACEMENTS
 * =================================================================
 */

#define ode_type       (cp_mem->cp_ode_type)
#define tn             (cp_mem->cp_tn)
#define zn             (cp_mem->cp_zn) 
#define tol_type       (cp_mem->cp_tol_type)
#define efun           (cp_mem->cp_efun)
#define e_data         (cp_mem->cp_e_data)
#define ewt            (cp_mem->cp_ewt)
#define gamma          (cp_mem->cp_gamma)
#define tempv          (cp_mem->cp_tempv)
#define uround         (cp_mem->cp_uround)

#define linit          (cp_mem->cp_linit)
#define lsetup         (cp_mem->cp_lsetup)
#define lsolve         (cp_mem->cp_lsolve)
#define lsetup_exists  (cp_mem->cp_lsetup_exists)

#define proj_enabled   (cp_mem->cp_proj_enabled)
#define proj_type      (cp_mem->cp_proj_type)
#define proj_norm      (cp_mem->cp_proj_norm)
#define cnstr_type     (cp_mem->cp_cnstr_type)
#define pfun           (cp_mem->cp_pfun)
#define p_data         (cp_mem->cp_p_data)
#define cfun           (cp_mem->cp_cfun)
#define c_data         (cp_mem->cp_c_data)

#define prjcoef        (cp_mem->cp_prjcoef)
#define yC             (cp_mem->cp_yC)
#define acorP          (cp_mem->cp_acorP)
#define ctemp          (cp_mem->cp_ctemp)
#define ctol           (cp_mem->cp_ctol)
#define tempvP1        (cp_mem->cp_tempvP1)
#define tempvP2        (cp_mem->cp_tempvP2)

#define linitP         (cp_mem->cp_linitP)
#define lsetupP        (cp_mem->cp_lsetupP)
#define lsolveP        (cp_mem->cp_lsolveP)
#define lsetupP_exists (cp_mem->cp_lsetupP_exists)

#define icprj_convcoef (cp_mem->cp_icprj_convcoef)
#define icprj_normtol  (cp_mem->cp_icprj_normtol)
#define icprj_maxrcvr  (cp_mem->cp_icprj_maxrcvr)
#define icprj_maxiter  (cp_mem->cp_icprj_maxiter)
#define yy0            (cp_mem->cp_yy0)
#define yp0            (cp_mem->cp_yp0)

/* 
 * =================================================================
 * EXPORTED FUNCTIONS
 * =================================================================
 */

/*
 * CPodeCalcIC 
 *
 * This calculates corrected initial conditions that are consistent 
 * with the invariant constraints and (for implicit-form ODEs) with 
 * the ODE system itself. It first projects the initial guess for 
 * the state vector (given by the user through CPodeInit or CPodeReInit) 
 * and then, if necessary, computes a state derivative vector as 
 * solution of F(t0, y0, y0') = 0.
 *
 */
int CPodeCalcIC(void *cpode_mem)
{
  CPodeMem cp_mem;
  int flag;

  /* Check if cpode_mem exists */
  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeCalcIC", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  /* Check if cpode_mem was allocated */
  if (cp_mem->cp_MallocDone == FALSE) {
    cpProcessError(cp_mem, CP_NO_MALLOC, "CPODES", "CPodeCalcIC", MSGCP_NO_MALLOC);
    return(CP_NO_MALLOC);
  }

  /* For explicit ODE, if projection is not enabled, there is nothing to do */
  if ( (ode_type == CP_EXPL) && (!proj_enabled) ) {
    cpProcessError(cp_mem, CP_WARNING, "CPODES", "CPodeCalcIC", MSGCP_IC_NO_WORK);
    return(CP_SUCCESS);
  }

  /* Perform initial setup */
  flag = cpicInitialSetup(cp_mem);
  if (flag != CP_SUCCESS) return(flag);

  /* Initialize various convergence test constants 
   *   icprj_convcoef - constant used to test weighted norms.
   *                    Here, we use a value 10 times smaller than the
   *                    one used in the projectin during integration.
   *   icprj_normtol  - stopping tolerance on max-norm of constraints
   *   icprj_maxrcvr  - maximum number of reductin in full Newton step
   *                    in attempting to correct recoverable constraint 
   *                    function errors.
   *   icprj_maxiter  - maximum number of nonlinear iterations
   */
  icprj_convcoef = PT1 * prjcoef;
  icprj_normtol = RPowerR(uround, HALF);
  icprj_maxrcvr = MAX_RECVR;
  icprj_maxiter = MAX_ITERS;

  /* Allocate space and initialize temporary vectors */
  yy0 = N_VClone(tempv);
  N_VScale(ONE, zn[0], yy0);
  if (ode_type == CP_IMPL) {
    yp0 = N_VClone(tempv);
    N_VScale(ONE, zn[1], yy0);
  }

  /* Compute y consistent with constraints */
  if (proj_enabled)
    flag = cpicDoProjection(cp_mem);
  
  /* Compute yp consistent with DE */
  if ( (flag == CP_SUCCESS) && (ode_type == CP_IMPL) )
    flag = cpicImplComputeYp(cp_mem);
  
  /* If successful, load yy0 and yp0 into Nordsieck array */
  if (flag == CP_SUCCESS) {
    N_VScale(ONE, yy0, zn[0]);
    if (ode_type == CP_IMPL)  N_VScale(ONE, yp0, zn[1]);
  }

  /* Free temporary space */
  N_VDestroy(yy0);
  if (ode_type == CP_IMPL)  N_VDestroy(yp0);

  /* On any type of failure, print message and return proper flag */
  if (flag != CP_SUCCESS) {
    cpicFailFlag(cp_mem, flag);
    return(flag);
  }

  /* Otherwise return success flag */
  return(CP_SUCCESS);

}

/* 
 * =================================================================
 *  PRIVATE FUNCTIONS
 * =================================================================
 */

/*  
 * cpicInitialSetup
 *
 * This routine performs input consistency checks for IC calculation.
 * If needed, it also checks the linear solver module and calls the
 * linear solver initialization routine.
 *
 * This function performs many of the same tests that cpInitialSetup
 * does at the first integration step.
 */
static int cpicInitialSetup(CPodeMem cp_mem)
{
  int flag;

  /* If projection is enabled, check if appropriate projection functions are available */
  if (proj_enabled) {

    switch (proj_type) {

    case CP_PROJ_USER:

      /* Check if user provided the projection function */
      if ( pfun == NULL ) {
        cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeCalcIC", MSGCP_NO_PFUN);
        return(CP_ILL_INPUT);
      }

      break;

    case CP_PROJ_INTERNAL:

      /* Check if user provided the constraint function */
      if (cfun == NULL) {
        cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeCalcIC", MSGCP_NO_CFUN);
        return(CP_ILL_INPUT);
      }
      /* Check that lsolveP exists */ 
      if ( lsolveP == NULL) {
        cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeCalcIC", MSGCP_PLSOLVE_NULL);
        return(CP_ILL_INPUT);
      }
      /* Call linitP if it exists */
      if ( linitP != NULL ) {
        flag = linitP(cp_mem);
        if (flag != 0) {
          cpProcessError(cp_mem, CP_PLINIT_FAIL, "CPODES", "CPodeCalcIC", MSGCP_PLINIT_FAIL);
          return(CP_PLINIT_FAIL);
        }
      }

      break;

    } 

  }

  /* For implicit ODE, we will always need a linear solver */
  if (ode_type == CP_IMPL) {
    /* Check that lsolve exists */ 
    if (lsolve == NULL) {
      cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeCalcIC", MSGCP_LSOLVE_NULL);
      return(CP_ILL_INPUT);
    }
    /* Call linit if it exists */
    if ( linit != NULL ) {
      flag = linit(cp_mem);
      if (flag != 0) {
        cpProcessError(cp_mem, CP_LINIT_FAIL, "CPODES", "CPodeCalcIC", MSGCP_LINIT_FAIL);
        return(CP_LINIT_FAIL);
      }
    }
  }

  /* Did the user provide efun? */
  if (tol_type != CP_WF) {
    efun = cpEwtSet;
    e_data = (void *)cp_mem;
  } else {
    if (efun == NULL) {
      cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeCalcIC", MSGCP_NO_EFUN);
      return(CP_ILL_INPUT);
    }
  }

  return(CP_SUCCESS);
}

/* 
 * cpicDoProjection
 *
 * This is the top-level function handling the projection onto the
 * invariant manifold. It either calls the user-supplied projection
 * function or, depending on the type of constraints, one of the
 * internal projection functions (cpicProjLinear or cpicProjNonlinear).
 *
 * For user supplied projection function, use tempv to store the
 * corection due to projection, acorP (tempv is not touched until it 
 * is potentially used in cpCompleteStep).
 *
 * For the internal projection function, space for acorP was allocated 
 * in CPodeProjInit.
 */
static int cpicDoProjection(CPodeMem cp_mem)
{
  int flag = CP_SUCCESS, retval;

  switch (proj_type) {

  case CP_PROJ_INTERNAL:

    /* Evaluate constraints at initial time and with the provided yy0 */
    retval = cfun(tn, yy0, ctemp, c_data);
    if (retval < 0) return(CP_CNSTRFUNC_FAIL);
    if (retval > 0) return(CP_FIRST_CNSTRFUNC_ERR);

    /* Perform projection step 
     * On a successful return, yy0 was updated. */
    if (cnstr_type == CP_CNSTR_NONLIN) flag = cpicProjNonlinear(cp_mem);
    else                               flag = cpicProjLinear(cp_mem);

    break;

  case CP_PROJ_USER:

    acorP = tempv;
    
    /* Call the user projection function (with err=NULL) */
    retval = pfun(tn, yy0, acorP, icprj_convcoef, NULL, p_data);
    if (retval != 0) return(CP_PROJFUNC_FAIL);

    /* Update yy0 */
    N_VLinearSum(ONE, yy0, ONE, acorP, yy0);

    break;

  }

  return(flag);

}

/*
 * cpicProjLinear
 *
 * This function performs the projection onto a linear invariant manifold
 */
static int cpicProjLinear(CPodeMem cp_mem)
{
  int retval;

  /* Call the lsetupP function to evaluate and factorize the 
   * Jacobian of constraints */
  retval = lsetupP(cp_mem, yy0, ctemp, tempvP1, tempvP2, tempv);
  if (retval != 0) return(CP_PLSETUP_FAIL);

  /* Call lsolveP (rhs is ctemp; solution in acorP) */
  retval = lsolveP(cp_mem, ctemp, acorP, yy0, ctemp, tempvP1, tempv);
  if (retval != 0) return(CP_PLSOLVE_FAIL);

  /* Update yy0 */
  N_VLinearSum(ONE, yy0, -ONE, acorP, yy0);

  return(CP_SUCCESS);
}

/*
 * cpicProjNonlinear
 *
 * This function performs the projection onto a nonlinear invariant manifold.
 *
 * The convergence tests are:
 *  ||c||_WL2 < icprj_convcoef AND ||c||_max < icprj_normtol
 *       OR
 *  ||p||_WRMS < icprj_convcoef 
 *
 * Default values:
 *   icprj_convcoef = 0.1 * prjcoef = 0.1 * 0.1
 *   icprj_normtol  = sqrt(uround)
 *
 * In other words, we stop when the constraints are within the user-specified
 * tolerances (but, to cover the case where ctol was too large, we also impose
 * that the maximum constraint violation is small enough). Otherwise, we declare 
 * convergence when the current correction becomes small enough (within the
 * integration tolerances).
 *
 * Convergence failure is declared on any unrecoverable function failure, or
 * if it is not possible to correct a recoverable function error, or if the
 * maximum number of iterations was reached with up-to-date Jacobian.
 *
 */
static int cpicProjNonlinear(CPodeMem cp_mem)
{
  N_Vector yy0_new;
  realtype pnorm, pnorm_p;
  realtype cnorm, cmax;
  realtype ccon, pcon, crate, ratio;
  booleantype callSetup, jacCurrent, cOK;
  int retval, m, ircvr, flag;

  /* Evaluate ewt at yy0 */
  flag = efun(yy0, ewt, e_data);
  if (flag != 0) return(CP_ILL_INPUT);
  
  /* Rename yC to yy0_new */
  yy0_new = yC;

  /* Initializations */
  m = 0;
  callSetup = TRUE;
  crate = ONE;
  pnorm = ZERO;
  pnorm_p = ZERO;

  /* Looping point for iterations */
  loop {

    /* 
     * 1. Evaluate Jacobian (if requested)
     */

    jacCurrent = FALSE;    
    if (lsetupP_exists && callSetup) {
      retval = lsetupP(cp_mem, yy0, ctemp, tempvP1, tempvP2, tempv);
      if (retval != 0) return(CP_PLSETUP_FAIL);
      jacCurrent = TRUE;
      callSetup = FALSE;
      crate = ONE;
      m = 0;
    }

#ifdef CPODES_DEBUG
    printf("Iteration %d\n",m+1);
    printf("  Jacobian current: %d\n",jacCurrent);
#endif

    /* 
     * 2. Solve for Newton step 
     */

    /* compute Newton step and load it into acorP */
    retval = lsolveP(cp_mem, ctemp, acorP, yy0, ctemp, tempvP1, tempv);
    /* If lsolveP failed unrecoverably, return */ 
    if (retval < 0) return(CP_PLSOLVE_FAIL);
    /* If lsolveP had a recoverable error, with up-to-date Jacobian, return */
    if (retval > 0) {
      if (!lsetupP_exists || jacCurrent) return(CP_NO_RECOVERY);
      callSetup = TRUE;
      continue;
    }

    /*
     * 3. Apply Newton step
     */

    /* Compute step length = ||acorP|| */
    pnorm = N_VWrmsNorm(acorP, ewt);
    ratio = ONE;
    
#ifdef CPODES_DEBUG
    printf("  Norm of full Newton step: %lg\n", pnorm);
#endif
    
    /* Attempt (at most icprj_maxrcvr times) to evaluate 
       constraint function at the new iterate */
    cOK = FALSE;
    for (ircvr=0; ircvr<icprj_maxrcvr; ircvr++) {
      
      /* compute new iterate */
      N_VLinearSum(ONE, yy0, -ONE, acorP, yy0_new);
      
      /* evaluate constraints at yy0_new */
      retval = cfun(tn, yy0_new, ctemp, c_data);
      
      /* if successful, accept current acorP */
      if (retval == 0) {cOK = TRUE; break;}
      
      /* if function failed unrecoverably, give up */
      if (retval < 0) return(CP_CNSTRFUNC_FAIL);
      
      /* function had a recoverable error; cut step in half */
      pnorm *= HALF;
      ratio *= HALF;
      N_VScale(HALF, acorP, acorP);
      
    }

    /* If cfunc failed recoverably MAX_RECVR times, give up */
    if (!cOK) return(CP_REPTD_CNSTRFUNC_ERR);

#ifdef CPODES_DEBUG
    printf("  Full Newton step cut by ratio: %lg\n", ratio);
#endif

    /*
     * 4. Evaluate various stoping test quantities
     */
    
    /* Weighted L-2 norm of constraints */
    cnorm = N_VWL2Norm(ctemp, ctol);
    ccon = cnorm / icprj_convcoef;

    /* Max-norm of constraints */
    cmax = N_VMaxNorm(ctemp);

    /* Estimated convergence rate */
    if (m > 0) crate = MAX(PRJ_CRDOWN * crate, pnorm/pnorm_p);

    /* Convergence test based on pnorm and crate */
    pcon = pnorm * MIN(ONE, crate) / icprj_convcoef;

#ifdef CPODES_DEBUG
    printf("  Stop test quantities:\n");
    printf("    crate: %lg\n", crate);
    printf("    cnorm: %lg\n", cnorm);
    printf("    ccon:  %lg  <? %lg\n", ccon, ONE);
    printf("    cmax:  %lg  <? %lg\n", cmax, icprj_normtol);
    printf("    pcon:  %lg  <? %lg\n", pcon, ONE);
#endif

    /*
     * 4. Test for convergence
     */

    /* If converged, load new solution and return */
    if ( (pcon <= ONE) || (ccon <= ONE && cmax <= icprj_normtol) ) {
      N_VScale(ONE, yy0_new, yy0);
      return(CP_SUCCESS);
    }
    
    /* Increment iteration counter */
    m++;

#ifdef CPODES_DEBUG
    printf("    m = %d ==? %d = maxiter\n",m,icprj_maxiter);
    if (m >=2)
      printf("    pnorm = %lg  >?  %lg = PRJ_RDIV*pnorm_p\n", pnorm, PRJ_RDIV*pnorm_p);
#endif

    /* Check if we reached max. iters. or if iters. seem to be diverging */
    if ((m == icprj_maxiter) || ((m >= 2) && (pnorm > PRJ_RDIV*pnorm_p))) {

      /* If the Jacobian is up-to-date, give up */
      if (!lsetupP_exists || jacCurrent) return(CP_PROJ_FAILURE);

      /* Otherwise, attempt to recover by re-evaluating the Jacobian */
      callSetup = TRUE;
      retval = cfun(tn, yy0, ctemp, c_data);
      continue;

    }

    /* Save norm of correction, update solution and ewt, and loop again */
    pnorm_p = pnorm;
    N_VScale(ONE, yy0_new, yy0);
    flag = efun(yy0, ewt, e_data);
    if (flag != 0) return(CP_ILL_INPUT);

  }

  return(CP_SUCCESS);
}


/*
 * cpicImplComputeYp
 *
 * For implicit-form ODEs, this function computes y' values consistent 
 * with the ODE, given y values consistent with the constraints.
 * In other words, it solves F(t0,y0,y0') = 0 for y0' using a Newton
 * iteration.
 */
static int cpicImplComputeYp(CPodeMem cp_mem)
{
  /*
  int retval;
  */

  /* Evaluate residual at initial time, using the (projected) yy0 
     and the initial guess yp0 */

  /*
  retval = fi(tn, yy0, yp0, ftemp, f_data);
  if (retval < 0) return(CP_ODEFUNC_FAIL);
  if (retval > 0) return(CP_FIRST_ODEFUNC_ERR);
  */

  /* Set gamma */

  /*
  gamma = ZERO;
  */
  

  return(CP_SUCCESS);
}

/*
 * cpicFailFlag
 *
 * This function is called if the calculation of consistent IC failed.
 * Depending on the flag value, it calls the error handler and returns
 * the proper flag.
 */
static void cpicFailFlag(CPodeMem cp_mem, int flag)
{
  switch(flag) {
  case CP_FIRST_CNSTRFUNC_ERR:
    cpProcessError(cp_mem, CP_FIRST_CNSTRFUNC_ERR, "CPODES", "CPodeCalcIC", MSGCP_IC_CNSTRFUNC_FIRST);
    break;
  case CP_CNSTRFUNC_FAIL:
    cpProcessError(cp_mem, CP_CNSTRFUNC_FAIL, "CPODES", "CPodeCalcIC", MSGCP_IC_CNSTRFUNC_FAILED);
    break;
  case CP_REPTD_CNSTRFUNC_ERR:
    cpProcessError(cp_mem, CP_REPTD_CNSTRFUNC_ERR, "CPODES", "CPodeCalcIC", MSGCP_IC_CNSTRFUNC_REPTD);
    break;
  case CP_PROJFUNC_FAIL:
    cpProcessError(cp_mem, CP_PROJFUNC_FAIL, "CPODES", "CPodeCalcIC", MSGCP_IC_PROJFUNC_FAILED);
    break;
  case CP_PLSETUP_FAIL:
    cpProcessError(cp_mem, CP_PLSETUP_FAIL, "CPODES", "CPodeCalcIC", MSGCP_IC_PLSETUP_FAILED);
    break;
  case CP_PLSOLVE_FAIL:
    cpProcessError(cp_mem, CP_PLSOLVE_FAIL, "CPODES", "CPodeCalcIC", MSGCP_IC_PLSOLVE_FAILED);
    break;
  case CP_NO_RECOVERY:
    cpProcessError(cp_mem, CP_NO_RECOVERY, "CPODES", "CPodeCalcIC", MSGCP_IC_NO_RECOVERY);
    break;
  case CP_PROJ_FAILURE:
    cpProcessError(cp_mem, CP_PROJ_FAILURE, "CPODES", "CPodeCalcIC", MSGCP_IC_PROJ_FAILS);
    break;
  case CP_ILL_INPUT:
    if (tol_type == CP_WF) cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeCalcIC", MSGCP_IC_EWT_FAIL);
    else                   cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeCalcIC", MSGCP_IC_EWT_BAD);
    break;
  case CP_FIRST_ODEFUNC_ERR:
    cpProcessError(cp_mem, CP_FIRST_ODEFUNC_ERR, "CPODES", "CPodeCalcIC", MSGCP_IC_ODEFUNC_FIRST);
    break;
  case CP_ODEFUNC_FAIL:
    cpProcessError(cp_mem, CP_ODEFUNC_FAIL, "CPODES", "CPodeCalcIC", MSGCP_IC_ODEFUNC_FAILED);
    break;

  }
}

