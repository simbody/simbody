/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2006/11/30 21:11:29 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban  @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementatino for the internal projection step.
 * -----------------------------------------------------------------
 */

/*
 * =================================================================
 * Import Header Files
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include "cpodes_private.h"
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

/* Maximum step size decrease on projection failure */
#define ETAPR  RCONST(0.25)

/*
 * =================================================================
 * Private Function Prototypes
 * =================================================================
 */

static int cpProjLinear(CPodeMem);
static int cpProjNonlinear(CPodeMem);
static int cpProjNonlinearIteration(CPodeMem cp_mem);

/* 
 * =================================================================
 * Readibility Constants
 * =================================================================
 */

#define tn             (cp_mem->cp_tn)
#define zn             (cp_mem->cp_zn) 
#define y              (cp_mem->cp_y)
#define nst            (cp_mem->cp_nst)
#define ewt            (cp_mem->cp_ewt)
#define acor           (cp_mem->cp_acor)
#define acnrm          (cp_mem->cp_acnrm)
#define ftemp          (cp_mem->cp_ftemp)
#define tempv          (cp_mem->cp_tempv)
#define q              (cp_mem->cp_q)
#define h              (cp_mem->cp_h)
#define next_h         (cp_mem->cp_next_h)
#define hmin           (cp_mem->cp_hmin)
#define hscale         (cp_mem->cp_hscale)
#define eta            (cp_mem->cp_eta)
#define etamax         (cp_mem->cp_etamax)
#define nscon          (cp_mem->cp_nscon)

#define proj_type      (cp_mem->cp_proj_type)
#define proj_norm      (cp_mem->cp_proj_norm)
#define cnstr_type     (cp_mem->cp_cnstr_type)
#define pfun           (cp_mem->cp_pfun)
#define p_data         (cp_mem->cp_p_data)
#define cfun           (cp_mem->cp_cfun)
#define c_data         (cp_mem->cp_c_data)

#define yC             (cp_mem->cp_yC)
#define acorP          (cp_mem->cp_acorP)
#define errP           (cp_mem->cp_errP)
#define ctol           (cp_mem->cp_ctol)
#define ctemp          (cp_mem->cp_ctemp)
#define tempvP1        (cp_mem->cp_tempvP1)
#define tempvP2        (cp_mem->cp_tempvP2)

#define lsetupP        (cp_mem->cp_lsetupP)
#define lsolveP        (cp_mem->cp_lsolveP)
#define lmultP         (cp_mem->cp_lmultP)
#define lsetupP_exists (cp_mem->cp_lsetupP_exists)
#define lsetupP_freq   (cp_mem->cp_lsetupP_freq)

#define prjcoef        (cp_mem->cp_prjcoef)
#define crateP         (cp_mem->cp_crateP)
#define nstlsetP       (cp_mem->cp_nstlsetP)
#define maxcorP        (cp_mem->cp_maxcorP)

#define project_err    (cp_mem->cp_project_err)
#define test_cnstr     (cp_mem->cp_test_cnstr)
#define first_proj     (cp_mem->cp_first_proj)
#define applyProj      (cp_mem->cp_applyProj)

#define nproj          (cp_mem->cp_nproj)
#define nprf           (cp_mem->cp_nprf)
#define nce            (cp_mem->cp_nce)
#define nsetupsP       (cp_mem->cp_nsetupsP)
#define maxnpf         (cp_mem->cp_maxnpf)

/*
 * =================================================================
 * Main interface function
 * =================================================================
 */

/* 
 * cpDoProjection
 *
 * For user supplied projection function, use ftemp as temporary storage
 * for the current error estimate (acor) and use tempv to store the
 * accumulated corection due to projection, acorP (tempv is not touched
 * until it is potentially used in cpCompleteStep).
 *
 * For the internal projection function, both ftemp and tempv are needed
 * for temporary storage, so space for acorP was allocated.
 */

int cpDoProjection(CPodeMem cp_mem, realtype saved_t, int *npfPtr)
{
  int flag, retval;
  realtype cnorm;

  switch (proj_type) {

  case CP_PROJ_INTERNAL:

    /* Evaluate constraints at current time and with the corrected y */
    retval = cfun(tn, y, ctemp, c_data);
    nce++;
    if (retval < 0) {flag = CP_CNSTRFUNC_FAIL; break;}
    if (retval > 0) {flag = CNSTRFUNC_RECVR; break;}

    /*
     * If activated, evaluate WL2 norm of constraint violation.
     * If the constraint violation is small enough, return. 
     */
    if (test_cnstr) {
      cnorm = N_VWL2Norm(ctemp, ctol);
      cnorm /= prjcoef;

#ifdef CPODES_DEBUG
      printf("      Constraint violation norm = %lg\n",cnorm);
#endif

      if (cnorm <= ONE) {
        applyProj = FALSE;
        return(CP_SUCCESS);
      }
    }

#ifdef CPODES_DEBUG
    else {
      printf("      No constraint testing\n");
    }
#endif

    /* Perform projection step 
     * On a successful return, the projection correction is available in acorP.
     * Also, if projection of the error estimate was enabled, the new error
     * estimate is available in errP and acnrm contains ||errP||_WRMS.
     */
    nproj++;
    if (cnstr_type == CP_CNSTR_NONLIN) flag = cpProjNonlinear(cp_mem);
    else                               flag = cpProjLinear(cp_mem);

    break;

  case CP_PROJ_USER:

#ifdef CPODES_DEBUG
    printf("      User-defined projection\n");
#endif

    /* Use ftemp to store errP and tempv to store acorP 
     * (recall that in this case we did not allocate memory
     * errP and acorP).
     */
    errP = ftemp;
    acorP = tempv;
    
    /* Copy acor into errP */
    N_VScale(ONE, acor, errP);

    /* Call the user projection function */
    retval = pfun(tn, y, acorP, prjcoef, errP, p_data);
    nproj++;
    if (retval < 0) {flag = CP_PROJFUNC_FAIL; break;}
    if (retval > 0) {flag = PROJFUNC_RECVR; break;}

    /* Recompute acnrm to be used in error test */
    acnrm = N_VWrmsNorm(errP, ewt);

    flag = CP_SUCCESS;

    break;

  }

#ifdef CPODES_DEBUG
  printf("      acnrm = %lg\n",acnrm);
#endif

  /* This is not the first projection anymore */
  first_proj = FALSE;

  /* If the projection was successful, return now. */
  if (flag == CP_SUCCESS) {
    applyProj = TRUE;
    return(CP_SUCCESS);
  }

  /* The projection failed. Increment nprf and restore zn */
  nprf++;
  cpRestore(cp_mem, saved_t);

  /* Return if lsetupP, lsolveP, cfun, or pfun failed unrecoverably */
  if (flag == CP_PLSETUP_FAIL)   return(CP_PLSETUP_FAIL);
  if (flag == CP_PLSOLVE_FAIL)   return(CP_PLSOLVE_FAIL);
  if (flag == CP_CNSTRFUNC_FAIL) return(CP_CNSTRFUNC_FAIL);
  if (flag == CP_PROJFUNC_FAIL)  return(CP_PROJFUNC_FAIL);

  /*  At this point, flag = CONV_FAIL or CNSTRFUNC_RECVR or PRJFUNC_RECVR; increment npf */
  (*npfPtr)++;
  etamax = ONE;
  
  /* If we had maxnpf failures or |h| = hmin, 
     return CP_PROJ_FAILURE or CP_REPTD_CNSTRFUNC_ERR or CP_REPTD_PROJFUNC_ERR. */
  if ((ABS(h) <= hmin*ONEPSM) || (*npfPtr == maxnpf)) {
    if (flag == CONV_FAIL)       return(CP_PROJ_FAILURE);
    if (flag == CNSTRFUNC_RECVR) return(CP_REPTD_CNSTRFUNC_ERR);    
    if (flag == PROJFUNC_RECVR)  return(CP_REPTD_PROJFUNC_ERR);    
  }

  /* Reduce step size; return to reattempt the step */
  eta = MAX(ETAPR, hmin / ABS(h));
  cpRescale(cp_mem);

  return(PREDICT_AGAIN);

}

/*
 * cpProjLinear
 *
 * This is the internal CPODES implementation of the projection function
 * for linear constraints. On entry, ctemp contains c(t,y).
 * 
 * Possible return values:
 *
 *   CP_SUCCESS         ---> continue with error test
 *
 *   CP_CNSTRFUNC_FAIL  -+  
 *   CP_PLSETUP_FAIL     |-> halt the integration 
 *   CP_PLSOLVE_FAIL    -+
 *
 *   CONV_FAIL          -+
 *   CNSTRFUNC_RECVR    -+-> predict again or stop if too many
 */

static int cpProjLinear(CPodeMem cp_mem)
{
  int retval;

#ifdef CPODES_DEBUG
  printf("      Internal linear projection\n");
#endif

  /* Call the lsetupP function to evaluate and factorize the 
   * Jacobian of constraints ONLY the first time we get here */
  if (first_proj) {
    retval = lsetupP(cp_mem, y, ctemp, tempvP1, tempvP2, tempv);
    nsetupsP++;
    if (retval < 0) return(CP_PLSETUP_FAIL);
    if (retval > 0) return(CONV_FAIL);
  }

  /* Call lsolveP (rhs is ctemp; solution in acorP) */
  retval = lsolveP(cp_mem, ctemp, acorP, y, ctemp, tempvP1, tempv);
  if (retval < 0) return(CP_PLSOLVE_FAIL);
  if (retval > 0) return(CONV_FAIL);

  /* Fix the sign of acorP */
  N_VScale(-ONE, acorP, acorP);

  /* If activated, compute error estimate of the projected method. */
  if (project_err) {

    /* Find G*acor and load it in tempvP2 */
    lmultP(cp_mem, acor, tempvP2);

    /* Call lsolveP (rhs is tempvP2; solution in errP) */
    retval = lsolveP(cp_mem, tempvP2, errP, y, ctemp, tempvP1, tempv);
    if (retval < 0) return(CP_PLSOLVE_FAIL);
    if (retval > 0) return(CONV_FAIL);

    /* Update error estimate: errP now contains the new estimated errors. */
    N_VLinearSum(ONE, acor, -ONE, errP, errP);

    /* Recompute acnrm to be used in error test */
    acnrm = N_VWrmsNorm(errP, ewt);

  }

  return(CP_SUCCESS);
}

/*
 * cpProjNonlinear
 *
 * This is the internal CPODES implementation of the projection function
 * for nonlinear constraints. On entry, ctemp contains c(t,y).
 * 
 * Possible return values:
 *
 *   CP_SUCCESS         ---> continue with error test
 *
 *   CP_CNSTRFUNC_FAIL  -+  
 *   CP_PLSETUP_FAIL     |-> halt the integration 
 *   CP_PLSOLVE_FAIL    -+
 *
 *   CONV_FAIL          -+
 *   CNSTRFUNC_RECVR    -+-> predict again or stop if too many
 */

static int cpProjNonlinear(CPodeMem cp_mem)
{
  booleantype callSetup;
  int retval;
  
#ifdef CPODES_DEBUG
  printf("      Internal nonlinear projection\n");
#endif

  /* If the linear solver provides a setup function, call it if:
   *   - it was never called before (first_proj = TRUE), or
   *   - enough steps passed from last evaluation */
  if (lsetupP_exists) {
    callSetup =  (first_proj) || (nst >= nstlsetP + lsetupP_freq);
  } else {
    callSetup = FALSE;
    crateP = ONE;
  }

  /* Save the corrected y, in case we need a second pass through the loop */
  N_VScale(ONE, y, yC);

  /* Begin the main loop. This loop is traversed at most twice. 
   * The second pass only occurs when the first pass had a recoverable
   * failure with old Jacobian data. */
  loop {

    /* If needed, call setup function */
    if (callSetup) {
      retval = lsetupP(cp_mem, y, ctemp, tempvP1, tempvP2, tempv);

#ifdef CPODES_DEBUG
      printf("         Linear solver setup return value = %d\n",retval);
#endif

      nsetupsP++;
      nstlsetP = nst;
      crateP = ONE;
      if (retval < 0) return(CP_PLSETUP_FAIL);
      if (retval > 0) return(CONV_FAIL);
    }

    /* Do the iteration */

#ifdef CPODES_DEBUG
    printf("         NonlinearIteration\n");
#endif

    retval = cpProjNonlinearIteration(cp_mem);

#ifdef CPODES_DEBUG
    printf("         NonlinearIteration return value = %d\n",retval);        
#endif

    /* On a recoverable failure, if setup was not called,
     * reload the corrected y, recompute constraints, and reattempt loop */
    if ( (retval > 0) && lsetupP_exists && !callSetup ) {
      N_VScale(ONE, yC, y);
      retval = cfun(tn, y, ctemp, c_data);
      nce++;
      if (retval < 0) return(CP_CNSTRFUNC_FAIL);
      if (retval > 0) return(CNSTRFUNC_RECVR);
      callSetup = TRUE;
      continue;
    }

    /* Break from loop */
    break;

  }

  return(retval);
}

/*
 * cpProjNonlinearIteration
 *
 * This routine performs the actual nonlinear iteration for the
 * projection step.
 *
 * Possible return values:
 *
 *   CP_SUCCESS  - the iteration was successful.
 *
 *   CP_PLSOLVE_FAIL - the linear solver solution failed unrecoverably.
 *   CP_CNSTRFUNC_FAIL - the constraint function failed unrecoverably.
 *
 *   CNSTRFUNC_RECVR - the constraint function failed recoverably.
 *   CONV_FAIL - a recoverable error occurred.
 */

static int cpProjNonlinearIteration(CPodeMem cp_mem)
{
  int m, retval;
  realtype del, delp, dcon;
  N_Vector corP;

  /* Initialize counter */
  m = 0;

  /* Initialize delp to avoid compiler warning message */
  del = delp = ZERO;

  /* Rename ftemp as corP */
  corP = ftemp;

  /* Set acorP to zero */
  N_VConst(ZERO, acorP);
  
  /* Set errP to acor */
  if (project_err) N_VScale(ONE, acor, errP);
 
  /* Looping point for iterations */
  loop {

#ifdef CPODES_DEBUG
    printf("            Iteration # %d\n",m);
#endif

    /* Call lsolveP (rhs is ctemp; solution in corP) */
    retval = lsolveP(cp_mem, ctemp, corP, y, ctemp, tempvP1, tempv);

#ifdef CPODES_DEBUG
    printf("            Linear solver solve return value = %d\n",retval);
#endif

    if (retval < 0) return(CP_PLSOLVE_FAIL);
    if (retval > 0) return(CONV_FAIL);

    /* Add correction to acorP and y and get its WRMS norm */
    N_VLinearSum(ONE, acorP, -ONE, corP, acorP);
    N_VLinearSum(ONE, y, -ONE, corP, y);
    del = N_VWrmsNorm(corP, ewt);

#ifdef CPODES_DEBUG
    printf("            Norm of correction:  del = %lg\n", del);
#endif

    /* If activated, find correction to the error estimate */
    if (project_err) {
      /* Compute rhs as G*errP and load it in tempvP2 */
      lmultP(cp_mem, errP, tempvP2);

      /* Call lsolveP (rhs is tempvP2; solution in corP) */
      retval = lsolveP(cp_mem, tempvP2, corP, y, ctemp, tempvP1, tempv);
      if (retval < 0) return(CP_PLSOLVE_FAIL);
      if (retval > 0) return(CONV_FAIL);

      /* Add correction to errP */
      N_VLinearSum(ONE, errP, -ONE, corP, errP);
    }

    /* Test for convergence.  If m > 0, an estimate of the convergence
       rate constant is stored in crateP, and used in the test.        */
    if (m > 0) crateP = MAX(PRJ_CRDOWN * crateP, del/delp);
    dcon = del * MIN(ONE, crateP) / prjcoef;

#ifdef CPODES_DEBUG
    printf("            Convergence test  dcon = %lg\n", dcon);
#endif

    if (dcon <= ONE) {

      if (project_err) acnrm = N_VWrmsNorm(errP, ewt);
      
      return(CP_SUCCESS);
    }
    m++;

    /* Stop at maxcorP iterations or if iter. seems to be diverging. */
    if ((m == maxcorP) || ((m >= 2) && (del > PRJ_RDIV*delp)))  return(CONV_FAIL);
    
    /* Save norm of correction, evaluate cfun, and loop again */
    delp = del;
    retval = cfun(tn, y, ctemp, c_data);
    nce++;
    if (retval < 0) return(CP_CNSTRFUNC_FAIL);
    if (retval > 0) return(CNSTRFUNC_RECVR);

  }

}
