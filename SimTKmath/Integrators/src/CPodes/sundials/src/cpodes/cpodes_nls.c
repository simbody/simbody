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
 * This is the implementatino for the nonlinear solvers for CPODES.
 * -----------------------------------------------------------------
 */

/*
 * NOTE: cpNlsFunctionalImpl is never used
 * (the combination CP_IMPL/CP_FUNCTIONAL is disallowed in CPodeCreate)
 */


/*
 * =================================================================
 * Import Header Files
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "cpodes_private.h"
#include <sundials/sundials_math.h>

/*
 * =================================================================
 * Private Function Prototypes
 * =================================================================
 */

static int cpNlsFunctionalExpl(CPodeMem cp_mem);
static int cpNlsNewtonExpl(CPodeMem cp_mem, int nflag);
static int cpNewtonIterationExpl(CPodeMem cp_mem);

static int cpNlsFunctionalImpl(CPodeMem cp_mem);
static int cpNlsNewtonImpl(CPodeMem cp_mem, int nflag);
static int cpNewtonIterationImpl(CPodeMem cp_mem);

/*
 * =================================================================
 * Readibility Constants
 * =================================================================
 */

#define ode_type       (cp_mem->cp_ode_type)
#define nls_type       (cp_mem->cp_nls_type)

#define fi             (cp_mem->cp_fi)
#define fe             (cp_mem->cp_fe)
#define f_data         (cp_mem->cp_f_data)

#define zn             (cp_mem->cp_zn)
#define ewt            (cp_mem->cp_ewt)
#define y              (cp_mem->cp_y)
#define yp             (cp_mem->cp_yp)
#define maxcor         (cp_mem->cp_maxcor)
#define acor           (cp_mem->cp_acor)
#define tempv          (cp_mem->cp_tempv)
#define ftemp          (cp_mem->cp_ftemp)
#define q              (cp_mem->cp_q)
#define h              (cp_mem->cp_h)
#define next_h         (cp_mem->cp_next_h)
#define hmin           (cp_mem->cp_hmin)
#define hscale         (cp_mem->cp_hscale)
#define eta            (cp_mem->cp_eta)
#define etamax         (cp_mem->cp_etamax)
#define tn             (cp_mem->cp_tn)
#define tq             (cp_mem->cp_tq)
#define rl1            (cp_mem->cp_rl1)
#define gamma          (cp_mem->cp_gamma)
#define gammap         (cp_mem->cp_gammap)
#define gamrat         (cp_mem->cp_gamrat)
#define crate          (cp_mem->cp_crate)
#define acnrm          (cp_mem->cp_acnrm)
#define maxncf         (cp_mem->cp_maxncf)
#define mnewt          (cp_mem->cp_mnewt)
#define nst            (cp_mem->cp_nst)
#define nfe            (cp_mem->cp_nfe)
#define ncfn           (cp_mem->cp_ncfn)
#define nni            (cp_mem->cp_nni)
#define nsetups        (cp_mem->cp_nsetups)
#define nscon          (cp_mem->cp_nscon)

#define lsetup         (cp_mem->cp_lsetup)
#define lsolve         (cp_mem->cp_lsolve)

#define jcur           (cp_mem->cp_jcur)
#define nstlset        (cp_mem->cp_nstlset)

#define lsetup_exists  (cp_mem->cp_lsetup_exists)

/*
 * =================================================================
 * Main interface function
 * =================================================================
 */

/*
 * cpNls
 *
 * This routine attempts to solve the nonlinear system associated
 * with a single implicit step of the linear multistep method.
 * Depending on nls_type, it calls cpNlsFunctional or cpNlsNewton
 * to do the work.
 *
 * The input arguments to cpNls are as follows:
 *  nflag is one of FIRST_CALL, PREV_CONV_FAIL, or PREV_ERR_FAIL.
 *  saved_t is the old value of tn (in case we need to restore)
 *  ncfPtr is a pointer accumulating nls failures at current step.
 *
 * If the nonlinear solver is successful, cpNls retuns CP_SUCCESS.
 *
 * If the nonlinear solver fails, the following actions are taken:
 *
 *  (1) ncfn and ncf=*ncfPtr are incremented and Nordsieck array zn
 *      is restored.
 *
 *  (2) If the solution of the nonlinear system failed due to an
 *      unrecoverable failure by lsetup, lsolve, or fun, we return
 *      an appropriate failure flag (CP_LSETUP_FAIL, CP_LSOLVE_FAIL,
 *      or CP_ODEFUNC_FAIL) which will tell cpStep to halt.
 *
 *  (3) Otherwise, a recoverable failure occurred when solving the
 *      nonlinear system (flag == CONV_FAIL or ODEFUNC_RECVR).
 *      In this case, if ncf is now equal to maxncf or |h| = hmin,
 *      we return the value CP_CONV_FAILURE (if flag=CONV_FAIL) or
 *      CP_REPTD_ODEFUNC_ERR (if flag=ODEFUNC_RECVR).
 *      If not, we return the value PREDICT_AGAIN, telling cpStep to
 *      reattempt the step (with nflag = PREV_CONV_FAIL).
 */

int cpNls(CPodeMem cp_mem, int nflag, realtype saved_t, int *ncfPtr)
{
  int flag = CP_SUCCESS;

  switch(nls_type) {

  case CP_FUNCTIONAL:
    switch(ode_type) {
    case CP_EXPL:
      flag = cpNlsFunctionalExpl(cp_mem);
      break;
    case CP_IMPL:
      flag = cpNlsFunctionalImpl(cp_mem);
      break;
    }
    break;

  case CP_NEWTON:
    switch(ode_type) {
    case CP_EXPL:
      flag = cpNlsNewtonExpl(cp_mem, nflag);
      break;
    case CP_IMPL:
      flag = cpNlsNewtonImpl(cp_mem, nflag);
      break;
    }
    break;

  }

#ifdef CPODES_DEBUG
  printf("      acnrm = %lg\n",acnrm);
#endif

  /* Return now if the nonlinear solver was successful */
  if (flag == CP_SUCCESS) return(CP_SUCCESS);

  /* The nonlinear solver failed. Increment ncfn and restore zn */
  ncfn++;
  cpRestore(cp_mem, saved_t);

  /* Return if lsetup, lsolve, or rhs failed unrecoverably */
  if (flag == CP_LSETUP_FAIL)  return(CP_LSETUP_FAIL);
  if (flag == CP_LSOLVE_FAIL)  return(CP_LSOLVE_FAIL);
  if (flag == CP_ODEFUNC_FAIL) return(CP_ODEFUNC_FAIL);

  /* At this point, flag = CONV_FAIL or ODEFUNC_RECVR; increment ncf */
  (*ncfPtr)++;
  etamax = ONE;

  /* If we had maxncf failures or |h| = hmin,
     return CP_CONV_FAILURE or CP_REPTD_ODEFUNC_ERR. */
  if ((ABS(h) <= hmin*ONEPSM) || (*ncfPtr == maxncf)) {
    if (flag == CONV_FAIL)     return(CP_CONV_FAILURE);
    if (flag == ODEFUNC_RECVR) return(CP_REPTD_ODEFUNC_ERR);
  }

  /* Reduce step size; return to reattempt the step */
  eta = MAX(ETACF, hmin / ABS(h));
  cpRescale(cp_mem);

  return(PREDICT_AGAIN);
}

/*
 * =================================================================
 * Functions for explicit-form ODE
 * =================================================================
 */

/*
 * cpNlsFunctionalExpl
 *
 * This routine attempts to solve the nonlinear system for an ODE in
 * explicit form using functional iteration (no matrices involved).
 *
 * Possible return values are:
 *
 *   CV_SUCCESS      --->  continue with error test
 *
 *   CV_ODEFUNC_FAIL --->  halt the integration
 *
 *   CONV_FAIL       -+
 *   ODEFUNC_RECVR   -+->  predict again or stop if too many
 *
 */

static int cpNlsFunctionalExpl(CPodeMem cp_mem)
{
  int retval, m;
  realtype del, delp, dcon;

#ifdef CPODES_DEBUG
  printf("      Functional, explicit ODE\n");
#endif

  /* Initialize counter */
  crate = ONE;
  m = 0;

  /* Evaluate f at predicted y */
  retval = fe(tn, zn[0], tempv, f_data);
  nfe++;
  if (retval < 0) return(CP_ODEFUNC_FAIL);
  if (retval > 0) return(ODEFUNC_RECVR);

  /* Initialize accumulated correction to 0 */
  N_VConst(ZERO, acor);

  /* Initialize delp to avoid compiler warning message */
  del = delp = ZERO;

  /* Loop until convergence; accumulate corrections in acor */
  loop {

    nni++;

    /* Correct y directly from the last f value */
    N_VLinearSum(h, tempv, -ONE, zn[1], tempv);
    N_VScale(rl1, tempv, tempv);
    N_VLinearSum(ONE, zn[0], ONE, tempv, y);

    /* Get WRMS norm of current correction to use in convergence test */
    N_VLinearSum(ONE, tempv, -ONE, acor, acor);
    del = N_VWrmsNorm(acor, ewt);

    /* Update accumulated correction */
    N_VScale(ONE, tempv, acor);

    /* Test for convergence.  If m > 0, an estimate of the convergence
       rate constant is stored in crate, and used in the test.        */
    if (m > 0) crate = MAX(NLS_CRDOWN * crate, del / delp);
    dcon = del * MIN(ONE, crate) / tq[4];
    if (dcon <= ONE) {
      acnrm = (m == 0) ? del : N_VWrmsNorm(acor, ewt);
      return(CP_SUCCESS);  /* Convergence achieved */
    }

    /* Stop at maxcor iterations or if iter. seems to be diverging */
    m++;
    if ((m==maxcor) || ((m >= 2) && (del > NLS_RDIV * delp))) return(CONV_FAIL);

    /* Save norm of correction, evaluate f, and loop again */
    delp = del;

    /* Evaluate f again */
    retval = fe(tn, y, tempv, f_data);
    nfe++;
    if (retval < 0) return(CP_ODEFUNC_FAIL);
    if (retval > 0) return(ODEFUNC_RECVR);

  }
}

/*
 * cpNlsNewtonExpl
 *
 * This routine handles the Newton iteration for an ODE in explicit form.
 * It calls lsetup if indicated, calls cpNewtonIterationExpl to perform
 * the actual Newton iteration, and retries a failed attempt at Newton
 * iteration if that is indicated.
 *
 * Possible return values:
 *
 *   CP_SUCCESS       ---> continue with error test
 *
 *   CP_ODEFUNC_FAIL  -+
 *   CP_LSETUP_FAIL    |-> halt the integration
 *   CP_LSOLVE_FAIL   -+
 *
 *   CONV_FAIL        -+
 *   ODEFUNC_RECVR    -+-> predict again or stop if too many
 *
 */

static int cpNlsNewtonExpl(CPodeMem cp_mem, int nflag)
{
  N_Vector vtemp1, vtemp2, vtemp3;
  int convfail, retval, ier;
  booleantype callSetup;

#ifdef CPODES_DEBUG
  printf("      Newton, explicit ODE\n");
#endif

  vtemp1 = acor;  /* rename acor as vtemp1 for readability  */
  vtemp2 = y;     /* rename y as vtemp2 for readability     */
  vtemp3 = tempv; /* rename tempv as vtemp3 for readability */

  /* Set flag convfail, input to lsetup for its evaluation decision */
  convfail = ((nflag == FIRST_CALL) || (nflag == PREV_ERR_FAIL)) ?
    CP_NO_FAILURES : CP_FAIL_OTHER;

  /* Decide whether or not to call setup routine (if one exists) */
  if (lsetup_exists) {
    callSetup = (nflag == PREV_CONV_FAIL) || (nflag == PREV_ERR_FAIL) ||
      (nst == 0) || (nst >= nstlset + NLS_MSBLS) || (ABS(gamrat-ONE) > DGMAX);
  } else {
    crate = ONE;
    callSetup = FALSE;
  }

  /* Begin the main loop. This loop is traversed at most twice.
   * The second pass only occurs when the first pass had a recoverable
   * failure with old Jacobian data.
   */
  loop {

    /* Evaluate function at predicted y */
    retval = fe(tn, zn[0], ftemp, f_data);
    nfe++;
    if (retval < 0) return(CP_ODEFUNC_FAIL);
    if (retval > 0) return(ODEFUNC_RECVR);

    /* If needed, call setup function (ypP = NULL in this case) */
    if (callSetup) {
      ier = lsetup(cp_mem, convfail, zn[0], NULL, ftemp, &jcur,
                   vtemp1, vtemp2, vtemp3);

#ifdef CPODES_DEBUG
      printf("         Linear solver setup return value = %d\n",ier);
#endif

      nsetups++;
      callSetup = FALSE;
      gamrat = crate = ONE;
      gammap = gamma;
      nstlset = nst;
      if (ier < 0) return(CP_LSETUP_FAIL);
      if (ier > 0) return(CONV_FAIL);
    }

    /* Set acor to zero and load prediction into y vector */
    N_VConst(ZERO, acor);
    N_VScale(ONE, zn[0], y);

    /* Do the Newton iteration */

#ifdef CPODES_DEBUG
    printf("         NonlinearIteration\n");
#endif

    ier = cpNewtonIterationExpl(cp_mem);

#ifdef CPODES_DEBUG
    printf("         NonlinearIteration return value = %d\n",ier);
#endif

    /* If there is a convergence failure and the Jacobian-related
       data appears not to be current, loop again with a call to lsetup
       in which convfail=CP_FAIL_BAD_J.  Otherwise return.                 */

    if ( (ier > 0) && (lsetup_exists) && (!jcur) ) {
      callSetup = TRUE;
      convfail = CP_FAIL_BAD_J;
      continue;
    }

    return(ier);

  }
}

/*
 * cpNewtonIterationExpl
 *
 * This routine performs the actual Newton iteration for an ODE system
 * in explicit form.
 *
 * Possible return values:
 *
 *   CP_SUCCESS  - the Newton iteration was successful.
 *
 *   CP_LSOLVE_FAIL - the linear solver solution failed unrecoverably.
 *   CP_ODEFUNC_FAIL - the ODE function failed unrecoverably.
 *
 *   ODEFUNC_RECVR - the ODE function failed recoverably.
 *   CONV_FAIL - a recoverable error occurred.
 */

static int cpNewtonIterationExpl(CPodeMem cp_mem)
{
  int m, retval;
  realtype del, delp, dcon;
  N_Vector b;

  mnewt = m = 0;

  /* Initialize delp to avoid compiler warning message */
  del = delp = ZERO;

  /* Looping point for Newton iteration */
  loop {

#ifdef CPODES_DEBUG
    printf("            Iteration # %d\n",m);
#endif
#ifdef CPODES_DEBUG_SERIAL
    printf("                zn[0] = "); N_VPrint_Serial(zn[0]);
    printf("                zn[1] = "); N_VPrint_Serial(zn[1]);
    printf("                ewt   = "); N_VPrint_Serial(ewt);
    printf("                acor  = "); N_VPrint_Serial(acor);
    printf("                ftemp = "); N_VPrint_Serial(ftemp);
#endif

    /* Evaluate the residual of the nonlinear system*/
    N_VLinearSum(rl1, zn[1], ONE, acor, tempv);
    N_VLinearSum(gamma, ftemp, -ONE, tempv, tempv);

    /* Call the lsolve function */
    b = tempv;

#ifdef CPODES_DEBUG
    printf("            Linear solver solve\n");
#endif
#ifdef CPODES_DEBUG_SERIAL
    printf("                rhs   = "); N_VPrint_Serial(b);
#endif

    retval = lsolve(cp_mem, b, ewt, y, NULL, ftemp);

#ifdef CPODES_DEBUG_SERIAL
    printf("                sol   = "); N_VPrint_Serial(b);
#endif
#ifdef CPODES_DEBUG
    printf("            Linear solver solve return value = %d\n",retval);
#endif

    nni++;
    if (retval < 0) return(CP_LSOLVE_FAIL);
    if (retval > 0) return(CONV_FAIL);

    /* Get WRMS norm of correction; add correction to acor and y */
    del = N_VWrmsNorm(b, ewt);

#ifdef CPODES_DEBUG
    printf("            Norm of correction:  del = %lg\n", del);
#endif

    N_VLinearSum(ONE, acor, ONE, b, acor);
    N_VLinearSum(ONE, zn[0], ONE, acor, y);

    /* Test for convergence.  If m > 0, an estimate of the convergence
       rate constant is stored in crate, and used in the test.        */
    if (m > 0) {
      crate = MAX(NLS_CRDOWN * crate, del/delp);
    }
    dcon = del * MIN(ONE, crate) / tq[4];

#ifdef CPODES_DEBUG
    printf("            Convergence test  dcon = %lg\n", dcon);
#endif

    if (dcon <= ONE) {
      acnrm = (m==0) ? del : N_VWrmsNorm(acor, ewt);

#ifdef CPODES_DEBUG_SERIAL
      printf("            acor = "); N_VPrint_Serial(acor);
#endif
#ifdef CPODES_DEBUG
      printf("            Accumulated correction norm = %lg\n", acnrm);
#endif

      jcur = FALSE;
      return(CP_SUCCESS); /* Nonlinear system was solved successfully */
    }

    mnewt = ++m;

    /* Stop at maxcor iterations or if iter. seems to be diverging. */
    if ((m == maxcor) || ((m >= 2) && (del > NLS_RDIV*delp))) return(CONV_FAIL);

    /* Save norm of correction, evaluate f, and loop again */
    delp = del;
    retval = fe(tn, y, ftemp, f_data);
    nfe++;
    if (retval < 0) return(CP_ODEFUNC_FAIL);
    if (retval > 0) return(ODEFUNC_RECVR);

  } /* end loop */
}


/*
 * =================================================================
 * Functions for implicit-form ODE
 * =================================================================
 */

/*
 * cpNlsFunctionalImpl
 *
 * This routine attempts to solve the nonlinear system for an ODE in
 * implicit form using functional iteration (no matrices involved).
 *
 * The system to be solved for y is written as
 *
 *    y = y - gamma * F( y, yP' + (y - yP)/gamma )
 *
 * At each iterate, y' is obtained as
 *
 *    y' = yP' + 1/gamma * ( y - yP )
 *
 * where yP and yP' are the predicted values for y and y'.
 *
 * The iterate updates are thus coded as:
 *    y' <- yP' + 1/gamma * ( y - yP )
 *    y  <- yP - gamma * F(y,y')
 *
 * Possible return values are:
 *
 *   CP_SUCCESS      --->  continue with error test
 *
 *   CP_ODEFUNC_FAIL --->  halt the integration
 *
 *   CONV_FAIL       -+
 *   ODEFUNC_RECVR   -+->  predict again or stop if too many
 *
 */

static int cpNlsFunctionalImpl(CPodeMem cp_mem)
{
  int retval, m;
  realtype del, delp, dcon;

#ifdef CPODES_DEBUG
  printf("      Functional, implicit ODE\n");
#endif

  /* Initialize counter */
  crate = ONE;
  m = 0;

  /* Evaluate f at predicted y and y' */
  N_VScale(ONE, zn[0], y);
  N_VScale(ONE/h, zn[1], yp);
  retval = fi(tn, y, yp, ftemp, f_data);
  nfe++;
  if (retval < 0) return(CP_ODEFUNC_FAIL);
  if (retval > 0) return(ODEFUNC_RECVR);

  /* Initialize accumulated correction to 0 */
  N_VConst(ZERO, acor);

  /* Initialize delp to avoid compiler warning message */
  del = delp = ZERO;

  /* Loop until convergence; accumulate corrections in acor */
  loop {

    nni++;

    /* Correct y directly from the last f value ( y <- y - gamma*f )  */
    N_VScale(-gamma, ftemp, tempv);
    N_VLinearSum(ONE, y, ONE, tempv, y);

    /* Get WRMS norm of current correction to use in convergence test */
    del = N_VWrmsNorm(tempv, ewt);

    /* Update accumulated correction */
    N_VLinearSum(ONE, acor, ONE, tempv, acor);

    /* Test for convergence.  If m > 0, an estimate of the convergence
       rate constant is stored in crate, and used in the test.        */
    if (m > 0) crate = MAX(NLS_CRDOWN * crate, del / delp);
    dcon = del * MIN(ONE, crate) / tq[4];
    if (dcon <= ONE) {
      acnrm = (m == 0) ? del : N_VWrmsNorm(acor, ewt);
      return(CP_SUCCESS);  /* Convergence achieved */
    }

    /* Stop at maxcor iterations or if iter. seems to be diverging */
    m++;
    if ((m==maxcor) || ((m >= 2) && (del > NLS_RDIV * delp))) return(CONV_FAIL);

    /* Save norm of correction */
    delp = del;

    /* Update y' and evaluate f again */
    N_VLinearSum(ONE, y, -ONE, zn[0], yp);
    N_VLinearSum(ONE/gamma, yp, ONE/h, zn[1], yp);
    retval = fi(tn, y, yp, ftemp, f_data);
    nfe++;
    if (retval < 0) return(CP_ODEFUNC_FAIL);
    if (retval > 0) return(ODEFUNC_RECVR);

  }
}

/*
 * cpNlsNewtonImpl
 *
 * This routine handles the Newton iteration for an ODE in implicit form.
 * It calls lsetup if indicated, calls cpNewtonIterationImpl to perform
 * the actual Newton iteration, and retries a failed attempt at Newton
 * iteration if that is indicated.
 *
 * Possible return values:
 *
 *   CP_SUCCESS       ---> continue with error test
 *
 *   CP_ODEFUNC_FAIL  -+
 *   CP_LSETUP_FAIL    |-> halt the integration
 *   CP_LSOLVE_FAIL   -+
 *
 *   CONV_FAIL        -+
 *   ODEFUNC_RECVR    -+-> predict again or stop if too many
 *
 */

static int cpNlsNewtonImpl(CPodeMem cp_mem, int nflag)
{
  N_Vector vtemp1, vtemp2, vtemp3;
  int convfail, retval;
  booleantype callSetup;

#ifdef CPODES_DEBUG
  printf("      Newton, implicit ODE\n");
#endif

  /* In this situation we do not pass convfail information to
     the linear sovler setup function */
  convfail = 0;

  /* Vectors used as temporary space in lsetup */
  vtemp1 = acor;
  vtemp2 = y;
  vtemp3 = tempv;

  /* If the linear solver provides a setup function, call it if:
   *   - we are at the first step (nst == 0), or
   *   - gamma changed significantly, or
   *   - enough steps passed from last evaluation
   */
  if (lsetup_exists) {
    callSetup = (nst == 0) || (ABS(gamrat-ONE) > DGMAX) || (nst >= nstlset + NLS_MSBLS);
  } else {
    callSetup = FALSE;
    crate = ONE;
  }


  /* Begin the main loop. This loop is traversed at most twice.
   * The second pass only occurs when the first pass had a recoverable
   * failure with old Jacobian data.
   */
  loop {

    /* Evaluate predicted y' from Nordsieck array */
    N_VScale(ONE/h, zn[1], yp);

    /* Evaluate residual at predicted y and y' */
    retval = fi(tn, zn[0], yp, ftemp, f_data);
    nfe++;
    if (retval < 0) return(CP_ODEFUNC_FAIL);
    if (retval > 0) return(ODEFUNC_RECVR);

    /* If needed, call setup function */
    if (callSetup) {
      retval = lsetup(cp_mem, 0, zn[0], yp, ftemp, &jcur,
                      vtemp1, vtemp2, vtemp3);
      nsetups++;
      nstlset = nst;
      crate = ONE;
      gammap = gamma;
      gamrat = ONE;
      if (retval < 0) return(CP_LSETUP_FAIL);
      if (retval > 0) return(CONV_FAIL);
    }

    /* Set acor to zero and load prediction into y vector */
    N_VConst(ZERO, acor);
    N_VScale(ONE, zn[0], y);

    /* Do the Newton iteration */
    retval = cpNewtonIterationImpl(cp_mem);

    /* On a recoverable failure, if setup was not called, reattempt loop */
    if ( (retval > 0) && lsetup_exists && !callSetup ) {
      callSetup = TRUE;
      continue;
    }

    /* Return (success or unrecoverable failure) */
    return(retval);

  }
}

/*
 * cpNewtonIterationImpl
 *
 * This routine performs actual the Newton iteration for an ODE system
 * given in implicit form.
 *
 * Possible return values:
 *
 *   CP_SUCCESS  - the Newton iteration was successful.
 *
 *   CP_LSOLVE_FAIL - the linear solver solution failed unrecoverably.
 *   CP_ODEFUNC_FAIL - the ODE function failed unrecoverably.
 *
 *   ODEFUNC_RECVR - the ODE function failed recoverably.
 *   CONV_FAIL - a recoverable error occurred.
 */

static int cpNewtonIterationImpl(CPodeMem cp_mem)
{
  int m, retval;
  realtype del, delp, dcon;

  /* Initialize counters */
  mnewt = m = 0;

  /* Initialize delp to avoid compiler warning message */
  del = delp = ZERO;

  /* Looping point for Newton iteration */
  loop {

    /* Evaluate the residual of the nonlinear system*/
    N_VScale(-gamma, ftemp, tempv);

    /* Call the lsolve function (solution in tempv) */
    retval = lsolve(cp_mem, tempv, ewt, y, yp, ftemp);
    nni++;
    if (retval < 0) return(CP_LSOLVE_FAIL);
    if (retval > 0) return(CONV_FAIL);

    /* Get WRMS norm of correction and add correction to acor, y, and yp */
    del = N_VWrmsNorm(tempv, ewt);
    N_VLinearSum(ONE, acor, ONE,       tempv, acor);
    N_VLinearSum(ONE, y,    ONE,       tempv, y);
    N_VLinearSum(ONE, yp,   ONE/gamma, tempv, yp);

    /* Test for convergence.  If m > 0, an estimate of the convergence
       rate constant is stored in crate, and used in the test.        */
    if (m > 0) crate = MAX(NLS_CRDOWN * crate, del/delp);
    dcon = del * MIN(ONE, crate) / tq[4];
    if (dcon <= ONE) {
      acnrm = (m==0) ? del : N_VWrmsNorm(acor, ewt);
      return(CP_SUCCESS);
    }

    mnewt = ++m;

    /* Stop at maxcor iterations or if iter. seems to be diverging. */
    if ((m == maxcor) || ((m >= 2) && (del > NLS_RDIV*delp)))  return(CONV_FAIL);

    /* Save norm of correction, evaluate f, and loop again */
    delp = del;
    retval = fi(tn, y, yp, ftemp, f_data);
    nfe++;
    if (retval < 0) return(CP_ODEFUNC_FAIL);
    if (retval > 0) return(ODEFUNC_RECVR);

  } /* end loop */
}

