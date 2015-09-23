/*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2007/04/06 20:33:24 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for the optional input and output
 * functions for the CPODES solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "cpodes_private.h"

#define lrw   (cp_mem->cp_lrw)
#define liw   (cp_mem->cp_liw)
#define lrw1  (cp_mem->cp_lrw1)
#define liw1  (cp_mem->cp_liw1)
#define lrw1Q (cp_mem->cp_lrw1Q)
#define liw1Q (cp_mem->cp_liw1Q)

/*
 * =================================================================
 * CPODES optional input functions
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Optional inputs to the main solver
 * -----------------------------------------------------------------
 */

/*
 * CPodeSetErrHandlerFn
 *
 * Specifies the error handler function
 */

int CPodeSetErrHandlerFn(void *cpode_mem, CPErrHandlerFn ehfun, void *eh_data)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeSetErrHandlerFn", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  cp_mem->cp_ehfun = ehfun;
  cp_mem->cp_eh_data = eh_data;

  return(CP_SUCCESS);
}

/*
 * CPodeSetErrFile
 *
 * Specifies the FILE pointer for output (NULL means no messages)
 */

int CPodeSetErrFile(void *cpode_mem, FILE *errfp)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeSetErrFile", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  cp_mem->cp_errfp = errfp;

  return(CP_SUCCESS);
}

/*
 * CPodeSetMaxOrd
 *
 * Specifies the maximum method order
 */

int CPodeSetMaxOrd(void *cpode_mem, int maxord)
{
  CPodeMem cp_mem;
  int qmax_alloc;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeSetMaxOrd", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (maxord <= 0) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeSetMaxOrd", MSGCP_NEG_MAXORD);
    return(CP_ILL_INPUT);
  }

  /* Cannot increase maximum order beyond the value that
     was used when allocating memory */
  qmax_alloc = cp_mem->cp_qmax_alloc;

  if (maxord > qmax_alloc) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeSetMaxOrd", MSGCP_BAD_MAXORD);
    return(CP_ILL_INPUT);
  }

  cp_mem->cp_qmax = maxord;

  return(CP_SUCCESS);
}

/*
 * CPodeSetMaxNumSteps
 *
 * Specifies the maximum number of integration steps
 */

int CPodeSetMaxNumSteps(void *cpode_mem, long int mxsteps)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeSetMaxNumSteps", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (mxsteps < 0) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeSetMaxNumSteps", MSGCP_NEG_MXSTEPS);
    return(CP_ILL_INPUT);
  }

  /* Passing 0 sets the default */
  if (mxsteps == 0)
    cp_mem->cp_mxstep = MXSTEP_DEFAULT;
  else
    cp_mem->cp_mxstep = mxsteps;

  return(CP_SUCCESS);
}

/*
 * CPodeSetMaxHnilWarns
 *
 * Specifies the maximum number of warnings for small h
 */

int CPodeSetMaxHnilWarns(void *cpode_mem, int mxhnil)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeSetMaxHnilWarns", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  cp_mem->cp_mxhnil = mxhnil;

  return(CP_SUCCESS);
}

/*
 *CPodeSetStabLimDet
 *
 * Turns on/off the stability limit detection algorithm
 */

int CPodeSetStabLimDet(void *cpode_mem, booleantype sldet)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeSetStabLimDet", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if( sldet && (cp_mem->cp_lmm_type != CP_BDF) ) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeSetStabLimDet", MSGCP_SET_SLDET);
    return(CP_ILL_INPUT);
  }

  cp_mem->cp_sldeton = sldet;

  return(CP_SUCCESS);
}

/*
 * CPodeSetInitStep
 *
 * Specifies the initial step size
 */

int CPodeSetInitStep(void *cpode_mem, realtype hin)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeSetInitStep", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  cp_mem->cp_hin = hin;

  return(CP_SUCCESS);
}

/*
 * CPodeSetMinStep
 *
 * Specifies the minimum step size
 */

int CPodeSetMinStep(void *cpode_mem, realtype hmin)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeSetMinStep", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (hmin<0) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeSetMinStep", MSGCP_NEG_HMIN);
    return(CP_ILL_INPUT);
  }

  /* Passing 0 sets hmin = zero */
  if (hmin == ZERO) {
    cp_mem->cp_hmin = HMIN_DEFAULT;
    return(CP_SUCCESS);
  }

  if (hmin * cp_mem->cp_hmax_inv > ONE) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeSetMinStep", MSGCP_BAD_HMIN_HMAX);
    return(CP_ILL_INPUT);
  }

  cp_mem->cp_hmin = hmin;

  return(CP_SUCCESS);
}

/*
 * CPodeSetMaxStep
 *
 * Specifies the maximum step size
 */

int CPodeSetMaxStep(void *cpode_mem, realtype hmax)
{
  realtype hmax_inv;
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeSetMaxStep", MSGCP_NO_MEM);
    return (CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (hmax < 0) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeSetMaxStep", MSGCP_NEG_HMAX);
    return(CP_ILL_INPUT);
  }

  /* Passing 0 sets hmax = infinity */
  if (hmax == ZERO) {
    cp_mem->cp_hmax_inv = HMAX_INV_DEFAULT;
    return(CP_SUCCESS);
  }

  hmax_inv = ONE/hmax;
  if (hmax_inv * cp_mem->cp_hmin > ONE) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeSetMaxStep", MSGCP_BAD_HMIN_HMAX);
    return(CP_ILL_INPUT);
  }

  cp_mem->cp_hmax_inv = hmax_inv;

  return(CP_SUCCESS);
}

/*
 * CPodeSetStopTime
 *
 * Specifies the time beyond which the integration is not to
 * proceed
 */

int CPodeSetStopTime(void *cpode_mem, realtype tstop)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeSetStopTime", MSGCP_NO_MEM);
    return (CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  cp_mem->cp_tstop = tstop;
  cp_mem->cp_tstopset = TRUE;

  return(CP_SUCCESS);
}

/*
 * CPodeSetMaxErrTestFails
 *
 * Specifies the maximum number of error test failures during one
 * step try.
 */

int CPodeSetMaxErrTestFails(void *cpode_mem, int maxnef)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeSetMaxErrTestFails", MSGCP_NO_MEM);
    return (CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  cp_mem->cp_maxnef = maxnef;

  return(CP_SUCCESS);
}

/*
 * CPodeSetMaxConvFails
 *
 * Specifies the maximum number of nonlinear convergence failures
 * during one step try.
 */

int CPodeSetMaxConvFails(void *cpode_mem, int maxncf)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeSetMaxConvFails", MSGCP_NO_MEM);
    return (CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  cp_mem->cp_maxncf = maxncf;

  return(CP_SUCCESS);
}

/*
 * CPodeSetMaxNonlinIters
 *
 * Specifies the maximum number of nonlinear iterations during
 * one solve.
 */

int CPodeSetMaxNonlinIters(void *cpode_mem, int maxcor)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeSetMaxNonlinIters", MSGCP_NO_MEM);
    return (CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  cp_mem->cp_maxcor = maxcor;

  return(CP_SUCCESS);
}

/*
 * CPodeSetNonlinConvCoef
 *
 * Specifies the coeficient in the nonlinear solver convergence
 * test
 */

int CPodeSetNonlinConvCoef(void *cpode_mem, realtype nlscoef)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeSetNonlinConvCoef", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  cp_mem->cp_nlscoef = nlscoef;

  return(CP_SUCCESS);
}

/*
 * CPodeSetTolerances
 *
 * Changes the integration tolerances between calls to CPode().
 * Here, only CP_SS or CP_SV are allowed.
 */

int CPodeSetTolerances(void *cpode_mem,
                       int tol_type, realtype reltol, void *abstol)
{
  CPodeMem cp_mem;
  booleantype neg_abstol;

  /* Check CPODES memory */
  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeSetTolerances", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  /* Check if cpode_mem was allocated */
  if (cp_mem->cp_MallocDone == FALSE) {
    cpProcessError(cp_mem, CP_NO_MALLOC, "CPODES", "CPodeSetTolerances", MSGCP_NO_MALLOC);
    return(CP_NO_MALLOC);
  }

  /* Check inputs for legal values */
  if ( (tol_type != CP_SS) && (tol_type != CP_SV) ) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeSetTolerances", MSGCP_BAD_ITOL);
    return(CP_ILL_INPUT);
  }
  if (abstol == NULL) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeSetTolerances", MSGCP_NULL_ABSTOL);
    return(CP_ILL_INPUT);
  }

  /* Check positivity of tolerances */
  if (reltol < ZERO) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeSetTolerances", MSGCP_BAD_RELTOL);
    return(CP_ILL_INPUT);
  }

  if (tol_type == CP_SS)
    neg_abstol = (*((realtype *)abstol) < ZERO);
  else
    neg_abstol = (N_VMin((N_Vector)abstol) < ZERO);

  if (neg_abstol) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeSetTolerances", MSGCP_BAD_ABSTOL);
    return(CP_ILL_INPUT);
  }

  /* Copy tolerances into memory */
  cp_mem->cp_tol_type = tol_type;
  cp_mem->cp_reltol   = reltol;

  if (tol_type == CP_SS) {
    cp_mem->cp_Sabstol = *((realtype *)abstol);
  } else {
    if ( !(cp_mem->cp_VabstolMallocDone) ) {
      cp_mem->cp_Vabstol = N_VClone(cp_mem->cp_ewt);
      lrw += lrw1;
      liw += liw1;
      cp_mem->cp_VabstolMallocDone = TRUE;
    }
    N_VScale(ONE, (N_Vector)abstol, cp_mem->cp_Vabstol);
  }

  return(CP_SUCCESS);
}

/*
 * CPodeSetEwtFn
 *
 * Specifies the user-provided function efun of type EwtSet
 * and a data pointer for efun.
 */

int CPodeSetEwtFn(void *cpode_mem, CPEwtFn efun, void *e_data)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeSetEwtFn", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  cp_mem->cp_efun   = efun;
  cp_mem->cp_e_data = e_data;

  return(CP_SUCCESS);
}

/*
 * CPodeSetRootDirection
 *
 * Specifies the direction of zero-crossings to be monitored.
 * The default is to monitor both crossings.
 */

int CPodeSetRootDirection(void *cpode_mem, int *rootdir)
{
  CPodeMem cp_mem;
  int i, nrt;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeSetRootDirection", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }

  cp_mem = (CPodeMem) cpode_mem;

  nrt = cp_mem->cp_nrtfn;
  if (nrt==0) {
    cpProcessError(NULL, CP_ILL_INPUT, "CPODES", "CPodeSetRootDirection", MSGCP_NO_ROOT);
    return(CP_ILL_INPUT);
  }

  for(i=0; i<nrt; i++) cp_mem->cp_rootdir[i] = rootdir[i];

  return(CP_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Optional inputs for projection step
 * -----------------------------------------------------------------
 */

/*
 *
 * CPodeSetProjUpdateErrEst
 */

int CPodeSetProjUpdateErrEst(void *cpode_mem, booleantype proj_err)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeSetProjUpdateErrEst", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  cp_mem->cp_project_err = proj_err;

  return(CP_SUCCESS);
}

/*
 * CPodeSetProjFrequency
 */

int CPodeSetProjFrequency(void *cpode_mem, long int proj_freq)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeSetProjFrequency", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (proj_freq < 0) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeSetProjFrequency", MSGCP_BAD_FREQ);
    return(CP_ILL_INPUT);
  }

  cp_mem->cp_proj_freq = proj_freq;

  return(CP_SUCCESS);
}

/*
 * CPodeSetProjTestCnstr
 */

int CPodeSetProjTestCnstr(void *cpode_mem, booleantype test_cnstr)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeSetProjTestCnstr", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  cp_mem->cp_test_cnstr = test_cnstr;

  return(CP_SUCCESS);
}

/*
 * CPodeSetProjLsetupFreq
 *
 * Spcifies the frequency of constraint Jacobian evaluations
 * (actually, the maximum number of steps allowed between calls
 * to the linear solver's setup function).
 */

int CPodeSetProjLsetupFreq(void *cpode_mem, long int lset_freq)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeSetProjLsetupFreq", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lset_freq <= 0) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeSetProjLsetupFreq", MSGCP_BAD_LSFREQ);
    return(CP_ILL_INPUT);
  }

  cp_mem->cp_lsetupP_freq = lset_freq;

  return(CP_SUCCESS);
}

/*
 * CPodeSetProjNonlinConvCoef
 *
 * Specifies the coeficient in the projection nonlinear solver convergence
 * test
 */

int CPodeSetProjNonlinConvCoef(void *cpode_mem, realtype prjcoef)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeSetProjNonlinConvCoef", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  cp_mem->cp_prjcoef = prjcoef;

  return(CP_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Optional inputs for quadrature integration
 * -----------------------------------------------------------------
 */

int CPodeSetQuadErrCon(void *cpode_mem, booleantype errconQ,
                       int tol_typeQ, realtype reltolQ, void *abstolQ)
{
  CPodeMem cp_mem;
  booleantype neg_abstol;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeSetQuadErrCon", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }

  cp_mem = (CPodeMem) cpode_mem;

  cp_mem->cp_errconQ = errconQ;

  /* Ckeck if quadrature was initialized? */

  if (cp_mem->cp_quadMallocDone == FALSE) {
    cpProcessError(cp_mem, CP_NO_QUAD, "CPODES", "CPodeSetQuadErrCon", MSGCP_NO_QUAD);
    return(CP_NO_QUAD);
  }

  /* Check inputs */

  if(errconQ == FALSE) {
    if (cp_mem->cp_VabstolQMallocDone) {
      N_VDestroy(cp_mem->cp_VabstolQ);
      lrw -= lrw1Q;
      liw -= liw1Q;
      cp_mem->cp_VabstolQMallocDone = FALSE;
    }
    return(CP_SUCCESS);
  }

  if ((tol_typeQ != CP_SS) && (tol_typeQ != CP_SV)) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeSetQuadErrCon", MSGCP_BAD_ITOLQ);
    return(CP_ILL_INPUT);
  }

  if (abstolQ == NULL) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeSetQuadErrCon", MSGCP_NULL_ABSTOLQ);
    return(CP_ILL_INPUT);
  }

  if (reltolQ < ZERO) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeSetQuadErrCon", MSGCP_BAD_RELTOLQ);
    return(CP_ILL_INPUT);
  }

  if (tol_typeQ == CP_SS)
    neg_abstol = (*((realtype *)abstolQ) < ZERO);
  else
    neg_abstol = (N_VMin((N_Vector)abstolQ) < ZERO);

  if (neg_abstol) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeSetQuadErrCon", MSGCP_BAD_ABSTOLQ);
    return(CP_ILL_INPUT);
  }

  /* See if we need to free or allocate memory */

  if ( (tol_typeQ != CP_SV) && (cp_mem->cp_VabstolQMallocDone) ) {
    N_VDestroy(cp_mem->cp_VabstolQ);
    lrw -= lrw1Q;
    liw -= liw1Q;
    cp_mem->cp_VabstolQMallocDone = FALSE;
  }

  if ( (tol_typeQ == CP_SV) && !(cp_mem->cp_VabstolQMallocDone) ) {
    cp_mem->cp_VabstolQ = N_VClone(cp_mem->cp_tempvQ);
    lrw += lrw1Q;
    liw += liw1Q;
    cp_mem->cp_VabstolQMallocDone = TRUE;
  }

  /* Copy tolerances into memory */

  cp_mem->cp_tol_typeQ = tol_typeQ;
  cp_mem->cp_reltolQ   = reltolQ;

  if (tol_typeQ == CP_SS)
    cp_mem->cp_SabstolQ = *((realtype *)abstolQ);
  else
    N_VScale(ONE, (N_Vector)abstolQ, cp_mem->cp_VabstolQ);

  return(CP_SUCCESS);
}

/*
 * =================================================================
 * CPODE optional output functions
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Optional outputs for IC calculation
 * -----------------------------------------------------------------
 */

int CPodeGetConsistentIC(void *cpode_mem, N_Vector yy0, N_Vector yp0)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeGetConsistentIC", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (cp_mem->cp_next_q != 0) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeGetConsistentIC", MSGCP_TOO_LATE);
    return(CP_ILL_INPUT);
  }

  if (yy0 != NULL) N_VScale(ONE, cp_mem->cp_zn[0], yy0);
  if (cp_mem->cp_ode_type == CP_IMPL)
    if (yp0 != NULL) N_VScale(ONE, cp_mem->cp_zn[1], yp0);

  return(CP_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Optional outputs for integration
 * -----------------------------------------------------------------
 */

/*
 * CPodeGetNumSteps
 *
 * Returns the current number of integration steps
 */

int CPodeGetNumSteps(void *cpode_mem, long int *nsteps)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeGetNumSteps", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  *nsteps = cp_mem->cp_nst;

  return(CP_SUCCESS);
}

/*
 * CPodeGetNumFctEvals
 *
 * Returns the current number of calls to f
 */

int CPodeGetNumFctEvals(void *cpode_mem, long int *nfevals)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeGetNumFctEvals", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  *nfevals = cp_mem->cp_nfe;

  return(CP_SUCCESS);
}

/*
 * CPodeGetNumLinSolvSetups
 *
 * Returns the current number of calls to the linear solver setup routine
 */

int CPodeGetNumLinSolvSetups(void *cpode_mem, long int *nlinsetups)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeGetNumLinSolvSetups", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  *nlinsetups = cp_mem->cp_nsetups;

  return(CP_SUCCESS);
}

/*
 * CPodeGetNumErrTestFails
 *
 * Returns the current number of error test failures
 */

int CPodeGetNumErrTestFails(void *cpode_mem, long int *netfails)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeGetNumErrTestFails", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  *netfails = cp_mem->cp_netf;

  return(CP_SUCCESS);
}

/*
 * CPodeGetLastOrder
 *
 * Returns the order on the last successful step
 */

int CPodeGetLastOrder(void *cpode_mem, int *qlast)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeGetLastOrder", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  *qlast = cp_mem->cp_qu;

  return(CP_SUCCESS);
}

/*
 * CPodeGetCurrentOrder
 *
 * Returns the order to be attempted on the next step
 */

int CPodeGetCurrentOrder(void *cpode_mem, int *qcur)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeGetCurrentOrder", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  *qcur = cp_mem->cp_next_q;

  return(CP_SUCCESS);
}

/*
 * CPodeGetNumStabLimOrderReds
 *
 * Returns the number of order reductions triggered by the stability
 * limit detection algorithm
 */

int CPodeGetNumStabLimOrderReds(void *cpode_mem, long int *nslred)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeGetNumStabLimOrderReds", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (cp_mem->cp_sldeton==FALSE)
    *nslred = 0;
  else
    *nslred = cp_mem->cp_nor;

  return(CP_SUCCESS);
}

/*
 * CPodeGetActualInitStep
 *
 * Returns the step size used on the first step
 */

int CPodeGetActualInitStep(void *cpode_mem, realtype *hinused)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeGetActualInitStep", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  *hinused = cp_mem->cp_h0u;

  return(CP_SUCCESS);
}

/*
 * CPodeGetLastStep
 *
 * Returns the step size used on the last successful step
 */

int CPodeGetLastStep(void *cpode_mem, realtype *hlast)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeGetLastStep", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  *hlast = cp_mem->cp_hu;

  return(CP_SUCCESS);
}

/*
 * CPodeGetCurrentStep
 *
 * Returns the step size to be attempted on the next step
 */

int CPodeGetCurrentStep(void *cpode_mem, realtype *hcur)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeGetCurrentStep", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  *hcur = cp_mem->cp_next_h;

  return(CP_SUCCESS);
}

/*
 * CPodeGetCurrentTime
 *
 * Returns the current value of the independent variable
 */

int CPodeGetCurrentTime(void *cpode_mem, realtype *tcur)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeGetCurrentTime", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  *tcur = cp_mem->cp_tn;

  return(CP_SUCCESS);
}

/*
 * CPodeGetTolScaleFactor
 *
 * Returns a suggested factor for scaling tolerances
 */

int CPodeGetTolScaleFactor(void *cpode_mem, realtype *tolsfact)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeGetTolScaleFactor", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  *tolsfact = cp_mem->cp_tolsf;

  return(CP_SUCCESS);
}

/*
 * CPodeGetErrWeights
 *
 * This routine returns the current weight vector.
 */

int CPodeGetErrWeights(void *cpode_mem, N_Vector eweight)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeGetErrWeights", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  N_VScale(ONE, cp_mem->cp_ewt, eweight);

  return(CP_SUCCESS);
}

/*
 * CPodeGetEstLocalErrors
 *
 * Returns an estimate of the local error
 */

int CPodeGetEstLocalErrors(void *cpode_mem, N_Vector ele)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeGetEstLocalErrors", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  N_VScale(ONE, cp_mem->cp_acor, ele);

  return(CP_SUCCESS);
}

/*
 * CPodeGetWorkSpace
 *
 * Returns integrator work space requirements
 */

int CPodeGetWorkSpace(void *cpode_mem, long int *lenrw, long int *leniw)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeGetWorkSpace", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  *leniw = liw;
  *lenrw = lrw;

  return(CP_SUCCESS);
}

/*
 * CPodeGetIntegratorStats
 *
 * Returns integrator statistics
 */

int CPodeGetIntegratorStats(void *cpode_mem, long int *nsteps, long int *nfevals,
                            long int *nlinsetups, long int *netfails, int *qlast,
                            int *qcur, realtype *hinused, realtype *hlast,
                            realtype *hcur, realtype *tcur)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeGetIntegratorStats", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  *nsteps = cp_mem->cp_nst;
  *nfevals = cp_mem->cp_nfe;
  *nlinsetups = cp_mem->cp_nsetups;
  *netfails = cp_mem->cp_netf;
  *qlast = cp_mem->cp_qu;
  *qcur = cp_mem->cp_next_q;
  *hinused = cp_mem->cp_h0u;
  *hlast = cp_mem->cp_hu;
  *hcur = cp_mem->cp_next_h;
  *tcur = cp_mem->cp_tn;

  return(CP_SUCCESS);
}

/*
 * CPodeGetNumGEvals
 *
 * Returns the current number of calls to g (for rootfinding)
 */

int CPodeGetNumGEvals(void *cpode_mem, long int *ngevals)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeGetNumGEvals", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  *ngevals = cp_mem->cp_nge;

  return(CP_SUCCESS);
}

/*
 * CPodeGetRootInfo
 *
 * Returns pointer to array rootsfound showing roots found
 */

int CPodeGetRootInfo(void *cpode_mem, int *rootsfound)
{
  CPodeMem cp_mem;
  int i, nrt;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeGetRootInfo", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  nrt = cp_mem->cp_nrtfn;

  for (i=0; i<nrt; i++) rootsfound[i] = cp_mem->cp_iroots[i];

  return(CP_SUCCESS);
}

/*
 * CPodeGetRootWindow
 *
 * Returns the window (tLo,tHi] within which root(s) were found. The result
 * is meaningless unless this is called right after a root-found return.
 */

int CPodeGetRootWindow(void *cpode_mem, realtype *tLo, realtype *tHi)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeGetRootWindow", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  *tLo = cp_mem->cp_tlo;
  *tHi = cp_mem->cp_thi;

  return(CP_SUCCESS);
}

/*
 * CPodeGetNumNonlinSolvIters
 *
 * Returns the current number of iterations in the nonlinear solver
 */

int CPodeGetNumNonlinSolvIters(void *cpode_mem, long int *nniters)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeGetNumNonlinSolvIters", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  *nniters = cp_mem->cp_nni;

  return(CP_SUCCESS);
}

/*
 * CPodeGetNumNonlinSolvConvFails
 *
 * Returns the current number of convergence failures in the
 * nonlinear solver
 */

int CPodeGetNumNonlinSolvConvFails(void *cpode_mem, long int *nncfails)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeGetNumNonlinSolvConvFails", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  *nncfails = cp_mem->cp_ncfn;

  return(CP_SUCCESS);
}

/*
 * CPodeGetNonlinSolvStats
 *
 * Returns nonlinear solver statistics
 */

int CPodeGetNonlinSolvStats(void *cpode_mem, long int *nniters,
                            long int *nncfails)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeGetNonlinSolvStats", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  *nniters = cp_mem->cp_nni;
  *nncfails = cp_mem->cp_ncfn;

  return(CP_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Optional outputs for projection step
 * -----------------------------------------------------------------
 */

int CPodeGetProjNumProj(void *cpode_mem, long int *nproj)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeGetNonlinSolvStats", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  *nproj    = cp_mem->cp_nproj;

  return(CP_SUCCESS);
}

int CPodeGetProjNumCnstrEvals(void *cpode_mem, long int *nce)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeGetNonlinSolvStats", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  *nce      = cp_mem->cp_nce;

  return(CP_SUCCESS);
}

int CPodeGetProjNumLinSolvSetups(void *cpode_mem, long int *nsetupsP)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeGetNonlinSolvStats", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  *nsetupsP = cp_mem->cp_nsetupsP;

  return(CP_SUCCESS);
}

int CPodeGetProjNumFailures(void *cpode_mem, long int *nprf)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeGetNonlinSolvStats", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  *nprf = cp_mem->cp_nprf;

  return(CP_SUCCESS);
}

int CPodeGetProjStats(void *cpode_mem, long int *nproj, long int *nce,
                      long int *nsetupsP, long int *nprf)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeGetNonlinSolvStats", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  *nproj    = cp_mem->cp_nproj;
  *nce      = cp_mem->cp_nce;
  *nsetupsP = cp_mem->cp_nsetupsP;
  *nprf     = cp_mem->cp_nprf;

  return(CP_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Optional outputs for quadrature integration
 * -----------------------------------------------------------------
 */

int CPodeGetQuadNumFunEvals(void *cpode_mem, long int *nqevals)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeGetQuadNumFunEvals", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }

  cp_mem = (CPodeMem) cpode_mem;

  if (cp_mem->cp_quadr==FALSE) {
    cpProcessError(cp_mem, CP_NO_QUAD, "CPODES", "CPodeGetQuadNumFunEvals", MSGCP_NO_QUAD);
    return(CP_NO_QUAD);
  }

  *nqevals = cp_mem->cp_nqe;

  return(CP_SUCCESS);
}

int CPodeGetQuadErrWeights(void *cpode_mem, N_Vector eQweight)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeGetQuadErrWeights", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }

  cp_mem = (CPodeMem) cpode_mem;

  if (cp_mem->cp_quadr==FALSE) {
    cpProcessError(cp_mem, CP_NO_QUAD, "CPODES", "CPodeGetQuadErrWeights", MSGCP_NO_QUAD);
    return(CP_NO_QUAD);
  }

  if(cp_mem->cp_errconQ) N_VScale(ONE, cp_mem->cp_ewtQ, eQweight);

  return(CP_SUCCESS);
}


/*-----------------------------------------------------------------*/

char *CPodeGetReturnFlagName(int flag)
{
  char *name;

  name = (char *)malloc(24*sizeof(char));

  switch(flag) {
  case CP_SUCCESS:
    sprintf(name,"CP_SUCCESS");
    break;
  case CP_TSTOP_RETURN:
    sprintf(name,"CP_TSTOP_RETURN");
    break;
  case CP_ROOT_RETURN:
    sprintf(name,"CP_ROOT_RETURN");
    break;
  case CP_TOO_MUCH_WORK:
    sprintf(name,"CP_TOO_MUCH_WORK");
    break;
  case CP_TOO_MUCH_ACC:
    sprintf(name,"CP_TOO_MUCH_ACC");
    break;
  case CP_ERR_FAILURE:
    sprintf(name,"CP_ERR_FAILURE");
    break;
  case CP_CONV_FAILURE:
    sprintf(name,"CP_CONV_FAILURE");
    break;
  case CP_LINIT_FAIL:
    sprintf(name,"CP_LINIT_FAIL");
    break;
  case CP_LSETUP_FAIL:
    sprintf(name,"CP_LSETUP_FAIL");
    break;
  case CP_LSOLVE_FAIL:
    sprintf(name,"CP_LSOLVE_FAIL");
    break;
  case CP_ODEFUNC_FAIL:
    sprintf(name,"CP_ODEFUNC_FAIL");
    break;
  case CP_FIRST_ODEFUNC_ERR:
    sprintf(name,"CP_FIRST_ODEFUNC_ERR");
    break;
  case CP_REPTD_ODEFUNC_ERR:
    sprintf(name,"CP_REPTD_ODEFUNC_ERR");
    break;
  case CP_UNREC_ODEFUNC_ERR:
    sprintf(name,"CP_UNREC_ODEFUNC_ERR");
    break;
  case CP_RTFUNC_FAIL:
    sprintf(name,"CP_RTFUNC_FAIL");
    break;
  case CP_MEM_FAIL:
    sprintf(name,"CP_MEM_FAIL");
    break;
  case CP_MEM_NULL:
    sprintf(name,"CP_MEM_NULL");
    break;
  case CP_ILL_INPUT:
    sprintf(name,"CP_ILL_INPUT");
    break;
  case CP_NO_MALLOC:
    sprintf(name,"CP_NO_MALLOC");
    break;
  case CP_BAD_K:
    sprintf(name,"CP_BAD_K");
    break;
  case CP_BAD_T:
    sprintf(name,"CP_BAD_T");
    break;
  case CP_BAD_DKY:
    sprintf(name,"CP_BAD_DKY");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}

