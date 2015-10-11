/*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2007/01/29 17:35:23 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This file contains implementations of the banded difference
 * quotient Jacobian-based preconditioner and solver routines for
 * use with the CPSPILS linear solvers.
 * -----------------------------------------------------------------
 */

/*
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include <cpodes/cpodes_sptfqmr.h>
#include <cpodes/cpodes_spbcgs.h>
#include <cpodes/cpodes_spgmr.h>

#include <sundials/sundials_math.h>

#include "cpodes_bandpre_impl.h"
#include "cpodes_private.h"

/*
 * =================================================================
 * FUNCTION SPECIFIC CONSTANTS
 * =================================================================
 */

#define MIN_INC_MULT RCONST(1000.0)

/*
 * =================================================================
 * PROTOTYPES FOR PRIVATE FUNCTIONS
 * =================================================================
 */

/* cpBandPrecSetupExpl and cpBandPrecSetupImpl */

static int cpBandPrecSetupExpl(realtype t, N_Vector y, N_Vector fy,
                               booleantype jok, booleantype *jcurPtr,
                               realtype gamma, void *bp_data,
                               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int cpBandPrecSetupImpl(realtype t, N_Vector y, N_Vector yp, N_Vector r,
                               realtype gamma, void *bp_data,
                               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* cpBandPrecSolveExpl and cpBandPrecSolveImpl */

static int cpBandPrecSolveExpl(realtype t, N_Vector y, N_Vector fy,
                               N_Vector b, N_Vector x,
                               realtype gamma, realtype delta,
                               int lr, void *bp_data, N_Vector tmp);

static int cpBandPrecSolveImpl(realtype t, N_Vector y, N_Vector yp, N_Vector r,
                               N_Vector b, N_Vector x,
                               realtype gamma, realtype delta,
                               void *bp_data, N_Vector tmp);

/* Difference quotient Jacobian calculation routines */

static int cpBandPDQJacExpl(CPBandPrecData pdata,
                            realtype t, N_Vector y, N_Vector fy,
                            N_Vector ftemp, N_Vector ytemp);

static int cpBandPDQJacImpl(CPBandPrecData pdata,
                            realtype t, realtype gamma,
                            N_Vector y, N_Vector yp, N_Vector r,
                            N_Vector ftemp, N_Vector ytemp, N_Vector yptemp);

/* Redability replacements */

#define ode_type (cp_mem->cp_ode_type)
#define vec_tmpl (cp_mem->cp_tempv)

/*
 * =================================================================
 * EXPORTED FUNCTIONS
 * =================================================================
 */

void *CPBandPrecAlloc(void *cpode_mem, int N, int mu, int ml)
{
  CPodeMem cp_mem;
  CPBandPrecData pdata;
  int mup, mlp, storagemu;

  if (cpode_mem == NULL) {
    cpProcessError(NULL, 0, "CPBANDPRE", "CPBandPrecAlloc", MSGBP_CPMEM_NULL);
    return(NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  /* Test if the NVECTOR package is compatible with the BAND preconditioner */
  if(vec_tmpl->ops->nvgetarraypointer == NULL) {
    cpProcessError(cp_mem, 0, "CPBANDPRE", "CPBandPrecAlloc", MSGBP_BAD_NVECTOR);
    return(NULL);
  }

  pdata = NULL;
  pdata = (CPBandPrecData) malloc(sizeof *pdata);  /* Allocate data memory */
  if (pdata == NULL) {
    cpProcessError(cp_mem, 0, "CPBANDPRE", "CPBandPrecAlloc", MSGBP_MEM_FAIL);
    return(NULL);
  }

  /* Load pointers and bandwidths into pdata block. */
  pdata->cpode_mem = cpode_mem;
  pdata->N = N;
  pdata->mu = mup = MIN(N-1, MAX(0,mu));
  pdata->ml = mlp = MIN(N-1, MAX(0,ml));

  /* Initialize nfeBP counter */
  pdata->nfeBP = 0;

  /* Allocate memory for saved banded Jacobian approximation. */
  pdata->savedJ = NULL;
  pdata->savedJ = NewBandMat(N, mup, mlp, mup);
  if (pdata->savedJ == NULL) {
    free(pdata); pdata = NULL;
    cpProcessError(cp_mem, 0, "CPBANDPRE", "CPBandPrecAlloc", MSGBP_MEM_FAIL);
    return(NULL);
  }

  /* Allocate memory for banded preconditioner. */
  storagemu = MIN(N-1, mup+mlp);
  pdata->savedP = NULL;
  pdata->savedP = NewBandMat(N, mup, mlp, storagemu);
  if (pdata->savedP == NULL) {
    DestroyMat(pdata->savedJ);
    free(pdata); pdata = NULL;
    cpProcessError(cp_mem, 0, "CPBANDPRE", "CPBandPrecAlloc", MSGBP_MEM_FAIL);
    return(NULL);
  }

  /* Allocate memory for pivot array. */
  pdata->pivots = NULL;
  pdata->pivots = NewIntArray(N);
  if (pdata->savedJ == NULL) {
    DestroyMat(pdata->savedP);
    DestroyMat(pdata->savedJ);
    free(pdata); pdata = NULL;
    cpProcessError(cp_mem, 0, "CPBANDPRE", "CPBandPrecAlloc", MSGBP_MEM_FAIL);
    return(NULL);
  }

  return((void *) pdata);
}

int CPBPSptfqmr(void *cpode_mem, int pretype, int maxl, void *p_data)
{
  CPodeMem cp_mem;
  int flag;

  flag = CPSptfqmr(cpode_mem, pretype, maxl);
  if(flag != CPSPILS_SUCCESS) return(flag);

  cp_mem = (CPodeMem) cpode_mem;

  if ( p_data == NULL ) {
    cpProcessError(cp_mem, CPBANDPRE_PDATA_NULL, "CPBANDPRE", "CPBPSptfqmr", MSGBP_PDATA_NULL);
    return(CPBANDPRE_PDATA_NULL);
  }

  switch (ode_type) {
  case CP_EXPL:
    flag = CPSpilsSetPreconditioner(cpode_mem, (void *)cpBandPrecSetupExpl, (void *)cpBandPrecSolveExpl, p_data);
    break;
  case CP_IMPL:
    flag = CPSpilsSetPreconditioner(cpode_mem, (void *)cpBandPrecSetupImpl, (void *)cpBandPrecSolveImpl, p_data);
    break;
  }
  if(flag != CPSPILS_SUCCESS) return(flag);

  return(CPSPILS_SUCCESS);
}

int CPBPSpbcg(void *cpode_mem, int pretype, int maxl, void *p_data)
{
  CPodeMem cp_mem;
  int flag;

  flag = CPSpbcg(cpode_mem, pretype, maxl);
  if(flag != CPSPILS_SUCCESS) return(flag);

  cp_mem = (CPodeMem) cpode_mem;

  if ( p_data == NULL ) {
    cpProcessError(cp_mem, CPBANDPRE_PDATA_NULL, "CPBANDPRE", "CPBPSpbcg", MSGBP_PDATA_NULL);
    return(CPBANDPRE_PDATA_NULL);
  }

  switch (ode_type) {
  case CP_EXPL:
    flag = CPSpilsSetPreconditioner(cpode_mem, (void *)cpBandPrecSetupExpl, (void *)cpBandPrecSolveExpl, p_data);
    break;
  case CP_IMPL:
    flag = CPSpilsSetPreconditioner(cpode_mem, (void *)cpBandPrecSetupImpl, (void *)cpBandPrecSolveImpl, p_data);
    break;
  }
  if(flag != CPSPILS_SUCCESS) return(flag);

  return(CPSPILS_SUCCESS);
}

int CPBPSpgmr(void *cpode_mem, int pretype, int maxl, void *p_data)
{
  CPodeMem cp_mem;
  int flag;

  flag = CPSpgmr(cpode_mem, pretype, maxl);
  if(flag != CPSPILS_SUCCESS) return(flag);

  cp_mem = (CPodeMem) cpode_mem;

  if ( p_data == NULL ) {
    cpProcessError(cp_mem, CPBANDPRE_PDATA_NULL, "CPBANDPRE", "CPBPSpgmr", MSGBP_PDATA_NULL);
    return(CPBANDPRE_PDATA_NULL);
  }

  switch (ode_type) {
  case CP_EXPL:
    flag = CPSpilsSetPreconditioner(cpode_mem, (void *)cpBandPrecSetupExpl, (void *)cpBandPrecSolveExpl, p_data);
    break;
  case CP_IMPL:
    flag = CPSpilsSetPreconditioner(cpode_mem, (void *)cpBandPrecSetupImpl, (void *)cpBandPrecSolveImpl, p_data);
    break;
  }
  if(flag != CPSPILS_SUCCESS) return(flag);

  return(CPSPILS_SUCCESS);
}

void CPBandPrecFree(void **bp_data)
{
  CPBandPrecData pdata;

  if (*bp_data == NULL) return;

  pdata = (CPBandPrecData) (*bp_data);
  DestroyMat(pdata->savedJ);
  DestroyMat(pdata->savedP);
  DestroyArray(pdata->pivots);

  free(*bp_data);
  *bp_data = NULL;

}

int CPBandPrecGetWorkSpace(void *bp_data, long int *lenrwBP, long int *leniwBP)
{
  CPBandPrecData pdata;
  int N, ml, mu, smu;

  if ( bp_data == NULL ) {
    cpProcessError(NULL, CPBANDPRE_PDATA_NULL, "CPBANDPRE", "CPBandPrecGetWorkSpace", MSGBP_PDATA_NULL);
    return(CPBANDPRE_PDATA_NULL);
  }

  pdata = (CPBandPrecData) bp_data;

  N   = pdata->N;
  mu  = pdata->mu;
  ml  = pdata->ml;
  smu = MIN( N-1, mu + ml);

  *leniwBP = pdata->N;
  *lenrwBP = N * ( 2*ml + smu + mu + 2 );

  return(CPBANDPRE_SUCCESS);
}

int CPBandPrecGetNumFctEvals(void *bp_data, long int *nfevalsBP)
{
  CPBandPrecData pdata;

  if (bp_data == NULL) {
    cpProcessError(NULL, CPBANDPRE_PDATA_NULL, "CPBANDPRE", "CPBandPrecGetNumFctEvals", MSGBP_PDATA_NULL);
    return(CPBANDPRE_PDATA_NULL);
  }

  pdata = (CPBandPrecData) bp_data;

  *nfevalsBP = pdata->nfeBP;

  return(CPBANDPRE_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CPBandPrecGetReturnFlagName
 * -----------------------------------------------------------------
 */

char *CPBandPrecGetReturnFlagName(int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case CPBANDPRE_SUCCESS:
    sprintf(name,"CPBANDPRE_SUCCESS");
    break;
  case CPBANDPRE_PDATA_NULL:
    sprintf(name,"CPBANDPRE_PDATA_NULL");
    break;
  case CPBANDPRE_FUNC_UNRECVR:
    sprintf(name,"CPBANDPRE_FUNC_UNRECVR");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}

/*
 * =================================================================
 *  PRIVATE FUNCTIONS
 * =================================================================
 */

/* Readability Replacements */

#define N      (pdata->N)
#define mu     (pdata->mu)
#define ml     (pdata->ml)
#define pivots (pdata->pivots)
#define savedJ (pdata->savedJ)
#define savedP (pdata->savedP)
#define nfeBP  (pdata->nfeBP)

/*
 * -----------------------------------------------------------------
 * cpBandPrecSetupExpl
 * -----------------------------------------------------------------
 * Together cpBandPrecSetupExpl and cpBandPrecSolveExpl use a banded
 * difference quotient Jacobian to create a preconditioner.
 * cpBandPrecSetupExpl calculates a new J, if necessary, then
 * calculates P = I - gamma*J, and does an LU factorization of P.
 *
 * The parameters of cpBandPrecSetupExpl are as follows:
 *
 * t       is the current value of the independent variable.
 *
 * y       is the current value of the dependent variable vector,
 *         namely the predicted value of y(t).
 *
 * fy      is the vector f(t,y).
 *
 * jok     is an input flag indicating whether Jacobian-related
 *         data needs to be recomputed, as follows:
 *           jok == FALSE means recompute Jacobian-related data
 *                  from scratch.
 *           jok == TRUE means that Jacobian data from the
 *                  previous cpBandPrecSetupExpl call will be reused
 *                  (with the current value of gamma).
 *         A cpBandPrecSetupExpl call with jok == TRUE should only
 *         occur after a call with jok == FALSE.
 *
 * *jcurPtr is a pointer to an output integer flag which is
 *          set by cpBandPrecsetupExpl as follows:
 *            *jcurPtr = TRUE if Jacobian data was recomputed.
 *            *jcurPtr = FALSE if Jacobian data was not recomputed,
 *                       but saved data was reused.
 *
 * gamma   is the scalar appearing in the Newton matrix.
 *
 * bp_data is a pointer to preconditoner data - the same as the
 *         bp_data parameter passed to CPSp*.
 *
 * tmp1, tmp2, and tmp3 are pointers to memory allocated
 *           for vectors of length N for work space. This
 *           routine uses only tmp1 and tmp2.
 *
 * The value to be returned by the cpBandPrecSetupExpl function is
 *   0  if successful, or
 *   1  if the band factorization failed.
 * -----------------------------------------------------------------
 */

static int cpBandPrecSetupExpl(realtype t, N_Vector y, N_Vector fy,
                               booleantype jok, booleantype *jcurPtr,
                               realtype gamma, void *bp_data,
                               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  CPBandPrecData pdata;
  CPodeMem cp_mem;
  int ier, retval;

  /* Assume matrix and pivots have already been allocated. */
  pdata = (CPBandPrecData) bp_data;

  cp_mem = (CPodeMem) pdata->cpode_mem;

  if (jok) {

    /* If jok = TRUE, use saved copy of J. */
    *jcurPtr = FALSE;
    BandCopy(savedJ, savedP, mu, ml);

  } else {

    /* If jok = FALSE, call cpBandPDQJac for new J value. */
    *jcurPtr = TRUE;
    BandZero(savedJ);

    retval = cpBandPDQJacExpl(pdata, t, y, fy, tmp1, tmp2);
    if (retval < 0) {
      cpProcessError(cp_mem, CPBANDPRE_FUNC_UNRECVR, "CPBANDPRE", "cpBandPrecSetupExpl", MSGBP_FUNC_FAILED);
      return(-1);
    }
    if (retval > 0) {
      return(1);
    }

    BandCopy(savedJ, savedP, mu, ml);

  }

  /* Scale and add I to get savedP = I - gamma*J. */
  BandScale(-gamma, savedP);
  BandAddI(savedP);

  /* Do LU factorization of matrix. */
  ier = BandGBTRF(savedP, pivots);

  /* Return 0 if the LU was complete; otherwise return 1. */
  if (ier > 0) return(1);
  return(0);
}

/*
 * -----------------------------------------------------------------
 * cpBandPrecSolveExpl
 * -----------------------------------------------------------------
 * cpBandPrecSolveExpl solves a linear system P x = b, where P is the
 * matrix computed by cpBandPrecSetupExpl.
 *
 * The parameters of cpBandPrecSolveExpl used here are as follows:
 *
 * b       is the right-hand side vector of the linear system.
 *
 * bp_data is a pointer to preconditioner data - the same as the
 *         bp_data parameter passed to CPSp*.
 *
 * x       is the output vector computed by cpBandPrecSolveExpl.
 *
 * The value returned by the cpBandPrecSolveExpl function is always 0,
 * indicating success.
 * -----------------------------------------------------------------
 */

static int cpBandPrecSolveExpl(realtype t, N_Vector y, N_Vector fy,
                               N_Vector b, N_Vector x,
                               realtype gamma, realtype delta,
                               int lr, void *bp_data, N_Vector tmp)
{
  CPBandPrecData pdata;
  realtype *xd;

  /* Assume matrix and pivots have already been allocated. */
  pdata = (CPBandPrecData) bp_data;

  /* Copy b to x. */
  N_VScale(ONE, b, x);

  /* Do band backsolve on the vector z. */
  xd = N_VGetArrayPointer(x);

  BandGBTRS(savedP, pivots, xd);

  return(0);
}

/*
 * -----------------------------------------------------------------
 * cpBandPrecSetupImpl
 * -----------------------------------------------------------------
 * Together cpBandPrecSetupImpl and cpBandPrecSolveImpl use a banded
 * difference quotient Jacobian to create a preconditioner.
 * cpBandPrecSetupImpl calculates a new J = dF/dy + gamma*dF/dy'
 * and does an LU factorization of P.
 * -----------------------------------------------------------------
 */

static int cpBandPrecSetupImpl(realtype t, N_Vector y, N_Vector yp, N_Vector r,
                               realtype gamma, void *bp_data,
                               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  CPBandPrecData pdata;
  CPodeMem cp_mem;
  int ier, retval;

  /* Assume matrix and pivots have already been allocated. */
  pdata = (CPBandPrecData) bp_data;

  cp_mem = (CPodeMem) pdata->cpode_mem;

  BandZero(savedJ);
  retval = cpBandPDQJacImpl(pdata, t, gamma, y, yp, r, tmp1, tmp2, tmp3);
  if (retval < 0) {
    cpProcessError(cp_mem, CPBANDPRE_FUNC_UNRECVR, "CPBANDPRE", "cpBandPrecSetupImpl", MSGBP_FUNC_FAILED);
    return(-1);
  }
  if (retval > 0) {
    return(1);
  }

  /* Do LU factorization of matrix. */
  ier = BandGBTRF(savedP, pivots);

  /* Return 0 if the LU was complete; otherwise return 1. */
  if (ier > 0) return(1);
  return(0);
}

/*
 * -----------------------------------------------------------------
 * cpBandPrecSolveImpl
 * -----------------------------------------------------------------
 * cpBandPrecSolveImpl solves a linear system P x = b, where P is the
 * matrix computed by cpBandPrecSetupImpl.
 *
 * The parameters of cpBandPrecSolveImpl used here are as follows:
 *
 * b       is the right-hand side vector of the linear system.
 *
 * bp_data is a pointer to preconditioner data - the same as the
 *         bp_data parameter passed to CPSp*.
 *
 * x       is the output vector computed by cpBandPrecSolveImpl.
 *
 * The value returned by the cpBandPrecSolveImpl function is always 0,
 * indicating success.
 * -----------------------------------------------------------------
 */

static int cpBandPrecSolveImpl(realtype t, N_Vector y, N_Vector yp, N_Vector r,
                               N_Vector b, N_Vector x,
                               realtype gamma, realtype delta,
                               void *bp_data, N_Vector tmp)
{
  CPBandPrecData pdata;
  realtype *xd;

  /* Assume matrix and pivots have already been allocated. */
  pdata = (CPBandPrecData) bp_data;

  /* Copy b to x. */
  N_VScale(ONE, b, x);

  /* Do band backsolve on the vector z. */
  xd = N_VGetArrayPointer(x);

  BandGBTRS(savedP, pivots, xd);

  return(0);
}

/*
 * =================================================================
 * DQ LOCAL JACOBIAN APROXIMATIONS
 * =================================================================
 */

#define ewt    (cp_mem->cp_ewt)
#define uround (cp_mem->cp_uround)
#define h      (cp_mem->cp_h)
#define fe     (cp_mem->cp_fe)
#define fi     (cp_mem->cp_fi)
#define f_data (cp_mem->cp_f_data)

/*
 * -----------------------------------------------------------------
 * cpBandPDQJacExpl
 * -----------------------------------------------------------------
 * This routine generates a banded difference quotient approximation to
 * the Jacobian of f(t,y). It assumes that a band matrix of type
 * BandMat is stored column-wise, and that elements within each column
 * are contiguous. This makes it possible to get the address of a column
 * of J via the macro BAND_COL and to write a simple for loop to set
 * each of the elements of a column in succession.
 * -----------------------------------------------------------------
 */

static int cpBandPDQJacExpl(CPBandPrecData pdata,
                            realtype t, N_Vector y, N_Vector fy,
                            N_Vector ftemp, N_Vector ytemp)
{
  CPodeMem cp_mem;
  realtype fnorm, minInc, inc, inc_inv, srur;
  int group, i, j, width, ngroups, i1, i2;
  realtype *col_j, *ewt_data, *fy_data, *ftemp_data, *y_data, *ytemp_data;
  int retval;

  cp_mem = (CPodeMem) pdata->cpode_mem;

  /* Obtain pointers to the data for ewt, fy, ftemp, y, ytemp. */
  ewt_data   = N_VGetArrayPointer(ewt);
  fy_data    = N_VGetArrayPointer(fy);
  ftemp_data = N_VGetArrayPointer(ftemp);
  y_data     = N_VGetArrayPointer(y);
  ytemp_data = N_VGetArrayPointer(ytemp);

  /* Load ytemp with y = predicted y vector. */
  N_VScale(ONE, y, ytemp);

  /* Set minimum increment based on uround and norm of f. */
  srur = RSqrt(uround);
  fnorm = N_VWrmsNorm(fy, ewt);
  minInc = (fnorm != ZERO) ?
           (MIN_INC_MULT * ABS(h) * uround * N * fnorm) : ONE;

  /* Set bandwidth and number of column groups for band differencing. */
  width = ml + mu + 1;
  ngroups = MIN(width, N);

  /* Loop over column groups. */
  for (group = 1; group <= ngroups; group++) {

    /* Increment all y_j in group. */
    for(j = group-1; j < N; j += width) {
      inc = MAX(srur*ABS(y_data[j]), minInc/ewt_data[j]);
      ytemp_data[j] += inc;
    }

    /* Evaluate f with incremented y. */

    retval = fe(t, ytemp, ftemp, f_data);
    nfeBP++;
    if (retval != 0) return(retval);

    /* Restore ytemp, then form and load difference quotients. */
    for (j = group-1; j < N; j += width) {
      ytemp_data[j] = y_data[j];
      col_j = BAND_COL(savedJ,j);
      inc = MAX(srur*ABS(y_data[j]), minInc/ewt_data[j]);
      inc_inv = ONE/inc;
      i1 = MAX(0, j-mu);
      i2 = MIN(j+ml, N-1);
      for (i=i1; i <= i2; i++)
        BAND_COL_ELEM(col_j,i,j) = inc_inv * (ftemp_data[i] - fy_data[i]);
    }
  }

  return(0);
}

/*
 * -----------------------------------------------------------------
 * cpBandPDQJacImpl
 * -----------------------------------------------------------------
 * This routine generates a banded difference quotient approximation
 * to the Jacobian dF/dy + gamma*dF/dy' and loads it directly into
 * savedP. It assumes that a band matrix of type BandMat is stored
 * column-wise, and that elements within each column are contiguous.
 * This makes it possible to get the address of a column of J via
 * the macro BAND_COL and to write a simple for loop to set each of
 * the elements of a column in succession.
 * -----------------------------------------------------------------
 */

static int cpBandPDQJacImpl(CPBandPrecData pdata,
                            realtype t, realtype gamma,
                            N_Vector y, N_Vector yp, N_Vector r,
                            N_Vector ftemp, N_Vector ytemp, N_Vector yptemp)

{
  CPodeMem cp_mem;
  realtype inc, inc_inv, yj, ypj, srur, ewtj;
  realtype *y_data, *yp_data, *ewt_data;
  realtype *ytemp_data, *yptemp_data, *ftemp_data, *r_data, *col_j;
  int i, j, i1, i2, width, group, ngroups;
  int retval = 0;

  cp_mem = (CPodeMem) pdata->cpode_mem;

  /* Obtain pointers to the data for all vectors used.  */
  ewt_data    = N_VGetArrayPointer(ewt);
  r_data      = N_VGetArrayPointer(r);
  y_data      = N_VGetArrayPointer(y);
  yp_data     = N_VGetArrayPointer(yp);
  ftemp_data  = N_VGetArrayPointer(ftemp);
  ytemp_data  = N_VGetArrayPointer(ytemp);
  yptemp_data = N_VGetArrayPointer(yptemp);

  /* Initialize ytemp and yptemp. */
  N_VScale(ONE, y, ytemp);
  N_VScale(ONE, yp, yptemp);

  /* Compute miscellaneous values for the Jacobian computation. */
  srur = RSqrt(uround);
  width = ml + mu + 1;
  ngroups = MIN(width, N);

  /* Loop over column groups. */
  for (group=1; group <= ngroups; group++) {

    /* Increment all y[j] and yp[j] for j in this group. */
    for (j=group-1; j<N; j+=width) {
        yj = y_data[j];
        ypj = yp_data[j];
        ewtj = ewt_data[j];

        /* Set increment inc to yj based on sqrt(uround)*abs(yj), with
           adjustments using ypj and ewtj if this is small, and a further
           adjustment to give it the same sign as h*ypj. */
        inc = MAX( srur * MAX( ABS(yj), ABS(h*ypj) ) , ONE/ewtj );

        if (h*ypj < ZERO) inc = -inc;
        inc = (yj + inc) - yj;

        /* Increment yj and ypj. */
        ytemp_data[j]  += gamma*inc;
        yptemp_data[j] += inc;
    }

    /* Call ODE fct. with incremented arguments. */
    retval = fi(t, ytemp, yptemp, ftemp, f_data);
    nfeBP++;
    if (retval != 0) break;

    /* Loop over the indices j in this group again. */
    for (j=group-1; j<N; j+=width) {

      /* Reset ytemp and yptemp components that were perturbed. */
      yj = ytemp_data[j]  = y_data[j];
      ypj = yptemp_data[j] = yp_data[j];
      col_j = BAND_COL(savedP, j);
      ewtj = ewt_data[j];

      /* Set increment inc exactly as above. */
      inc = MAX( srur * MAX( ABS(yj), ABS(h*ypj) ) , ONE/ewtj );
      if (h*ypj < ZERO) inc = -inc;
      inc = (yj + inc) - yj;

      /* Load the difference quotient Jacobian elements for column j. */
      inc_inv = ONE/inc;
      i1 = MAX(0, j-mu);
      i2 = MIN(j+ml,N-1);
      for (i=i1; i<=i2; i++)
        BAND_COL_ELEM(col_j,i,j) = inc_inv*(ftemp_data[i]-r_data[i]);
    }

  }

  return(retval);
}
