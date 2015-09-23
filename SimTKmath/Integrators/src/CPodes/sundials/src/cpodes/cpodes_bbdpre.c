/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2006/11/24 19:09:18 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This file contains implementations of routines for a
 * band-block-diagonal preconditioner, i.e. a block-diagonal
 * matrix with banded blocks, for use with CPODES, a CPSPILS linear
 * solver, and the parallel implementation of NVECTOR.
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

#include "cpodes_bbdpre_impl.h"
#include "cpodes_private.h"

#include <sundials/sundials_math.h>

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

/* cpBBDPrecSetupExpl and cpBBDPrecSetupImpl */

static int cpBBDPrecSetupExpl(realtype t, N_Vector y, N_Vector fy,
                              booleantype jok, booleantype *jcurPtr,
                              realtype gamma, void *bbd_data,
                              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int cpBBDPrecSetupImpl(realtype t, N_Vector y, N_Vector yp, N_Vector r,
                              realtype gamma, void *bbd_data,
                              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* cpBBDPrecSolveExpl and cpBBDPrecSolveImpl */

static int cpBBDPrecSolveExpl(realtype t, N_Vector y, N_Vector fy,
                              N_Vector b, N_Vector x,
                              realtype gamma, realtype delta,
                              int lr, void *bbd_data, N_Vector tmp);

static int cpBBDPrecSolveImpl(realtype t, N_Vector y, N_Vector yp, N_Vector r,
                              N_Vector b, N_Vector x,
                              realtype gamma, realtype delta,
                              void *bbd_data, N_Vector tmp);

/* Difference quotient Jacobian calculation routines */

static int cpBBDDQJacExpl(CPBBDPrecData pdata, realtype t,
                          N_Vector y, N_Vector gy,
                          N_Vector ytemp, N_Vector gtemp);

static int cpBBDDQJacImpl(CPBBDPrecData pdata, realtype tt, realtype gamma,
                          N_Vector yy, N_Vector yp, N_Vector gref,
                          N_Vector ytemp, N_Vector yptemp, N_Vector gtemp);

/* Redability replacements */

#define ode_type (cp_mem->cp_ode_type)
#define uround   (cp_mem->cp_uround)
#define vec_tmpl (cp_mem->cp_tempv)

/*
 * =================================================================
 * EXPORTED FUNCTIONS
 * =================================================================
 */

void *CPBBDPrecAlloc(void *cpode_mem, int Nlocal,
                     int mudq, int mldq, int mukeep, int mlkeep,
                     realtype dqrely,
                     void *gloc, CPBBDCommFn cfn)
{
  CPodeMem cp_mem;
  CPBBDPrecData pdata;
  N_Vector tmp4;
  int muk, mlk, storage_mu;

  if (cpode_mem == NULL) {
    cpProcessError(NULL, 0, "CPBBDPRE", "CPBBDPrecAlloc", MSGBBDP_CPMEM_NULL);
    return(NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  /* Test if the NVECTOR package is compatible with the BLOCK BAND preconditioner */
  if(vec_tmpl->ops->nvgetarraypointer == NULL) {
    cpProcessError(cp_mem, 0, "CPBBDPRE", "CPBBDPrecAlloc", MSGBBDP_BAD_NVECTOR);
    return(NULL);
  }

  /* Allocate data memory */
  pdata = NULL;
  pdata = (CPBBDPrecData) malloc(sizeof *pdata);
  if (pdata == NULL) {
    cpProcessError(cp_mem, 0, "CPBBDPRE", "CPBBDPrecAlloc", MSGBBDP_MEM_FAIL);
    return(NULL);
  }

  /* Set pointers to gloc and cfn; load half-bandwidths */

  pdata->cpode_mem = cpode_mem;

  switch (ode_type) {
  case CP_EXPL:
    pdata->glocE = (CPBBDLocalRhsFn) gloc;
    pdata->glocI = NULL;
    break;
  case CP_IMPL:
    pdata->glocI = (CPBBDLocalResFn) gloc;
    pdata->glocE = NULL;
    break;
  }

  pdata->cfn = cfn;

  pdata->mudq = MIN(Nlocal-1, MAX(0,mudq));
  pdata->mldq = MIN(Nlocal-1, MAX(0,mldq));

  muk = MIN(Nlocal-1, MAX(0,mukeep));
  mlk = MIN(Nlocal-1, MAX(0,mlkeep));
  pdata->mukeep = muk;
  pdata->mlkeep = mlk;

  /* Allocate memory for saved Jacobian */
  pdata->savedJ = NewBandMat(Nlocal, muk, mlk, muk);
  if (pdata->savedJ == NULL) {
    free(pdata); pdata = NULL;
    cpProcessError(cp_mem, 0, "CPBBDPRE", "CPBBDPrecAlloc", MSGBBDP_MEM_FAIL);
    return(NULL);
  }

  /* Allocate memory for preconditioner matrix */
  storage_mu = MIN(Nlocal-1, muk + mlk);
  pdata->savedP = NULL;
  pdata->savedP = NewBandMat(Nlocal, muk, mlk, storage_mu);
  if (pdata->savedP == NULL) {
    DestroyMat(pdata->savedJ);
    free(pdata); pdata = NULL;
    cpProcessError(cp_mem, 0, "CPBBDPRE", "CPBBDPrecAlloc", MSGBBDP_MEM_FAIL);
    return(NULL);
  }
  /* Allocate memory for pivots */
  pdata->pivots = NULL;
  pdata->pivots = NewIntArray(Nlocal);
  if (pdata->savedJ == NULL) {
    DestroyMat(pdata->savedP);
    DestroyMat(pdata->savedJ);
    free(pdata); pdata = NULL;
    cpProcessError(cp_mem, 0, "CPBBDPRE", "CPBBDPrecAlloc", MSGBBDP_MEM_FAIL);
    return(NULL);
  }
  /* Allocate tmp4 for use by cpBBDDQJacImpl */
  tmp4 = NULL;
  tmp4 = N_VClone(vec_tmpl);
  if (tmp4 == NULL){
    DestroyMat(pdata->savedP);
    DestroyMat(pdata->savedJ);
    DestroyArray(pdata->pivots);
    free(pdata); pdata = NULL;
    cpProcessError(cp_mem, 0, "CPBBDPRE", "CPBBDPrecAlloc", MSGBBDP_MEM_FAIL);
    return(NULL);
  }
  pdata->tmp4 = tmp4;

  /* Set pdata->dqrely based on input dqrely (0 implies default). */
  pdata->dqrely = (dqrely > ZERO) ? dqrely : RSqrt(uround);

  /* Store Nlocal to be used in cpBBDPrecSetupExpl and cpBBDPrecSetupImpl */
  pdata->n_local = Nlocal;

  /* Set work space sizes and initialize nge */
  pdata->rpwsize = Nlocal*(muk + 2*mlk + storage_mu + 2);
  pdata->ipwsize = Nlocal;
  pdata->nge = 0;

  return((void *)pdata);
}

int CPBBDSptfqmr(void *cpode_mem, int pretype, int maxl, void *bbd_data)
{
  CPodeMem cp_mem;
  int flag;

  flag = CPSptfqmr(cpode_mem, pretype, maxl);
  if(flag != CPSPILS_SUCCESS) return(flag);

  cp_mem = (CPodeMem) cpode_mem;

  if (bbd_data == NULL) {
    cpProcessError(cp_mem, CPBBDPRE_PDATA_NULL, "CPBBDPRE", "CPBBDSptfqmr", MSGBBDP_PDATA_NULL);
    return(CPBBDPRE_PDATA_NULL);
  }

  switch (ode_type) {
  case CP_EXPL:
    flag = CPSpilsSetPreconditioner(cpode_mem, (void *)cpBBDPrecSetupExpl, (void *)cpBBDPrecSolveExpl, bbd_data);
    break;
  case CP_IMPL:
    flag = CPSpilsSetPreconditioner(cpode_mem, (void *)cpBBDPrecSetupImpl, (void *)cpBBDPrecSolveImpl, bbd_data);
    break;
  }
  if(flag != CPSPILS_SUCCESS) return(flag);

  return(CPSPILS_SUCCESS);
}

int CPBBDSpbcg(void *cpode_mem, int pretype, int maxl, void *bbd_data)
{
  CPodeMem cp_mem;
  int flag;

  flag = CPSpbcg(cpode_mem, pretype, maxl);
  if(flag != CPSPILS_SUCCESS) return(flag);

  cp_mem = (CPodeMem) cpode_mem;

  if (bbd_data == NULL) {
    cpProcessError(cp_mem, CPBBDPRE_PDATA_NULL, "CPBBDPRE", "CPBBDSpbcg", MSGBBDP_PDATA_NULL);
    return(CPBBDPRE_PDATA_NULL);
  }

  switch (ode_type) {
  case CP_EXPL:
    flag = CPSpilsSetPreconditioner(cpode_mem, (void *)cpBBDPrecSetupExpl, (void *)cpBBDPrecSolveExpl, bbd_data);
    break;
  case CP_IMPL:
    flag = CPSpilsSetPreconditioner(cpode_mem, (void *)cpBBDPrecSetupImpl, (void *)cpBBDPrecSolveImpl, bbd_data);
    break;
  }
  if(flag != CPSPILS_SUCCESS) return(flag);

  return(CPSPILS_SUCCESS);
}

int CPBBDSpgmr(void *cpode_mem, int pretype, int maxl, void *bbd_data)
{
  CPodeMem cp_mem;
  int flag;

  flag = CPSpgmr(cpode_mem, pretype, maxl);
  if(flag != CPSPILS_SUCCESS) return(flag);

  cp_mem = (CPodeMem) cpode_mem;

  if (bbd_data == NULL) {
    cpProcessError(cp_mem, CPBBDPRE_PDATA_NULL, "CPBBDPRE", "CPBBDSpgmr", MSGBBDP_PDATA_NULL);
    return(CPBBDPRE_PDATA_NULL);
  }

  switch (ode_type) {
  case CP_EXPL:
    flag = CPSpilsSetPreconditioner(cpode_mem, (void *)cpBBDPrecSetupExpl, (void *)cpBBDPrecSolveExpl, bbd_data);
    break;
  case CP_IMPL:
    flag = CPSpilsSetPreconditioner(cpode_mem, (void *)cpBBDPrecSetupImpl, (void *)cpBBDPrecSolveImpl, bbd_data);
    break;
  }
  if(flag != CPSPILS_SUCCESS) return(flag);

  return(CPSPILS_SUCCESS);
}

int CPBBDPrecReInit(void *bbd_data, int mudq, int mldq,
                    realtype dqrely,
                    void *gloc, CPBBDCommFn cfn)
{
  CPBBDPrecData pdata;
  CPodeMem cp_mem;
  int Nlocal;

  if (bbd_data == NULL) {
    cpProcessError(NULL, CPBBDPRE_PDATA_NULL, "CPBBDPRE", "CPBBDPrecReInit", MSGBBDP_PDATA_NULL);
    return(CPBBDPRE_PDATA_NULL);
  }
  pdata  = (CPBBDPrecData) bbd_data;
  cp_mem = (CPodeMem) pdata->cpode_mem;

  /* Set pointers to gloc and cfn; load half-bandwidths */
  switch (ode_type) {
  case CP_EXPL:
    pdata->glocE = (CPBBDLocalRhsFn) gloc;
    pdata->glocI = NULL;
    break;
  case CP_IMPL:
    pdata->glocI = (CPBBDLocalResFn) gloc;
    pdata->glocE = NULL;
    break;
  }
  pdata->cfn = cfn;
  Nlocal = pdata->n_local;
  pdata->mudq = MIN(Nlocal-1, MAX(0,mudq));
  pdata->mldq = MIN(Nlocal-1, MAX(0,mldq));

  /* Set pdata->dqrely based on input dqrely (0 implies default). */
  pdata->dqrely = (dqrely > ZERO) ? dqrely : RSqrt(uround);

  /* Re-initialize nge */
  pdata->nge = 0;

  return(CPBBDPRE_SUCCESS);
}

void CPBBDPrecFree(void **bbd_data)
{
  CPBBDPrecData pdata;

  if (*bbd_data == NULL) return;

  pdata = (CPBBDPrecData) (*bbd_data);
  DestroyMat(pdata->savedJ);
  DestroyMat(pdata->savedP);
  DestroyArray(pdata->pivots);
  N_VDestroy(pdata->tmp4);

  free(*bbd_data);
  *bbd_data = NULL;

}

int CPBBDPrecGetWorkSpace(void *bbd_data, long int *lenrwBBDP, long int *leniwBBDP)
{
  CPBBDPrecData pdata;

  if (bbd_data == NULL) {
    cpProcessError(NULL, CPBBDPRE_PDATA_NULL, "CPBBDPRE", "CPBBDPrecGetWorkSpace", MSGBBDP_PDATA_NULL);
    return(CPBBDPRE_PDATA_NULL);
  }

  pdata = (CPBBDPrecData) bbd_data;

  *lenrwBBDP = pdata->rpwsize;
  *leniwBBDP = pdata->ipwsize;

  return(CPBBDPRE_SUCCESS);
}

int CPBBDPrecGetNumGfnEvals(void *bbd_data, long int *ngevalsBBDP)
{
  CPBBDPrecData pdata;

  if (bbd_data == NULL) {
    cpProcessError(NULL, CPBBDPRE_PDATA_NULL, "CPBBDPRE", "CPBBDPrecGetNumGfnEvals", MSGBBDP_PDATA_NULL);
    return(CPBBDPRE_PDATA_NULL);
  }

  pdata = (CPBBDPrecData) bbd_data;

  *ngevalsBBDP = pdata->nge;

  return(CPBBDPRE_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CPBBDPrecGetReturnFlagName
 * -----------------------------------------------------------------
 */

char *CPBBDPrecGetReturnFlagName(int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case CPBBDPRE_SUCCESS:
    sprintf(name,"CPBBDPRE_SUCCESS");
    break;
  case CPBBDPRE_PDATA_NULL:
    sprintf(name,"CPBBDPRE_PDATA_NULL");
    break;
  case CPBBDPRE_FUNC_UNRECVR:
    sprintf(name,"CPBBDPRE_FUNC_UNRECVR");
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

#define Nlocal (pdata->n_local)
#define mudq   (pdata->mudq)
#define mldq   (pdata->mldq)
#define mukeep (pdata->mukeep)
#define mlkeep (pdata->mlkeep)
#define dqrely (pdata->dqrely)
#define glocE  (pdata->glocE)
#define glocI  (pdata->glocI)
#define cfn    (pdata->cfn)
#define savedJ (pdata->savedJ)
#define savedP (pdata->savedP)
#define pivots (pdata->pivots)
#define nge    (pdata->nge)

/*
 * -----------------------------------------------------------------
 * Function: cpBBDPrecSetupExpl
 * -----------------------------------------------------------------
 * cpBBDPrecSetupExpl generates and factors a banded block of the
 * preconditioner matrix on each processor, via calls to the
 * user-supplied glocE and cfn functions. It uses difference
 * quotient approximations to the Jacobian elements.
 *
 * cpBBDPrecSetupExpl calculates a new J,if necessary, then
 * calculates P = I - gamma*J, and does an LU factorization of P.
 *
 * The parameters of cpBBDPrecSetupExpl used here are as follows:
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
 *           jok == TRUE  means that Jacobian data from the
 *                  previous cpBBDPrecSetupExpl call can be reused
 *                  (with the current value of gamma).
 *         A cpBBDPrecSetupExpl call with jok == TRUE should only
 *         occur after a call with jok == FALSE.
 *
 * jcurPtr is a pointer to an output integer flag which is
 *         set by cpBBDPrecSetupExpl as follows:
 *           *jcurPtr = TRUE if Jacobian data was recomputed.
 *           *jcurPtr = FALSE if Jacobian data was not recomputed,
 *                      but saved data was reused.
 *
 * gamma   is the scalar appearing in the Newton matrix.
 *
 * bbd_data  is a pointer to user data returned by CPBBDPrecAlloc.
 *
 * tmp1, tmp2, and tmp3 are pointers to memory allocated for
 *           NVectors which are be used by cpBBDPrecSetupExpl
 *           as temporary storage or work space.
 *
 * Return value:
 * The value returned by this cpBBDPrecSetupExpl function is the int
 *   0  if successful,
 *   1  for a recoverable error (step will be retried).
 * -----------------------------------------------------------------
 */

static int cpBBDPrecSetupExpl(realtype t, N_Vector y, N_Vector fy,
                              booleantype jok, booleantype *jcurPtr,
                              realtype gamma, void *bbd_data,
                              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int ier;
  CPBBDPrecData pdata;
  CPodeMem cp_mem;
  int retval;

  pdata = (CPBBDPrecData) bbd_data;

  cp_mem = (CPodeMem) pdata->cpode_mem;

  if (jok) {

    /* If jok = TRUE, use saved copy of J */
    *jcurPtr = FALSE;
    BandCopy(savedJ, savedP, mukeep, mlkeep);

  } else {

    /* Otherwise call cpBBDDQJacExpl for new J value */
    *jcurPtr = TRUE;
    BandZero(savedJ);

    retval = cpBBDDQJacExpl(pdata, t, y, tmp1, tmp2, tmp3);
    if (retval < 0) {
      cpProcessError(cp_mem, CPBBDPRE_FUNC_UNRECVR, "CPBBDPRE", "cpBBDPrecSetup", MSGBBDP_FUNC_FAILED);
      return(-1);
    }
    if (retval > 0) {
      return(1);
    }

    BandCopy(savedJ, savedP, mukeep, mlkeep);

  }

  /* Scale and add I to get P = I - gamma*J */
  BandScale(-gamma, savedP);
  BandAddI(savedP);

  /* Do LU factorization of P in place */
  ier = BandGBTRF(savedP, pivots);

  /* Return 0 if the LU was complete; otherwise return 1 */
  if (ier > 0) return(1);
  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function: cpBBDPrecSolveExpl
 * -----------------------------------------------------------------
 * cpBBDPrecSolveExpl solves a linear system P z = r, with the
 * band-block-diagonal preconditioner matrix P generated and
 * factored by cpBBDPrecSetupExpl.
 *
 * The parameters of cpBBDPrecSolveExpl used here are as follows:
 *
 *   r - right-hand side vector of the linear system.
 *
 *   bbd_data - pointer to the preconditioner data returned by
 *              CPBBDPrecAlloc.
 *
 *   z - output vector computed by cpBBDPrecSolveExpl.
 *
 * The value returned by the cpBBDPrecSolveExpl function is always
 * 0, indicating success.
 * -----------------------------------------------------------------
 */

static int cpBBDPrecSolveExpl(realtype t, N_Vector y, N_Vector fy,
                              N_Vector b, N_Vector x,
                              realtype gamma, realtype delta,
                              int lr, void *bbd_data, N_Vector tmp)
{
  CPBBDPrecData pdata;
  realtype *xd;

  pdata = (CPBBDPrecData) bbd_data;

  /* Copy b to x, then do backsolve and return */
  N_VScale(ONE, b, x);

  xd = N_VGetArrayPointer(x);

  BandGBTRS(savedP, pivots, xd);

  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : cpBBDPrecSetupImpl
 * -----------------------------------------------------------------
 * cpBBDPrecSetupImpl generates a band-block-diagonal preconditioner
 * matrix, where the local block (on this processor) is a band
 * matrix. Each local block is computed by a difference quotient
 * scheme via calls to the user-supplied routines glocal, gcomm.
 * After generating the block in the band matrix savedP, this routine
 * does an LU factorization in place in savedP.
 *
 * The cpBBDPrecSetupImpl parameters used here are as follows:
 *
 * t is the current value of the independent variable t.
 *
 * yy is the current value of the dependent variable vector,
 *    namely the predicted value of y(t).
 *
 * yp is the current value of the derivative vector y',
 *    namely the predicted value of y'(t).
 *
 * gamma is the scalar in the system Jacobian, proportional to h.
 *
 * bbd_data is a pointer to user preconditioner data returned by
 *    CPBBDPrecAlloc.
 *
 * tmp1, tmp2, tmp3 are pointers to vectors of type N_Vector,
 *    used for temporary storage or work space.
 *
 * Return value:
 * The value returned by cpBBDPrecSetupImpl function is a int
 * flag indicating whether it was successful. This value is
 *    0    if successful,
 *  > 0    for a recoverable error (step will be retried), or
 *  < 0    for a nonrecoverable error (step fails).
 * -----------------------------------------------------------------
 */

static int cpBBDPrecSetupImpl(realtype t, N_Vector y, N_Vector yp, N_Vector r,
                              realtype gamma, void *bbd_data,
                              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int ier, retval;
  CPBBDPrecData pdata;
  CPodeMem cp_mem;

  pdata =(CPBBDPrecData) bbd_data;

  cp_mem = (CPodeMem) pdata->cpode_mem;

  /* Call cpBBDDQJacImpl for a new Jacobian calculation and store in savedP. */
  BandZero(savedP);
  retval = cpBBDDQJacImpl(pdata, t, gamma, y, yp,
                          tmp1, tmp2, tmp3, pdata->tmp4);
  if (retval < 0) {
    cpProcessError(cp_mem, CPBBDPRE_FUNC_UNRECVR, "CPBBDPRE", "cpBBDPrecSetupImpl", MSGBBDP_FUNC_FAILED);
    return(-1);
  }
  if (retval > 0) {
    return(+1);
  }

  /* Do LU factorization of preconditioner block in place (in savedP). */
  ier = BandGBTRF(savedP, pivots);

  /* Return 0 if the LU was complete, or +1 otherwise. */
  if (ier > 0) return(+1);
  return(0);
}


/*
 * -----------------------------------------------------------------
 * Function: cpBBDPrecSolveImpl
 * -----------------------------------------------------------------
 * The function cpBBDPrecSolve computes a solution to the linear
 * system P x = b, where P is the left preconditioner defined by
 * the routine cpBBDPrecSetupImpl.
 *
 * The cpBBDPrecSolveImpl parameters used here are as follows:
 *
 * b is the input right-hand side vector r.
 *
 * x is the computed solution vector.
 *
 * bbd_data is a pointer to user preconditioner data returned
 *     by CPBBDPrecAlloc.
 *
 * The arguments t, y, yp, r, gamma, delta, and tmp are NOT used.
 *
 * cpBBDPrecSolveImpl always returns 0, indicating success.
 * -----------------------------------------------------------------
 */

static int cpBBDPrecSolveImpl(realtype t, N_Vector y, N_Vector yp, N_Vector r,
                              N_Vector b, N_Vector x,
                              realtype gamma, realtype delta,
                              void *bbd_data, N_Vector tmp)
{
  CPBBDPrecData pdata;
  realtype *xd;

  pdata = (CPBBDPrecData) bbd_data;

  /* Copy b to x, do the backsolve, and return. */
  N_VScale(ONE, b, x);

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
#define h      (cp_mem->cp_h)
#define f_data (cp_mem->cp_f_data)

/*
 * -----------------------------------------------------------------
 * Function: cpBBDDQJacExpl
 * -----------------------------------------------------------------
 * This routine generates a banded difference quotient approximation
 * to the local block of the Jacobian of g(t,y). It assumes that a
 * band matrix of type BandMat is stored columnwise, and that elements
 * within each column are contiguous. All matrix elements are generated
 * as difference quotients, by way of calls to the user routine glocE.
 * By virtue of the band structure, the number of these calls is
 * bandwidth + 1, where bandwidth = mldq + mudq + 1.
 * But the band matrix kept has bandwidth = mlkeep + mukeep + 1.
 * This routine also assumes that the local elements of a vector are
 * stored contiguously.
 * -----------------------------------------------------------------
 */

static int cpBBDDQJacExpl(CPBBDPrecData pdata, realtype t,
                          N_Vector y, N_Vector gy,
                          N_Vector ytemp, N_Vector gtemp)
{
  CPodeMem cp_mem;
  realtype gnorm, minInc, inc, inc_inv;
  int group, i, j, width, ngroups, i1, i2;
  realtype *y_data, *ewt_data, *gy_data, *gtemp_data, *ytemp_data, *col_j;
  int retval;

  cp_mem = (CPodeMem) pdata->cpode_mem;

  /* Load ytemp with y = predicted solution vector */
  N_VScale(ONE, y, ytemp);

  /* Call cfn and glocE to get base value of g(t,y) */
  if (cfn != NULL) {
    retval = cfn(Nlocal, t, y, NULL, f_data);
    if (retval != 0) return(retval);
  }

  retval = glocE(Nlocal, t, ytemp, gy, f_data);
  nge++;
  if (retval != 0) return(retval);

  /* Obtain pointers to the data for various vectors */
  y_data     =  N_VGetArrayPointer(y);
  gy_data    =  N_VGetArrayPointer(gy);
  ewt_data   =  N_VGetArrayPointer(ewt);
  ytemp_data =  N_VGetArrayPointer(ytemp);
  gtemp_data =  N_VGetArrayPointer(gtemp);

  /* Set minimum increment based on uround and norm of g */
  gnorm = N_VWrmsNorm(gy, ewt);
  minInc = (gnorm != ZERO) ?
           (MIN_INC_MULT * ABS(h) * uround * Nlocal * gnorm) : ONE;

  /* Set bandwidth and number of column groups for band differencing */
  width = mldq + mudq + 1;
  ngroups = MIN(width, Nlocal);

  /* Loop over groups */
  for (group=1; group <= ngroups; group++) {

    /* Increment all y_j in group */
    for(j=group-1; j < Nlocal; j+=width) {
      inc = MAX(dqrely*ABS(y_data[j]), minInc/ewt_data[j]);
      ytemp_data[j] += inc;
    }

    /* Evaluate g with incremented y */
    retval = glocE(Nlocal, t, ytemp, gtemp, f_data);
    nge++;
    if (retval != 0) return(retval);

    /* Restore ytemp, then form and load difference quotients */
    for (j=group-1; j < Nlocal; j+=width) {
      ytemp_data[j] = y_data[j];
      col_j = BAND_COL(savedJ,j);
      inc = MAX(dqrely*ABS(y_data[j]), minInc/ewt_data[j]);
      inc_inv = ONE/inc;
      i1 = MAX(0, j-mukeep);
      i2 = MIN(j+mlkeep, Nlocal-1);
      for (i=i1; i <= i2; i++)
        BAND_COL_ELEM(col_j,i,j) =
          inc_inv * (gtemp_data[i] - gy_data[i]);
    }
  }

  return(0);
}


/*
 * -----------------------------------------------------------------
 * cpBBDDQJacImpl
 * -----------------------------------------------------------------
 * This routine generates a banded difference quotient approximation
 * to the local block of the Jacobian of G(t,y,y'). It assumes that
 * a band matrix of type BandMat is stored column-wise, and that
 * elements within each column are contiguous.
 *
 * All matrix elements are generated as difference quotients, by way
 * of calls to the user routine glocI. By virtue of the band
 * structure, the number of these calls is bandwidth + 1, where
 * bandwidth = mldq + mudq + 1. But the band matrix kept has
 * bandwidth = mlkeep + mukeep + 1. This routine also assumes that
 * the local elements of a vector are stored contiguously.
 *
 * Return values are: 0 (success), > 0 (recoverable error),
 * or < 0 (nonrecoverable error).
 * -----------------------------------------------------------------
 */

static int cpBBDDQJacImpl(CPBBDPrecData pdata, realtype t, realtype gamma,
                          N_Vector y, N_Vector yp, N_Vector gref,
                          N_Vector ytemp, N_Vector yptemp, N_Vector gtemp)
{
  CPodeMem cp_mem;
  realtype inc, inc_inv;
  int  retval;
  int group, i, j, width, ngroups, i1, i2;
  realtype *ydata, *ypdata, *ytempdata, *yptempdata, *grefdata, *gtempdata;
  realtype *ewtdata;
  realtype *col_j, yj, ypj, ewtj;

  cp_mem = (CPodeMem) pdata->cpode_mem;

  /* Initialize ytemp and yptemp. */
  N_VScale(ONE, y, ytemp);
  N_VScale(ONE, yp, yptemp);

  /* Obtain pointers as required to the data array of vectors. */
  ydata     = N_VGetArrayPointer(y);
  ypdata    = N_VGetArrayPointer(yp);
  gtempdata = N_VGetArrayPointer(gtemp);
  ewtdata   = N_VGetArrayPointer(ewt);
  ytempdata = N_VGetArrayPointer(ytemp);
  yptempdata= N_VGetArrayPointer(yptemp);
  grefdata = N_VGetArrayPointer(gref);

  /* Call cfn and glocI to get base value of G(t,y,y'). */
  if (cfn != NULL) {
    retval = cfn(Nlocal, t, y, yp, f_data);
    if (retval != 0) return(retval);
  }

  retval = glocI(Nlocal, t, y, yp, gref, f_data);
  nge++;
  if (retval != 0) return(retval);

  /* Set bandwidth and number of column groups for band differencing. */
  width = mldq + mudq + 1;
  ngroups = MIN(width, Nlocal);

  /* Loop over groups. */
  for(group = 1; group <= ngroups; group++) {

    /* Loop over the components in this group. */
    for(j = group-1; j < Nlocal; j += width) {
      yj = ydata[j];
      ypj = ypdata[j];
      ewtj = ewtdata[j];

      /* Set increment inc to yj based on rel_yy*abs(yj), with
         adjustments using ypj and ewtj if this is small, and a further
         adjustment to give it the same sign as hh*ypj. */
      inc = dqrely*MAX(ABS(yj), MAX( ABS(h*ypj), ONE/ewtj));
      if (h*ypj < ZERO) inc = -inc;
      inc = (yj + inc) - yj;

      /* Increment yj and ypj. */
      ytempdata[j]  += gamma*inc;
      yptempdata[j] += inc;

    }

    /* Evaluate G with incremented y and yp arguments. */
    retval = glocI(Nlocal, t, ytemp, yptemp, gtemp, f_data);
    nge++;
    if (retval != 0) return(retval);

    /* Loop over components of the group again; restore ytemp and yptemp. */
    for(j = group-1; j < Nlocal; j += width) {
      yj  = ytempdata[j]  = ydata[j];
      ypj = yptempdata[j] = ypdata[j];
      ewtj = ewtdata[j];

      /* Set increment inc as before .*/
      inc = dqrely*MAX(ABS(yj), MAX( ABS(h*ypj), ONE/ewtj));
      if (h*ypj < ZERO) inc = -inc;
      inc = (yj + inc) - yj;

      /* Form difference quotients and load into savedP. */
      inc_inv = ONE/inc;
      col_j = BAND_COL(savedP,j);
      i1 = MAX(0, j-mukeep);
      i2 = MIN(j+mlkeep, Nlocal-1);
      for(i = i1; i <= i2; i++)
        BAND_COL_ELEM(col_j,i,j) = inc_inv * (gtempdata[i] - grefdata[i]);
    }
  }

  return(0);
}
