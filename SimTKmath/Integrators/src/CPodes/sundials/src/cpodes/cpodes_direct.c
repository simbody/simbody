/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006/11/22 00:12:48 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for the CPDIRECT linear solvers
 * -----------------------------------------------------------------
 */

/*
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include "cpodes_direct_impl.h"
#include "cpodes_private.h"

#include <sundials/sundials_math.h>

/*
 * =================================================================
 * FUNCTION SPECIFIC CONSTANTS
 * =================================================================
 */

/* Constant for DQ Jacobian approximation */
#define MIN_INC_MULT RCONST(1000.0)

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/*
 * =================================================================
 * READIBILITY REPLACEMENTS
 * =================================================================
 */

#define ode_type       (cp_mem->cp_ode_type)
#define fe             (cp_mem->cp_fe)
#define fi             (cp_mem->cp_fi)
#define f_data         (cp_mem->cp_f_data)
#define uround         (cp_mem->cp_uround)
#define nst            (cp_mem->cp_nst)
#define tn             (cp_mem->cp_tn)
#define h              (cp_mem->cp_h)
#define gamma          (cp_mem->cp_gamma)
#define gammap         (cp_mem->cp_gammap)
#define gamrat         (cp_mem->cp_gamrat)
#define ewt            (cp_mem->cp_ewt)
#define lmem           (cp_mem->cp_lmem)

#define cfun           (cp_mem->cp_cfun)
#define c_data         (cp_mem->cp_c_data)
#define lmemP          (cp_mem->cp_lmemP)

#define mtype          (cpdls_mem->d_type)
#define n              (cpdls_mem->d_n)
#define ml             (cpdls_mem->d_ml)
#define mu             (cpdls_mem->d_mu)
#define smu            (cpdls_mem->d_smu)
#define djacE          (cpdls_mem->d_djacE)
#define djacI          (cpdls_mem->d_djacI)
#define bjacE          (cpdls_mem->d_bjacE)
#define bjacI          (cpdls_mem->d_bjacI)
#define M              (cpdls_mem->d_M)
#define savedJ         (cpdls_mem->d_savedJ)
#define pivots         (cpdls_mem->d_pivots)
#define nstlj          (cpdls_mem->d_nstlj)
#define nje            (cpdls_mem->d_nje)
#define nfeDQ          (cpdls_mem->d_nfeDQ)
#define J_data         (cpdls_mem->d_J_data)
#define last_flag      (cpdls_mem->d_last_flag)

#define djacP          (cpdlsP_mem->d_jacP)
#define JP_data        (cpdlsP_mem->d_JP_data)
#define njeP           (cpdlsP_mem->d_njeP)
#define nceDQ          (cpdlsP_mem->d_nceDQ)

/*
 * =================================================================
 * EXPORTED FUNCTIONS FOR IMPLICIT INTEGRATION
 * =================================================================
 */

/*
 * CPDlsSetJacFn specifies the (dense or band) Jacobian function.
 */
int CPDlsSetJacFn(void *cpode_mem, void *jac, void *jac_data)
{
  CPodeMem cp_mem;
  CPDlsMem cpdls_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPDIRECT_MEM_NULL, "CPDIRECT", "CPDlsSetJacFn", MSGD_CPMEM_NULL);
    return(CPDIRECT_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(cp_mem, CPDIRECT_LMEM_NULL, "CPDIRECT", "CPDlsSetJacFn", MSGD_LMEM_NULL);
    return(CPDIRECT_LMEM_NULL);
  }
  cpdls_mem = (CPDlsMem) lmem;

  switch (ode_type) {

  case CP_EXPL:
    if (mtype == SUNDIALS_DENSE)
      djacE = (CPDlsDenseJacExplFn) jac;
    else if (mtype == SUNDIALS_BAND)
      bjacE = (CPDlsBandJacExplFn) jac;
    break;

  case CP_IMPL:
    if (mtype == SUNDIALS_DENSE)
      djacI = (CPDlsDenseJacImplFn) jac;
    else if (mtype == SUNDIALS_BAND)
      bjacI = (CPDlsBandJacImplFn) jac;
    break;

  }

  J_data = jac_data;

  return(CPDIRECT_SUCCESS);
}

/*
 * CPDlsGetWorkSpace returns the length of workspace allocated for the
 * CPDIRECT linear solver.
 */
int CPDlsGetWorkSpace(void *cpode_mem, long int *lenrwLS, long int *leniwLS)
{
  CPodeMem cp_mem;
  CPDlsMem cpdls_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPDIRECT_MEM_NULL, "CPDIRECT", "CPDlsGetWorkSpace", MSGD_CPMEM_NULL);
    return(CPDIRECT_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(cp_mem, CPDIRECT_LMEM_NULL, "CPDIRECT", "CPDlsGetWorkSpace", MSGD_LMEM_NULL);
    return(CPDIRECT_LMEM_NULL);
  }
  cpdls_mem = (CPDlsMem) lmem;

  switch (ode_type) {

  case CP_EXPL:

    if (mtype == SUNDIALS_DENSE) {
      *lenrwLS = 2*n*n;
      *leniwLS = n;
    } else if (mtype == SUNDIALS_BAND) {
      *lenrwLS = n*(smu + mu + 2*ml + 2);
      *leniwLS = n;
    }

    break;

  case CP_IMPL:

    if (mtype == SUNDIALS_DENSE) {
      *lenrwLS = n*n;
      *leniwLS = n;
    } else if (mtype == SUNDIALS_BAND) {
      *lenrwLS = n*(smu + ml + 2);
      *leniwLS = n;
    }

    break;

  }

  return(CPDIRECT_SUCCESS);
}

/*
 * CPDlsGetNumJacEvals returns the number of Jacobian evaluations.
 */
int CPDlsGetNumJacEvals(void *cpode_mem, long int *njevals)
{
  CPodeMem cp_mem;
  CPDlsMem cpdls_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPDIRECT_MEM_NULL, "CPDIRECT", "CPDlsGetNumJacEvals", MSGD_CPMEM_NULL);
    return(CPDIRECT_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(cp_mem, CPDIRECT_LMEM_NULL, "CPDIRECT", "CPDlsGetNumJacEvals", MSGD_LMEM_NULL);
    return(CPDIRECT_LMEM_NULL);
  }
  cpdls_mem = (CPDlsMem) lmem;

  *njevals = nje;

  return(CPDIRECT_SUCCESS);
}

/*
 * CPDlsGetNumFctEvals returns the number of calls to the ODE function
 * needed for the DQ Jacobian approximation.
 */
int CPDlsGetNumFctEvals(void *cpode_mem, long int *nfevalsLS)
{
  CPodeMem cp_mem;
  CPDlsMem cpdls_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPDIRECT_MEM_NULL, "CPDIRECT", "CPDlsGetNumRhsEvals", MSGD_CPMEM_NULL);
    return(CPDIRECT_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(cp_mem, CPDIRECT_LMEM_NULL, "CPDIRECT", "CPDlsGetNumRhsEvals", MSGD_LMEM_NULL);
    return(CPDIRECT_LMEM_NULL);
  }
  cpdls_mem = (CPDlsMem) lmem;

  *nfevalsLS = nfeDQ;

  return(CPDIRECT_SUCCESS);
}

/*
 * CPDlsGetReturnFlagName returns the name associated with a CPDIRECT
 * return value.
 */
char *CPDlsGetReturnFlagName(int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case CPDIRECT_SUCCESS:
    sprintf(name,"CPDIRECT_SUCCESS");
    break;
  case CPDIRECT_MEM_NULL:
    sprintf(name,"CPDIRECT_MEM_NULL");
    break;
  case CPDIRECT_LMEM_NULL:
    sprintf(name,"CPDIRECT_LMEM_NULL");
    break;
  case CPDIRECT_ILL_INPUT:
    sprintf(name,"CPDIRECT_ILL_INPUT");
    break;
  case CPDIRECT_MEM_FAIL:
    sprintf(name,"CPDIRECT_MEM_FAIL");
    break;
  case CPDIRECT_JACFUNC_UNRECVR:
    sprintf(name,"CPDIRECT_JACFUNC_UNRECVR");
    break;
  case CPDIRECT_JACFUNC_RECVR:
    sprintf(name,"CPDIRECT_JACFUNC_RECVR");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}

/*
 * CPDlsGetLastFlag returns the last flag set in a CPDIRECT function.
 */
int CPDlsGetLastFlag(void *cpode_mem, int *flag)
{
  CPodeMem cp_mem;
  CPDlsMem cpdls_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPDIRECT_MEM_NULL, "CPDIRECT", "CPDlsGetLastFlag", MSGD_CPMEM_NULL);
    return(CPDIRECT_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(cp_mem, CPDIRECT_LMEM_NULL, "CPDIRECT", "CPDlsGetLastFlag", MSGD_LMEM_NULL);
    return(CPDIRECT_LMEM_NULL);
  }
  cpdls_mem = (CPDlsMem) lmem;

  *flag = last_flag;

  return(CPDIRECT_SUCCESS);
}

/*
 * =================================================================
 * EXPORTED FUNCTIONS FOR PROJECTION
 * =================================================================
 */

/*
 * CPDlsProjSetJacFn specifies the constraint Jacobian function.
 */
int CPDlsProjSetJacFn(void *cpode_mem, void *jacP, void *jacP_data)
{
  CPodeMem cp_mem;
  CPDlsProjMem cpdlsP_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPDIRECT_MEM_NULL, "CPDIRECT", "CPDlsProjSetJacFn", MSGD_CPMEM_NULL);
    return(CPDIRECT_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmemP == NULL) {
    cpProcessError(cp_mem, CPDIRECT_LMEM_NULL, "CPDIRECT", "CPDlsProjSetJacFn", MSGD_LMEM_NULL);
    return(CPDIRECT_LMEM_NULL);
  }
  cpdlsP_mem = (CPDlsProjMem) lmemP;

  djacP = (CPDlsDenseProjJacFn) jacP;
  JP_data = jacP_data;

  return(CPDIRECT_SUCCESS);
}

/*
 * CPDlsProjGetNumJacEvals returns the number of constraint Jacobian evaluations
 */
int CPDlsProjGetNumJacEvals(void *cpode_mem, long int *njPevals)
{
  CPodeMem cp_mem;
  CPDlsProjMem cpdlsP_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPDIRECT_MEM_NULL, "CPDIRECT", "CPDlsProjGetNumJacEvals", MSGD_CPMEM_NULL);
    return(CPDIRECT_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmemP == NULL) {
    cpProcessError(cp_mem, CPDIRECT_LMEM_NULL, "CPDIRECT", "CPDlsProjGetNumJacEvals", MSGD_LMEM_NULL);
    return(CPDIRECT_LMEM_NULL);
  }
  cpdlsP_mem = (CPDlsProjMem) lmemP;

  *njPevals = njeP;

  return(CPDIRECT_SUCCESS);
}

/*
 * CPDlsProjGetNumFctEvals returns the number of constraint function
 * evaluations for computing the DQ constraint Jacobian.
 */
int CPDlsProjGetNumFctEvals(void *cpode_mem, long int *ncevalsLS)
{
  CPodeMem cp_mem;
  CPDlsProjMem cpdlsP_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPDIRECT_MEM_NULL, "CPDIRECT", "CPDlsProjGetNumFctEvals", MSGD_CPMEM_NULL);
    return(CPDIRECT_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmemP == NULL) {
    cpProcessError(cp_mem, CPDIRECT_LMEM_NULL, "CPDIRECT", "CPDlsProjGetNumFctEvals", MSGD_LMEM_NULL);
    return(CPDIRECT_LMEM_NULL);
  }
  cpdlsP_mem = (CPDlsProjMem) lmemP;

  *ncevalsLS = nceDQ;

  return(CPDIRECT_SUCCESS);
}

/*
 * =================================================================
 * DENSE DQ JACOBIAN APPROXIMATIONS FOR IMPLICIT INTEGRATION
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * cpDlsDenseDQJacExpl
 * -----------------------------------------------------------------
 * This routine generates a dense difference quotient approximation to
 * the Jacobian of f(t,y). It assumes that a dense matrix of type
 * DlsMat is stored column-wise, and that elements within each column
 * are contiguous. The address of the jth column of J is obtained via
 * the macro DENSE_COL and this pointer is associated with an N_Vector
 * using the N_VGetArrayPointer/N_VSetArrayPointer functions.
 * Finally, the actual computation of the jth column of the Jacobian is
 * done with a call to N_VLinearSum.
 * -----------------------------------------------------------------
 */

int cpDlsDenseDQJacExpl(int N, realtype t,
                        N_Vector y, N_Vector fy,
                        DlsMat Jac, void *jac_data,
                        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype fnorm, minInc, inc, inc_inv, yjsaved, srur;
  realtype *tmp2_data, *y_data, *ewt_data;
  N_Vector ftemp, jthCol;
  int j;
  int retval = 0;

  CPodeMem cp_mem;
  CPDlsMem cpdls_mem;

  /* jac_data points to cpode_mem */
  cp_mem = (CPodeMem) jac_data;
  cpdls_mem = (CPDlsMem) lmem;

  /* Save pointer to the array in tmp2 */
  tmp2_data = N_VGetArrayPointer(tmp2);

  /* Rename work vectors for readibility */
  ftemp = tmp1;
  jthCol = tmp2;

  /* Obtain pointers to the data for ewt, y */
  ewt_data = N_VGetArrayPointer(ewt);
  y_data   = N_VGetArrayPointer(y);

  /* Set minimum increment based on uround and norm of f */
  srur = RSqrt(uround);
  fnorm = N_VWrmsNorm(fy, ewt);
  minInc = (fnorm != ZERO) ?
           (MIN_INC_MULT * ABS(h) * uround * N * fnorm) : ONE;

  for (j = 0; j < N; j++) {

    /* Generate the jth col of J(tn,y) */

    N_VSetArrayPointer(DENSE_COL(Jac,j), jthCol);

    yjsaved = y_data[j];
    inc = MAX(srur*ABS(yjsaved), minInc/ewt_data[j]);
    y_data[j] += inc;

    retval = fe(t, y, ftemp, f_data);
    nfeDQ++;
    if (retval != 0) break;

    y_data[j] = yjsaved;

    inc_inv = ONE/inc;
    N_VLinearSum(inc_inv, ftemp, -inc_inv, fy, jthCol);

    DENSE_COL(Jac,j) = N_VGetArrayPointer(jthCol);
  }

  /* Restore original array pointer in tmp2 */
  N_VSetArrayPointer(tmp2_data, tmp2);

  return(retval);
}

/*
 * -----------------------------------------------------------------
 * cpDlsDenseDQJacImpl
 * -----------------------------------------------------------------
 * This routine generates a dense difference quotient approximation to
 * the Jacobian F_y' + gamma*F_y. It assumes that a dense matrix of type
 * DlsMat is stored column-wise, and that elements within each column
 * are contiguous. The address of the jth column of J is obtained via
 * the macro DENSE_COL and this pointer is associated with an N_Vector
 * using the N_VGetArrayPointer/N_VSetArrayPointer functions.
 * Finally, the actual computation of the jth column of the Jacobian is
 * done with a call to N_VLinearSum.
 * -----------------------------------------------------------------
 */
int cpDlsDenseDQJacImpl(int N, realtype t, realtype gm,
                        N_Vector y, N_Vector yp, N_Vector r,
                        DlsMat Jac, void *jac_data,
                        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype inc, inc_inv, yj, ypj, srur;
  realtype *tmp2_data, *y_data, *yp_data, *ewt_data;
  N_Vector ftemp, jthCol;
  int j;
  int retval = 0;

  CPodeMem cp_mem;
  CPDlsMem  cpdls_mem;

  /* jac_data points to cpode_mem */
  cp_mem = (CPodeMem) jac_data;
  cpdls_mem = (CPDlsMem) lmem;

  /* Save pointer to the array in tmp2 */
  tmp2_data = N_VGetArrayPointer(tmp2);

  /* Rename work vectors for readibility */
  ftemp = tmp1;
  jthCol = tmp2;

  /* Obtain pointers to the data for ewt, y, and yp */
  ewt_data = N_VGetArrayPointer(ewt);
  y_data   = N_VGetArrayPointer(y);
  yp_data  = N_VGetArrayPointer(yp);

  /* Set minimum increment based on uround and norm of f */
  srur = RSqrt(uround);

  /* Generate each column of the Jacobian M = F_y' + gamma * F_y
     as delta(F)/delta(y_j). */
  for (j = 0; j < N; j++) {

    /* Set data address of jthCol, and save y_j and yp_j values. */
    N_VSetArrayPointer(DENSE_COL(Jac,j), jthCol);
    yj = y_data[j];
    ypj = yp_data[j];

    /* Set increment inc to y_j based on sqrt(uround)*abs(y_j), with
    adjustments using yp_j and ewt_j if this is small, and a further
    adjustment to give it the same sign as h*yp_j. */
    inc = MAX( srur * MAX( ABS(yj), ABS(h*ypj) ) , ONE/ewt_data[j] );
    if (h*ypj < ZERO) inc = -inc;
    inc = (yj + inc) - yj;

    /* Increment y_j and yp_j, call res, and break on error return. */
    y_data[j]  += gamma*inc;
    yp_data[j] += inc;
    retval = fi(t, y, yp, ftemp, f_data);
    nfeDQ++;
    if (retval != 0) break;

    /* Generate the jth col of J(tn,y) */
    inc_inv = ONE/inc;
    N_VLinearSum(inc_inv, ftemp, -inc_inv, r, jthCol);

    DENSE_COL(Jac,j) = N_VGetArrayPointer(jthCol);

    /* Reset y_j, yp_j */
    y_data[j] = yj;
    yp_data[j] = ypj;
  }

  /* Restore original array pointer in tmp2 */
  N_VSetArrayPointer(tmp2_data, tmp2);

  return(retval);
}

/*
 * =================================================================
 *  BAND DQ JACOBIAN APPROXIMATIONS FOR IMPLICIT INTEGRATION
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * cpDlsBandDQJacExpl
 * -----------------------------------------------------------------
 * This routine generates a banded difference quotient approximation to
 * the Jacobian of f(t,y).  It assumes that a band matrix of type
 * DlsMat is stored column-wise, and that elements within each column
 * are contiguous. This makes it possible to get the address of a column
 * of J via the macro BAND_COL and to write a simple for loop to set
 * each of the elements of a column in succession.
 * -----------------------------------------------------------------
 */

int cpDlsBandDQJacExpl(int N, int mupper, int mlower,
                       realtype t, N_Vector y, N_Vector fy,
                       DlsMat Jac, void *jac_data,
                       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  N_Vector ftemp, ytemp;
  realtype fnorm, minInc, inc, inc_inv, srur;
  realtype *col_j, *ewt_data, *fy_data, *ftemp_data, *y_data, *ytemp_data;
  int group, i, j, width, ngroups, i1, i2;
  int retval = 0;

  CPodeMem cp_mem;
  CPDlsMem cpdls_mem;

  /* jac_dat points to cpode_mem */
  cp_mem = (CPodeMem) jac_data;
  cpdls_mem = (CPDlsMem) lmem;

  /* Rename work vectors for use as temporary values of y and f */
  ftemp = tmp1;
  ytemp = tmp2;

  /* Obtain pointers to the data for ewt, fy, ftemp, y, ytemp */
  ewt_data   = N_VGetArrayPointer(ewt);
  fy_data    = N_VGetArrayPointer(fy);
  ftemp_data = N_VGetArrayPointer(ftemp);
  y_data     = N_VGetArrayPointer(y);
  ytemp_data = N_VGetArrayPointer(ytemp);

  /* Load ytemp with y = predicted y vector */
  N_VScale(ONE, y, ytemp);

  /* Set minimum increment based on uround and norm of f */
  srur = RSqrt(uround);
  fnorm = N_VWrmsNorm(fy, ewt);
  minInc = (fnorm != ZERO) ?
           (MIN_INC_MULT * ABS(h) * uround * N * fnorm) : ONE;

  /* Set bandwidth and number of column groups for band differencing */
  width = mlower + mupper + 1;
  ngroups = MIN(width, N);

  /* Loop over column groups. */
  for (group=1; group <= ngroups; group++) {

    /* Increment all y_j in group */
    for(j=group-1; j < N; j+=width) {
      inc = MAX(srur*ABS(y_data[j]), minInc/ewt_data[j]);
      ytemp_data[j] += inc;
    }

    /* Evaluate f with incremented y */

    retval = fe(tn, ytemp, ftemp, f_data);
    nfeDQ++;
    if (retval != 0) break;

    /* Restore ytemp, then form and load difference quotients */
    for (j=group-1; j < N; j+=width) {
      ytemp_data[j] = y_data[j];
      col_j = BAND_COL(Jac,j);
      inc = MAX(srur*ABS(y_data[j]), minInc/ewt_data[j]);
      inc_inv = ONE/inc;
      i1 = MAX(0, j-mupper);
      i2 = MIN(j+mlower, N-1);
      for (i=i1; i <= i2; i++)
        BAND_COL_ELEM(col_j,i,j) = inc_inv * (ftemp_data[i] - fy_data[i]);
    }
  }

  return(retval);
}

/*
 * -----------------------------------------------------------------
 * cpDlsBandDQJacImpl
 * -----------------------------------------------------------------
 */

int cpDlsBandDQJacImpl(int N, int mupper, int mlower,
                       realtype t, realtype gm,
                       N_Vector y, N_Vector yp, N_Vector r,
                       DlsMat Jac, void *jac_data,
                       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  N_Vector ftemp, ytemp, yptemp;
  realtype inc, inc_inv, yj, ypj, srur, ewtj;
  realtype *y_data, *yp_data, *ewt_data;
  realtype *ytemp_data, *yptemp_data, *ftemp_data, *r_data, *col_j;
  int i, j, i1, i2, width, group, ngroups;
  int retval = 0;

  CPodeMem cp_mem;
  CPDlsMem cpdls_mem;

  /* jac_data points to cpode_mem */
  cp_mem = (CPodeMem) jac_data;
  cpdls_mem = (CPDlsMem) lmem;

  ftemp = tmp1; /* Rename work vector for use as the perturbed residual. */
  ytemp = tmp2; /* Rename work vector for use as a temporary for yy. */
  yptemp= tmp3; /* Rename work vector for use as a temporary for yp. */

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
  width = mlower + mupper + 1;
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
    retval = fi(tn, ytemp, yptemp, ftemp, f_data);
    nfeDQ++;
    if (retval != 0) break;

    /* Loop over the indices j in this group again. */
    for (j=group-1; j<N; j+=width) {

      /* Reset ytemp and yptemp components that were perturbed. */
      yj = ytemp_data[j]  = y_data[j];
      ypj = yptemp_data[j] = yp_data[j];
      col_j = BAND_COL(Jac, j);
      ewtj = ewt_data[j];

      /* Set increment inc exactly as above. */
      inc = MAX( srur * MAX( ABS(yj), ABS(h*ypj) ) , ONE/ewtj );
      if (h*ypj < ZERO) inc = -inc;
      inc = (yj + inc) - yj;

      /* Load the difference quotient Jacobian elements for column j. */
      inc_inv = ONE/inc;
      i1 = MAX(0, j-mupper);
      i2 = MIN(j+mlower,N-1);
      for (i=i1; i<=i2; i++)
        BAND_COL_ELEM(col_j,i,j) = inc_inv*(ftemp_data[i]-r_data[i]);
    }

  }

  return(retval);
}

/*
 * =================================================================
 *  DENSE DQ JACOBIAN APPROXIMATIONS FOR PROJECTION
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * cpDlsDenseProjDQJac
 * -----------------------------------------------------------------
 * This routine generates a dense difference quotient approximation
 * to the transpose of the Jacobian of c(t,y). It loads it into a
 * dense matrix of type DlsMat stored column-wise with elements
 * within each column contiguous. The address of the jth column of
 * J is obtained via the macro DENSE_COL and this pointer is
 * associated with an N_Vector using the N_VGetArrayPointer and
 * N_VSetArrayPointer functions.
 * Finally, the actual computation of the jth column of the Jacobian
 * transposed is done with a call to N_VLinearSum.
 * -----------------------------------------------------------------
 */
int cpDlsDenseProjDQJac(int Nc, int Ny, realtype t,
                        N_Vector y, N_Vector cy,
                        DlsMat Jac, void *jac_data,
                        N_Vector c_tmp1, N_Vector c_tmp2)
{
  realtype inc, inc_inv, yj, srur;
  realtype *y_data, *ewt_data, *jthCol_data;
  N_Vector ctemp, jthCol;
  int i, j;
  int retval = 0;

  CPodeMem cp_mem;
  CPDlsProjMem cpdlsP_mem;

  /* jac_data points to cpode_mem */
  cp_mem = (CPodeMem) jac_data;
  cpdlsP_mem = (CPDlsProjMem) lmemP;

  /* Rename work vectors for readibility */
  ctemp  = c_tmp1;
  jthCol = c_tmp2;

  /* Obtain pointers to the data for ewt and y */
  ewt_data = N_VGetArrayPointer(ewt);
  y_data   = N_VGetArrayPointer(y);

  /* Obtain pointer to the data for jthCol */
  jthCol_data = N_VGetArrayPointer(jthCol);

  /* Set minimum increment based on uround and norm of f */
  srur = RSqrt(uround);

  /* Generate each column of the Jacobian G = dc/dy as delta(c)/delta(y_j). */
  for (j = 0; j < Ny; j++) {

    /* Save the y_j values. */
    yj = y_data[j];

    /* Set increment inc to y_j based on sqrt(uround)*abs(y_j),
       with an adjustment using ewt_j if this is small */
    inc = MAX( srur * ABS(yj) , ONE/ewt_data[j] );
    inc = (yj + inc) - yj;

    /* Increment y_j, call cfun, and break on error return. */
    y_data[j]  += inc;
    retval = cfun(t, y, ctemp, c_data);
    nceDQ++;
    if (retval != 0) break;

    /* Generate the jth col of G(tn,y) */
    inc_inv = ONE/inc;
    N_VLinearSum(inc_inv, ctemp, -inc_inv, cy, jthCol);

    /* Copy the j-th column of G into the j-th row of Jac */
    for (i = 0; i < Nc ; i++) {
      DENSE_ELEM(Jac,j,i) = jthCol_data[i];
    }

    /* Reset y_j */
    y_data[j] = yj;
  }

  return(retval);

}


