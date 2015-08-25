/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2006/11/22 00:12:48 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for the CPBAND linear solver.
 * -----------------------------------------------------------------
 */

/* 
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include <cpodes/cpodes_band.h>
#include "cpodes_direct_impl.h"
#include "cpodes_private.h"

#include <sundials/sundials_math.h>

/* 
 * =================================================================
 * PROTOTYPES FOR PRIVATE FUNCTIONS
 * =================================================================
 */

/* CPBAND linit, lsetup, lsolve, and lfree routines */
static int cpBandInit(CPodeMem cp_mem);
static int cpBandSetup(CPodeMem cp_mem, int convfail, 
                       N_Vector yP, N_Vector ypP, N_Vector fctP,
                       booleantype *jcurPtr,
                       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int cpBandSolve(CPodeMem cp_mem, N_Vector b, N_Vector weight,
                       N_Vector yC, N_Vector ypC, N_Vector fctC);
static void cpBandFree(CPodeMem cp_mem);

/*
 * =================================================================
 * READIBILITY REPLACEMENTS
 * =================================================================
 */

#define ode_type       (cp_mem->cp_ode_type)
#define lmm_type       (cp_mem->cp_lmm_type)
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
#define tempv          (cp_mem->cp_tempv)

#define linit          (cp_mem->cp_linit)
#define lsetup         (cp_mem->cp_lsetup)
#define lsolve         (cp_mem->cp_lsolve)
#define lfree          (cp_mem->cp_lfree)
#define lmem           (cp_mem->cp_lmem)
#define lsetup_exists  (cp_mem->cp_lsetup_exists)

#define mtype          (cpdls_mem->d_type)
#define n              (cpdls_mem->d_n)
#define jacE           (cpdls_mem->d_bjacE)
#define jacI           (cpdls_mem->d_bjacI)
#define M              (cpdls_mem->d_M)
#define mu             (cpdls_mem->d_mu)
#define ml             (cpdls_mem->d_ml)
#define smu            (cpdls_mem->d_smu)
#define pivots         (cpdls_mem->d_pivots)
#define savedJ         (cpdls_mem->d_savedJ)
#define nstlj          (cpdls_mem->d_nstlj)
#define nje            (cpdls_mem->d_nje)
#define nfeDQ          (cpdls_mem->d_nfeDQ)
#define J_data         (cpdls_mem->d_J_data)
#define last_flag      (cpdls_mem->d_last_flag)

/* 
 * =================================================================
 * EXPORTED FUNCTIONS
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * CPBand
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the band linear solver module.  CPBand first calls
 * the existing lfree routine if this is not NULL.  It then sets the
 * cp_linit, cp_lsetup, cp_lsolve, and cp_lfree fields in (*cpode_mem)
 * to be CPBandInit, CPBandSetup, CPBandSolve, and CPBandFree,
 * respectively.  It allocates memory for a structure of type
 * CPDlsMemRec and sets the cp_lmem field in (*cpode_mem) to the
 * address of this structure.  It sets lsetup_exists in (*cpode_mem) to be
 * TRUE, d_mu to be mupper, d_ml to be mlower, and initializes the d_bjacE
 * and d_bjacI to NULL.
 * Finally, it allocates memory for M, savedJ, and pivot.  The CPBand
 * return value is SUCCESS = 0, LMEM_FAIL = -1, or LIN_ILL_INPUT = -2.
 *
 * NOTE: The band linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, CPBand will first 
 *       test for compatible a compatible N_Vector internal
 *       representation by checking that the function 
 *       N_VGetArrayPointer exists.
 * -----------------------------------------------------------------
 */
                  
int CPBand(void *cpode_mem, int N, int mupper, int mlower)
{
  CPodeMem cp_mem;
  CPDlsMem cpdls_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPDIRECT_MEM_NULL, "CPBAND", "CPBand", MSGD_CPMEM_NULL);
    return(CPDIRECT_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  /* Test if the NVECTOR package is compatible with the BAND solver */
  if (tempv->ops->nvgetarraypointer == NULL) {
    cpProcessError(cp_mem, CPDIRECT_ILL_INPUT, "CPBAND", "CPBand", MSGD_BAD_NVECTOR);
    return(CPDIRECT_ILL_INPUT);
  }

  if (lfree != NULL) lfree(cp_mem);

  /* Set four main function fields in cp_mem */  
  linit  = cpBandInit;
  lsetup = cpBandSetup;
  lsolve = cpBandSolve;
  lfree  = cpBandFree;
  
  /* Get memory for CPDlsMemRec */
  cpdls_mem = NULL;
  cpdls_mem = (CPDlsMem) malloc(sizeof(CPDlsMemRec));
  if (cpdls_mem == NULL) {
    cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPBAND", "CPBand", MSGD_MEM_FAIL);
    return(CPDIRECT_MEM_FAIL);
  }

  /* Set matrix type */
  mtype = SUNDIALS_BAND;

  /* Set default Jacobian routine and Jacobian data */
  jacE = NULL;
  jacI = NULL;
  J_data = NULL;

  last_flag = CPDIRECT_SUCCESS;
  lsetup_exists = TRUE;
  
  /* Load problem dimension */
  n = N;

  /* Load half-bandwiths in cpdls_mem */
  ml = mlower;
  mu = mupper;

  /* Test ml and mu for legality */
  if ((ml < 0) || (mu < 0) || (ml >= N) || (mu >= N)) {
    cpProcessError(cp_mem, CPDIRECT_ILL_INPUT, "CPBAND", "CPBand", MSGD_BAD_SIZES);
    free(cpdls_mem);
    return(CPDIRECT_ILL_INPUT);
  }

  /* Set extended upper half-bandwidth for M (required for pivoting) */
  smu = MIN(N-1, mu + ml);

  /* Allocate memory for M, savedJ, and pivot arrays */
  M = NULL;
  pivots = NULL;
  savedJ = NULL;

  M = NewBandMat(N, mu, ml, smu);
  if (M == NULL) {
    cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPBAND", "CPBand", MSGD_MEM_FAIL);
    free(cpdls_mem);
    return(CPDIRECT_MEM_FAIL);
  }  
  pivots = NewIntArray(N);
  if (pivots == NULL) {
    cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPBAND", "CPBand", MSGD_MEM_FAIL);
    DestroyMat(M);
    free(cpdls_mem);
    return(CPDIRECT_MEM_FAIL);
  }
  if (ode_type == CP_EXPL) {
    savedJ = NewBandMat(N, mu, ml, mu);
    if (savedJ == NULL) {
      cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPBAND", "CPBand", MSGD_MEM_FAIL);
      DestroyMat(M);
      DestroyArray(pivots);
      free(cpdls_mem);
      return(CPDIRECT_MEM_FAIL);
    }
  }

  /* Attach linear solver memory to integrator memory */
  lmem = cpdls_mem;

  return(CPDIRECT_SUCCESS);
}

/* 
 * =================================================================
 *  PRIVATE FUNCTIONS FOR IMPLICIT INTEGRATION
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * cpBandInit
 * -----------------------------------------------------------------
 * This routine does remaining initializations specific to the band
 * linear solver.
 * -----------------------------------------------------------------
 */

static int cpBandInit(CPodeMem cp_mem)
{
  CPDlsMem cpdls_mem;

  cpdls_mem = (CPDlsMem) lmem;

  nje   = 0;
  nfeDQ = 0;
  nstlj = 0;

  if (ode_type == CP_EXPL && jacE == NULL) {
    jacE = cpDlsBandDQJacExpl;
    J_data = cp_mem;
  } 
  
  if (ode_type == CP_IMPL && jacI == NULL) {
    jacI = cpDlsBandDQJacImpl;
    J_data = cp_mem;
  }

  last_flag = CPDIRECT_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * cpBandSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the band linear solver.
 * It makes a decision whether or not to call the Jacobian evaluation
 * routine based on various state variables, and if not it uses the 
 * saved copy (for explicit ODE only). In any case, it constructs 
 * the Newton matrix M = I - gamma*J or M = F_y' - gamma*F_y, updates 
 * counters, and calls the band LU factorization routine.
 * -----------------------------------------------------------------
 */

static int cpBandSetup(CPodeMem cp_mem, int convfail, 
                       N_Vector yP, N_Vector ypP, N_Vector fctP,
                        booleantype *jcurPtr,
                        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  CPDlsMem cpdls_mem;
  booleantype jbad, jok;
  realtype dgamma;
  int ier, retval;

  cpdls_mem = (CPDlsMem) lmem;

  switch (ode_type) {

  case CP_EXPL:

    /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
    dgamma = ABS((gamma/gammap) - ONE);
    jbad = (nst == 0) || (nst > nstlj + CPD_MSBJ) ||
      ((convfail == CP_FAIL_BAD_J) && (dgamma < CPD_DGMAX)) ||
      (convfail == CP_FAIL_OTHER);
    jok = !jbad;
    
    if (jok) {
      
      /* If jok = TRUE, use saved copy of J */
      *jcurPtr = FALSE;
      BandCopy(savedJ, M, mu, ml);
      
    } else {
      
      /* If jok = FALSE, call jac routine for new J value */
      nje++;
      nstlj = nst;
      *jcurPtr = TRUE;

      BandZero(M);
      retval = jacE(n, mu, ml, tn, yP, fctP, M, J_data, tmp1, tmp2, tmp3);
      if (retval == 0) {
        BandCopy(M, savedJ, mu, ml);
      }else if (retval < 0) {
        cpProcessError(cp_mem, CPDIRECT_JACFUNC_UNRECVR, "CPBAND", "CPBandSetup", MSGD_JACFUNC_FAILED);
        last_flag = CPDIRECT_JACFUNC_UNRECVR;
        return(-1);
      }else if (retval > 0) {
        last_flag = CPDIRECT_JACFUNC_RECVR;
        return(1);
      }

    }
  
    /* Scale and add I to get M = I - gamma*J */
    BandScale(-gamma, M);
    BandAddI(M);

    break;

  case CP_IMPL:

    BandZero(M);
    retval = jacI(n, mu, ml, tn, gamma, yP, ypP, fctP, M, J_data, tmp1, tmp2, tmp3);
    if (retval < 0) {
      cpProcessError(cp_mem, CPDIRECT_JACFUNC_UNRECVR, "CPBAND", "CPBandSetup", MSGD_JACFUNC_FAILED);
      last_flag = CPDIRECT_JACFUNC_UNRECVR;
      return(-1);
    } else if (retval > 0) {
      last_flag = CPDIRECT_JACFUNC_RECVR;
      return(+1);
    }

    break;

  }

  /* Do LU factorization of M */
  ier = BandGBTRF(M, pivots);

  /* Return 0 if the LU was complete; otherwise return 1 */
  last_flag = ier;
  if (ier>0) return(1);
  return(0);
}

/*
 * -----------------------------------------------------------------
 * cpBandSolve
 * -----------------------------------------------------------------
 * This routine handles the solve operation for the band linear solver
 * by calling the band backsolve routine.  The return value is 0.
 * -----------------------------------------------------------------
 */

static int cpBandSolve(CPodeMem cp_mem, N_Vector b, N_Vector weight,
                       N_Vector yC, N_Vector ypC, N_Vector fctC)
{
  CPDlsMem cpdls_mem;
  realtype *bd;

  cpdls_mem = (CPDlsMem) lmem;

  bd = N_VGetArrayPointer(b);

  BandGBTRS(M, pivots, bd);

  /* For BDF, scale the correction to account for change in gamma */
  if ((lmm_type == CP_BDF) && (gamrat != ONE)) {
    N_VScale(TWO/(ONE + gamrat), b, b);
  }

  last_flag = CPDIRECT_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * cpBandFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the band linear solver.
 * -----------------------------------------------------------------
 */

static void cpBandFree(CPodeMem cp_mem)
{
  CPDlsMem cpdls_mem;

  cpdls_mem = (CPDlsMem) lmem;

  DestroyMat(M);
  DestroyArray(pivots);
  if (ode_type == CP_EXPL) DestroyMat(savedJ);
  free(cpdls_mem); 
  cpdls_mem = NULL;
}

