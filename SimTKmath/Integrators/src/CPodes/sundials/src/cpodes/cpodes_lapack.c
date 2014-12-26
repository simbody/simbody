/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2006/12/01 22:48:57 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for a CPODES dense linear solver
 * using BLAS and LAPACK functions.
 * -----------------------------------------------------------------
 */

/*
 * NOTE: the only operations (all O(n)) that do not use Blas/Lapack 
 *       functions are:
 *
 *   - matrix plus identity (I-gamma*J)
 *     (in lsetup -> use functions from sundials_lapack)
 *   - diagonal matrix times vector (D^(-1)*x) 
 *     (in lsolveP for QR, QRP, and SC -> hardcoded)
 *   - permutation of a diagonal matrix (P*D*P^T) 
 *     (in lsolveP for LU; P uses pivots from dgetrv -> hardcoded)
 *   - permutation matrix times vector (P^T*x)
 *     (in lsolveP for LU; P uses pivots from dgetr -> hardcoded)
 *   - permutation matrix times vector (P^T*b)
 *     (in lsolveP for QRP; P uses pivots from dgeqp3 -> hardcoded)
 */

/*
 * TODO:
 *   
 *   cplLUcomputeKD
 *   cplQRcomputeKD
 *   cplSCcomputeKD
 */

/* 
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include <cpodes/cpodes_lapack.h>
#include "cpodes_direct_impl.h"
#include "cpodes_private.h"

#include <sundials/sundials_math.h>

/* 
 * =================================================================
 * PROTOTYPES FOR PRIVATE FUNCTIONS
 * =================================================================
 */

/* CPLAPACK DENSE linit, lsetup, lsolve, and lfree routines */ 
static int cpLapackDenseInit(CPodeMem cp_mem);
static int cpLapackDenseSetup(CPodeMem cp_mem, int convfail, 
                              N_Vector yP, N_Vector ypP, N_Vector fctP, 
                              booleantype *jcurPtr,
                              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int cpLapackDenseSolve(CPodeMem cp_mem, N_Vector b, N_Vector weight,
                              N_Vector yC, N_Vector ypC, N_Vector fctC);
static void cpLapackDenseFree(CPodeMem cp_mem);

/* CPLAPACK BAND linit, lsetup, lsolve, and lfree routines */ 
static int cpLapackBandInit(CPodeMem cp_mem);
static int cpLapackBandSetup(CPodeMem cp_mem, int convfail, 
                             N_Vector yP, N_Vector ypP, N_Vector fctP, 
                             booleantype *jcurPtr,
                             N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int cpLapackBandSolve(CPodeMem cp_mem, N_Vector b, N_Vector weight,
                             N_Vector yC, N_Vector ypC, N_Vector fctC);
static void cpLapackBandFree(CPodeMem cp_mem);


/* CPLAPACK DENSE linitP, lsetupP, lsolveP, lmultP, and lfreeP routines */
static int cpLapackDenseProjInit(CPodeMem cp_mem);
static int cpLapackDenseProjSetup(CPodeMem cp_mem, N_Vector y, N_Vector cy,
                                  N_Vector c_tmp1, N_Vector c_tmp2, N_Vector s_tmp1);
static int cpLapackDenseProjSolve(CPodeMem cp_mem, N_Vector b, N_Vector x,
                                  N_Vector y, N_Vector cy,
                                  N_Vector c_tmp1, N_Vector s_tmp1);
static void cpLapackDenseProjMult(CPodeMem cp_mem, N_Vector x, N_Vector Gx);
static void cpLapackDenseProjFree(CPodeMem cp_mem);

/* Private functions for LU, QR, and SC projection */
static void cplLUcomputeKD(CPodeMem cp_mem, N_Vector d);
static void cplQRcomputeKD(CPodeMem cp_mem, N_Vector d);
static void cplSCcomputeKD(CPodeMem cp_mem, N_Vector d);

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

#define linit          (cp_mem->cp_linit)
#define lsetup         (cp_mem->cp_lsetup)
#define lsolve         (cp_mem->cp_lsolve)
#define lfree          (cp_mem->cp_lfree)
#define lmem           (cp_mem->cp_lmem)
#define tempv          (cp_mem->cp_tempv)
#define lsetup_exists  (cp_mem->cp_lsetup_exists)

#define pnorm          (cp_mem->cp_proj_norm)

#define linitP         (cp_mem->cp_linitP)
#define lsetupP        (cp_mem->cp_lsetupP)
#define lsolveP        (cp_mem->cp_lsolveP)
#define lmultP         (cp_mem->cp_lmultP)
#define lfreeP         (cp_mem->cp_lfreeP)
#define lmemP          (cp_mem->cp_lmemP)
#define lsetupP_exists (cp_mem->cp_lsetupP_exists)

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

#define ny             (cpdlsP_mem->d_ny)
#define nc             (cpdlsP_mem->d_nc)
#define nr             (cpdlsP_mem->d_nr)
#define djacP          (cpdlsP_mem->d_jacP)
#define JP_data        (cpdlsP_mem->d_JP_data)
#define ftype          (cpdlsP_mem->d_ftype)
#define G              (cpdlsP_mem->d_G)
#define savedG         (cpdlsP_mem->d_savedG)
#define K              (cpdlsP_mem->d_K)
#define pivotsP        (cpdlsP_mem->d_pivotsP)
#define beta           (cpdlsP_mem->d_beta)
#define wrk            (cpdlsP_mem->d_wrk)
#define len_wrk        (cpdlsP_mem->d_len_wrk)
#define nstljP         (cpdlsP_mem->d_nstljP)
#define njeP           (cpdlsP_mem->d_njeP)
#define nceDQ          (cpdlsP_mem->d_nceDQ)

/* 
 * =================================================================
 * EXPORTED FUNCTIONS
 * =================================================================
 */
              
/*
 * -----------------------------------------------------------------
 * CPLapackDense
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the linear solver module.  CPLapackDense first
 * calls the existing lfree routine if this is not NULL.  Then it sets
 * the cp_linit, cp_lsetup, cp_lsolve, cp_lfree fields in (*cpode_mem)
 * to be cpLapackDenseInit, cpLapackDenseSetup, cpLapackDenseSolve, 
 * and cpLapackDenseFree, respectively.  It allocates memory for a 
 * structure of type CPDlsMemRec and sets the cp_lmem field in 
 * (*cpode_mem) to the address of this structure.  It sets lsetup_exists 
 * in (*cpode_mem) to TRUE, and the d_jac field to the default 
 * cpDlsDenseDQJac. Finally, it allocates memory for M, pivots, and 
 * (if needed) savedJ.
 * The return value is SUCCESS = 0, or LMEM_FAIL = -1.
 *
 * NOTE: The dense linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, CPLapackDense will first 
 *       test for a compatible N_Vector internal representation 
 *       by checking that N_VGetArrayPointer and N_VSetArrayPointer 
 *       exist.
 * -----------------------------------------------------------------
 */
int CPLapackDense(void *cpode_mem, int N)
{
  CPodeMem cp_mem;
  CPDlsMem cpdls_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPDIRECT_MEM_NULL, "CPLAPACK", "CPLapackDense", MSGD_CPMEM_NULL);
    return(CPDIRECT_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  /* Test if the NVECTOR package is compatible with the LAPACK solver */
  if (tempv->ops->nvgetarraypointer == NULL ||
      tempv->ops->nvsetarraypointer == NULL) {
    cpProcessError(cp_mem, CPDIRECT_ILL_INPUT, "CPLAPACK", "CPLapackDense", MSGD_BAD_NVECTOR);
    return(CPDIRECT_ILL_INPUT);
  }

  if (lfree !=NULL) lfree(cp_mem);

  /* Set four main function fields in cp_mem */
  linit  = cpLapackDenseInit;
  lsetup = cpLapackDenseSetup;
  lsolve = cpLapackDenseSolve;
  lfree  = cpLapackDenseFree;

  /* Get memory for CPDlsMemRec */
  cpdls_mem = NULL;
  cpdls_mem = (CPDlsMem) malloc(sizeof(CPDlsMemRec));
  if (cpdls_mem == NULL) {
    cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPLAPACK", "CPLapackDense", MSGD_MEM_FAIL);
    return(CPDIRECT_MEM_FAIL);
  }

  /* Set matrix type */
  mtype = SUNDIALS_DENSE;

  /* Set default Jacobian routine and Jacobian data */
  djacE = NULL;
  djacI = NULL;
  J_data = NULL;

  last_flag = CPDIRECT_SUCCESS;
  lsetup_exists = TRUE;

  /* Set problem dimension */
  n = N;

  /* Allocate memory for M, pivot array, and (if needed) savedJ */
  M = NULL;
  pivots = NULL;
  savedJ = NULL;

  M = NewDenseMat(N, N);
  if (M == NULL) {
    cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPLAPACK", "CPLapackDense", MSGD_MEM_FAIL);
    free(cpdls_mem);
    return(CPDIRECT_MEM_FAIL);
  }
  pivots = NewIntArray(N);
  if (pivots == NULL) {
    cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPLAPACK", "CPLapackDense", MSGD_MEM_FAIL);
    DestroyMat(M);
    free(cpdls_mem);
    return(CPDIRECT_MEM_FAIL);
  }
  if (ode_type == CP_EXPL) {
    savedJ = NewDenseMat(N, N);
    if (savedJ == NULL) {
      cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPLAPACK", "CPLapackDense", MSGD_MEM_FAIL);
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
 * -----------------------------------------------------------------
 * CPLapackBand
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the band linear solver module. It first calls
 * the existing lfree routine if this is not NULL.  It then sets the
 * cp_linit, cp_lsetup, cp_lsolve, and cp_lfree fields in (*cpode_mem)
 * to be cpLapackBandInit, cpLapackBandSetup, cpLapackBandSolve, 
 * and cpLapackBandFree, respectively.  It allocates memory for a 
 * structure of type CPLapackBandMemRec and sets the cp_lmem field in 
 * (*cpode_mem) to the address of this structure.  It sets lsetup_exists 
 * in (*cpode_mem) to be TRUE, mu to be mupper, ml to be mlower, and 
 * the jacE and jacI field to NULL.
 * Finally, it allocates memory for M, pivots, and (if needed) savedJ.  
 * The CPLapackBand return value is CPDIRECT_SUCCESS = 0, 
 * CPDIRECT_MEM_FAIL = -1, or CPDIRECT_ILL_INPUT = -2.
 *
 * NOTE: The CPLAPACK linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, CPLapackBand will first 
 *       test for compatible a compatible N_Vector internal
 *       representation by checking that the function 
 *       N_VGetArrayPointer exists.
 * -----------------------------------------------------------------
 */                  
int CPLapackBand(void *cpode_mem, int N, int mupper, int mlower)
{
  CPodeMem cp_mem;
  CPDlsMem cpdls_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPDIRECT_MEM_NULL, "CPLAPACK", "CPLapackBand", MSGD_CPMEM_NULL);
    return(CPDIRECT_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  /* Test if the NVECTOR package is compatible with the BAND solver */
  if (tempv->ops->nvgetarraypointer == NULL) {
    cpProcessError(cp_mem, CPDIRECT_ILL_INPUT, "CPLAPACK", "CPLapackBand", MSGD_BAD_NVECTOR);
    return(CPDIRECT_ILL_INPUT);
  }

  if (lfree != NULL) lfree(cp_mem);

  /* Set four main function fields in cp_mem */  
  linit  = cpLapackBandInit;
  lsetup = cpLapackBandSetup;
  lsolve = cpLapackBandSolve;
  lfree  = cpLapackBandFree;
  
  /* Get memory for CPDlsMemRec */
  cpdls_mem = NULL;
  cpdls_mem = (CPDlsMem) malloc(sizeof(CPDlsMemRec));
  if (cpdls_mem == NULL) {
    cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPLAPACK", "CPLapackBand", MSGD_MEM_FAIL);
    return(CPDIRECT_MEM_FAIL);
  }

  /* Set matrix type */
  mtype = SUNDIALS_BAND;

  /* Set default Jacobian routine and Jacobian data */
  bjacE = NULL;
  bjacI = NULL;
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
    cpProcessError(cp_mem, CPDIRECT_ILL_INPUT, "CPLAPACK", "CPLapackBand", MSGD_BAD_SIZES);
    free(cpdls_mem);
    return(CPDIRECT_ILL_INPUT);
  }

  /* Set extended upper half-bandwith for M (required for pivoting) */
  smu = MIN(N-1, mu + ml);

  /* Allocate memory for M, savedJ, and pivot arrays */
  M = NULL;
  pivots = NULL;
  savedJ = NULL;

  M = NewBandMat(N, mu, ml, smu);
  if (M == NULL) {
    cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPLAPACK", "CPLapackBand", MSGD_MEM_FAIL);
    free(cpdls_mem);
    return(CPDIRECT_MEM_FAIL);
  }  
  pivots = NewIntArray(N);
  if (pivots == NULL) {
    cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPLAPACK", "CPLapackBand", MSGD_MEM_FAIL);
    DestroyMat(M);
    free(cpdls_mem);
    return(CPDIRECT_MEM_FAIL);
  }
  if (ode_type == CP_EXPL) {
    savedJ = NewBandMat(N, mu, ml, smu);
    if (savedJ == NULL) {
      cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPLAPACK", "CPLapackBand", MSGD_MEM_FAIL);
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
 * -----------------------------------------------------------------
 * CPLapackDenseProj
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the dense linear solver module for the projection.
 *
 * NOTE: The dense linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, CPLapackDenseProj will
 *       first test for compatible a compatible N_Vector internal
 *       representation by checking that N_VGetArrayPointer and
 *       N_VSetArrayPointer exist.
 * -----------------------------------------------------------------
 */
int CPLapackDenseProj(void *cpode_mem, int Nc, int Ny, int fact_type)
{
  CPodeMem cp_mem;
  CPDlsProjMem cpdlsP_mem;
  int ier, mone = -1;
  realtype tmp;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPDIRECT_MEM_NULL, "CPLAPACK", "CPLapackDenseProj", MSGD_CPMEM_NULL);
    return(CPDIRECT_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  /* Test if the NVECTOR package is compatible with the LAPACK solver */
  if (tempv->ops->nvgetarraypointer == NULL ||
      tempv->ops->nvsetarraypointer == NULL) {
    cpProcessError(cp_mem, CPDIRECT_ILL_INPUT, "CPLAPACK", "CPLapackDenseProj", MSGD_BAD_NVECTOR);
    return(CPDIRECT_ILL_INPUT);
  }

  /* Check if fact_type has a legal value */
  if ( (fact_type != CPDIRECT_LU) && 
       (fact_type != CPDIRECT_QR) && 
       (fact_type != CPDIRECT_SC) &&
       (fact_type != CPDIRECT_QRP) ) {
    cpProcessError(cp_mem, CPDIRECT_ILL_INPUT, "CPLAPACK", "CPLapackDenseProj", MSGD_BAD_FACT);
    return(CPDIRECT_ILL_INPUT);
  }

  if (lfreeP !=NULL) lfreeP(cp_mem);

  /* Set the five function fields in cp_mem */
  linitP  = cpLapackDenseProjInit;
  lsetupP = cpLapackDenseProjSetup;
  lsolveP = cpLapackDenseProjSolve;
  lmultP  = cpLapackDenseProjMult;
  lfreeP  = cpLapackDenseProjFree;

  /* Get memory for CPDlsProjMemRec */
  cpdlsP_mem = NULL;
  cpdlsP_mem = (CPDlsProjMem) malloc(sizeof(CPDlsProjMemRec));
  if (cpdlsP_mem == NULL) {
    cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPLAPACK", "CPLapackDenseProj", MSGD_MEM_FAIL);
    return(CPDIRECT_MEM_FAIL);
  }

  /* Set default Jacobian routine and Jacobian data */
  djacP = NULL;
  JP_data = NULL;

  lsetupP_exists = TRUE;

  /* Allocate memory */
  G = NULL;
  pivotsP = NULL;
  K = NULL;
  beta = NULL;
  wrk = NULL;
     
  /* Allocate memory for G */
  G = NewDenseMat(Ny, Nc);
  if (G == NULL) {
    cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPLAPACK", "CPLapackDenseProj", MSGD_MEM_FAIL);
    free(cpdlsP_mem);
    return(CPDIRECT_MEM_FAIL);
  }
  savedG = NewDenseMat(Ny, Nc);
  if (savedG == NULL) {
    cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPLAPACK", "CPLapackDenseProj", MSGD_MEM_FAIL);
    DestroyMat(G);
    free(cpdlsP_mem);
    return(CPDIRECT_MEM_FAIL);
  }

  /* Allocate additional work space, depending on factorization */
  switch(fact_type) {

  case CPDIRECT_LU:
    /* Allocate space for pivotsP and K */
    pivotsP = NewIntArray(Nc);
    if (pivotsP == NULL) {
      cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPLAPACK", "CPLapackDenseProj", MSGD_MEM_FAIL);
      DestroyMat(savedG);
      DestroyMat(G);
      free(cpdlsP_mem);
      return(CPDIRECT_MEM_FAIL);
    }
    K = NewDenseMat(Ny-Nc, Ny-Nc);
    if (K == NULL) {
      cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPLAPACK", "CPLapackDenseProj", MSGD_MEM_FAIL);
      DestroyArray(pivotsP);
      DestroyMat(savedG);
      DestroyMat(G);
      free(cpdlsP_mem);
      return(CPDIRECT_MEM_FAIL);
    }
    break;

  case CPDIRECT_QR:
    /* Allocate space for beta */
    beta = NewRealArray(Nc);
    if (beta == NULL) {
      cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPLAPACK", "CPLapackDenseProj", MSGD_MEM_FAIL);
      DestroyMat(savedG);
      DestroyMat(G);
      free(cpdlsP_mem);
      return(CPDIRECT_MEM_FAIL);      
    }
    /* Find optimal length of work array */
    dgeqrf_f77(&Ny, &Nc, G->data, &Ny, beta, &tmp, &mone, &ier);
    /* Allocate space for wrk */
    len_wrk = (int)tmp;
    wrk = NewRealArray(len_wrk);
    if (wrk == NULL) {
      cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPLAPACK", "CPLapackDenseProj", MSGD_MEM_FAIL);
      DestroyArray(beta);
      DestroyMat(savedG);
      DestroyMat(G);
      free(cpdlsP_mem);
      return(CPDIRECT_MEM_FAIL);      
    }
    /* If projecting in WRMS norm, allocate space for K=Q^T*D^(-1)*Q */
    if (pnorm == CP_PROJ_ERRNORM) {
      K = NewDenseMat(Nc, Nc);
      if (K == NULL) {
        cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPLAPACK", "CPLapackDenseProj", MSGD_MEM_FAIL);
        DestroyArray(wrk);
        DestroyArray(beta);
        DestroyMat(savedG);
        DestroyMat(G);
        free(cpdlsP_mem);
        return(CPDIRECT_MEM_FAIL);
      }
    }
    break;

  case CPDIRECT_QRP:
    /* Allocate space for pivotsP */
    pivotsP = NewIntArray(Nc);
    if (pivotsP == NULL) {
      cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPLAPACK", "CPLapackDenseProj", MSGD_MEM_FAIL);
      DestroyMat(savedG);
      DestroyMat(G);
      free(cpdlsP_mem);
      return(CPDIRECT_MEM_FAIL);
    }
    /* Allocate space for beta */
    beta = NewRealArray(Nc);
    if (beta == NULL) {
      cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPLAPACK", "CPLapackDenseProj", MSGD_MEM_FAIL);
      DestroyArray(pivotsP);
      DestroyMat(savedG);
      DestroyMat(G);
      free(cpdlsP_mem);
      return(CPDIRECT_MEM_FAIL);      
    }
    /* Find optimal length of work array */
    dgeqp3_f77(&Ny, &Nc, G->data, &Ny, pivotsP, beta, &tmp, &mone, &ier);
    /* Allocate space for wrk */
    len_wrk = (int)tmp;
    wrk = NewRealArray(len_wrk);
    if (wrk == NULL) {
      cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPLAPACK", "CPLapackDenseProj", MSGD_MEM_FAIL);
      DestroyArray(beta);
      DestroyArray(pivotsP);
      DestroyMat(savedG);
      DestroyMat(G);
      free(cpdlsP_mem);
      return(CPDIRECT_MEM_FAIL);      
    }
    /* If projecting in WRMS norm, allocate space for K=Q^T*D^(-1)*Q */
    if (pnorm == CP_PROJ_ERRNORM) {
      K = NewDenseMat(Nc, Nc);
      if (K == NULL) {
        cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPLAPACK", "CPLapackDenseProj", MSGD_MEM_FAIL);
        DestroyArray(wrk);
        DestroyArray(beta);
        DestroyArray(pivotsP);
        DestroyMat(savedG);
        DestroyMat(G);
        free(cpdlsP_mem);
        return(CPDIRECT_MEM_FAIL);
      }
    }
    break;

  case CPDIRECT_SC:
    /* Allocate space for K = G * D^(-1) * G^T */
    K = NewDenseMat(Nc, Nc);
    if (K == NULL) {
      cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPLAPACK", "CPLapackDenseProj", MSGD_MEM_FAIL);
      DestroyMat(savedG);
      DestroyMat(G);
      free(cpdlsP_mem);
      return(CPDIRECT_MEM_FAIL);
    }

    break;

  }

  /* Copy inputs into memory */
  nc    = Nc;        /* number of constraints */
  ny    = Ny;        /* number of states      */
  ftype = fact_type; /* factorization type    */

  /* Attach linear solver memory to integrator memory */
  lmemP = cpdlsP_mem;

  return(CPDIRECT_SUCCESS);
}

/* 
 * =================================================================
 *  PRIVATE FUNCTIONS FOR IMPLICIT INTEGRATION WITH DENSE JACOBIANS
 * =================================================================
 */

/*
 * cpLapackDenseInit does remaining initializations specific to the dense
 * linear solver.
 */
static int cpLapackDenseInit(CPodeMem cp_mem)
{
  CPDlsMem cpdls_mem;

  cpdls_mem = (CPDlsMem) lmem;
  
  nje   = 0;
  nfeDQ = 0;
  nstlj = 0;
  
  if (ode_type == CP_EXPL && djacE == NULL) {
    djacE = cpDlsDenseDQJacExpl;
    J_data = cp_mem;
  } 
  
  if (ode_type == CP_IMPL && djacI == NULL) {
    djacI = cpDlsDenseDQJacImpl;
    J_data = cp_mem;
  }

  last_flag = CPDIRECT_SUCCESS;
  return(0);
}

/*
 * cpLapackDenseSetup does the setup operations for the dense linear solver.
 * It makes a decision whether or not to call the Jacobian evaluation
 * routine based on various state variables, and if not it uses the 
 * saved copy (for explicit ODE only). In any case, it constructs 
 * the Newton matrix M = I - gamma*J or M = F_y' - gamma*F_y, updates 
 * counters, and calls the dense LU factorization routine.
 */
static int cpLapackDenseSetup(CPodeMem cp_mem, int convfail,
                              N_Vector yP, N_Vector ypP, N_Vector fctP,
                              booleantype *jcurPtr,
                              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  CPDlsMem cpdls_mem;
  realtype dgamma, fact;
  booleantype jbad, jok;
  int ier, retval, one = 1;

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
      dcopy_f77(&(savedJ->ldata), savedJ->data, &one, M->data, &one);

    } else {

      /* If jok = FALSE, call jac routine for new J value */
      nje++;
      nstlj = nst;
      *jcurPtr = TRUE;

      retval = djacE(n, tn, yP, fctP, M, J_data, tmp1, tmp2, tmp3);
      if (retval == 0) {
        dcopy_f77(&(M->ldata), M->data, &one, savedJ->data, &one);
      } else if (retval < 0) {
        cpProcessError(cp_mem, CPDIRECT_JACFUNC_UNRECVR, "CPLAPACK", "cpLapackDenseSetup", MSGD_JACFUNC_FAILED);
        last_flag = CPDIRECT_JACFUNC_UNRECVR;
        return(-1);
      } else if (retval > 0) {
        last_flag = CPDIRECT_JACFUNC_RECVR;
        return(1);
      }

    }

    /* Scale J by - gamma */
    fact = -gamma;
    dscal_f77(&(M->ldata), &fact, M->data, &one);

    /* Add identity to get M = I - gamma*J*/
    LapackDenseAddI(M);

    break;

  case CP_IMPL:

    /* Call Jacobian function */
    nje++;
    retval = djacI(n, tn, gamma, yP, ypP, fctP, M, J_data, tmp1, tmp2, tmp3);
    if (retval == 0) {
      break;
    } else if (retval < 0) {
      cpProcessError(cp_mem, CPDIRECT_JACFUNC_UNRECVR, "CPLAPACK", "cpLapackDenseSetup", MSGD_JACFUNC_FAILED);
      last_flag = CPDIRECT_JACFUNC_UNRECVR;
      return(-1);
    } else if (retval > 0) {
      last_flag = CPDIRECT_JACFUNC_RECVR;
      return(1);
    }
  
    break;

  }

  /* Do LU factorization of M */
  dgetrf_f77(&n, &n, M->data, &(M->ldim), pivots, &ier);

  /* Return 0 if the LU was complete; otherwise return 1 */
  last_flag = ier;
  if (ier > 0) return(1);
  return(0);
}

/*
 * cpLapackDenseSolve handles the solve operation for the dense linear solver
 * by calling the dense backsolve routine.
 */
static int cpLapackDenseSolve(CPodeMem cp_mem, N_Vector b, N_Vector weight,
                              N_Vector yC, N_Vector ypC, N_Vector fctC)
{
  CPDlsMem cpdls_mem;
  realtype *bd, fact;
  int ier, one = 1;

  cpdls_mem = (CPDlsMem) lmem;
  
  bd = N_VGetArrayPointer(b);

  dgetrs_f77("N", &n, &one, M->data, &(M->ldim), pivots, bd, &n, &ier, 1); 
  if (ier > 0) return(1);

  /* For BDF, scale the correction to account for change in gamma */
  if ((lmm_type == CP_BDF) && (gamrat != ONE)) {
    fact = TWO/(ONE + gamrat);
    dscal_f77(&n, &fact, bd, &one); 
  }
  
  last_flag = CPDIRECT_SUCCESS;
  return(0);
}

/*
 * cpLapackDenseFree frees memory specific to the dense linear solver.
 */
static void cpLapackDenseFree(CPodeMem cp_mem)
{
  CPDlsMem  cpdls_mem;

  cpdls_mem = (CPDlsMem) lmem;
  
  DestroyMat(M);
  DestroyArray(pivots);
  if (ode_type == CP_EXPL) DestroyMat(savedJ);
  free(cpdls_mem); 
  cpdls_mem = NULL;
}

/* 
 * =================================================================
 *  PRIVATE FUNCTIONS FOR IMPLICIT INTEGRATION WITH BAND JACOBIANS
 * =================================================================
 */

/*
 * cpLapackBandInit does remaining initializations specific to the band
 * linear solver.
 */
static int cpLapackBandInit(CPodeMem cp_mem)
{
  CPDlsMem cpdls_mem;

  cpdls_mem = (CPDlsMem) lmem;

  nje   = 0;
  nfeDQ = 0;
  nstlj = 0;

  if (ode_type == CP_EXPL && bjacE == NULL) {
    bjacE = cpDlsBandDQJacExpl;
    J_data = cp_mem;
  } 
  
  if (ode_type == CP_IMPL && bjacI == NULL) {
    bjacI = cpDlsBandDQJacImpl;
    J_data = cp_mem;
  }

  last_flag = CPDIRECT_SUCCESS;
  return(0);
}

/*
 * cpLapackBandSetup does the setup operations for the band linear solver.
 * It makes a decision whether or not to call the Jacobian evaluation
 * routine based on various state variables, and if not it uses the 
 * saved copy (for explicit ODE only). In any case, it constructs 
 * the Newton matrix M = I - gamma*J or M = F_y' - gamma*F_y, updates 
 * counters, and calls the band LU factorization routine.
 */
static int cpLapackBandSetup(CPodeMem cp_mem, int convfail, 
                             N_Vector yP, N_Vector ypP, N_Vector fctP, 
                             booleantype *jcurPtr,
                             N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  CPDlsMem cpdls_mem;
  realtype dgamma, fact;
  booleantype jbad, jok;
  int ier, retval, one = 1;

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
      dcopy_f77(&(savedJ->ldata), savedJ->data, &one, M->data, &one);
      
    } else {
      
      /* If jok = FALSE, call jac routine for new J value */
      nje++;
      nstlj = nst;
      *jcurPtr = TRUE;

      retval = bjacE(n, mu, ml, tn, yP, fctP, M, J_data, tmp1, tmp2, tmp3);
      if (retval == 0) {
        dcopy_f77(&(M->ldata), M->data, &one, savedJ->data, &one);
      } else if (retval < 0) {
        cpProcessError(cp_mem, CPDIRECT_JACFUNC_UNRECVR, "CPLAPACK", "cpLapackBandSetup", MSGD_JACFUNC_FAILED);
        last_flag = CPDIRECT_JACFUNC_UNRECVR;
        return(-1);
      } else if (retval > 0) {
        last_flag = CPDIRECT_JACFUNC_RECVR;
        return(1);
      }

    }
  
    /* Scale J by - gamma */
    fact = -gamma;
    dscal_f77(&(M->ldata), &fact, M->data, &one);

    /* Add identity to get M = I - gamma*J*/
    LapackBandAddI(M);

    break;

  case CP_IMPL:

    /* Call Jacobian function */
    nje++;
    retval = bjacI(n, mu, ml, tn, gamma, yP, ypP, fctP, M, J_data, tmp1, tmp2, tmp3);
    if (retval == 0) {
      break;
    } else if (retval < 0) {
      cpProcessError(cp_mem, CPDIRECT_JACFUNC_UNRECVR, "CPLAPACK", "cpLapackBandSetup", MSGD_JACFUNC_FAILED);
      last_flag = CPDIRECT_JACFUNC_UNRECVR;
      return(-1);
    } else if (retval > 0) {
      last_flag = CPDIRECT_JACFUNC_RECVR;
      return(+1);
    }

    break;

  }

  /* Do LU factorization of M */
  dgbtrf_f77(&n, &n, &ml, &mu, M->data, &(M->ldim), pivots, &ier);

  /* Return 0 if the LU was complete; otherwise return 1 */
  last_flag = ier;
  if (ier > 0) return(1);
  return(0);

}

/*
 * cpLapackBandSolve handles the solve operation for the band linear solver
 * by calling the band backsolve routine.
 */
static int cpLapackBandSolve(CPodeMem cp_mem, N_Vector b, N_Vector weight,
                             N_Vector yC, N_Vector ypC, N_Vector fctC)
{
  CPDlsMem cpdls_mem;
  realtype *bd, fact;
  int ier, one = 1;

  cpdls_mem = (CPDlsMem) lmem;

  bd = N_VGetArrayPointer(b);

  dgbtrs_f77("N", &n, &ml, &mu, &one, M->data, &(M->ldim), pivots, bd, &n, &ier, 1);
  if (ier > 0) return(1);

  /* For BDF, scale the correction to account for change in gamma */
  if ((lmm_type == CP_BDF) && (gamrat != ONE)) {
    fact = TWO/(ONE + gamrat);
    dscal_f77(&n, &fact, bd, &one); 
  }

  last_flag = CPDIRECT_SUCCESS;
  return(0);
}

/*
 * cpLapackBandFree frees memory specific to the band linear solver.
 */
static void cpLapackBandFree(CPodeMem cp_mem)
{
  CPDlsMem  cpdls_mem;

  cpdls_mem = (CPDlsMem) lmem;
  
  DestroyMat(M);
  DestroyArray(pivots);
  if (ode_type == CP_EXPL) DestroyMat(savedJ);
  free(cpdls_mem); 
  cpdls_mem = NULL;
}

/* 
 * =================================================================
 *  PRIVATE FUNCTIONS FOR PROJECTION WITH DENSE JACOBIANS
 * =================================================================
 */

/*
 * cpLapackDenseProjInit does remaining initializations specific to
 * the dense linear solver.
 */
static int cpLapackDenseProjInit(CPodeMem cp_mem)
{
  CPDlsProjMem cpdlsP_mem;

  cpdlsP_mem = (CPDlsProjMem) lmemP;
  
  njeP   = 0;
  nceDQ  = 0;
  nstljP = 0;
  
  if (djacP == NULL) {
    djacP = cpDlsDenseProjDQJac;
    JP_data = cp_mem;
  }  

  return(0);
}

/*
 * cpLapackDenseProjSetup does the setup operations for the dense 
 * linear solver.
 * It calls the Jacobian evaluation routine and, depending on ftype,
 * it performs various factorizations.
 */
static int cpLapackDenseProjSetup(CPodeMem cp_mem, N_Vector y, N_Vector cy,
                                  N_Vector c_tmp1, N_Vector c_tmp2, N_Vector s_tmp1)
{
  int ier;
  CPDlsProjMem cpdlsP_mem;
  realtype *col_i, rim1, ri;
  int i, j, nd, one = 1;
  int retval;

  realtype coef_1 = ONE, coef_0 = ZERO;

  cpdlsP_mem = (CPDlsProjMem) lmemP;

  nd = ny-nc;

  /* Call Jacobian function (G will contain the Jacobian transposed) */
  retval = djacP(nc, ny, tn, y, cy, G, JP_data, c_tmp1, c_tmp2);
  njeP++;
  if (retval < 0) {
    cpProcessError(cp_mem, CPDIRECT_JACFUNC_UNRECVR, "CPLAPACK", "cpLapackDenseProjSetup", MSGD_JACFUNC_FAILED);
    return(-1);
  } else if (retval > 0) {
    return(1);
  }

  /* Save Jacobian before factorization for possible use by lmultP */
  dcopy_f77(&(G->ldata), G->data, &one, savedG->data, &one);

  /* Factorize G, depending on ftype */
  switch (ftype) {

  case CPDIRECT_LU:

    /* 
     * LU factorization of G^T with partial pivoting
     *    P*G^T = | U1^T | * L^T
     *            | U2^T |
     * After factorization, P is encoded in pivotsP and
     * G^T is overwritten with U1 (nc by nc unit upper triangular), 
     * U2 ( nc by ny-nc rectangular), and L (nc by nc lower triangular).
     * Return ier if factorization failed. 
     */
    dgetrf_f77(&ny, &nc, G->data, &ny, pivotsP, &ier);
    if (ier != 0) return(ier);

    /* 
     * Build S = U1^{-1} * U2 (in place, S overwrites U2) 
     * For each row j of G, j = nc,...,ny-1, perform
     * a backward substitution (row version).
     * After this step, G^T contains U1, S, and L.
     */
    for (j=nc; j<ny; j++)
      dtrsv_f77("L", "T", "U", &nc, G->data, &ny, (G->data + j), &ny, 1, 1, 1);

    /*   
     * Build K = D1 + S^T * D2 * S 
     * S^T is stored in g_mat[nc...ny-1][0...nc]
     * Compute and store only the lower triangular part of K.
     */
    if (pnorm == CP_PROJ_L2NORM) {
      dsyrk_f77("L", "N", &nd, &nc, &coef_1, (G->data + nc), &ny, &coef_0, K->data, &nd, 1, 1);
      LapackDenseAddI(K);
    } else {
      cplLUcomputeKD(cp_mem, s_tmp1);
    }

    /*
     * Perform Cholesky decomposition of K: K = C*C^T
     * After factorization, the lower triangular part of K contains C.
     * Return ier if factorization failed. 
     */
    dpotrf_f77("L", &nd, K->data, &nd, &ier, 1);
    if (ier != 0) return(ier);

    break;

  case CPDIRECT_QR:

    /* 
     * QR factorization of G^T: G^T = Q*R
     * After factorization, the upper triangular part of G^T 
     * contains the matrix R. The lower trapezoidal part of
     * G^T, together with the array beta, encodes the orthonormal
     * columns of Q as elementary reflectors.
     */
    dgeqrf_f77(&ny, &nc, G->data, &ny, beta, wrk, &len_wrk, &ier);
    if (ier != 0) return(ier);

    /* If projecting in WRMS norm */
    if (pnorm == CP_PROJ_ERRNORM) {
      /* Build K = Q^T * D^(-1) * Q */
      cplQRcomputeKD(cp_mem, s_tmp1);
      /* Perform Cholesky decomposition of K */
      dpotrf_f77("L", &nc, K->data, &nc, &ier, 1);
      if (ier != 0) return(ier);
    }

    break;

  case CPDIRECT_QRP:

    /* 
     * QR with pivoting factorization of G^T: G^T * P = Q * R.
     * After factorization, the upper triangular part of G^T 
     * contains the matrix R. The lower trapezoidal part of
     * G^T, together with the array beta, encodes the orthonormal
     * columns of Q as elementary reflectors.
     * The pivots are stored in 'pivotsP'.
     */
    for (i=0; i<nc; i++) pivotsP[i] = 0;
    dgeqp3_f77(&ny, &nc, G->data, &ny, pivotsP, beta, wrk, &len_wrk, &ier);
    if (ier != 0) return(ier);

    /*
     * Determine the number of independent constraints.
     * After the QR factorization, the diagonal elements of R should 
     * be in decreasing order of their absolute values.
     */
    rim1 = ABS(G->data[0]);
    for (i=1, nr=1; i<nc; i++, nr++) {
      col_i = G->cols[i];
      ri = ABS(col_i[i]);
      if (ri < 100*uround) break;
      if (ri/rim1 < RPowerR(uround, THREE/FOUR)) break;
    }

    /* If projecting in WRMS norm */
    if (pnorm == CP_PROJ_ERRNORM) {
      /* Build K = Q^T * D^(-1) * Q */
      cplQRcomputeKD(cp_mem, s_tmp1);
      /* Perform Cholesky decomposition of K */
      dpotrf_f77("L", &nc, K->data, &nc, &ier, 1);
      if (ier != 0) return(ier);
    }

    break;

  case CPDIRECT_SC:

    /* 
     * Build K = G*D^(-1)*G^T
     * G^T is stored in g_mat[0...ny-1][0...nc]
     * Compute and store only the lower triangular part of K.
     */
    if (pnorm == CP_PROJ_L2NORM) {
      dsyrk_f77("L", "T", &nc, &ny, &coef_1, G->data, &ny, &coef_0, K->data, &nc, 1, 1);
    } else {
      cplSCcomputeKD(cp_mem, s_tmp1);
    }

    /* 
     * Perform Cholesky decomposition of K: K = C*C^T
     * After factorization, the lower triangular part of K contains C.
     * Return 1 if factorization failed. 
     */
    dpotrf_f77("L", &nc, K->data, &nc, &ier, 1);
    if (ier != 0) return(ier);

    break;

  }

  return(0);
}

/*
 * cpLapackDenseProjSolve handles the solve operation for the dense linear 
 * linear solver.
 */
static int cpLapackDenseProjSolve(CPodeMem cp_mem, N_Vector b, N_Vector x,
                                  N_Vector y, N_Vector cy,
                                  N_Vector c_tmp1, N_Vector s_tmp1)
{
  CPDlsProjMem cpdlsP_mem;
  realtype *bd, *xd;
  realtype  *ewt_data, *d_data, *da_data, tmp;
  int nd, i, j, k, pk, ier, one = 1;
  realtype coef_1 = ONE, coef_0 = ZERO, coef_m1 = -ONE;

  cpdlsP_mem = (CPDlsProjMem) lmemP;

  ewt_data = N_VGetArrayPointer(ewt);
  bd = N_VGetArrayPointer(b);
  xd = N_VGetArrayPointer(x);
  d_data = N_VGetArrayPointer(s_tmp1);

  nd = ny - nc;

  /* Solve the linear system, depending on ftype */
  switch (ftype) {

  case CPDIRECT_LU:

    /* Solve L*U1*alpha = bd
     *   (a) solve L*beta = bd using fwd. subst.
     *   (b) solve U1*alpha = beta using bckwd. subst
     * where L^T and U1^T are stored in G[0...nc-1][0...nc-1].
     * beta and then alpha overwrite bd.
     */
    dtrsv_f77("U", "T", "N", &nc, G->data, &ny, bd, &one, 1, 1, 1);
    dtrsv_f77("L", "T", "U", &nc, G->data, &ny, bd, &one, 1, 1, 1);

    /* Compute S^T * (D1 * alpha)
     * alpha is stored in bd.
     * S^T is stored in g_mat[nc...ny-1][0...nc-1].
     * Store result in x2 = x[nc...ny-1].
     */
    if (pnorm == CP_PROJ_ERRNORM) {
      
      /* Load squared error weights into d */
      for (k=0; k<ny; k++) d_data[k] = ewt_data[k] * ewt_data[k];
      /* Permute elements of d, based on pivotsP. 
       * Note that pivot information from dgetrf is 1-based.
       * Swap d[k] and d[pivotsP[k]]. */
      for (k=0; k<nc; k++) {
        pk = pivotsP[k];
        if (pk != k) {
          tmp = d_data[k];
          d_data[k]  = d_data[pk];
          d_data[pk] = tmp;
        }
      }
      /* Compute D1*alpha and store it into da_data */
      da_data = N_VGetArrayPointer(c_tmp1);
      for(k=0; k<nc; k++) da_data[k] = d_data[k] * bd[k];
      /* Compute S^T * D1 * alpha = S^T * da */
      dgemv_f77("N", &ny, &nc, &coef_1, (G->data + nc), &ny, da_data, &one, &coef_0, (xd + nc), &one, 1);
      
    } else {
      
      /* Compute S^T * alpha */
      dgemv_f77("N", &nd, &nc, &coef_1, (G->data + nc), &ny, bd, &one, &coef_0, (xd + nc), &one, 1);

    }

    /* Solve K*x2 = S^T*D1*alpha, using the Cholesky decomposition available in K.
     * S^T*D1*alpha is stored in x2 = x[nc...ny-1].
     */
    dpotrs_f77("L", &nd, &one, K->data, &nd, (xd + nc), &nd, &ier, 1);
    if (ier != 0) return(ier);

    /* Compute x1 = alpha - S*x2 
     * alpha is stored in bd.
     * x2 is stored in x[nc...ny-1].
     * S^T is stored in g_mat[nc...ny-1][0...nc-1].
     * Store result in x1 = x[0...nc-1].
     */
    dcopy_f77(&nc, bd, &one, xd, &one);
    dgemv_f77("T", &nd, &nc, &coef_m1, (G->data + nc), &ny, (xd + nc), &one, &coef_1, xd, &one, 1);

    /* Compute P^T * x, where P is encoded into pivotsP.
     * Note that pivot information from dgetrf is 1-based.
     * Store result in x.
     */
    for (k=nc-1; k>=0; k--) {
      pk = pivotsP[k]-1;
      if(pk != k) {
        tmp = xd[k];
        xd[k] = xd[pk];
        xd[pk] = tmp;
      }
    }

    break;

  case CPDIRECT_QR:

    /* 
     * Solve R^T*alpha = bd using fwd. subst. (row version)
     * The upper triangular matrix R is stored in g_mat[0...nc-1][0...nc-1]
     * alpha overwrites bd.
     */
    dtrsv_f77("U", "T", "N", &nc, G->data, &ny, bd, &one, 1, 1, 1);

    /* If projecting in WRMS norm, solve K*beta = alpha */
    if (pnorm == CP_PROJ_ERRNORM) {
      dpotrs_f77("L", &nc, &one, K->data, &nc, bd, &nc, &ier, 1);
      if (ier != 0) return(ier);
    }

    /* 
     * Compute x = Q1*alpha 
     *
     * Since we cannot really use the "thin" QR decomposition, we
     * first need to initialize xd = [alpha; 0].
     */
    for (k=0; k<nc; k++)  xd[k] = bd[k];
    for (k=nc; k<ny; k++) xd[k] = ZERO;
    dormqr_f77("L", "N", &ny, &one, &nc, G->data, &ny, beta, xd, &ny, wrk, &len_wrk, &ier, 1, 1);
    if (ier != 0) return(ier);

    /* If projecting in WRMS norm, scale x by D^(-1) */
    if (pnorm == CP_PROJ_ERRNORM) {
      for (i=0; i<ny; i++)
        xd[i] /= ewt_data[i]*ewt_data[i];
    }

    break;    


  case CPDIRECT_QRP:

    /* Compute P^T * b, where P is encoded into pivotsP.
     * If pivotsP[j] = k, then, for the factorization, the j-th column of G^T*P 
     * was the k-th column of G^T.
     * Therefore, to compute P^T*b, we must move the k-th element of b to the
     * j-th position, for j=1,2,... This is a forward permutation.
     * Note that pivot information from dgeqp3 is 1-based.
     * Store result in b.
     */
    for (i=1; i<=nc; i++) pivotsP[i-1] = -pivotsP[i-1];

    for (i=1; i<=nc; i++) {
      
      if (pivotsP[i-1] > 0) continue;

      j = i;
      pivotsP[j-1] = -pivotsP[j-1];
      pk = pivotsP[j-1];

      while (pivotsP[pk-1] < 0) {

        tmp = bd[j-1];
        bd[j-1] = bd[pk-1];
        bd[pk-1] = tmp;

        pivotsP[pk-1] = -pivotsP[pk-1];
        j = pk;
        pk = pivotsP[pk-1];
      }
      
    }

    /* 
     * Solve R11^T * alpha = P^T * bd using fwd. subst. (row version)
     * The upper triangular matrix R is stored in g_mat[0...nr-1][0...nr-1]
     * P^T * bd is available in bd.
     * We only consider the first nr components in P^T*bd.
     * alpha overwrites bd.
     */
    dtrsv_f77("U", "T", "N", &nr, G->data, &ny, bd, &one, 1, 1, 1);

    /* If projecting in WRMS norm, solve K*beta = alpha */
    if (pnorm == CP_PROJ_ERRNORM) {
      dpotrs_f77("L", &nc, &one, K->data, &nc, bd, &nc, &ier, 1);
      if (ier != 0) return(ier);
    }

    /* 
     * Compute x = Q1*alpha 
     *
     * Since we cannot really use the "thin" QR decomposition, we
     * first need to initialize xd = [alpha; 0].
     */
    for (k=0; k<nr; k++)  xd[k] = bd[k];
    for (k=nr; k<ny; k++) xd[k] = ZERO;
    dormqr_f77("L", "N", &ny, &one, &nc, G->data, &ny, beta, xd, &ny, wrk, &len_wrk, &ier, 1, 1);
    if (ier != 0) return(ier);

    /* If projecting in WRMS norm, scale x by D^(-1) */
    if (pnorm == CP_PROJ_ERRNORM) {
      for (i=0; i<ny; i++)
        xd[i] /= ewt_data[i]*ewt_data[i];
    }

    break;    


  case CPDIRECT_SC:

    /* 
     * Solve K*xi = bd, using the Cholesky decomposition available in K.
     * xi overwrites bd.
     */
    dpotrs_f77("L", &nc, &one, K->data, &nc, bd, &nc, &ier, 1);
    if (ier != 0) return(ier);

    /* Compute x = G^T * xi
     * G^T is stored in g_mat[0...ny-1][0...nc-1]
     * xi is available in bd.
     */
    dgemv_f77("N", &ny, &nc, &coef_1, G->data, &ny, bd, &one, &coef_0, xd, &one, 1);

    /* If projecting in WRMS norm, scale x by D^(-1) */
    if (pnorm == CP_PROJ_ERRNORM) {
      for (i=0; i<ny; i++)
        xd[i] /= ewt_data[i]*ewt_data[i];
    }

    break;

  }

  return(0);
}

/*
 * cpDenseProjMult computes the Jacobian-vector product used a saved 
 * Jacobian copy.
 */
static void cpLapackDenseProjMult(CPodeMem cp_mem, N_Vector x, N_Vector Gx)
{
  CPDlsProjMem cpdlsP_mem;
  realtype coef_1 = ONE, coef_0 = ZERO;
  realtype *xd, *Gxd;
  int one = 1;

  cpdlsP_mem = (CPDlsProjMem) lmemP;

  xd = N_VGetArrayPointer(x);
  Gxd = N_VGetArrayPointer(Gx);

  dgemv_f77("T", &ny, &nc, &coef_1, savedG->data, &ny, xd, &one, &coef_0, Gxd, &one, 1);
}

/*
 * cpLapackDenseProjFree frees memory specific to the dense linear solver.
 */
static void cpLapackDenseProjFree(CPodeMem cp_mem)
{
  CPDlsProjMem cpdlsP_mem;

  cpdlsP_mem = (CPDlsProjMem) lmemP;
  
  DestroyMat(G);
  DestroyMat(savedG);
  switch (ftype) {
  case CPDIRECT_LU:
    DestroyArray(pivotsP);
    DestroyMat(K);
    break;
  case CPDIRECT_QR:
    DestroyArray(wrk);
    DestroyArray(beta);
    if (pnorm == CP_PROJ_ERRNORM) DestroyMat(K);
    break;
  case CPDIRECT_QRP:
    DestroyArray(wrk);
    DestroyArray(beta);
    DestroyArray(pivotsP);
    if (pnorm == CP_PROJ_ERRNORM) DestroyMat(K);
    break;
  case CPDIRECT_SC:
    DestroyMat(K);
    break;
  }

  free(cpdlsP_mem); 
  cpdlsP_mem = NULL;
}

/*
 * -----------------------------------------------------------------
 * Private functions for LU-, QR-, and SC-based projection
 * -----------------------------------------------------------------
 */

/*
 * Compute the lower triangle of K = D1 + S^T*D2*S,
 * D = diag(D1, D2) = P*W*P^T, W is a diagonal matrix
 * containing the squared error weights, and P is the 
 * permutation matrix encoded into pivotsP.
 * D1 has length nc and D2 has length (ny-nc).
 */
static void cplLUcomputeKD(CPodeMem cp_mem, N_Vector d)
{
  /* RADU:: implement this ... */
}

static void cplQRcomputeKD(CPodeMem cp_mem, N_Vector d)
{
  /* RADU:: implement this ... */
}

static void cplSCcomputeKD(CPodeMem cp_mem, N_Vector d)
{
  /* RADU:: implement this ... */
}

