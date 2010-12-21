/*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2007/10/26 21:48:38 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for the CPDENSE linear solver.
 * -----------------------------------------------------------------
 */


/*
 * TODO:
 *
 *   cpdQRcomputeKD
 *   cpdSCcomputeKD
 */


/* 
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include <cpodes/cpodes_dense.h>
#include "cpodes_direct_impl.h"
#include "cpodes_private.h"

#include <sundials/sundials_math.h>

/* 
 * =================================================================
 * PROTOTYPES FOR PRIVATE FUNCTIONS
 * =================================================================
 */

/* CPDENSE linit, lsetup, lsolve, and lfree routines */
 
static int cpDenseInit(CPodeMem cp_mem);
static int cpDenseSetup(CPodeMem cp_mem, int convfail, 
                        N_Vector yP, N_Vector ypP, N_Vector fctP, 
                        booleantype *jcurPtr,
                        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int cpDenseSolve(CPodeMem cp_mem, N_Vector b, N_Vector weight,
                        N_Vector yC, N_Vector ypC, N_Vector fctC);
static void cpDenseFree(CPodeMem cp_mem);

/* CPDENSE linitP, lsetupP, lsolveP, lmultP, and lfreeP routines */
static int cpDenseProjInit(CPodeMem cp_mem);
static int cpDenseProjSetup(CPodeMem cp_mem, N_Vector y, N_Vector cy,
                            N_Vector c_tmp1, N_Vector c_tmp2, N_Vector s_tmp1);
static int cpDenseProjSolve(CPodeMem cp_mem, N_Vector b, N_Vector x,
                            N_Vector y, N_Vector cy,
                            N_Vector c_tmp1, N_Vector s_tmp1);
static void cpDenseProjMult(CPodeMem cp_mem, N_Vector x, N_Vector Gx);
static void cpDenseProjFree(CPodeMem cp_mem);

/* Private functions for LU, QR, and SC projection */
static void cpdLUcomputeKI(CPodeMem cp_mem);
static void cpdLUcomputeKD(CPodeMem cp_mem, N_Vector d);
static void cpdQRcomputeKD(CPodeMem cp_mem, N_Vector d);
static void cpdSCcomputeKI(CPodeMem cp_mem);
static void cpdSCcomputeKD(CPodeMem cp_mem, N_Vector d);

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
#define jacE           (cpdls_mem->d_djacE)
#define jacI           (cpdls_mem->d_djacI)
#define M              (cpdls_mem->d_M)
#define savedJ         (cpdls_mem->d_savedJ)
#define pivots         (cpdls_mem->d_pivots)
#define nstlj          (cpdls_mem->d_nstlj)
#define nje            (cpdls_mem->d_nje)
#define nfeDQ          (cpdls_mem->d_nfeDQ)
#define J_data         (cpdls_mem->d_J_data)
#define last_flag      (cpdls_mem->d_last_flag)

#define nc             (cpdlsP_mem->d_nc)
#define ny             (cpdlsP_mem->d_ny)
#define jacP           (cpdlsP_mem->d_jacP)
#define JP_data        (cpdlsP_mem->d_JP_data)
#define ftype          (cpdlsP_mem->d_ftype)
#define G              (cpdlsP_mem->d_G)
#define savedG         (cpdlsP_mem->d_savedG)
#define K              (cpdlsP_mem->d_K)
#define pivotsP        (cpdlsP_mem->d_pivotsP)
#define beta           (cpdlsP_mem->d_beta)
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
 * CPDense
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the dense linear solver module.  CPDense first
 * calls the existing lfree routine if this is not NULL.  Then it sets
 * the cp_linit, cp_lsetup, cp_lsolve, cp_lfree fields in (*cpode_mem)
 * to be cpDenseInit, cpDenseSetup, cpDenseSolve, and cpDenseFree,
 * respectively.  It allocates memory for a structure of type
 * CPDlsMemRec and sets the cp_lmem field in (*cpode_mem) to the
 * address of this structure.  It sets lsetup_exists in (*cpode_mem) to
 * TRUE, and the d_jac field to the default cpDlsDenseDQJac.
 * Finally, it allocates memory for M, pivots, and (if needed) savedJ.
 * The return value is SUCCESS = 0, or LMEM_FAIL = -1.
 *
 * NOTE: The dense linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, CPDense will first 
 *       test for compatible a compatible N_Vector internal
 *       representation by checking that N_VGetArrayPointer and
 *       N_VSetArrayPointer exist.
 * -----------------------------------------------------------------
 */

int CPDense(void *cpode_mem, int N)
{
  CPodeMem cp_mem;
  CPDlsMem cpdls_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPDIRECT_MEM_NULL, "CPDENSE", "CPDense", MSGD_CPMEM_NULL);
    return(CPDIRECT_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  /* Test if the NVECTOR package is compatible with the DENSE solver */
  if (tempv->ops->nvgetarraypointer == NULL ||
      tempv->ops->nvsetarraypointer == NULL) {
    cpProcessError(cp_mem, CPDIRECT_ILL_INPUT, "CPDENSE", "CPDense", MSGD_BAD_NVECTOR);
    return(CPDIRECT_ILL_INPUT);
  }

  if (lfree !=NULL) lfree(cp_mem);

  /* Set four main function fields in cp_mem */
  linit  = cpDenseInit;
  lsetup = cpDenseSetup;
  lsolve = cpDenseSolve;
  lfree  = cpDenseFree;

  /* Get memory for CPDlsMemRec */
  cpdls_mem = NULL;
  cpdls_mem = (CPDlsMem) malloc(sizeof(CPDlsMemRec));
  if (cpdls_mem == NULL) {
    cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPDENSE", "CPDense", MSGD_MEM_FAIL);
    return(CPDIRECT_MEM_FAIL);
  }

  /* Set matrix type */
  mtype = SUNDIALS_DENSE;

  /* Set default Jacobian routine and Jacobian data */
  jacE = NULL;
  jacI = NULL;
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
    cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPDENSE", "CPDense", MSGD_MEM_FAIL);
    free(cpdls_mem);
    return(CPDIRECT_MEM_FAIL);
  }
  pivots = NewIntArray(N);
  if (pivots == NULL) {
    cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPDENSE", "CPDense", MSGD_MEM_FAIL);
    DestroyMat(M);
    free(cpdls_mem);
    return(CPDIRECT_MEM_FAIL);
  }
  if (ode_type == CP_EXPL) {
    savedJ = NewDenseMat(N, N);
    if (savedJ == NULL) {
      cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPDENSE", "CPDense", MSGD_MEM_FAIL);
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
 * CPDenseProj
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the dense linear solver module for the projection.
 *
 * NOTE: The dense linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, CPDense will first 
 *       test for compatible a compatible N_Vector internal
 *       representation by checking that N_VGetArrayPointer and
 *       N_VSetArrayPointer exist.
 * -----------------------------------------------------------------
 */

int CPDenseProj(void *cpode_mem, int Nc, int Ny, int fact_type)
{
  CPodeMem cp_mem;
  CPDlsProjMem cpdlsP_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPDIRECT_MEM_NULL, "CPDENSE", "CPDenseProj", MSGD_CPMEM_NULL);
    return(CPDIRECT_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  /* Test if the NVECTOR package is compatible with the DENSE solver */
  if (tempv->ops->nvgetarraypointer == NULL ||
      tempv->ops->nvsetarraypointer == NULL) {
    cpProcessError(cp_mem, CPDIRECT_ILL_INPUT, "CPDENSE", "CPDenseProj", MSGD_BAD_NVECTOR);
    return(CPDIRECT_ILL_INPUT);
  }

  /* Check if fact_type has a legal value */
  if ( (fact_type != CPDIRECT_LU) && (fact_type != CPDIRECT_QR) && (fact_type != CPDIRECT_SC) ) {
    cpProcessError(cp_mem, CPDIRECT_ILL_INPUT, "CPDENSE", "CPDenseProj", MSGD_BAD_FACT);
    return(CPDIRECT_ILL_INPUT);
  }

  if (lfreeP !=NULL) lfreeP(cp_mem);

  /* Set the five function fields in cp_mem */
  linitP  = cpDenseProjInit;
  lsetupP = cpDenseProjSetup;
  lsolveP = cpDenseProjSolve;
  lmultP  = cpDenseProjMult;
  lfreeP  = cpDenseProjFree;

  /* Get memory for CPDlsProjMemRec */
  cpdlsP_mem = NULL;
  cpdlsP_mem = (CPDlsProjMem) malloc(sizeof(CPDlsProjMemRec));
  if (cpdlsP_mem == NULL) {
    cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPDENSE", "CPDenseProj", MSGD_MEM_FAIL);
    return(CPDIRECT_MEM_FAIL);
  }

  lsetupP_exists = TRUE;

  /* Initialize all internal pointers to NULL */
  G = NULL;
  K = NULL;
  pivotsP = NULL;
  beta = NULL;

  /* Allocate memory for G and other work space */
  G = NewDenseMat(Ny, Nc);
  if (G == NULL) {
    cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPDENSE", "CPDenseProj", MSGD_MEM_FAIL);
    free(cpdlsP_mem);
    return(CPDIRECT_MEM_FAIL);
  }
  savedG = NewDenseMat(Ny, Nc);
  if (savedG == NULL) {
    cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPDENSE", "CPDenseProj", MSGD_MEM_FAIL);
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
      cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPDENSE", "CPDenseProj", MSGD_MEM_FAIL);
      DestroyMat(savedG);
      DestroyMat(G);
      free(cpdlsP_mem);
      return(CPDIRECT_MEM_FAIL);
    }
    K = NewDenseMat(Ny-Nc, Ny-Nc);
    if (K == NULL) {
      cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPDENSE", "CPDenseProj", MSGD_MEM_FAIL);
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
      cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPDENSE", "CPDenseProj", MSGD_MEM_FAIL);
      DestroyMat(savedG);
      DestroyMat(G);
      free(cpdlsP_mem);
      return(CPDIRECT_MEM_FAIL);      
    }
    /* If projecting in WRMS norm, allocate space for K=Q^T*D^(-1)*Q */
    if (pnorm == CP_PROJ_ERRNORM) {
      K = NewDenseMat(Nc, Nc);
      if (K == NULL) {
        cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPDENSE", "CPDenseProj", MSGD_MEM_FAIL);
        DestroyArray(beta);
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
      cpProcessError(cp_mem, CPDIRECT_MEM_FAIL, "CPDENSE", "CPDenseProj", MSGD_MEM_FAIL);
      DestroyMat(savedG);
      DestroyMat(G);
      free(cpdlsP_mem);
      return(CPDIRECT_MEM_FAIL);
    }

    break;

  }

  /* Set default Jacobian routine and Jacobian data */
  jacP = NULL;
  JP_data = NULL;

  lsetupP_exists = TRUE;

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
 *  PRIVATE FUNCTIONS FOR IMPLICIT INTEGRATION
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * cpDenseInit
 * -----------------------------------------------------------------
 * This routine does remaining initializations specific to the dense
 * linear solver.
 * -----------------------------------------------------------------
 */

static int cpDenseInit(CPodeMem cp_mem)
{
  CPDlsMem cpdls_mem;

  cpdls_mem = (CPDlsMem) lmem;
  
  nje   = 0;
  nfeDQ = 0;
  nstlj = 0;
  
  if (ode_type == CP_EXPL && jacE == NULL) {
    jacE = cpDlsDenseDQJacExpl;
    J_data = cp_mem;
  } 
  
  if (ode_type == CP_IMPL && jacI == NULL) {
    jacI = cpDlsDenseDQJacImpl;
    J_data = cp_mem;
  }

  last_flag = CPDIRECT_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * cpDenseSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the dense linear solver.
 * It makes a decision whether or not to call the Jacobian evaluation
 * routine based on various state variables, and if not it uses the 
 * saved copy (for explicit ODE only). In any case, it constructs 
 * the Newton matrix M = I - gamma*J or M = F_y' - gamma*F_y, updates 
 * counters, and calls the dense LU factorization routine.
 * -----------------------------------------------------------------
 */

static int cpDenseSetup(CPodeMem cp_mem, int convfail,
                        N_Vector yP, N_Vector ypP, N_Vector fctP,
                        booleantype *jcurPtr,
                        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  booleantype jbad, jok;
  realtype dgamma;
  long int ier;
  CPDlsMem cpdls_mem;
  int retval;

  cpdls_mem = (CPDlsMem) lmem;

  switch (ode_type) {

  case CP_EXPL:

    /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
    dgamma = ABS((gamma/gammap) - ONE);
    jbad = (nst == 0) || (nst > nstlj + CPD_MSBJ) ||
      ((convfail == CP_FAIL_BAD_J) && (dgamma < CPD_DGMAX)) ||
      (convfail == CP_FAIL_OTHER);
    jok = !jbad;
    
    /* Test if it is enough to use a saved Jacobian copy */
    if (jok) {
      *jcurPtr = FALSE;
      DenseCopy(savedJ, M);
    } else {
      nstlj = nst;
      *jcurPtr = TRUE;
      DenseZero(M);
      retval = jacE(n, tn, yP, fctP, M, J_data, tmp1, tmp2, tmp3);
      nje++;
      if (retval < 0) {
        cpProcessError(cp_mem, CPDIRECT_JACFUNC_UNRECVR, "CPDENSE", "cpDenseSetup", MSGD_JACFUNC_FAILED);
        last_flag = CPDIRECT_JACFUNC_UNRECVR;
        return(-1);
      }
      if (retval > 0) {
        last_flag = CPDIRECT_JACFUNC_RECVR;
        return(1);
      }
      DenseCopy(M, savedJ);
    }
  
    /* Scale and add I to get M = I - gamma*J */
    DenseScale(-gamma, M);
    DenseAddI(M);

    break;

  case CP_IMPL:

    /* Initialize Jacobian to 0 and call Jacobian function */
    DenseZero(M);
    retval = jacI(n, tn, gamma, yP, ypP, fctP, M, J_data, tmp1, tmp2, tmp3);
    nje++;
    if (retval < 0) {
      cpProcessError(cp_mem, CPDIRECT_JACFUNC_UNRECVR, "CPDENSE", "cpDenseSetup", MSGD_JACFUNC_FAILED);
      last_flag = CPDIRECT_JACFUNC_UNRECVR;
      return(-1);
    }
    if (retval > 0) {
      last_flag = CPDIRECT_JACFUNC_RECVR;
      return(1);
    }
  
    break;

  }

  /* Do LU factorization of M */
  ier = DenseGETRF(M, pivots); 

  /* Return 0 if the LU was complete; otherwise return 1 */
  last_flag = ier;
  if (ier > 0) return(1);
  return(0);
}

/*
 * -----------------------------------------------------------------
 * cpDenseSolve
 * -----------------------------------------------------------------
 * This routine handles the solve operation for the dense linear solver
 * by calling the dense backsolve routine.  The returned value is 0.
 * -----------------------------------------------------------------
 */

static int cpDenseSolve(CPodeMem cp_mem, N_Vector b, N_Vector weight,
                        N_Vector yC, N_Vector ypC, N_Vector fctC)
{
  CPDlsMem cpdls_mem;
  realtype *bd;

  cpdls_mem = (CPDlsMem) lmem;
  
  bd = N_VGetArrayPointer(b);

  DenseGETRS(M, pivots, bd);

  /* For BDF, scale the correction to account for change in gamma */
  if ((lmm_type == CP_BDF) && (gamrat != ONE)) {
    N_VScale(TWO/(ONE + gamrat), b, b);
  }
  
  last_flag = CPDIRECT_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * cpDenseFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the dense linear solver.
 * -----------------------------------------------------------------
 */

static void cpDenseFree(CPodeMem cp_mem)
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
 *  PRIVATE FUNCTIONS FOR PROJECTION
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * cpDenseProjInit
 * -----------------------------------------------------------------
 * This routine does remaining initializations specific to the dense
 * linear solver.
 * -----------------------------------------------------------------
 */

static int cpDenseProjInit(CPodeMem cp_mem)
{
  CPDlsProjMem cpdlsP_mem;

  cpdlsP_mem = (CPDlsProjMem) lmemP;
  
  njeP   = 0;
  nceDQ  = 0;
  nstljP = 0;
  
  if (jacP == NULL) {
    jacP = cpDlsDenseProjDQJac;
    JP_data = cp_mem;
  }  

  return(0);
}

/*
 * -----------------------------------------------------------------
 * cpDenseProjSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the dense linear solver.
 * It calls the Jacobian evaluation routine and, depending on ftype,
 * it performs various factorizations.
 * -----------------------------------------------------------------
 */

static int cpDenseProjSetup(CPodeMem cp_mem, N_Vector y, N_Vector cy,
                            N_Vector c_tmp1, N_Vector c_tmp2, N_Vector s_tmp1)
{
  long int ier;
  CPDlsProjMem cpdlsP_mem;
  realtype **g_mat, *col_i, *s_tmpd;
  long int i, j, k;
  int retval;

  cpdlsP_mem = (CPDlsProjMem) lmemP;

  g_mat = G->cols;

  /* 
   * Initialize Jacobian matrix to 0 and call Jacobian function 
   * G will contain the Jacobian transposed. 
   */
  DenseZero(G);
  retval = jacP(nc, ny, tn, y, cy, G, JP_data, c_tmp1, c_tmp2);
  njeP++;
  if (retval < 0) {
    cpProcessError(cp_mem, CPDIRECT_JACFUNC_UNRECVR, "CPDENSE", "cpDenseProjSetup", MSGD_JACFUNC_FAILED);
    return(-1);
  } else if (retval > 0) {
    return(1);
  }

  /* Save Jacobian before factorization for possible use by lmultP */
  DenseCopy(G, savedG);

  /* Factorize G, depending on ftype */
  switch (ftype) {

  case CPDIRECT_LU:

    /* 
     * LU factorization of G^T
     *      
     *    P*G^T =  | U1^T | * L^T     
     *             | U2^T |
     *
     * After factorization, P is encoded in pivotsP and
     * G^T is overwritten with U1 (nc by nc unit upper triangular), 
     * U2 ( nc by ny-nc rectangular), and L (nc by nc lower triangular).
     *
     * Return 1 if factorization failed. 
     */
    ier = DenseGETRF(G, pivotsP); 
    if (ier > 0) return(1);

    /* 
     * Build S = U1^{-1} * U2 (in place, S overwrites U2) 
     * For each row j of G, j = nc,...,ny-1, perform
     * a backward substitution (row version).
     *
     * After this step, G^T contains U1, S, and L.
     */
    for (j=nc; j<ny; j++) {
      for (i=nc-2; i>=0; i--) {
        col_i = g_mat[i];
        for (k=i+1; k<nc; k++) g_mat[i][j] -= col_i[k]*g_mat[k][j];
      }      
    }

    /* 
     * Build K = D1 + S^T * D2 * S 
     * Compute and store only the lower triangular part of K.
     */
    if (pnorm == CP_PROJ_L2NORM) cpdLUcomputeKI(cp_mem);
    else                         cpdLUcomputeKD(cp_mem, s_tmp1);

    /* 
     * Perform Cholesky decomposition of K (in place, gaxpy version)
     * 
     *     K = C*C^T
     *
     * After factorization, the lower triangular part of K contains C.
     *
     * Return 1 if factorization failed. 
     */
    ier = DensePOTRF(K);
    if (ier > 0) return(1);

    break;

  case CPDIRECT_QR:

    /* 
     * Thin QR factorization of G^T
     *
     *   G^T = Q * R
     *
     * After factorization, the upper trianguler part of G^T 
     * contains the matrix R. The lower trapezoidal part of
     * G^T, together with the array beta, encodes the orthonormal
     * columns of Q as elementary reflectors.
     */

    /* Use s_tmp1 as workspace */
    s_tmpd = N_VGetArrayPointer(s_tmp1);
    ier = DenseGEQRF(G, beta, s_tmpd);

    /* If projecting in WRMS norm */
    if (pnorm == CP_PROJ_ERRNORM) {
      /* Build K = Q^T * D^(-1) * Q */
      cpdQRcomputeKD(cp_mem, s_tmp1);
      /* Perform Cholesky decomposition of K */
      ier = DensePOTRF(K);
      if (ier > 0) return(1);
    }

    break;

  case CPDIRECT_SC:

    /* 
     * Build K = G*D^(-1)*G^T
     * Compute and store only the lower triangular part of K.
     */
    if (pnorm == CP_PROJ_L2NORM) cpdSCcomputeKI(cp_mem);
    else                         cpdSCcomputeKD(cp_mem, s_tmp1);

    /* 
     * Perform Cholesky decomposition of K (in place, gaxpy version)
     * 
     *     K = C*C^T
     *
     * After factorization, the lower triangular part of K contains C.
     *
     * Return 1 if factorization failed. 
     */
    ier = DensePOTRF(K);
    if (ier > 0) return(1);

    break;

  }

  return(0);
}

/*
 * -----------------------------------------------------------------
 * cpDenseProjSolve
 * -----------------------------------------------------------------
 * This routine handles the solve operation for the dense linear 
 * solver. The returned value is 0.
 * -----------------------------------------------------------------
 */

static int cpDenseProjSolve(CPodeMem cp_mem, N_Vector b, N_Vector x,
                            N_Vector y, N_Vector cy,
                            N_Vector c_tmp1, N_Vector s_tmp1)
{
  CPDlsProjMem cpdlsP_mem;
  realtype **g_mat, *bd, *xd, *col_i, *s_tmpd;
  realtype  *ewt_data, *d_data, *da_data, tmp;
  long int nd, i, k, pk;

  cpdlsP_mem = (CPDlsProjMem) lmemP;

  g_mat = G->cols;

  ewt_data = N_VGetArrayPointer(ewt);
  bd = N_VGetArrayPointer(b);
  xd = N_VGetArrayPointer(x);
  d_data = N_VGetArrayPointer(s_tmp1);

  nd = ny - nc;

  /* Solve the linear system, depending on ftype */
  switch (ftype) {

  case CPDIRECT_LU:

    /* 
     * Solve L*U1*alpha = bd
     *   (a) solve L*beta = bd using fwd. subst. (row version)
     *   (b) solve U1*alpha = beta using bckwd. subst (row version) 
     * where L^T and U1^T are stored in G[0...nc-1][0...nc-1].
     * beta and then alpha overwrite bd.
     */
    bd[0] /= g_mat[0][0];
    for (i=1; i<nc; i++) {
      col_i = g_mat[i];
      for (k=0; k<i; k++) bd[i] -= col_i[k]*bd[k];
      bd[i] /= col_i[i];
    }
    for (i=nc-2; i>=0; i--) {
      col_i = g_mat[i];
      for (k=i+1; k<nc; k++) bd[i] -= col_i[k]*bd[k];
    }  

    /* 
     * Compute S^T * (D1 * alpha)
     * alpha is stored in bd.
     * S^T is stored in g_mat[nc...ny-1][0...nc-1].
     * Store result in x2 = x[nc...ny-1].
     */
    if (pnorm == CP_PROJ_ERRNORM) {

      /* Load squared error weights into d */
      for (k=0; k<ny; k++) d_data[k] = ewt_data[k] * ewt_data[k];
      /* Permute elements of d, based on pivotsP. Swap d[k] and d[pivotsP[k]]. */
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
      for(i=0; i<nd; i++) {
        xd[nc+i] = ZERO;
        for(k=0; k<nc; k++) xd[nc+i] += g_mat[k][nc+i]*da_data[k];
      }

    } else {

      /* Compute S^T * alpha */
      for(i=0; i<nd; i++) {
        xd[nc+i] = ZERO;
        for(k=0; k<nc; k++) xd[nc+i] += g_mat[k][nc+i]*bd[k];
      }

    }

    /* 
     * Solve K*x2 = S^T*D1*alpha, using the Cholesky decomposition available in K.
     * S^T*D1*alpha is stored in x2 = x[nc...ny-1].
     */
    DensePOTRS(K, &xd[nc]);

    /* 
     * Compute x1 = alpha - S*x2 
     * alpha is stored in bd.
     * x2 is stored in x[nc...ny-1].
     * S^T is stored in g_mat[nc...ny-1][0...nc-1].
     * Store result in x1 = x[0...nc-1].
     */
    for (i=0; i<nc; i++) {
      xd[i] = bd[i];
      col_i = g_mat[i];
      for (k=nc; k<ny; k++) xd[i] -= col_i[k]*xd[k];
    }

    /* 
     * Compute P^T * x, where P is encoded into pivotsP.
     * Store result in x.
     */
    for (k=nc-1; k>=0; k--) {
      pk = pivotsP[k];
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
     * alpha overwrites bd.
     */
    bd[0] /= g_mat[0][0];
    for (i=1; i<nc; i++) {
      col_i = g_mat[i];
      for (k=0; k<i; k++) bd[i] -= bd[k]*col_i[k];
      bd[i] /= col_i[i];
    }

    /* If projecting in WRMS norm, solve K*beta = alpha */
    if (pnorm == CP_PROJ_ERRNORM) DensePOTRS(K, bd);

    /* Compute x = Q*alpha */
    s_tmpd = N_VGetArrayPointer(s_tmp1);
    DenseORMQR(G, beta, bd, xd, s_tmpd); 

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
    DensePOTRS(K, bd);

    /* Compute x = G^T * xi
     * G^T is stored in g_mat[0...ny-1][0...nc-1]
     */
    for(i=0; i<ny; i++) {
      xd[i] = ZERO;
      for(k=0; k<nc; k++) xd[i] += g_mat[k][i]*bd[k];
    }

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
 * -----------------------------------------------------------------
 * cpDenseProjMult
 * -----------------------------------------------------------------
 * This routine computes the Jacobian-vector product used a saved 
 * Jacobian copy.
 * -----------------------------------------------------------------
 */

static void cpDenseProjMult(CPodeMem cp_mem, N_Vector x, N_Vector Gx)
{
  CPDlsProjMem cpdlsP_mem;
  realtype **g_mat, *col_j, *x_d, *Gx_d;
  long int j, k;

  cpdlsP_mem = (CPDlsProjMem) lmemP;

  g_mat = savedG->cols;
  x_d   = N_VGetArrayPointer(x);
  Gx_d  = N_VGetArrayPointer(Gx);

  for (j=0; j<nc; j++) {
    col_j = g_mat[j];
    Gx_d[j] = ZERO;
    for (k=0; k<ny; k++) Gx_d[j] += col_j[k]*x_d[k];
  }
}

/*
 * -----------------------------------------------------------------
 * cpDenseProjFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the dense linear solver.
 * -----------------------------------------------------------------
 */

static void cpDenseProjFree(CPodeMem cp_mem)
{
  CPDlsProjMem  cpdlsP_mem;

  cpdlsP_mem = (CPDlsProjMem) lmemP;
  
  DestroyMat(G);
  DestroyMat(savedG);
  switch (ftype) {
  case CPDIRECT_LU:
    DestroyArray(pivotsP);
    DestroyMat(K);
    break;
  case CPDIRECT_QR:
    DestroyArray(beta);
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
 * Compute the lower triangle of K = I + S^T*S
 */
static void cpdLUcomputeKI(CPodeMem cp_mem)
{
  CPDlsProjMem cpdlsP_mem;
  realtype **g_mat, **k_mat, *k_col_j;
  long int nd, i, j, k;

  cpdlsP_mem = (CPDlsProjMem) lmemP;

  g_mat = G->cols;
  k_mat = K->cols;

  nd = ny-nc;

  /* Load K column by column */
  for (j=0; j<nd; j++) {

    k_col_j = k_mat[j];
    
    for (i=j; i<nd; i++) {
      k_col_j[i] = ZERO;
      for (k=0; k<nc; k++) k_col_j[i] += g_mat[k][nc+i]*g_mat[k][nc+j];
    }
    
    k_col_j[j] += ONE;

  }
  
}

/*
 * Compute the lower triangle of K = D1 + S^T*D2*S,
 * D = diag(D1, D2) = P*W*P^T, W is a diagonal matrix
 * containing the squared error weights, and P is the 
 * permutation matrix encoded into pivotsP.
 * D1 has length nc and D2 has length (ny-nc).
 */
static void cpdLUcomputeKD(CPodeMem cp_mem, N_Vector d)
{
  CPDlsProjMem cpdlsP_mem;
  realtype *d_data, *ewt_data, tmp;
  realtype **g_mat, **k_mat, *k_col_j;
  long int nd, i, j, k, pk;

  cpdlsP_mem = (CPDlsProjMem) lmemP;

  ewt_data = N_VGetArrayPointer(ewt);
  d_data   = N_VGetArrayPointer(d);

  g_mat = G->cols;
  k_mat = K->cols;

  nd = ny-nc;

  /* Load squared error weights into d */
  for (k=0; k<ny; k++) d_data[k] = ewt_data[k] * ewt_data[k];

  /* Permute elements of d, based on pivotsP. Swap d[k] and d[pivotsP[k]]. */
  for (k=0; k<nc; k++) {
    pk = pivotsP[k];
    if (pk != k) {
      tmp = d_data[k];
      d_data[k]  = d_data[pk];
      d_data[pk] = tmp;
    }
  }

  /* load K column by column */
  for (j=0; j<nd; j++) {

    k_col_j = k_mat[j];
 
    for (i=j; i<nd; i++) {
      k_col_j[i] = ZERO;
      for(k=0; k<nc; k++) 
        k_col_j[i] += g_mat[k][nc+i] * d_data[k]*d_data[k] * g_mat[k][nc+j];
    }

    k_col_j[j] += d_data[j]*d_data[j];

  }

}

static void cpdQRcomputeKD(CPodeMem cp_mem, N_Vector d)
{
  /* RADU:: implement this ... */
}


/*
 * Compute the lower triangle of K = G * G^T
 * G^T is available in g_mat[0...ny][0...nc].
 */
static void cpdSCcomputeKI(CPodeMem cp_mem)
{
  CPDlsProjMem cpdlsP_mem;
  realtype **g_mat, **k_mat, *k_col_j;
  long int i, j, k;

  cpdlsP_mem = (CPDlsProjMem) lmemP;

  g_mat = G->cols;
  k_mat = K->cols;

  /* Load K column by column */
  for (j=0; j<nc; j++) {
    k_col_j = k_mat[j];
    for (i=0; i<nc; i++) {
      k_col_j[i] = ZERO;
      for (k=0; k<ny; k++) k_col_j[i] += g_mat[i][k]*g_mat[j][k];
    }    
  }

}

static void cpdSCcomputeKD(CPodeMem cp_mem, N_Vector d)
{
  /* RADU:: implement this ... */
}


