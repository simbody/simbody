/*
 * -----------------------------------------------------------------
 * $Revision: 1.9 $
 * $Date: 2007/10/26 21:51:29 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban  @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for the main CPODES integrator.
 * It is independent of the CPODES linear solver in use.
 * -----------------------------------------------------------------
 *
 * CONTENTS
 *
 * Exported functions:
 *   CPodeCreate   - create CPODES solver object
 *   CPodeInit     - initialize solver
 *   CPodeReInit   - re-initialize solver
 *   CPodeProjInit - initialize internal projection algorithm
 *   CPodeProjDefine - initialize user-provided projection
 *   CPode         - main solver function
 *   CPodeGetDky   - dense output function
 *   CPodeFree     - memory deallocation
 *
 * Private functions
 *   - Memory allocation/deallocation and initialization functions
 *   - Initial step size evaluation functions
 *   - Main step function
 *   - LMM-related functions
 *   - Nonlinear and error test failure handlers
 *   - Succesful step completion functions
 *   - BDF stability limit detection functions
 *   - Projection functions
 *   - Internal error weight evaluation functions
 *   - Error reporting functions
 */

/* 
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "cpodes_private.h"
#include <sundials/sundials_math.h>

/* 
 * =================================================================
 * FUNCTION SPECIFIC CONSTANTS
 * =================================================================
 */

/*
 * Algorithmic constants
 * ---------------------
 */

/* constants used in estimating an initial step size */
#define HLB_FACTOR RCONST(100.0)
#define HUB_FACTOR RCONST(0.1)
#define H_BIAS     HALF
#define H_MAXITERS 4

/* constants used in taking a step */
#define THRESH RCONST(1.5)
#define ETAMX1 RCONST(10000.0) 
#define ETAMX2 RCONST(10.0)
#define ETAMX3 RCONST(10.0)
#define ETAMXF RCONST(0.2)
#define ETAMIN RCONST(0.1)
#define ADDON  RCONST(0.000001)
#define BIAS1  RCONST(6.0)
#define BIAS2  RCONST(6.0)
#define BIAS3  RCONST(10.0)

/* 
 * =================================================================
 * PROTOTYPES FOR PRIVATE FUNCTIONS
 * =================================================================
 */

/* Memory allocation/deallocation and initialization functions */
static booleantype cpCheckNvector(N_Vector tmpl);
static booleantype cpAllocVectors(CPodeMem cp_mem, N_Vector tmpl, int tol);
static void cpFreeVectors(CPodeMem cp_mem);
static booleantype cpProjAlloc(CPodeMem cp_mem, N_Vector c_tmpl, N_Vector s_tmpl);
static void cpProjFree(CPodeMem cp_mem);
static booleantype cpQuadAlloc(CPodeMem cp_mem, N_Vector q_tmpl);
static void cpQuadFree(CPodeMem cp_mem);

static int cpInitialSetup(CPodeMem cp_mem);

/* Initial step size evaluation functions */
static int cpHin(CPodeMem cp_mem, realtype tout);
static realtype cpUpperBoundH0(CPodeMem cp_mem, realtype tdist);
static int cpHinExpl(CPodeMem cp_mem, realtype hlb, realtype hub, int sign, realtype *h0);
static int cpYppNorm(CPodeMem cp_mem, realtype hg, realtype *yppnorm);
static int cpHinImpl(CPodeMem cp_mem, realtype tdist, realtype *h0);

/* Main step function */
static int cpStep(CPodeMem cp_mem);

/* Functions acting on Nordsieck history array */
static void cpPredict(CPodeMem cp_mem);
static void cpCorrect(CPodeMem cp_mem);

/* Quadrature correction function */
static int cpQuadNls(CPodeMem cp_mem, realtype saved_t, int *ncfPtr);

/* LMM-related functions */
static void cpAdjustParams(CPodeMem cp_mem);
static void cpAdjustOrder(CPodeMem cp_mem, int deltaq);
static void cpAdjustAdams(CPodeMem cp_mem, int deltaq);
static void cpAdjustBDF(CPodeMem cp_mem, int deltaq);
static void cpIncreaseBDF(CPodeMem cp_mem);
static void cpDecreaseBDF(CPodeMem cp_mem);
static void cpSet(CPodeMem cp_mem);
static void cpSetAdams(CPodeMem cp_mem);
static realtype cpAdamsStart(CPodeMem cp_mem, realtype m[]);
static void cpAdamsFinish(CPodeMem cp_mem, realtype m[], realtype M[], realtype hsum);
static realtype cpAltSum(int iend, realtype a[], int k);
static void cpSetBDF(CPodeMem cp_mem);
static void cpSetTqBDF(CPodeMem cp_mem, realtype hsum, realtype alpha0,
                       realtype alpha0_hat, realtype xi_inv, realtype xistar_inv);

/* Error test function */
static int cpDoErrorTest(CPodeMem cp_mem, realtype saved_t, realtype acor_norm,
                         int *nefPtr, realtype *dsmPtr);

/* Succesful step completion functions */
static void cpCompleteStep(CPodeMem cp_mem);
static void cpPrepareNextStep(CPodeMem cp_mem, realtype dsm);
static void cpSetEta(CPodeMem cp_mem);
static realtype cpComputeEtaqm1(CPodeMem cp_mem);
static realtype cpComputeEtaqp1(CPodeMem cp_mem);
static void cpChooseEta(CPodeMem cp_mem);
static int  cpHandleFailure(CPodeMem cp_mem,int flag);

/* BDF stability limit detection functions */
static void cpBDFStab(CPodeMem cp_mem);
static int cpSLdet(CPodeMem cp_mem);

/* Internal error weight evaluation functions */
static int cpEwtSetSS(CPodeMem cp_mem, N_Vector ycur, N_Vector weight);
static int cpEwtSetSV(CPodeMem cp_mem, N_Vector ycur, N_Vector weight);
static int cpQuadEwtSet(CPodeMem cp_mem, N_Vector qcur, N_Vector weightQ);
static int cpQuadEwtSetSS(CPodeMem cp_mem, N_Vector qcur, N_Vector weightQ);
static int cpQuadEwtSetSV(CPodeMem cp_mem, N_Vector qcur, N_Vector weightQ);

/* Function for combined norms */
static realtype cpQuadUpdateNorm(CPodeMem cp_mem, realtype old_nrm,
                                 N_Vector xQ, N_Vector wQ);
static realtype cpQuadUpdateDsm(CPodeMem cp_mem, realtype old_dsm, 
                                realtype dsmQ);

/* Error reporting functions */
static void cpErrHandler(int error_code, const char *module,
                         const char *function, char *msg, void *data);

/* 
 * =================================================================
 * EXPORTED FUNCTIONS
 * =================================================================
 */

/* 
 * CPodeCreate
 *
 * CPodeCreate creates an internal memory block for a problem to 
 * be solved by CPODES.
 * If successful, CPodeCreate returns a pointer to the problem memory. 
 * This pointer should be passed to CPodeInit.  
 * If an initialization error occurs, CPodeCreate prints an error 
 * message to standard err and returns NULL. 
 */

void *CPodeCreate(int ode_type, int lmm_type, int nls_type)
{
  int maxord;
  CPodeMem cp_mem;

  /* Test inputs */
  if ((ode_type != CP_EXPL) && (ode_type != CP_IMPL)) {
    cpProcessError(NULL, 0, "CPODES", "CPodeCreate", MSGCP_BAD_ODE);
    return(NULL);
  }
  if ((lmm_type != CP_ADAMS) && (lmm_type != CP_BDF)) {
    cpProcessError(NULL, 0, "CPODES", "CPodeCreate", MSGCP_BAD_LMM);
    return(NULL);
  }
  if ((nls_type != CP_FUNCTIONAL) && (nls_type != CP_NEWTON)) {
    cpProcessError(NULL, 0, "CPODES", "CPodeCreate", MSGCP_BAD_NLS);
    return(NULL);
  }
  if ((ode_type == CP_IMPL) && (nls_type == CP_FUNCTIONAL)) {
    cpProcessError(NULL, 0, "CPODES", "CPodeCreate", MSGCP_BAD_ODE_NLS);
    return(NULL);
  }

  /* Allocate space for solver object */
  cp_mem = NULL;
  cp_mem = (CPodeMem) malloc(sizeof(struct CPodeMemRec));
  if (cp_mem == NULL) {
    cpProcessError(NULL, 0, "CPODES", "CPodeCreate", MSGCP_CPMEM_FAIL);
    return(NULL);
  }

  maxord = (lmm_type == CP_ADAMS) ? ADAMS_Q_MAX : BDF_Q_MAX;

  /* Copy input parameters into cp_mem */
  cp_mem->cp_ode_type = ode_type;
  cp_mem->cp_lmm_type = lmm_type;
  cp_mem->cp_nls_type = nls_type;

  /* Set uround */
  cp_mem->cp_uround = UNIT_ROUNDOFF;

  /* Initialize required values */
  cp_mem->cp_fi       = NULL;
  cp_mem->cp_fe       = NULL;
  cp_mem->cp_f_data   = NULL;

  /* Set default values for integrator optional inputs */
  cp_mem->cp_efun     = cpEwtSet;
  cp_mem->cp_e_data   = (void *) cp_mem;
  cp_mem->cp_ehfun    = cpErrHandler;
  cp_mem->cp_eh_data  = (void *) cp_mem;
  cp_mem->cp_errfp    = stderr;
  cp_mem->cp_qmax     = maxord;
  cp_mem->cp_mxstep   = MXSTEP_DEFAULT;
  cp_mem->cp_mxhnil   = MXHNIL_DEFAULT;
  cp_mem->cp_sldeton  = FALSE;
  cp_mem->cp_hin      = ZERO;
  cp_mem->cp_hmin     = HMIN_DEFAULT;
  cp_mem->cp_hmax_inv = HMAX_INV_DEFAULT;
  cp_mem->cp_tstopset = FALSE;
  cp_mem->cp_maxcor   = NLS_MAXCOR;
  cp_mem->cp_maxnef   = MXNEF;
  cp_mem->cp_maxncf   = MXNCF;
  cp_mem->cp_nlscoef  = NLS_TEST_COEF;

  /* Set the linear solver addresses to NULL. */
  cp_mem->cp_linit  = NULL;
  cp_mem->cp_lsetup = NULL;
  cp_mem->cp_lsolve = NULL;
  cp_mem->cp_lfree  = NULL;
  cp_mem->cp_lmem   = NULL;
  cp_mem->cp_lsetup_exists = FALSE;

  /* Initialize projection variables */
  cp_mem->cp_proj_enabled  = FALSE;
  cp_mem->cp_proj_type     = CP_PROJ_INTERNAL;
  cp_mem->cp_cnstr_type    = CP_CNSTR_NONLIN;
  cp_mem->cp_cfun          = NULL;
  cp_mem->cp_c_data        = NULL;
  cp_mem->cp_pfun          = NULL;
  cp_mem->cp_p_data        = NULL;
  cp_mem->cp_prjcoef       = PRJ_TEST_COEF;
  cp_mem->cp_maxcorP       = PRJ_MAXCOR;
  cp_mem->cp_project_err   = TRUE;
  cp_mem->cp_test_cnstr    = FALSE;
  cp_mem->cp_proj_freq     = 1;
  cp_mem->cp_lsetupP_freq  = PRJ_MSBLS;
  cp_mem->cp_maxnpf        = MXNPF;

  /* Initialize quadrature variables */
  cp_mem->cp_quadr    = FALSE;
  cp_mem->cp_qfun     = NULL;
  cp_mem->cp_q_data   = NULL;
  cp_mem->cp_errconQ  = FALSE;

  /* Set the linear solver addresses to NULL. */
  cp_mem->cp_linitP  = NULL;
  cp_mem->cp_lsetupP = NULL;
  cp_mem->cp_lsolveP = NULL;
  cp_mem->cp_lfreeP  = NULL;
  cp_mem->cp_lmemP   = NULL;
  cp_mem->cp_lsetupP_exists = FALSE;

  /* Initialize root finding variables */
  cp_mem->cp_doRootfinding = FALSE;
  cp_mem->cp_nrtfn         = 0;
  cp_mem->cp_glo           = NULL;
  cp_mem->cp_ghi           = NULL;
  cp_mem->cp_grout         = NULL;
  cp_mem->cp_gactive       = NULL;
  cp_mem->cp_iroots        = NULL;
  cp_mem->cp_rootdir       = NULL;
  cp_mem->cp_gfun          = NULL;
  cp_mem->cp_g_data        = NULL;

  /* Set the saved value qmax_alloc */
  cp_mem->cp_qmax_alloc = maxord;
  cp_mem->cp_qmax_allocQ = maxord;

  /* Initialize lrw and liw */
  cp_mem->cp_lrw = 58 + 2*L_MAX + NUM_TESTS;
  cp_mem->cp_liw = 40;

  /* No mallocs have been done yet */
  cp_mem->cp_MallocDone         = FALSE;
  cp_mem->cp_VabstolMallocDone  = FALSE;
  cp_mem->cp_projMallocDone     = FALSE;
  cp_mem->cp_quadMallocDone     = FALSE;
  cp_mem->cp_VabstolQMallocDone = FALSE;
  cp_mem->cp_rootMallocDone     = FALSE;

  /* Return pointer to CPODES memory block */

  return((void *)cp_mem);
}

/*-----------------------------------------------------------------*/

#define ode_type (cp_mem->cp_ode_type)
#define lmm_type (cp_mem->cp_lmm_type)
#define nls_type (cp_mem->cp_nls_type)
#define lrw      (cp_mem->cp_lrw)
#define liw      (cp_mem->cp_liw)

/*-----------------------------------------------------------------*/

/*
 * CPodeInit
 * 
 * CPodeInit allocates and initializes memory for a problem. All 
 * problem inputs are checked for errors. If any error occurs during 
 * initialization, it is reported to the file whose file pointer is 
 * errfp and an error flag is returned. Otherwise, it returns CP_SUCCESS
 */

int CPodeInit(void *cpode_mem,
              void *fun, void *f_data,
              realtype t0, N_Vector y0, N_Vector yp0,
              int tol_type, realtype reltol, void *abstol)
{
  CPodeMem cp_mem;
  booleantype nvectorOK, allocOK, neg_abstol;
  long int lrw1, liw1;
  int i,k;

  /* Check cpode_mem */
  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeInit", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  /* Check for legal input parameters */
  if (fun == NULL) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeInit", MSGCP_NULL_F);
    return(CP_ILL_INPUT);
  }
  if (y0==NULL) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeInit", MSGCP_NULL_Y0);
    return(CP_ILL_INPUT);
  }
  if ( (ode_type==CP_IMPL) && (yp0==NULL) ) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeInit", MSGCP_NULL_YP0);
    return(CP_ILL_INPUT);
  }
  if ((tol_type != CP_SS) && (tol_type != CP_SV) && (tol_type != CP_WF)) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeInit", MSGCP_BAD_ITOL);
    return(CP_ILL_INPUT);
  }
  /* Test if all required vector operations are implemented */
  nvectorOK = cpCheckNvector(y0);
  if(!nvectorOK) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeInit", MSGCP_BAD_NVECTOR);
    return(CP_ILL_INPUT);
  }

  /* Test tolerances */
  if (tol_type != CP_WF) {

    if (abstol == NULL) {
      cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeInit", MSGCP_NULL_ABSTOL);
      return(CP_ILL_INPUT);
    }

    if (reltol < ZERO) {
      cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeInit", MSGCP_BAD_RELTOL);
      return(CP_ILL_INPUT);
    }

    if (tol_type == CP_SS)
      neg_abstol = (*((realtype *)abstol) < ZERO);
    else
      neg_abstol = (N_VMin((N_Vector)abstol) < ZERO);

    if (neg_abstol) {
      cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeInit", MSGCP_BAD_ABSTOL);
      return(CP_ILL_INPUT);
    }

  }

  /* Set space requirements for one N_Vector */
  if (y0->ops->nvspace != NULL) {
    N_VSpace(y0, &lrw1, &liw1);
  } else {
    lrw1 = 0;
    liw1 = 0;
  }
  cp_mem->cp_lrw1 = lrw1;
  cp_mem->cp_liw1 = liw1;

  /* Allocate the vectors (using y0 as a template) */
  allocOK = cpAllocVectors(cp_mem, y0, tol_type);
  if (!allocOK) {
    cpProcessError(cp_mem, CP_MEM_FAIL, "CPODES", "CPodeInit", MSGCP_MEM_FAIL);
    return(CP_MEM_FAIL);
  }

  /* 
   * All error checking is complete at this point 
   */

  /* Copy tolerances into memory */
  cp_mem->cp_tol_type = tol_type;
  cp_mem->cp_reltol   = reltol;      

  if (tol_type == CP_SS) {
    cp_mem->cp_Sabstol = *((realtype *)abstol);
  } else if (tol_type == CP_SV) {
    N_VScale(ONE, (N_Vector)abstol, cp_mem->cp_Vabstol);
  }

  /* Copy the input parameters into CPODES state */
  if (ode_type == CP_EXPL) cp_mem->cp_fe = (CPRhsFn) fun;
  else                     cp_mem->cp_fi = (CPResFn) fun;
  cp_mem->cp_f_data = f_data;
  cp_mem->cp_tn = t0;

  /* Set step parameters */
  cp_mem->cp_q      = 1;
  cp_mem->cp_L      = 2;
  cp_mem->cp_qwait  = cp_mem->cp_L;
  cp_mem->cp_etamax = ETAMX1;
  cp_mem->cp_qu     = 0;
  cp_mem->cp_hu     = ZERO;
  cp_mem->cp_tolsf  = ONE;

  /* Initialize the history array zn */
  N_VScale(ONE, y0, cp_mem->cp_zn[0]);
  if(ode_type==CP_IMPL) N_VScale(ONE, yp0, cp_mem->cp_zn[1]);

  /* Initialize all the counters */
  cp_mem->cp_nst     = 0;
  cp_mem->cp_nfe     = 0;
  cp_mem->cp_ncfn    = 0;
  cp_mem->cp_netf    = 0;
  cp_mem->cp_nni     = 0;
  cp_mem->cp_nsetups = 0;
  cp_mem->cp_nhnil   = 0;
  cp_mem->cp_nstlset = 0;
  cp_mem->cp_nscon   = 0;
  cp_mem->cp_nge     = 0;

  cp_mem->cp_irfnd   = 0;

  cp_mem->cp_nproj    = 0;
  cp_mem->cp_nprf     = 0;
  cp_mem->cp_nce      = 0;
  cp_mem->cp_nstlprj  = 0;
  cp_mem->cp_nsetupsP = 0;
  cp_mem->cp_first_proj = TRUE;

  /* Initialize other integrator optional outputs */
  cp_mem->cp_h0u     = ZERO;
  cp_mem->cp_next_h  = ZERO;
  cp_mem->cp_next_q  = 0;

  /* 
   * Initialize Stablilty Limit Detection data.
   * NOTE: We do this even if stab lim det was not turned on yet.
   *       This way, the user can turn it on at any later time.
   */
  cp_mem->cp_nor = 0;
  for (i = 1; i <= 5; i++)
    for (k = 1; k <= 3; k++) 
      cp_mem->cp_ssdat[i-1][k-1] = ZERO;

  /* Problem has been successfully initialized */
  cp_mem->cp_MallocDone = TRUE;

  return(CP_SUCCESS);
}

/*-----------------------------------------------------------------*/

#define lrw1 (cp_mem->cp_lrw1)
#define liw1 (cp_mem->cp_liw1)

/*-----------------------------------------------------------------*/

/*
 * CPodeReInit
 *
 * CPodeReInit re-initializes CPODE's memory for a problem, assuming
 * it has already been allocated in a prior CPodeInit call.
 * All problem specification inputs are checked for errors.
 * If any error occurs during initialization, it is reported to the
 * file whose file pointer is errfp.
 * The return value is CP_SUCCESS = 0 if no errors occurred, or
 * a negative value otherwise.
 */

int CPodeReInit(void *cpode_mem, 
                void *fun, void *f_data,
                realtype t0, N_Vector y0, N_Vector yp0,
                int tol_type, realtype reltol, void *abstol)
{
  CPodeMem cp_mem;
  booleantype neg_abstol;
  int i,k;
 
  /* Check cpode_mem */
  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeReInit", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  /* Check if cpode_mem was allocated */
  if (cp_mem->cp_MallocDone == FALSE) {
    cpProcessError(cp_mem, CP_NO_MALLOC, "CPODES", "CPodeReInit", MSGCP_NO_MALLOC);
    return(CP_NO_MALLOC);
  }

  /* Check for legal input parameters */
  if (fun == NULL) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeReInit", MSGCP_NULL_F);
    return(CP_ILL_INPUT);
  }
  if (y0 == NULL) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeReInit", MSGCP_NULL_Y0);
    return(CP_ILL_INPUT);
  }
  if ( (ode_type==CP_IMPL) && (yp0==NULL) ) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeReInit", MSGCP_NULL_YP0);
    return(CP_ILL_INPUT);
  }
  if ((tol_type != CP_SS) && (tol_type != CP_SV) && (tol_type != CP_WF)) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeReInit", MSGCP_BAD_ITOL);
    return(CP_ILL_INPUT);
  }

  /* Test tolerances */
  if (tol_type != CP_WF) {

    if (abstol == NULL) {
      cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeReInit", MSGCP_NULL_ABSTOL);
      return(CP_ILL_INPUT);
    }

    if (reltol < ZERO) {
      cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeReInit", MSGCP_BAD_RELTOL);
      return(CP_ILL_INPUT);
    }
    
    if (tol_type == CP_SS) {
      neg_abstol = (*((realtype *)abstol) < ZERO);
    } else {
      neg_abstol = (N_VMin((N_Vector)abstol) < ZERO);
    }
    
    if (neg_abstol) {
      cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeReInit", MSGCP_BAD_ABSTOL);
      return(CP_ILL_INPUT);
    }

  }

  /* 
   * All error checking is complete at this point 
   */  

  /* Copy tolerances into memory */
  cp_mem->cp_tol_type = tol_type;
  cp_mem->cp_reltol   = reltol;    

  if (tol_type == CP_SS) {
    cp_mem->cp_Sabstol = *((realtype *)abstol);
  } else if (tol_type == CP_SV) {
    if ( !(cp_mem->cp_VabstolMallocDone) ) {
      cp_mem->cp_Vabstol = N_VClone(y0);
      lrw += lrw1;
      liw += liw1;
      cp_mem->cp_VabstolMallocDone = TRUE;
    }
    N_VScale(ONE, (N_Vector)abstol, cp_mem->cp_Vabstol);
  }
  
  /* Copy the input parameters into CPODES state */
  if (ode_type == CP_EXPL) cp_mem->cp_fe = (CPRhsFn) fun;
  else                     cp_mem->cp_fi = (CPResFn) fun;
  cp_mem->cp_f_data = f_data;
  cp_mem->cp_tn = t0;
  
  /* Set step parameters */
  cp_mem->cp_q      = 1;
  cp_mem->cp_L      = 2;
  cp_mem->cp_qwait  = cp_mem->cp_L;
  cp_mem->cp_etamax = ETAMX1;
  cp_mem->cp_qu     = 0;
  cp_mem->cp_hu     = ZERO;
  cp_mem->cp_tolsf  = ONE;

  /* Initialize the history array zn */
  N_VScale(ONE, y0,  cp_mem->cp_zn[0]);
  if(ode_type==CP_IMPL) N_VScale(ONE, yp0, cp_mem->cp_zn[1]);
 
  /* Initialize all the counters */
  cp_mem->cp_nst     = 0;
  cp_mem->cp_nfe     = 0;
  cp_mem->cp_ncfn    = 0;
  cp_mem->cp_netf    = 0;
  cp_mem->cp_nni     = 0;
  cp_mem->cp_nsetups = 0;
  cp_mem->cp_nhnil   = 0;
  cp_mem->cp_nstlset = 0;
  cp_mem->cp_nscon   = 0;
  cp_mem->cp_nge     = 0;

  cp_mem->cp_irfnd   = 0;

  cp_mem->cp_nproj    = 0;
  cp_mem->cp_nprf     = 0;
  cp_mem->cp_nce      = 0;
  cp_mem->cp_nstlprj  = 0;
  cp_mem->cp_nsetupsP = 0;
  cp_mem->cp_first_proj = TRUE;

  /* Initialize other integrator optional outputs */
  cp_mem->cp_h0u     = ZERO;
  cp_mem->cp_next_h  = ZERO;
  cp_mem->cp_next_q  = 0;

  /* Initialize Stablilty Limit Detection data */
  cp_mem->cp_nor = 0;
  for (i = 1; i <= 5; i++)
    for (k = 1; k <= 3; k++) 
      cp_mem->cp_ssdat[i-1][k-1] = ZERO;
  
  /* Problem has been successfully re-initialized */
  return(CP_SUCCESS);
}

/*-----------------------------------------------------------------*/

/*
 * CPodeProjInit
 *
 */

int CPodeProjInit(void *cpode_mem, int proj_norm, 
                  int cnstr_type, CPCnstrFn cfun, void *c_data, 
                  N_Vector ctol)
{
  CPodeMem cp_mem;
  booleantype allocOK;

  /* Check cpode_mem */
  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeProjInit", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  /* Check if cpode_mem was allocated */
  if (cp_mem->cp_MallocDone == FALSE) {
    cpProcessError(cp_mem, CP_NO_MALLOC, "CPODES", "CPodeProjInit", MSGCP_NO_MALLOC);
    return(CP_NO_MALLOC);
  }

  /* Check for legal input parameters */
  if ((proj_norm != CP_PROJ_L2NORM) && (proj_norm != CP_PROJ_ERRNORM) ) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeProjInit", MSGCP_BAD_NORM);
    return(CP_ILL_INPUT);
  }
  if ((cnstr_type != CP_CNSTR_LIN) && (cnstr_type != CP_CNSTR_NONLIN) ) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeProjInit", MSGCP_BAD_CNSTR);
    return(CP_ILL_INPUT);
  }
  if (cfun == NULL) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeProjInit", MSGCP_NULL_C);
    return(CP_ILL_INPUT);
  }
  if (ctol == NULL) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeProjInit", MSGCP_CTOL_NULL);
    return(CP_ILL_INPUT);    
  }

  /* Allocate memory for vectors used in the projection algorithm,
     using ctol and tempv as templates */
  if (!cp_mem->cp_projMallocDone) {
    allocOK = cpProjAlloc(cp_mem, ctol, cp_mem->cp_tempv);
    if (!allocOK) {
      cpProcessError(cp_mem, CP_MEM_FAIL, "CPODES", "CPodeProjInit", MSGCP_MEM_FAIL);
      return(CP_MEM_FAIL);
    }
    cp_mem->cp_projMallocDone = TRUE;
  }

  /* Set variable values in CPODES memory block */
  cp_mem->cp_proj_norm  = proj_norm;
  cp_mem->cp_cnstr_type = cnstr_type;
  cp_mem->cp_cfun       = cfun;
  cp_mem->cp_c_data     = c_data;

  /* Copy 1/ctol into memory */
  N_VScale(ONE, ctol, cp_mem->cp_ctol);
  N_VInv(cp_mem->cp_ctol, cp_mem->cp_ctol);

  /* internal projection is now enabled */
  cp_mem->cp_proj_type = CP_PROJ_INTERNAL;
  cp_mem->cp_proj_enabled = TRUE;

  return(CP_SUCCESS);
}

/*
 * CPodeProjDefine
 */

int CPodeProjDefine(void *cpode_mem, CPProjFn pfun, void *p_data)
{
  CPodeMem cp_mem;

  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeProjDefine", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  cp_mem->cp_pfun   = pfun;
  cp_mem->cp_p_data = p_data;

  /* user-defined projection is now enabled */
  cp_mem->cp_proj_type = CP_PROJ_USER;
  cp_mem->cp_proj_enabled = TRUE;

  return(CP_SUCCESS);
}

/*-----------------------------------------------------------------*/

/*
 * CPodeQuadInit
 *
 * CPodeQuadInit allocates and initializes quadrature related 
 * memory for a problem. All problem specification inputs are 
 * checked for errors. If any error occurs during initialization, 
 * it is reported to the error handler function.
 * The return value is CP_SUCCESS = 0 if no errors occurred, or
 * a negative value otherwise.
 */

int CPodeQuadInit(void *cpode_mem, CPQuadFn qfun, void *q_data, N_Vector q0)
{
  CPodeMem cp_mem;
  booleantype allocOK;
  long int lrw1Q, liw1Q;

  /* Check cpode_mem */
  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeQuadInit", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  /* Set space requirements for one N_Vector */
  N_VSpace(q0, &lrw1Q, &liw1Q);
  cp_mem->cp_lrw1Q = lrw1Q;
  cp_mem->cp_liw1Q = liw1Q;

  /* Allocate the vectors (using q0 as a template) */
  allocOK = cpQuadAlloc(cp_mem, q0);
  if (!allocOK) {
    cpProcessError(cp_mem, CP_MEM_FAIL, "CPODES", "CPodeQuadInit", MSGCP_MEM_FAIL);
    return(CP_MEM_FAIL);
  }

  /* Initialize znQ[0] in the history array */
  N_VScale(ONE, q0, cp_mem->cp_znQ[0]);

  /* Copy the input parameters into CPODES state */
  cp_mem->cp_qfun   = qfun;
  cp_mem->cp_q_data = q_data;

  /* Initialize counters */
  cp_mem->cp_nqe   = 0;
  cp_mem->cp_netfQ = 0;

  /* Quadrature integration turned ON */
  cp_mem->cp_quadr = TRUE;
  cp_mem->cp_quadMallocDone = TRUE;

  /* Quadrature initialization was successfull */
  return(CP_SUCCESS);
}

/*-----------------------------------------------------------------*/

#define lrw1Q (cp_mem->cp_lrw1Q)
#define liw1Q (cp_mem->cp_liw1Q)

/*-----------------------------------------------------------------*/

/*
 * CPodeQuadReInit
 *
 * CPodeQuadReInit re-initializes CPODES's quadrature related memory 
 * for a problem, assuming it has already been allocated in prior 
 * calls to CPodeMalloc and CPodeQuadInit. 
 * All problem specification inputs are checked for errors.
 * If any error occurs during initialization, it is reported to the
 * error handler function.
 * The return value is CP_SUCCESS = 0 if no errors occurred, or
 * a negative value otherwise.
 */

int CPodeQuadReInit(void *cpode_mem, CPQuadFn qfun, void *q_data, N_Vector q0)
{
  CPodeMem cp_mem;

  /* Check cpode_mem */
  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeQuadReInit", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  /* Ckeck if quadrature was initialized? */
  if (cp_mem->cp_quadMallocDone == FALSE) {
    cpProcessError(cp_mem, CP_NO_QUAD, "CPODES", "CPodeQuadReInit", MSGCP_NO_QUAD);
    return(CP_NO_QUAD);
  }

  /* Initialize znQ[0] in the history array */
  N_VScale(ONE, q0, cp_mem->cp_znQ[0]);

  /* Copy the input parameters into CPODES state */
  cp_mem->cp_qfun   = qfun;
  cp_mem->cp_q_data = q_data;

  /* Initialize counters */
  cp_mem->cp_nqe   = 0;
  cp_mem->cp_netfQ = 0;

  /* Quadrature integration turned ON */
  cp_mem->cp_quadr = TRUE;

  /* Quadrature re-initialization was successfull */
  return(CP_SUCCESS);
}

/*
 * =================================================================
 * READIBILITY REPLACEMENTS
 * =================================================================
 */

#define fi             (cp_mem->cp_fi)
#define fe             (cp_mem->cp_fe)
#define f_data         (cp_mem->cp_f_data)
#define efun           (cp_mem->cp_efun)
#define e_data         (cp_mem->cp_e_data)
#define qmax           (cp_mem->cp_qmax)
#define mxstep         (cp_mem->cp_mxstep)
#define mxhnil         (cp_mem->cp_mxhnil)
#define sldeton        (cp_mem->cp_sldeton)
#define hin            (cp_mem->cp_hin)
#define hmin           (cp_mem->cp_hmin)
#define hmax_inv       (cp_mem->cp_hmax_inv)
#define istop          (cp_mem->cp_istop)
#define tstop          (cp_mem->cp_tstop)
#define tstopset       (cp_mem->cp_tstopset)
#define maxncf         (cp_mem->cp_maxncf)
#define maxnef         (cp_mem->cp_maxnef)
#define nlscoef        (cp_mem->cp_nlscoef)
#define tol_type       (cp_mem->cp_tol_type)
#define reltol         (cp_mem->cp_reltol)
#define Sabstol        (cp_mem->cp_Sabstol)
#define Vabstol        (cp_mem->cp_Vabstol)

#define uround         (cp_mem->cp_uround)
#define zn             (cp_mem->cp_zn)
#define ewt            (cp_mem->cp_ewt)
#define y              (cp_mem->cp_y)
#define yp             (cp_mem->cp_yp)
#define acor           (cp_mem->cp_acor)
#define tempv          (cp_mem->cp_tempv)
#define ftemp          (cp_mem->cp_ftemp)
#define q              (cp_mem->cp_q)
#define qprime         (cp_mem->cp_qprime)
#define next_q         (cp_mem->cp_next_q)
#define qwait          (cp_mem->cp_qwait)
#define L              (cp_mem->cp_L)
#define h              (cp_mem->cp_h)
#define hprime         (cp_mem->cp_hprime)
#define next_h         (cp_mem->cp_next_h)
#define eta            (cp_mem->cp_eta)
#define etaqm1         (cp_mem->cp_etaqm1)
#define etaq           (cp_mem->cp_etaq)
#define etaqp1         (cp_mem->cp_etaqp1)
#define nscon          (cp_mem->cp_nscon)
#define hscale         (cp_mem->cp_hscale)
#define tn             (cp_mem->cp_tn)
#define tau            (cp_mem->cp_tau)
#define tq             (cp_mem->cp_tq)
#define l              (cp_mem->cp_l)
#define p              (cp_mem->cp_p)
#define rl1            (cp_mem->cp_rl1)
#define gamma          (cp_mem->cp_gamma)
#define gammap         (cp_mem->cp_gammap)
#define gamrat         (cp_mem->cp_gamrat)
#define acnrm          (cp_mem->cp_acnrm)
#define etamax         (cp_mem->cp_etamax)
#define nst            (cp_mem->cp_nst)
#define nfe            (cp_mem->cp_nfe)
#define ncfn           (cp_mem->cp_ncfn)
#define netf           (cp_mem->cp_netf)
#define nhnil          (cp_mem->cp_nhnil)
#define qu             (cp_mem->cp_qu)
#define hu             (cp_mem->cp_hu)
#define saved_tq5      (cp_mem->cp_saved_tq5)
#define indx_acor      (cp_mem->cp_indx_acor)
#define tolsf          (cp_mem->cp_tolsf)
#define lsetup_exists  (cp_mem->cp_lsetup_exists)
#define nor            (cp_mem->cp_nor)
#define ssdat          (cp_mem->cp_ssdat)

#define proj_enabled   (cp_mem->cp_proj_enabled)
#define proj_freq      (cp_mem->cp_proj_freq)
#define proj_type      (cp_mem->cp_proj_type)
#define cnstr_type     (cp_mem->cp_cnstr_type)
#define cfun           (cp_mem->cp_cfun)
#define c_data         (cp_mem->cp_c_data)
#define pfun           (cp_mem->cp_pfun)
#define p_data         (cp_mem->cp_p_data)
#define yC             (cp_mem->cp_yC)
#define acorP          (cp_mem->cp_acorP)
#define errP           (cp_mem->cp_errP)
#define ctol           (cp_mem->cp_ctol)
#define ctemp          (cp_mem->cp_ctemp)
#define tempvP1        (cp_mem->cp_tempvP1)
#define tempvP2        (cp_mem->cp_tempvP2)
#define nstlprj        (cp_mem->cp_nstlprj)
#define applyProj      (cp_mem->cp_applyProj)

#define quadr          (cp_mem->cp_quadr)
#define qfun           (cp_mem->cp_qfun)
#define q_data         (cp_mem->cp_q_data)
#define errconQ        (cp_mem->cp_errconQ)
#define tol_typeQ      (cp_mem->cp_tol_typeQ)
#define reltolQ        (cp_mem->cp_reltolQ)
#define SabstolQ       (cp_mem->cp_SabstolQ)
#define VabstolQ       (cp_mem->cp_VabstolQ)
#define znQ            (cp_mem->cp_znQ)
#define ewtQ           (cp_mem->cp_ewtQ)
#define acorQ          (cp_mem->cp_acorQ)
#define yQ             (cp_mem->cp_yQ)
#define tempvQ         (cp_mem->cp_tempvQ)
#define acnrmQ         (cp_mem->cp_acnrmQ)
#define nqe            (cp_mem->cp_nqe)
#define netfQ          (cp_mem->cp_netfQ)
#define quadMallocDone (cp_mem->cp_quadMallocDone)

#define doRootfinding  (cp_mem->cp_doRootfinding)
#define tlo            (cp_mem->cp_tlo)
#define tretlast       (cp_mem->cp_tretlast)
#define toutc          (cp_mem->cp_toutc)
#define taskc          (cp_mem->cp_taskc)
#define irfnd          (cp_mem->cp_irfnd)

/*-----------------------------------------------------------------*/

/*
 * CPode
 *
 * This routine is the main driver of the CPODES package. 
 *
 * It integrates over a time interval defined by the user, by calling
 * cpStep to do internal time steps.
 *
 * The first time that CPode is called for a successfully initialized
 * problem, it computes a tentative initial step size h.
 *
 * CPode supports four modes, specified by mode: CP_NORMAL, CP_ONE_STEP,
 * CP_NORMAL_TSTOP, and CP_ONE_STEP_TSTOP.
 * In the CP_NORMAL mode, the solver steps until it reaches or passes tout
 * and then interpolates to obtain y(tout).
 * In the CP_ONE_STEP mode, it takes one internal step and returns.
 * CP_NORMAL_TSTOP and CP_ONE_STEP_TSTOP are similar to CP_NORMAL and CP_ONE_STEP,
 * respectively, but the integration never proceeds past tstop (which
 * must have been defined through a call to CPodeSetStopTime).
 */

int CPode(void *cpode_mem, realtype tout, realtype *tret,
          N_Vector yout, N_Vector ypout, int mode)
{
  CPodeMem cp_mem;
  long int nstloc;
  int retval, hflag, kflag, istate, ier, task, irfndp;
  realtype troundoff, rh, nrm;

  /*
   * -------------------------------------
   * 1. Check and process inputs
   * -------------------------------------
   */

  /* Check if cpode_mem exists */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPode", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  /* Check if cpode_mem was allocated */
  if (cp_mem->cp_MallocDone == FALSE) {
    cpProcessError(cp_mem, CP_NO_MALLOC, "CPODES", "CPode", MSGCP_NO_MALLOC);
    return(CP_NO_MALLOC);
  }
  
  /* Check for yout != NULL */
  if (yout == NULL) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPode", MSGCP_YOUT_NULL);
    return(CP_ILL_INPUT);
  }
  y = yout;

  /* Check for ypout != NULL */
  if (ypout == NULL ) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPode", MSGCP_YPOUT_NULL);
    return(CP_ILL_INPUT);
  }
  yp = ypout;

  /* Check for tret != NULL */
  if (tret == NULL) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPode", MSGCP_TRET_NULL);
    return(CP_ILL_INPUT);
  }

  /* Check for valid mode */
  if ((mode != CP_NORMAL)       && 
      (mode != CP_ONE_STEP)     &&
      (mode != CP_NORMAL_TSTOP) &&
      (mode != CP_ONE_STEP_TSTOP) ) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPode", MSGCP_BAD_MODE);
    return(CP_ILL_INPUT);
  }

  /* Split mode into task and istop */
  if ((mode == CP_NORMAL_TSTOP) || (mode == CP_ONE_STEP_TSTOP)) {
    if ( tstopset == FALSE ) {
      cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPode", MSGCP_NO_TSTOP);
      return(CP_ILL_INPUT);
    }
    istop = TRUE;
  } else {
    istop = FALSE;
  }
  if ((mode == CP_NORMAL) || (mode == CP_NORMAL_TSTOP)) {
    task = CP_NORMAL; toutc = tout;
  } else {
    task = CP_ONE_STEP;
  }
  taskc = task;

  /*
   * ----------------------------------------
   * 2. Initializations performed only at
   *    the first step (nst=0):
   *    - initial setup
   *    - compute initial step size
   *    - check for approach to tstop
   *    - check for approach to a root
   *    - scale zn[1] by h
   * ----------------------------------------
   */

  if (nst == 0) {

    /* Initial setup */
    ier = cpInitialSetup(cp_mem);
    if (ier!= CP_SUCCESS) return(ier);
    
    /* In CP_EXPL mode, call fun at (t0, y0) to set zn[1] = f'(t0,y0) */
    if (ode_type == CP_EXPL) {
      retval = fe(tn, zn[0], zn[1], f_data); 
      nfe++;
      if (retval < 0) {
        cpProcessError(cp_mem, CP_ODEFUNC_FAIL, "CPODES", "CPode", MSGCP_ODEFUNC_FAILED, tn);
        return(CP_ODEFUNC_FAIL);
      }
      if (retval > 0) {
        cpProcessError(cp_mem, CP_FIRST_ODEFUNC_ERR, "CPODES", "CPode", MSGCP_ODEFUNC_FIRST);
        return(CP_FIRST_ODEFUNC_ERR);
      }
    }

    /* If computing any quadratures, call qfun at (t0,y0) to set znQ[1] = q'(t0,y0) */
    if (quadr) {
      retval = qfun(tn, zn[0], znQ[1], q_data);
      nqe++;
      if (retval < 0) {
        cpProcessError(cp_mem, CP_QUADFUNC_FAIL, "CPODES", "CPode", MSGCP_QUADFUNC_FAILED, tn);
        return(CP_QUADFUNC_FAIL);
      }
      if (retval > 0) {
        cpProcessError(cp_mem, CP_FIRST_QUADFUNC_ERR, "CPODES", "CPode", MSGCP_QUADFUNC_FIRST);
        return(CP_FIRST_QUADFUNC_ERR);
      }
    }

    /* Set initial h (from H0 or cpHin). */
    h = hin;
    if ( (h != ZERO) && ((tout-tn)*h < ZERO) ) {
      cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPode", MSGCP_BAD_H0);
      return(CP_ILL_INPUT);
    }
    if (h == ZERO) {
      hflag = cpHin(cp_mem, tout);
      if (hflag != CP_SUCCESS) {
        istate = cpHandleFailure(cp_mem, hflag);
        return(istate);
      }
    }
    rh = ABS(h)*hmax_inv;
    if (rh > ONE) h /= rh;
    if (ABS(h) < hmin) h *= hmin/ABS(h);

    /* Check for approach to tstop */
    if (istop) {
      if ( (tstop - tn)*h < ZERO ) {
        cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPode", MSGCP_BAD_TSTOP, tn);
        return(CP_ILL_INPUT);
      }
      if ( (tn + h - tstop)*h > ZERO ) 
        h = (tstop - tn)*(ONE-FOUR*uround);
    }

    cp_mem->cp_h0u = h;

    /* Check for zeros of root function g at and near t0. */
    if (doRootfinding) {
      retval = cpRcheck1(cp_mem);
      if (retval == CP_RTFUNC_FAIL) {
        cpProcessError(cp_mem, CP_RTFUNC_FAIL, "CPODES", "cpRcheck1", MSGCP_RTFUNC_FAILED, tn);
        return(CP_RTFUNC_FAIL);
      }
    }

    /* Scale zn[1] and (if needed) znQ[1] by h */
    hscale = h; 
    hprime = h;
    N_VScale(h, zn[1], zn[1]);
    if (quadr) N_VScale(h, znQ[1], znQ[1]);

  }

  /*
   * ------------------------------------------------------
   * 3. At following steps, perform stop tests:
   *    - check for root in last step
   *    - check if we passed tstop
   *    - check if we passed tout (NORMAL mode)
   *    - check if current tn was returned (ONE_STEP mode)
   *    - check if we are close to tstop
   *      (adjust step size if needed)
   * -------------------------------------------------------
   */

  if (nst > 0) {

    /* Estimate an infinitesimal time interval to be used as
       a roundoff for time quantities (based on current time 
       and step size) */
    troundoff = FUZZ_FACTOR*uround*(ABS(tn) + ABS(h));

    /* First, check for a root in the last step taken, other than the
       last root found, if any.  If task = CP_ONE_STEP and y(tn) was not
       returned because of an intervening root, return y(tn) now.     */
    if (doRootfinding) {

      irfndp = irfnd;
      
      retval = cpRcheck2(cp_mem);

      if (retval == CP_RTFUNC_FAIL) {
        cpProcessError(cp_mem, CP_RTFUNC_FAIL, "CPODES", "cpRcheck2", MSGCP_RTFUNC_FAILED, tlo);
        return(CP_RTFUNC_FAIL);
      } else if (retval == RTFOUND) {
        tretlast = *tret = tlo;
        return(CP_ROOT_RETURN);
      }

      /* If tn is distinct from tretlast (within roundoff),
         check remaining interval for roots */
      if ( ABS(tn - tretlast) > troundoff ) {

        retval = cpRcheck3(cp_mem);

        if (retval == CP_SUCCESS) {     /* no root found */
          irfnd = 0;
          if ((irfndp == 1) && (task == CP_ONE_STEP)) {
            tretlast = *tret = tn;
            cpGetSolution(cp_mem, tn, yout, ypout);
            return(CP_SUCCESS);
          }
        } else if (retval == RTFOUND) {  /* a new root was found */
          irfnd = 1;
          tretlast = *tret = tlo;
          return(CP_ROOT_RETURN);
        } else if (retval == CP_RTFUNC_FAIL) {  /* g failed */
          cpProcessError(cp_mem, CP_RTFUNC_FAIL, "CPODES", "cpRcheck3", MSGCP_RTFUNC_FAILED, tlo);
          return(CP_RTFUNC_FAIL);
        }

      }

    } /* end of root stop check */

    /* Test for tn past tstop */
    if ( istop && ((tstop - tn)*h < ZERO) ) {
      cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPode", MSGCP_BAD_TSTOP, tn);
      return(CP_ILL_INPUT);
    }

    /* In CP_NORMAL mode, test if tout was reached */
    if ( (task == CP_NORMAL) && ((tn-tout)*h >= ZERO) ) {
      tretlast = *tret = tout;
      ier =  cpGetSolution(cp_mem, tout, yout, ypout);
      if (ier != CP_SUCCESS) {
        cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPode", MSGCP_BAD_TOUT, tout);
        return(CP_ILL_INPUT);
      }
      return(CP_SUCCESS);
    }

    /* In CP_ONE_STEP mode, test if tn was returned */
    if ( task == CP_ONE_STEP && ABS(tn - tretlast) > troundoff ) {
      tretlast = *tret = tn;
      cpGetSolution(cp_mem, tn, yout, ypout);
      return(CP_SUCCESS);
    }

    /* Test for tn at tstop or near tstop */
    if ( istop ) {

      if ( ABS(tn - tstop) <= troundoff) {
        ier =  cpGetSolution(cp_mem, tstop, yout, ypout);
        if (ier != CP_SUCCESS) {
          cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPode", MSGCP_BAD_TSTOP, tn);
          return(CP_ILL_INPUT);
        }
        tretlast = *tret = tstop;
        return(CP_TSTOP_RETURN);
      }
      
      /* If next step would overtake tstop, adjust stepsize */
      if ( (tn + hprime - tstop)*h > ZERO ) {
        hprime = (tstop - tn)*(ONE-FOUR*uround);
        eta = hprime/h;
      }

    }
    
  }

  /*
   * --------------------------------------------------
   * 4. Looping point for internal steps
   *
   *    4.1. check for errors (too many steps, too much
   *         accuracy requested, step size too small)
   *    4.2. take a new step (call cpStep)
   *    4.3. stop on error 
   *    4.4. perform stop tests:
   *         - check for root in last step
   *         - check if tout was passed
   *         - check if close to tstop
   *         - check if in ONE_STEP mode (must return)
   * --------------------------------------------------
   */

  nstloc = 0;
  loop {
   
    next_h = h;
    next_q = q;
    
    /* Reset and check ewt, ewtQ */
    if (nst > 0) {

      ier = efun(zn[0], ewt, e_data);
      if (ier != 0) {
        if (tol_type == CP_WF) cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPode", MSGCP_EWT_NOW_FAIL, tn);
        else                   cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPode", MSGCP_EWT_NOW_BAD, tn);
        istate = CP_ILL_INPUT;
        tretlast = *tret = tn;
        cpGetSolution(cp_mem, tn, yout, ypout);
        break;
      }

      if (quadr && errconQ) {
        ier = cpQuadEwtSet(cp_mem, znQ[0], ewtQ);
        if(ier != 0) {
          cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPode", MSGCP_EWTQ_NOW_BAD, tn);
          istate = CP_ILL_INPUT;
          tretlast = *tret = tn;
          N_VScale(ONE, zn[0], yout);
          break;
        }
      }

    }
    
    /* Check for too many steps */
    if (nstloc >= mxstep) {
      cpProcessError(cp_mem, CP_TOO_MUCH_WORK, "CPODES", "CPode", MSGCP_MAX_STEPS, tn);
      istate = CP_TOO_MUCH_WORK;
      tretlast = *tret = tn;
      cpGetSolution(cp_mem, tn, yout, ypout);
      break;
    }

    /* Check for too much accuracy requested */
    nrm = N_VWrmsNorm(zn[0], ewt);
    if (quadr && errconQ) nrm = cpQuadUpdateNorm(cp_mem, nrm, znQ[0], ewtQ); 
    tolsf = uround * nrm;
    if (tolsf > ONE) {
      cpProcessError(cp_mem, CP_TOO_MUCH_ACC, "CPODES", "CPode", MSGCP_TOO_MUCH_ACC, tn);
      istate = CP_TOO_MUCH_ACC;
      tretlast = *tret = tn;
      cpGetSolution(cp_mem, tn, yout, ypout);
      tolsf *= TWO;
      break;
    } else {
      tolsf = ONE;
    }

    /* Check for h below roundoff level in tn */
    if (tn + h == tn) {
      nhnil++;
      if (nhnil <= mxhnil) 
        cpProcessError(cp_mem, CP_WARNING, "CPODES", "CPode", MSGCP_HNIL, tn, h);
      if (nhnil == mxhnil) 
        cpProcessError(cp_mem, CP_WARNING, "CPODES", "CPode", MSGCP_HNIL_DONE);
    }

    /* Call cpStep to take a step */
    kflag = cpStep(cp_mem);

    /* Process failed step cases, and exit loop */
    if (kflag != CP_SUCCESS) {
      istate = cpHandleFailure(cp_mem, kflag);
      tretlast = *tret = tn;
      cpGetSolution(cp_mem, tn, yout, ypout);
      break;
    }
    
    nstloc++;

    /* Check for root in last step taken. */
    if (doRootfinding) {

      retval = cpRcheck3(cp_mem);

      if (retval == RTFOUND) {  /* A new root was found */
        irfnd = 1;
        istate = CP_ROOT_RETURN;
        tretlast = *tret = tlo;
        break;
      } else if (retval == CP_RTFUNC_FAIL) { /* g failed */
        cpProcessError(cp_mem, CP_RTFUNC_FAIL, "CPODES", "cpRcheck3", MSGCP_RTFUNC_FAILED, tlo);
        istate = CP_RTFUNC_FAIL;
        break;
      }

    }

    /* In NORMAL mode, check if tout reached */
    if ( (task == CP_NORMAL) &&  (tn-tout)*h >= ZERO ) {
      istate = CP_SUCCESS;
      tretlast = *tret = tout;
      (void) cpGetSolution(cp_mem, tout, yout, ypout);
      next_q = qprime;
      next_h = hprime;
      break;
    }

    /* Check if tn is at tstop or near tstop */
    if ( istop ) {

      troundoff = FUZZ_FACTOR*uround*(ABS(tn) + ABS(h));
      if ( ABS(tn - tstop) <= troundoff) {
        (void) cpGetSolution(cp_mem, tstop, yout, ypout);
        tretlast = *tret = tstop;
        istate = CP_TSTOP_RETURN;
        break;
      }

      if ( (tn + hprime - tstop)*h > ZERO ) {
        hprime = (tstop - tn)*(ONE-FOUR*uround);
        eta = hprime/h;
      }

    }

    /* In ONE_STEP mode, copy y and exit loop */
    if (task == CP_ONE_STEP) {
      istate = CP_SUCCESS;
      tretlast = *tret = tn;
      cpGetSolution(cp_mem, tn, yout, ypout);
      next_q = qprime;
      next_h = hprime;
      break;
    }

  } /* end looping for internal steps */

  return(istate);
}

/*-----------------------------------------------------------------*/

/*
 * CPodeGetDky
 *
 * This routine computes the k-th derivative of the interpolating
 * polynomial at the time t and stores the result in the vector dky.
 * The formula is:
 *         q 
 *  dky = SUM c(j,k) * (t - tn)^(j-k) * h^(-j) * zn[j] , 
 *        j=k 
 * where c(j,k) = j*(j-1)*...*(j-k+1), q is the current order, and
 * zn[j] is the j-th column of the Nordsieck history array.
 *
 */

int CPodeGetDky(void *cpode_mem, realtype t, int k, N_Vector dky)
{
  realtype s, c, r;
  realtype tfuzz, tp, tn1;
  int i, j;
  CPodeMem cp_mem;
  
  /* Check all inputs for legality */
 
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeGetDky", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (dky == NULL) {
    cpProcessError(cp_mem, CP_BAD_DKY, "CPODES", "CPodeGetDky", MSGCP_NULL_DKY);
    return(CP_BAD_DKY);
  }

  if ((k < 0) || (k > q)) {
    cpProcessError(cp_mem, CP_BAD_K, "CPODES", "CPodeGetDky", MSGCP_BAD_K);
    return(CP_BAD_K);
  }
  
  /* Allow for some slack */
  tfuzz = FUZZ_FACTOR * uround * (ABS(tn) + ABS(hu));
  if (hu < ZERO) tfuzz = -tfuzz;
  tp = tn - hu - tfuzz;
  tn1 = tn + tfuzz;
  if ((t-tp)*(t-tn1) > ZERO) {
    cpProcessError(cp_mem, CP_BAD_T, "CPODES", "CPodeGetDky", MSGCP_BAD_T, t, tn-hu, tn);
    return(CP_BAD_T);
  }

  /* Sum the differentiated interpolating polynomial */

  s = (t - tn) / h;
  for (j=q; j >= k; j--) {
    c = ONE;
    for (i=j; i >= j-k+1; i--) c *= i;
    if (j == q) {
      N_VScale(c, zn[q], dky);
    } else {
      N_VLinearSum(c, zn[j], s, dky, dky);
    }
  }
  if (k == 0) return(CP_SUCCESS);
  r = RPowerI(h,-k);
  N_VScale(r, dky, dky);
  return(CP_SUCCESS);
}

/* 
 * CPodeGetQuad
 *
 * This routine extracts quadrature solution into yQout.
 * This is just a wrapper that calls CPodeGetQuadDky with k=0                    
 */
 
int CPodeGetQuad(void *cpode_mem, realtype t, N_Vector yQout)
{
  return(CPodeGetQuadDky(cpode_mem,t,0,yQout));
}

/*
 * CPodeGetQuadDky
 *
 * CPodeQuadDky computes the kth derivative of the yQ function at
 * time t, where tn-hu <= t <= tn, tn denotes the current         
 * internal time reached, and hu is the last internal step size   
 * successfully used by the solver. The user may request 
 * k=0, 1, ..., qu, where qu is the current order. 
 * The derivative vector is returned in dky. This vector 
 * must be allocated by the caller. It is only legal to call this         
 * function after a successful return from CPode with quadrature
 * computation enabled.
 */

int CPodeGetQuadDky(void *cpode_mem, realtype t, int k, N_Vector dkyQ)
{ 
  realtype s, c, r;
  realtype tfuzz, tp, tn1;
  int i, j;
  CPodeMem cp_mem;
  
  /* Check all inputs for legality */
  
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeGetQuadDky", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;  

  if(quadr != TRUE) {
    cpProcessError(cp_mem, CP_NO_QUAD, "CPODES", "CPodeGetQuadDky", MSGCP_NO_QUAD);
    return(CP_NO_QUAD);
  }

  if (dkyQ == NULL) {
    cpProcessError(cp_mem, CP_BAD_DKY, "CPODES", "CPodeGetQuadDky", MSGCP_NULL_DKY);
    return(CP_BAD_DKY);
  }
  
  if ((k < 0) || (k > q)) {
    cpProcessError(cp_mem, CP_BAD_K, "CPODES", "CPodeGetQuadDky", MSGCP_BAD_K);
    return(CP_BAD_K);
  }
  
  /* Allow for some slack */
  tfuzz = FUZZ_FACTOR * uround * (ABS(tn) + ABS(hu));
  if (hu < ZERO) tfuzz = -tfuzz;
  tp = tn - hu - tfuzz;
  tn1 = tn + tfuzz;
  if ((t-tp)*(t-tn1) > ZERO) {
    cpProcessError(cp_mem, CP_BAD_T, "CPODES", "CPodeGetQuadDky", MSGCP_BAD_T);
    return(CP_BAD_T);
  }
  
  /* Sum the differentiated interpolating polynomial */
  
  s = (t - tn) / h;
  for (j=q; j >= k; j--) {
    c = ONE;
    for (i=j; i >= j-k+1; i--) c *= i;
    if (j == q) {
      N_VScale(c, znQ[q], dkyQ);
    } else {
      N_VLinearSum(c, znQ[j], s, dkyQ, dkyQ);
    }
  }
  if (k == 0) return(CP_SUCCESS);
  r = RPowerI(h,-k);
  N_VScale(r, dkyQ, dkyQ);
  return(CP_SUCCESS);
  
}

/*
 * CPodeFree
 *
 * This routine frees the problem memory allocated by CPodeInit.
 * Such memory includes all the vectors allocated by cpAllocVectors,
 * and the memory lmem for the linear solver (deallocated by a call
 * to lfree).
 */

void CPodeFree(void **cpode_mem)
{
  CPodeMem cp_mem;

  if (*cpode_mem == NULL) return;

  cp_mem = (CPodeMem) (*cpode_mem);
  
  cpFreeVectors(cp_mem);

  CPodeQuadFree(cp_mem);

  if (nls_type == CP_NEWTON && cp_mem->cp_lfree != NULL) 
    cp_mem->cp_lfree(cp_mem);

  if (cp_mem->cp_lfreeP != NULL)
    cp_mem->cp_lfreeP(cp_mem);

  if (cp_mem->cp_rootMallocDone) cpRootFree(cp_mem);

  if (cp_mem->cp_projMallocDone) cpProjFree(cp_mem);

  free(*cpode_mem);
  *cpode_mem = NULL;
}

/*
 * CPodeQuadFree
 *
 * CPodeQuadFree frees the problem memory in cpode_mem allocated
 * for quadrature integration. Its only argument is the pointer
 * cpode_mem returned by CPodeCreate. 
 */

void CPodeQuadFree(void *cpode_mem)
{
  CPodeMem cp_mem;
  
  if (cpode_mem == NULL) return;
  cp_mem = (CPodeMem) cpode_mem;

  if(quadMallocDone) {
    cpQuadFree(cp_mem);
    quadMallocDone = FALSE;
    quadr = FALSE;
  }
}


/* 
 * =================================================================
 *  PRIVATE FUNCTIONS
 * =================================================================
 */

/* 
 * -----------------------------------------------------------------
 * Memory allocation/deallocation and initialization functions
 * -----------------------------------------------------------------
 */

/*
 * cpCheckNvector
 * This routine checks if all required vector operations are present.
 * If any of them is missing it returns FALSE.
 */

static booleantype cpCheckNvector(N_Vector tmpl)
{
  if((tmpl->ops->nvclone     == NULL) ||
     (tmpl->ops->nvdestroy   == NULL) ||
     (tmpl->ops->nvlinearsum == NULL) ||
     (tmpl->ops->nvconst     == NULL) ||
     (tmpl->ops->nvprod      == NULL) ||
     (tmpl->ops->nvdiv       == NULL) ||
     (tmpl->ops->nvscale     == NULL) ||
     (tmpl->ops->nvabs       == NULL) ||
     (tmpl->ops->nvinv       == NULL) ||
     (tmpl->ops->nvaddconst  == NULL) ||
     (tmpl->ops->nvmaxnorm   == NULL) ||
     (tmpl->ops->nvwrmsnorm  == NULL) ||
     (tmpl->ops->nvmin       == NULL))
    return(FALSE);
  else
    return(TRUE);
}

/*
 * cpAllocVectors
 *
 * This routine allocates the CPODES vectors ewt, acor, tempv, ftemp,
 * and zn[0], ..., zn[maxord]. If tol_type=CP_SV, it also allocates 
 * space for Vabstol. If all memory allocations are successful, 
 * cpAllocVectors returns TRUE. Otherwise all allocated memory is 
 * freed and cpAllocVectors returns FALSE. This routine also sets the 
 * optional outputs lrw and liw, which are (respectively) the lengths 
 * of the real and integer work spaces allocated here.
 */

static booleantype cpAllocVectors(CPodeMem cp_mem, N_Vector tmpl, int tol)
{
  int i, j;

  /* Allocate ewt, acor, tempv, ftemp */
  
  ewt = N_VClone(tmpl);
  if (ewt == NULL) return(FALSE);

  acor = N_VClone(tmpl);
  if (acor == NULL) {
    N_VDestroy(ewt);
    return(FALSE);
  }

  tempv = N_VClone(tmpl);
  if (tempv == NULL) {
    N_VDestroy(ewt);
    N_VDestroy(acor);
    return(FALSE);
  }

  ftemp = N_VClone(tmpl);
  if (ftemp == NULL) {
    N_VDestroy(tempv);
    N_VDestroy(ewt);
    N_VDestroy(acor);
    return(FALSE);
  }

  /* Allocate zn[0] ... zn[qmax] */

  for (j=0; j <= qmax; j++) {
    zn[j] = N_VClone(tmpl);
    if (zn[j] == NULL) {
      N_VDestroy(ewt);
      N_VDestroy(acor);
      N_VDestroy(tempv);
      N_VDestroy(ftemp);
      for (i=0; i < j; i++) N_VDestroy(zn[i]);
      return(FALSE);
    }
  }

  /* Update solver workspace lengths  */
  lrw += (qmax + 5)*lrw1;
  liw += (qmax + 5)*liw1;

  if (tol == CP_SV) {
    Vabstol = N_VClone(tmpl);
    if (Vabstol == NULL) {
      N_VDestroy(ewt);
      N_VDestroy(acor);
      N_VDestroy(tempv);
      N_VDestroy(ftemp);
      for (i=0; i <= qmax; i++) N_VDestroy(zn[i]);
      return(FALSE);
    }
    lrw += lrw1;
    liw += liw1;
    cp_mem->cp_VabstolMallocDone = TRUE;
  }

  /* Store the value of qmax used here */
  cp_mem->cp_qmax_alloc = qmax;

  return(TRUE);
}

/*  
 * cpFreeVectors
 *
 * This routine frees the CPODES vectors allocated in cpAllocVectors.
 */

static void cpFreeVectors(CPodeMem cp_mem)
{
  int j, maxord;
  
  maxord = cp_mem->cp_qmax_alloc;

  N_VDestroy(ewt);
  N_VDestroy(acor);
  N_VDestroy(tempv);
  N_VDestroy(ftemp);
  for(j=0; j <= maxord; j++) N_VDestroy(zn[j]);

  lrw -= (maxord + 5)*lrw1;
  liw -= (maxord + 5)*liw1;

  if (cp_mem->cp_VabstolMallocDone) {
    N_VDestroy(Vabstol);
    lrw -= lrw1;
    liw -= liw1;
  }
}

/*
 * cpProjAlloc allocates memory for the internal projection functions.
 */

static booleantype cpProjAlloc(CPodeMem cp_mem, N_Vector c_tmpl, N_Vector s_tmpl)
{

  /* Vectors cloned from c_tmpl (length M) */

  ctol = N_VClone(c_tmpl);
  if (ctol == NULL) return(FALSE);

  ctemp = N_VClone(c_tmpl);
  if (ctemp == NULL) {
    N_VDestroy(ctol);
    return(FALSE);
  }

  tempvP1 = N_VClone(c_tmpl);
  if (tempvP1 == NULL) {
    N_VDestroy(ctol);
    N_VDestroy(ctemp);
    return(FALSE);
  }

  tempvP2 = N_VClone(c_tmpl);
  if (tempvP2 == NULL) {
    N_VDestroy(ctol);
    N_VDestroy(ctemp);
    N_VDestroy(tempvP1);
    return(FALSE);
  }

  /* Vectors cloned from s_tmpl (length N) */

  acorP = N_VClone(s_tmpl);
  if (acorP == NULL) {
    N_VDestroy(ctol);
    N_VDestroy(ctemp);
    N_VDestroy(tempvP1);
    N_VDestroy(tempvP2);
    return(FALSE);
  }

  yC = N_VClone(s_tmpl);
  if (yC == NULL) {
    N_VDestroy(ctol);
    N_VDestroy(ctemp);
    N_VDestroy(tempvP1);
    N_VDestroy(tempvP2);
    N_VDestroy(acorP);
    return(FALSE);
  }
  
  errP = N_VClone(s_tmpl);
  if (errP == NULL) {
    N_VDestroy(ctol);
    N_VDestroy(ctemp);
    N_VDestroy(tempvP1);
    N_VDestroy(tempvP2);
    N_VDestroy(acorP);
    N_VDestroy(yC);
    return(FALSE);
  }

  
  return(TRUE);
}

/*
 * cpProjFree frees the memory allocated in cpProjAlloc.
 */

static void cpProjFree(CPodeMem cp_mem)
{

  N_VDestroy(ctol);
  N_VDestroy(ctemp);
  N_VDestroy(tempvP1);
  N_VDestroy(tempvP2);
  N_VDestroy(acorP);
  N_VDestroy(yC);
  N_VDestroy(errP);

  lrw -= lrw1;
  liw -= liw1;

}

/*
 * cpQuadAlloc allocates memory for the quadrature integration
 *
 * NOTE: Space for ewtQ is allocated even when errconQ=FALSE, 
 * although in this case, ewtQ is never used. The reason for this
 * decision is to allow the user to re-initialize the quadrature
 * computation with errconQ=TRUE, after an initialization with
 * errconQ=FALSE, without new memory allocation within 
 * CPodeQuadReInit.
 */

static booleantype cpQuadAlloc(CPodeMem cp_mem, N_Vector q_tmpl)
{
  int i, j;

  /* Allocate ewtQ */
  ewtQ = N_VClone(q_tmpl);
  if (ewtQ == NULL) {
    return(FALSE);
  }
  
  /* Allocate acorQ */
  acorQ = N_VClone(q_tmpl);
  if (acorQ == NULL) {
    N_VDestroy(ewtQ);
    return(FALSE);
  }

  /* Allocate yQ */
  yQ = N_VClone(q_tmpl);
  if (yQ == NULL) {
    N_VDestroy(ewtQ);
    N_VDestroy(acorQ);
    return(FALSE);
  }

  /* Allocate tempvQ */
  tempvQ = N_VClone(q_tmpl);
  if (tempvQ == NULL) {
    N_VDestroy(ewtQ);
    N_VDestroy(acorQ);
    N_VDestroy(yQ);
    return(FALSE);
  }

  /* Allocate zQn[0] ... zQn[maxord] */

  for (j=0; j <= qmax; j++) {
    znQ[j] = N_VClone(q_tmpl);
    if (znQ[j] == NULL) {
      N_VDestroy(ewtQ);
      N_VDestroy(acorQ);
      N_VDestroy(yQ);
      N_VDestroy(tempvQ);
      for (i=0; i < j; i++) N_VDestroy(znQ[i]);
      return(FALSE);
    }
  }

  /* Store the value of qmax used here */
  cp_mem->cp_qmax_allocQ = qmax;

  /* Update solver workspace lengths */
  lrw += (qmax + 5)*lrw1Q;
  liw += (qmax + 5)*liw1Q;

  return(TRUE);
}

/*
 * cpQuadFree frees the memory allocated in cpQuadAlloc.
 */

static void cpQuadFree(CPodeMem cp_mem)
{
  int j, maxord;
  
  maxord = cp_mem->cp_qmax_allocQ;

  N_VDestroy(ewtQ);
  N_VDestroy(acorQ);
  N_VDestroy(yQ);
  N_VDestroy(tempvQ);
  
  for (j=0; j<=maxord; j++) N_VDestroy(znQ[j]);

  lrw -= (maxord + 5)*lrw1Q;
  liw -= (maxord + 5)*liw1Q;

  if (cp_mem->cp_VabstolQMallocDone) {
    N_VDestroy(VabstolQ);
    lrw -= lrw1Q;
    liw -= liw1Q;
  }

  cp_mem->cp_VabstolQMallocDone = FALSE;
}

/*  
 * cpInitialSetup
 *
 * This routine performs input consistency checks at the first step.
 * If needed, it also checks the linear solver module and calls the
 * linear solver initialization routine.
 */

static int cpInitialSetup(CPodeMem cp_mem)
{
  int ier;

  /* Did the user provide efun? */
  if (tol_type != CP_WF) {
    efun = cpEwtSet;
    e_data = (void *)cp_mem;
  } else {
    if (efun == NULL) {
      cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "cpInitialSetup", MSGCP_NO_EFUN);
      return(CP_ILL_INPUT);
    }
  }

  /* Evaluate error weights at initial time */
  ier = efun(zn[0], ewt, e_data);
  if (ier != 0) {
    if (tol_type == CP_WF) cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "cpInitialSetup", MSGCP_EWT_FAIL);
    else                   cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "cpInitialSetup", MSGCP_BAD_EWT);
    return(CP_ILL_INPUT);
  }
  
  /* Evaluate quadrature error weights */
  if (quadr && errconQ) {
    ier = cpQuadEwtSet(cp_mem, znQ[0], ewtQ);
    if (ier != 0) {
      cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPInitialSetup", MSGCP_BAD_EWTQ);
      return(CP_ILL_INPUT);
    }
  }

  if (!quadr) errconQ = FALSE;

  /* Check if lsolve function exists (if needed)
     and call linit function (if it exists) */
  if (nls_type == CP_NEWTON) {
    if (cp_mem->cp_lsolve == NULL) {
      cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "cpInitialSetup", MSGCP_LSOLVE_NULL);
      return(CP_ILL_INPUT);
    }
    if (cp_mem->cp_linit != NULL) {
      ier = cp_mem->cp_linit(cp_mem);
      if (ier != 0) {
        cpProcessError(cp_mem, CP_LINIT_FAIL, "CPODES", "cpInitialSetup", MSGCP_LINIT_FAIL);
        return(CP_LINIT_FAIL);
      }
    }
  }

  /* Tests related to coordinate projection */
  if (proj_enabled) {

    switch (proj_type) {

    case CP_PROJ_USER:

      /* Check if user provided the projection function */
      if ( pfun == NULL ) {
        cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "cpInitialSetup", MSGCP_NO_PFUN);
        return(CP_ILL_INPUT);
      }

      break;

    case CP_PROJ_INTERNAL:

      /* Check if user provided the constraint function */
      if (cfun == NULL) {
        cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "cpInitialSetup", MSGCP_NO_CFUN);
        return(CP_ILL_INPUT);
      }
      /* Check lsolveP exists */ 
      if ( cp_mem->cp_lsolveP == NULL) {
        cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "cpInitialSetup", MSGCP_PLSOLVE_NULL);
        return(CP_ILL_INPUT);
      }
      /* Call linitP if it exists */
      if ( cp_mem->cp_linitP != NULL ) {
        ier = cp_mem->cp_linitP(cp_mem);
        if (ier != 0) {
          cpProcessError(cp_mem, CP_LINIT_FAIL, "CPODES", "cpInitialSetup", MSGCP_PLINIT_FAIL);
          return(CP_PLINIT_FAIL);
        }
      }

      break;

    } 

  }

  return(CP_SUCCESS);
}

/* 
 * -----------------------------------------------------------------
 * Initial step size evaluation functions
 * -----------------------------------------------------------------
 */

/*
 * cpHin
 *
 * This routine computes a tentative initial step size h0. 
 * If tout is too close to tn (= t0), then cpHin returns CP_TOO_CLOSE
 * and h remains uninitialized. 
 *
 */

static int cpHin(CPodeMem cp_mem, realtype tout)
{
  int sign, retval;
  realtype hlb, hub, hg, h0;
  realtype tdiff, tdist, tround;

  /* If tout is too close to tn, give up */
  if ((tdiff = tout-tn) == ZERO) return(CP_TOO_CLOSE);
  sign = (tdiff > ZERO) ? 1 : -1;
  tdist = ABS(tdiff);
  tround = uround * MAX(ABS(tn), ABS(tout));
  if (tdist < TWO*tround) return(CP_TOO_CLOSE);
  
  /* Set lower and upper bounds on h0*/
  hlb = HLB_FACTOR * tround;
  hub = cpUpperBoundH0(cp_mem, tdist);

  /* Compute geometric mean of lower and upper bounds */
  hg  = RSqrt(hlb*hub);

  /* If the bounds cross each other, use hg as initial step size */
  if (hub < hlb) {
    if (sign == -1) h = -hg;
    else            h =  hg;
    return(CP_SUCCESS);
  }

  /* Estimate initial step size */
  if (ode_type == CP_EXPL) retval = cpHinExpl(cp_mem, hlb, hub, sign, &h0);
  else                     retval = cpHinImpl(cp_mem, tdist, &h0);

  /* Apply bounds and attach sign */
  if (h0 < hlb) h0 = hlb;
  if (h0 > hub) h0 = hub;
  if (sign == -1) h0 = -h0;
  h = h0;

  return(retval);
}

/*
 * cpUpperBoundH0
 *
 * This routine sets an upper bound on abs(h0) based on
 * tdist = tn - t0 and the values of y[i]/y'[i] and yQ[i]/yQ'[i].
 */

static realtype cpUpperBoundH0(CPodeMem cp_mem, realtype tdist)
{
  realtype hub_inv, hubQ_inv, hub;
  N_Vector temp1, temp2;
  N_Vector tempQ1, tempQ2;

  /* 
   * Bound based on |y0|/|y0'| -- allow at most an increase of
   * HUB_FACTOR in y0 (based on a forward Euler step). The weight 
   * factor is used as a safeguard against zero components in y0. 
   */

  temp1 = tempv;
  temp2 = acor;

  N_VAbs(zn[0], temp2);
  efun(zn[0], temp1, e_data);
  N_VInv(temp1, temp1);
  N_VLinearSum(HUB_FACTOR, temp2, ONE, temp1, temp1);

  N_VAbs(zn[1], temp2);

  N_VDiv(temp2, temp1, temp1);
  hub_inv = N_VMaxNorm(temp1);

  /* Bound based on |yQ|/|yQ'| */
  
  if (quadr && errconQ) {

    tempQ1 = tempvQ;
    tempQ2 = acorQ;

    N_VAbs(znQ[0], tempQ2);
    cpQuadEwtSet(cp_mem, znQ[0], tempQ1);
    N_VInv(tempQ1, tempQ1);
    N_VLinearSum(HUB_FACTOR, tempQ2, ONE, tempQ1, tempQ1);
    
    N_VAbs(znQ[1], tempQ2);
    
    N_VDiv(tempQ2, tempQ1, tempQ1);
    hubQ_inv = N_VMaxNorm(tempQ1);

    if (hubQ_inv > hub_inv) hub_inv = hubQ_inv;

  }

  /*
   * bound based on tdist -- allow at most a step of magnitude
   * HUB_FACTOR * tdist
   */

  hub = HUB_FACTOR*tdist;

  /* Use the smaler of the two */

  if (hub*hub_inv > ONE) hub = ONE/hub_inv;

  return(hub);
}

/*
 * cpHinExpl
 *
 * This routine computes a tentative initial step size h0. 
 * If the ODE function fails unrecoverably, cpHinExpl returns CP_ODEFUNC_FAIL.
 * If the ODE function fails recoverably too many times and recovery is
 * not possible, cpHinExpl returns CP_REPTD_ODEFUNC_ERR.
 * Otherwise, cpHinExpl sets h0 to the chosen value and returns CP_SUCCESS.
 *
 * The algorithm used seeks to find h0 as a solution of
 *       (WRMS norm of (h0^2 ypp / 2)) = 1, 
 * where ypp = estimated second derivative of y.
 *
 * We start with an initial estimate equal to the geometric mean of the
 * lower and upper bounds on the step size.
 *
 * Loop up to H_MAXITERS times to find h0.
 * Stop if new and previous values differ by a factor < 2.
 * Stop if hnew/hg > 2 after one iteration, as this probably means
 * that the ypp value is bad because of cancellation error.        
 *  
 * For each new proposed hg, we allow H_MAXITERS attempts to
 * resolve a possible recoverable failure from f() by reducing
 * the proposed stepsize by a factor of 0.2. If a legal stepsize
 * still cannot be found, fall back on a previous value if possible,
 * or else return CP_REPTD_ODEFUNC_ERR.
 *
 * Finally, we apply a bias (0.5).
 */

static int cpHinExpl(CPodeMem cp_mem, realtype hlb, realtype hub, int sign, realtype *h0)
{
  int retval, count1, count2;
  realtype hg, hgs, hs, hnew, hrat, yppnorm;
  booleantype hgOK, hnewOK;

  /* Initial estimate = geometric mean of computed bounds */
  hg = RSqrt(hlb*hub);

  /* Outer loop */

  hnewOK = FALSE;
  hs = hg;         /* safeguard against 'uninitialized variable' warning */

  for(count1 = 1; count1 <= H_MAXITERS; count1++) {

    /* Attempts to estimate ypp */

    hgOK = FALSE;

    for (count2 = 1; count2 <= H_MAXITERS; count2++) {
      hgs = hg*sign;
      retval = cpYppNorm(cp_mem, hgs, &yppnorm);
      /* If fun() or qfun() failed unrecoverably, give up */
      if (retval < 0) return(retval);
      /* If successful, we can use ypp */
      if (retval == CP_SUCCESS) {hgOK = TRUE; break;}
      /* fun() or qfun() failed recoverably; cut step size and test it again */
      hg *= PT2;
    }

    /* If fun() or qfun() failed recoverably H_MAXITERS times */

    if (!hgOK) {
      /* Exit if this is the first or second pass. No recovery possible */
      if (count1 <= 2) {
        if (retval == ODEFUNC_RECVR)  return(CP_REPTD_ODEFUNC_ERR);
        if (retval == QUADFUNC_RECVR) return(CP_REPTD_QUADFUNC_ERR);
      }
      /* We have a fall-back option. The value hs is a previous hnew which
         passed through f(). Use it and break */
      hnew = hs;
      break;
    }

    /* The proposed step size is feasible. Save it. */
    hs = hg;

    /* If the stopping criteria was met, or if this is the last pass, stop */
    if ( (hnewOK) || (count1 == H_MAXITERS))  {hnew = hg; break;}

    /* Propose new step size */
    hnew = (yppnorm*hub*hub > TWO) ? RSqrt(TWO/yppnorm) : RSqrt(hg*hub);
    hrat = hnew/hg;
    
    /* Accept hnew if it does not differ from hg by more than a factor of 2 */
    if ((hrat > HALF) && (hrat < TWO)) {
      hnewOK = TRUE;
    }

    /* After one pass, if ypp seems to be bad, use fall-back value. */
    if ((count1 > 1) && (hrat > TWO)) {
      hnew = hg;
      hnewOK = TRUE;
    }

    /* Send this value back through f() */
    hg = hnew;

  }

  /* Apply bias factor */
  *h0 = H_BIAS*hnew;

  return(CP_SUCCESS);
}

/*
 * cpYppNorm
 *
 * When the ODE is given in explicit form, this routine computes 
 * an estimate of the second derivative of y using a difference 
 * quotient, and returns its WRMS norm.
 */

static int cpYppNorm(CPodeMem cp_mem, realtype hg, realtype *yppnorm)
{
  int retval;

  /* y <- h*y'(t) + y(t) */
  N_VLinearSum(hg, zn[1], ONE, zn[0], y);

  /* tempv <- fun(t+h, h*y'(t)+y(t)) */
  retval = fe(tn+hg, y, tempv, f_data);
  nfe++;
  if (retval < 0) return(CP_ODEFUNC_FAIL);
  if (retval > 0) return(ODEFUNC_RECVR);

  /* tempvQ <- qfun(t+h, h*y'(t)+y(t)) */
  if (quadr && errconQ) {
    retval = qfun(tn+hg, y, tempvQ, q_data);
    nqe++;
    if (retval < 0) return(CP_QUADFUNC_FAIL);
    if (retval > 0) return(QUADFUNC_RECVR);
  }

  /* tempv <- fun(t+h, h*y'(t)+y(t)) - y'(t) */
  /* tempv <- ypp */
  N_VLinearSum(ONE, tempv, -ONE, zn[1], tempv);
  N_VScale(ONE/hg, tempv, tempv);

  /* tempvQ <- qfun(t+h, h*y'(t)+y(t)) - yQ'(t) */
  /* tempvQ <- yQdd */
  if (quadr && errconQ) {
    N_VLinearSum(ONE, tempvQ, -ONE, znQ[1], tempvQ);
    N_VScale(ONE/hg, tempvQ, tempvQ);
  }

  /* Estimate ||y''|| */
  *yppnorm = N_VWrmsNorm(tempv, ewt);
  if (quadr && errconQ) *yppnorm = cpQuadUpdateNorm(cp_mem, *yppnorm, tempvQ, ewtQ);

  return(CP_SUCCESS);
}

/*
 * cpHinImpl
 * 
 * For ODE systems given in implicit form, use
 *    h0 = min ( 0.001 * tdist , 0.5 / ||yp(t0)||_wrms , 0.5 / ||yQp(t0)||_wrms )
 */

static int cpHinImpl(CPodeMem cp_mem, realtype tdist, realtype *h0)
{
  realtype ypnorm;

  *h0 = PT001 * tdist;

  ypnorm = N_VWrmsNorm(zn[1], ewt);
  if (quadr && errconQ) ypnorm = cpQuadUpdateNorm(cp_mem, ypnorm, znQ[1], ewtQ);

  if (ypnorm > HALF/(*h0)) *h0 = HALF/ypnorm;

  return(CP_SUCCESS);
}

/* 
 * -----------------------------------------------------------------
 * Main step function
 * -----------------------------------------------------------------
 */

/* 
 * cpStep
 *
 * This routine performs one internal cpode step, from tn to tn + h.
 * It calls other routines to do all the work.
 *
 * The main operations done here are as follows:
 * - preliminary adjustments if a new step size was chosen;
 * - prediction of the Nordsieck history array zn at tn + h;
 * - setting of multistep method coefficients and test quantities;
 * - solution of the nonlinear system;
 * - testing the local error;
 * - updating zn and other state data if successful;
 * - resetting stepsize and order for the next step.
 * - if SLDET is on, check for stability, reduce order if necessary.
 * On a failure in the nonlinear system solution or error test, the
 * step may be reattempted, depending on the nature of the failure.
 */

static int cpStep(CPodeMem cp_mem)
{
  realtype saved_t, dsm, dsmQ;
  int ncf, npf, nef;
  int nflag, kflag, pflag, eflag, qflag;
  booleantype doProjection;

  saved_t = tn;
  ncf = npf = nef = 0;
  nflag = FIRST_CALL;

#ifdef CPODES_DEBUG
  printf("\nStep # %ld   t = %lg\n",nst,tn);
#endif

  /* If the step size was changed, adjust method parameters */
  if ((nst > 0) && (hprime != h)) {

#ifdef CPODES_DEBUG
    printf("   Adjust parameters\n");
    printf("      qprime = %d   hprime = %g\n", qprime, hprime);
#endif

    cpAdjustParams(cp_mem);

#ifdef CPODES_DEBUG
    printf("   Done  adjust parameters\n");
#endif

  }

  /* Looping point for attempts to take a step */

  loop {  

    /* Prediction */

    cpPredict(cp_mem);  
    cpSet(cp_mem);

#ifdef CPODES_DEBUG
    printf("   Predict\n");
#endif

    /* Solve the nonlinear system to correct the states */

#ifdef CPODES_DEBUG
    printf("   Nonlinear solver\n");
#endif

    kflag = cpNls(cp_mem, nflag, saved_t, &ncf);

#ifdef CPODES_DEBUG
    printf("   Nonlinear solver return flag = %d\n",kflag);
#endif

    if (kflag == PREDICT_AGAIN) {
      /* If we need to predict again, set nflag=PREV_CONV_FAIL */
      nflag = PREV_CONV_FAIL;
      continue;
    } else if (kflag != CP_SUCCESS) {
      /* Recovery is not possible. */
      return(kflag);
    }

    /* Nonlinear solve successful */
    ncf = 0;
    
    /* Do we need to perform projection? */

#ifdef CPODES_DEBUG
    printf("   Perform projection? ");
#endif

    doProjection = proj_enabled && (proj_freq > 0) && 
      ( (nst==0) || (nst >= nstlprj + proj_freq) );
    applyProj = FALSE;

#ifdef CPODES_DEBUG
    printf("%d\n",doProjection);
#endif

    /* If needed, perform the projection step */

    if ( doProjection ) {

#ifdef CPODE_DEBUG
      printf("   Projection\n");
#endif

      pflag = cpDoProjection(cp_mem, saved_t, &npf);

#ifdef CPODES_DEBUG
      printf("   Projection return flag = %d\n",pflag);
#endif

      if (pflag == PREDICT_AGAIN) {
        /* If we need to predict again, set nflag=PREV_PROJ_FAIL */
        nflag = PREV_PROJ_FAIL;
        continue;
      } else if (pflag != CP_SUCCESS) {
        /* Recovery is not possible. */
        return(pflag);
      }

      /* Projection successful */
      npf = 0;

    }

    /* Perform error test */

#ifdef CPODES_DEBUG
    printf("   Error test\n");
#endif

    eflag = cpDoErrorTest(cp_mem, saved_t, acnrm, &nef, &dsm);

#ifdef CPODES_DEBUG
    printf("   Error test return flag = %d\n",eflag);
#endif

    if (eflag == PREDICT_AGAIN) {
      /* If we need to predict again set nflag=PREV_ERR_FAIL */
      nflag = PREV_ERR_FAIL;
      continue;
    } else if (eflag != CP_SUCCESS) {
      /* Recovery is not possible */
      return(eflag);
    }

    /* State error test successful */
    nef = 0;

    /* Deal with quadrature variables if needed */

    if (quadr) {

      /* Correct quadrature variables */

#ifdef CPODES_DEBUG
      printf("   Correct quadratures\n");
#endif

      qflag = cpQuadNls(cp_mem, saved_t, &ncf);

#ifdef CPODES_DEBUG
      printf("   Correct quadratures return flag = %d\n",qflag);
#endif

      if (qflag == PREDICT_AGAIN) {
        /* If we need to predict again, set nflag=PREV_CONV_FAIL */
        nflag = PREV_CONV_FAIL;
        continue;
      } else if (qflag != CP_SUCCESS) {
        /* Recovery is not possible. */
        return(qflag);
      }

      /* Quadrature correction successful */
      ncf = 0;

      /* If needed, perform error test on quadrature variables */

      if (errconQ) {

#ifdef CPODES_DEBUG
        printf("   Quadratures error test\n");
#endif

        eflag = cpDoErrorTest(cp_mem, saved_t, acnrmQ, &nef, &dsmQ);

#ifdef CPODES_DEBUG
        printf("   Quadratures error test return flag = %d\n",eflag);
#endif

        if (eflag == PREDICT_AGAIN) {
          /* If we need to predict again set nflag=PREV_ERR_FAIL */
          nflag = PREV_ERR_FAIL;
          continue;
        } else if (eflag != CP_SUCCESS) {
          /* Recovery is not possible */
          return(eflag);
        }
        
        /* Quadrature error test successful */
        nef = 0;

        /* Update 'dsm' with 'dsmQ' (to be used in cpPrepareNextStep) */
        dsm = cpQuadUpdateDsm(cp_mem, dsm, dsmQ);

      }

    }

    /* Everything went fine. Make final updates and exit loop */
    nstlprj = nst;
    break;
    
  }

  /* Apply corrections to Nordsieck history arrays, 
   * update data, and 
   * consider change of step and/or order. */
  cpCorrect(cp_mem);
  cpCompleteStep(cp_mem); 
  cpPrepareNextStep(cp_mem, dsm); 

  /* If Stablilty Limit Detection is turned on 
   * check if order must be reduced */
  if (sldeton) cpBDFStab(cp_mem);

  etamax = (nst <= SMALL_NST) ? ETAMX2 : ETAMX3;

  /* Finally, we rescale the acor array to be 
   * the estimated local error vector. */
  N_VScale(ONE/tq[2], acor, acor);
  if (quadr) N_VScale(ONE/tq[2], acorQ, acorQ);

  return(CP_SUCCESS);
      
}

/*
 * cpGetSolution
 *
 * This routine evaluates y(t) and y'(t) as the value and derivative of 
 * the interpolating polynomial at the independent variable t, and stores
 * the results in the vectors yret and ypret.
 *
 * The formulas are:
 *                  q 
 *  yret = zn[0] + SUM (t - tn)^j * h^(-j) * zn[j] , 
 *                 j=1 
 *  and
 *           q 
 *  ypret = SUM j * (t - tn)^(j-1) * h^(-j) * zn[j] , 
 *          j=1
 *  
 * This function is called by CPode() with t = tout, t = tn, or t = tstop.
 *
 */

int cpGetSolution(void *cpode_mem, realtype t, N_Vector yret, N_Vector ypret)
{
  realtype s, c, d;
  realtype tfuzz, tp, tn1;
  int j;
  CPodeMem cp_mem;

  cp_mem = (CPodeMem) cpode_mem;

  /* Allow for some slack */
  tfuzz = FUZZ_FACTOR * uround * (ABS(tn) + ABS(hu));
  if (hu < ZERO) tfuzz = -tfuzz;
  tp = tn - hu - tfuzz;
  tn1 = tn + tfuzz;
  if ((t-tp)*(t-tn1) > ZERO) return(CP_BAD_T);

  /* Initialize yret = zn[0] and ypret = 0 */
  N_VScale(ONE, zn[0], yret);
  N_VConst(ZERO, ypret);

  /* Accumulate multiples of columns zn[j] into yret and ypret. */
  s = (t - tn) / h;
  c = s;
  d = ONE/h;
  for (j=1; j<=q; j++) {
    N_VLinearSum(ONE,  yret, c, zn[j],  yret);
    N_VLinearSum(ONE, ypret, d, zn[j], ypret);
    d = (j+1)*c/h;
    c *= s;
  }

  return(CP_SUCCESS);
}

/* 
 * -----------------------------------------------------------------
 * Functions acting on the Nordsieck history array
 * -----------------------------------------------------------------
 */

/*
 * cpPredict
 *
 * This routine advances tn by the tentative step size h, and computes
 * the predicted array z_n(0), which is overwritten on zn.
 * The prediction of zn is done by repeated additions.
 * In TSTOP mode, it is possible for tn + h to be past tstop by roundoff,
 * and in that case, we reset tn (after incrementing by h) to tstop.
 */

static void cpPredict(CPodeMem cp_mem)
{
  int j, k;
  
  tn += h;
  if (istop) {
    if ((tn - tstop)*h > ZERO) tn = tstop;
  }

  for (k = 1; k <= q; k++)
    for (j = q; j >= k; j--) 
      N_VLinearSum(ONE, zn[j-1], ONE, zn[j], zn[j-1]); 

  if (quadr) {
    for (k = 1; k <= q; k++)
      for (j = q; j >= k; j--) 
        N_VLinearSum(ONE, znQ[j-1], ONE, znQ[j], znQ[j-1]);
  }
  
}

/*
 * cpRescale
 *
 * This routine rescales the Nordsieck array by multiplying the
 * jth column zn[j] by eta^j, j = 1, ..., q.  Then the value of
 * h is rescaled by eta, and hscale is reset to h.
 *
 * This function is used in cpNls, cpDoProjection, and cpDoerrorTest
 */

void cpRescale(CPodeMem cp_mem)
{
  int j;
  realtype factor;
  
  factor = eta;
  for (j=1; j <= q; j++) {
    N_VScale(factor, zn[j], zn[j]);
    if (quadr) N_VScale(factor, znQ[j], znQ[j]);
    factor *= eta;
  }
  h = hscale * eta;
  next_h = h;
  hscale = h;
  nscon = 0;
}

/*
 * cpRestore
 *
 * This routine restores the value of tn to saved_t and undoes the
 * prediction.  After execution of cpRestore, the Nordsieck array zn has
 * the same values as before the call to cpPredict.
 *
 * This function is used in cpNls, cpDoProjection, and cpDoerrorTest
 */

void cpRestore(CPodeMem cp_mem, realtype saved_t)
{
  int j, k;
  
  tn = saved_t;

  for (k = 1; k <= q; k++)
    for (j = q; j >= k; j--)
      N_VLinearSum(ONE, zn[j-1], -ONE, zn[j], zn[j-1]);

  if (quadr) {
    for (k = 1; k <= q; k++)
      for (j = q; j >= k; j--)
        N_VLinearSum(ONE, znQ[j-1], -ONE, znQ[j], znQ[j-1]);
  }

}

/*
 * cpCorrect
 *
 * This routine applies the corrections to the zn array.
 * The correction to zn is done by repeated additions.
 */

static void cpCorrect(CPodeMem cp_mem)
{
  int j;
 
  for (j=0; j <= q; j++)
    N_VLinearSum(l[j], acor, ONE, zn[j], zn[j]);

  if (applyProj) {
    for (j=0; j <= q; j++)
      N_VLinearSum(p[j], acorP, ONE, zn[j], zn[j]); 
  }
  
  if (quadr) {
    for (j=0; j <= q; j++) 
      N_VLinearSum(l[j], acorQ, ONE, znQ[j], znQ[j]);
  }

}

/* 
 * -----------------------------------------------------------------
 * Quadrature correction function
 * -----------------------------------------------------------------
 */

/*
 * cpQuadNls
 *
 * This function attempts to correct the quadrature variables.
 *
 * The only way this function can fail is through the user's quad
 * function. If quad fails unrecoverably we return immediately. If
 * it fails with a recoverable error, we consider a new prediction.
 *
 * We count any failure here as a convergence failure (we do not
 * distinguish these from convergence failures in the nonlinear 
 * solver for the state correction).
 *
 * NOTE: This is technically not a nonlinear solver, but this name is 
 * used to reflect the similarity to cpNls.
 */

static int cpQuadNls(CPodeMem cp_mem, realtype saved_t, int *ncfPtr)
{
  int retval;

  /* Save quadrature correction in acorQ */
  retval = qfun(tn, y, acorQ, q_data);
  nqe++;

  if (retval == 0) {

    /* Find correction */
    N_VLinearSum(h, acorQ, -ONE, znQ[1], acorQ);
    N_VScale(rl1, acorQ, acorQ);

    /* Compute its WRMS norm */
    acnrmQ = N_VWrmsNorm(acorQ, ewtQ);

    /* Apply correction to quadrature variables */
    N_VLinearSum(ONE, znQ[0], ONE, acorQ, yQ);

    return(CP_SUCCESS);

  }

  /* Increment ncfn and restore zn */
  ncfn++;
  cpRestore(cp_mem, saved_t);

  if (retval < 0) return(CP_QUADFUNC_FAIL);
 
  /* At this point, qfun had a recoverable error; increment ncf */
  (*ncfPtr)++;
  etamax = ONE;

  /* If we had maxncf failures or |h| = hmin, 
     return CP_REPTD_QUADFUNC_ERR. */
  if ((ABS(h) <= hmin*ONEPSM) || (*ncfPtr == maxncf))
    return(CP_REPTD_QUADFUNC_ERR);    

  /* Reduce step size; return to reattempt the step */
  eta = MAX(ETACF, hmin / ABS(h));
  cpRescale(cp_mem);

  return(PREDICT_AGAIN);
}

/* 
 * -----------------------------------------------------------------
 * LMM-related functions
 * -----------------------------------------------------------------
 */

/*
 * cpAdjustParams
 *
 * This routine is called when a change in step size was decided upon,
 * and it handles the required adjustments to the history array zn.
 * If there is to be a change in order, we call cpAdjustOrder and reset
 * q, L = q+1, and qwait.  Then in any case, we call cpRescale, which
 * resets h and rescales the Nordsieck array.
 */

static void cpAdjustParams(CPodeMem cp_mem)
{
  if (qprime != q) {
    cpAdjustOrder(cp_mem, qprime-q);
    q = qprime;
    L = q+1;
    qwait = L;
  }
  cpRescale(cp_mem);
}

/*
 * cpAdjustOrder
 *
 * This routine is a high level routine which handles an order
 * change by an amount deltaq (= +1 or -1). If a decrease in order
 * is requested and q==2, then the routine returns immediately.
 * Otherwise cpAdjustAdams or cpAdjustBDF is called to handle the
 * order change (depending on the value of lmm_type).
 */

static void cpAdjustOrder(CPodeMem cp_mem, int deltaq)
{
  if ((q==2) && (deltaq != 1)) return;
  
  switch(lmm_type){
  case CP_ADAMS: 
    cpAdjustAdams(cp_mem, deltaq);
    break;
  case CP_BDF:   
    cpAdjustBDF(cp_mem, deltaq);
    break;
  }
}

/*
 * cpAdjustAdams
 *
 * This routine adjusts the history array on a change of order q by
 * deltaq, in the case that lmm_type == CP_ADAMS.
 */

static void cpAdjustAdams(CPodeMem cp_mem, int deltaq)
{
  int i, j;
  realtype xi, hsum;

  /* On an order increase, set new column of zn to zero and return */
  
  if (deltaq==1) {
    N_VConst(ZERO, zn[L]);
    if (quadr) N_VConst(ZERO, znQ[L]);
    return;
  }

  /*
   * On an order decrease, each zn[j] is adjusted by a multiple of zn[q].
   * The coeffs. in the adjustment are the coeffs. of the polynomial:
   *        x
   * q * INT { u * ( u + xi_1 ) * ... * ( u + xi_{q-2} ) } du 
   *        0
   * where xi_j = [t_n - t_(n-j)]/h => xi_0 = 0
   */

  for (i=0; i <= qmax; i++) l[i] = ZERO;
  l[1] = ONE;
  hsum = ZERO;
  for (j=1; j <= q-2; j++) {
    hsum += tau[j];
    xi = hsum / hscale;
    for (i=j+1; i >= 1; i--) l[i] = l[i]*xi + l[i-1];
  }
  
  for (j=1; j <= q-2; j++) l[j+1] = q * (l[j] / (j+1));
  
  for (j=2; j < q; j++)
    N_VLinearSum(-l[j], zn[q], ONE, zn[j], zn[j]);

  if (quadr) {
    for (j=2; j < q; j++)
      N_VLinearSum(-l[j], znQ[q], ONE, znQ[j], znQ[j]);
  }

}

/*
 * cpAdjustBDF
 *
 * This is a high level routine which handles adjustments to the
 * history array on a change of order by deltaq in the case that 
 * lmm_type == CP_BDF.  cpAdjustBDF calls cpIncreaseBDF if deltaq = +1 and 
 * cpDecreaseBDF if deltaq = -1 to do the actual work.
 */

static void cpAdjustBDF(CPodeMem cp_mem, int deltaq)
{
  switch(deltaq) {
  case 1 : 
    cpIncreaseBDF(cp_mem);
    return;
  case -1: 
    cpDecreaseBDF(cp_mem);
    return;
  }
}

/*
 * cpIncreaseBDF
 *
 * This routine adjusts the history array on an increase in the 
 * order q in the case that lmm_type == CP_BDF.  
 * A new column zn[q+1] is set equal to a multiple of the saved 
 * vector (= acor) in zn[indx_acor].  Then each zn[j] is adjusted by
 * a multiple of zn[q+1].  The coefficients in the adjustment are the 
 * coefficients of the polynomial x*x*(x+xi_1)*...*(x+xi_j),
 * where xi_j = [t_n - t_(n-j)]/h.
 */

static void cpIncreaseBDF(CPodeMem cp_mem)
{
  realtype alpha0, alpha1, prod, xi, xiold, hsum, A1;
  int i, j;
  
  for (i=0; i <= qmax; i++) l[i] = ZERO;
  l[2] = alpha1 = prod = xiold = ONE;
  alpha0 = -ONE;
  hsum = hscale;
  if (q > 1) {
    for (j=1; j < q; j++) {
      hsum += tau[j+1];
      xi = hsum / hscale;
      prod *= xi;
      alpha0 -= ONE / (j+1);
      alpha1 += ONE / xi;
      for (i=j+2; i >= 2; i--) l[i] = l[i]*xiold + l[i-1];
      xiold = xi;
    }
  }
  A1 = (-alpha0 - alpha1) / prod;

  /* 
   * zn[indx_acor] contains the value Delta_n = y_n - y_n(0) 
   * This value was stored there at the previous successful
   * step (in cpCompleteStep) 
   * 
   * A1 contains dbar = (1/xi* - 1/xi_q)/prod(xi_j)
   */

  N_VScale(A1, zn[indx_acor], zn[L]);
  for (j=2; j <= q; j++)
    N_VLinearSum(l[j], zn[L], ONE, zn[j], zn[j]);
  
  if (quadr) {
    N_VScale(A1, znQ[indx_acor], znQ[L]);
    for (j=2; j <= q; j++)
      N_VLinearSum(l[j], znQ[L], ONE, znQ[j], znQ[j]);
  }

}

/*
 * cpDecreaseBDF
 *
 * This routine adjusts the history array on a decrease in the 
 * order q in the case that lmm_type == CP_BDF.  
 * Each zn[j] is adjusted by a multiple of zn[q].  The coefficients
 * in the adjustment are the coefficients of the polynomial
 *   x*x*(x+xi_1)*...*(x+xi_j), where xi_j = [t_n - t_(n-j)]/h.
 */

static void cpDecreaseBDF(CPodeMem cp_mem)
{
  realtype hsum, xi;
  int i, j;
  
  for (i=0; i <= qmax; i++) l[i] = ZERO;
  l[2] = ONE;
  hsum = ZERO;
  for(j=1; j <= q-2; j++) {
    hsum += tau[j];
    xi = hsum /hscale;
    for (i=j+2; i >= 2; i--) l[i] = l[i]*xi + l[i-1];
  }
  
  for(j=2; j < q; j++)
    N_VLinearSum(-l[j], zn[q], ONE, zn[j], zn[j]);

  if (quadr) {
    for (j=2; j < q; j++)
      N_VLinearSum(-l[j], znQ[q], ONE, znQ[j], znQ[j]);
  }

}

/*
 * cpSet
 *
 * This routine is a high level routine which calls cpSetAdams or
 * cpSetBDF to set the coefficients l of the polynomial Lambda 
 * (used in applying the NLS correction to zn), the coefficients p 
 * of the polynomial Phi (used in applying the projection correction 
 * to zn), the test quantity array tq (used in the NLS convergence 
 * test and in the error test), and the related variables rl1, gamma, 
 * and gamrat.
 */

static void cpSet(CPodeMem cp_mem)
{
  switch(lmm_type) {
  case CP_ADAMS: 
    cpSetAdams(cp_mem);
    break;
  case CP_BDF:
    cpSetBDF(cp_mem);
    break;
  }
  rl1 = ONE / l[1];
  gamma = h * rl1;
  if (nst == 0) gammap = gamma;
  gamrat = (nst > 0) ? gamma / gammap : ONE;  /* protect x / x != 1.0 */
}

/*
 * cpSetAdams
 *
 * This routine handles the computation of l, p, and tq for the
 * case lmm_type == CP_ADAMS.
 *
 * The components of the array l are the coefficients of a
 * polynomial Lambda(x) = l_0 + l_1 x + ... + l_q x^q, that
 * satisfies the following q+1 conditions:
 *     Lambda(-1) = 0, Lambda(0) = 1
 *     Lambda'(-xi_i) = 0, i=1,2,...,q-1.
 * Here x = [t - t_n]  /h and xi_i = [t_n - t_(n-i)] / h.
 *
 * For the Adams case, the polynomial Phi(x) is identical to Lambda(x)
 * and therefore p_i = l_i, i=0,1,2,...,q.
 *
 * The array tq is set to test quantities used in the convergence
 * test, the error test, and the selection of h at a new order.
 */

static void cpSetAdams(CPodeMem cp_mem)
{
  realtype m[L_MAX], M[3], hsum;
  int i;

  if (q == 1) {

    l[0] = l[1] = tq[1] = tq[5] = ONE;
    tq[2] = TWO;
    tq[3] = TWELVE;
    tq[4] = nlscoef * tq[2];       /* = 0.1 * tq[2] */

  } else {
  
    hsum = cpAdamsStart(cp_mem, m);
    
    M[0] = cpAltSum(q-1, m, 1);
    M[1] = cpAltSum(q-1, m, 2);
    
    cpAdamsFinish(cp_mem, m, M, hsum);

  }

  for(i=0; i<=q; i++) p[i] = l[i];

  return;
}

/*
 * cpAdamsStart
 *
 * This routine generates in m[] the coefficients of the product
 * polynomial needed for the Adams l and tq coefficients for q > 1.
 */

static realtype cpAdamsStart(CPodeMem cp_mem, realtype m[])
{
  realtype hsum, xi_inv, sum;
  int i, j;
  
  hsum = h;
  m[0] = ONE;
  for (i=1; i <= q; i++) m[i] = ZERO;
  for (j=1; j < q; j++) {
    if ((j==q-1) && (qwait == 1)) {
      sum = cpAltSum(q-2, m, 2);
      tq[1] = m[q-2] / (q * sum);
    }
    xi_inv = h / hsum;
    for (i=j; i >= 1; i--) m[i] += m[i-1] * xi_inv;
    hsum += tau[j];
    /* The m[i] are coefficients of product(1 to j) (1 + x/xi_i) */
  }
  return(hsum);
}

/*
 * cpAdamsFinish
 *
 * This routine completes the calculation of the Adams l and tq.
 */

static void cpAdamsFinish(CPodeMem cp_mem, realtype m[], realtype M[], realtype hsum)
{
  int i;
  realtype M0_inv, xi, xi_inv;
  
  M0_inv = ONE / M[0];
  
  l[0] = ONE;
  for (i=1; i <= q; i++) l[i] = M0_inv * (m[i-1] / i);
  xi = hsum / h;
  xi_inv = ONE / xi;
  
  tq[2] = xi * M[0] / M[1];
  tq[5] = xi / l[q];

  if (qwait == 1) {
    for (i=q; i >= 1; i--) m[i] += m[i-1] * xi_inv;
    M[2] = cpAltSum(q, m, 2);
    tq[3] = L * M[0] / M[2];
  }

  tq[4] = nlscoef * tq[2];
}

/*  
 * cpAltSum
 *
 * cpAltSum returns the value of the alternating sum
 *   sum (i= 0 ... iend) [ (-1)^i * (a[i] / (i + k)) ].
 * If iend < 0 then cpAltSum returns 0.
 * This operation is needed to compute the integral, from -1 to 0,
 * of a polynomial x^(k-1) M(x) given the coefficients of M(x).
 */

static realtype cpAltSum(int iend, realtype a[], int k)
{
  int i, sign;
  realtype sum;
  
  if (iend < 0) return(ZERO);
  
  sum = ZERO;
  sign = 1;
  for (i=0; i <= iend; i++) {
    sum += sign * (a[i] / (i+k));
    sign = -sign;
  }
  return(sum);
}

/*
 * cpSetBDF
 *
 * This routine computes the coefficients l and tq in the case
 * lmm_type == CP_BDF.  cpSetBDF calls cpSetTqBDF to set the test
 * quantity array tq. 
 * 
 * The components of the array l are the coefficients of a
 * polynomial Lambda(x) = l_0 + l_1 x + ... + l_q x^q, given by
 *                                 q-1
 * Lambda(x) = (1 + x / xi*_q) * PRODUCT (1 + x / xi_i)
 *                                 i=1 
 * 
 * The components of the array p are the coefficients of a
 * polynomial Phi(x) = p_0 + p_1 x + ... + p_q x^q, given by
 *            q
 * Phi(x) = PRODUCT (1 + x / xi_i)
 *           i=1 
 *
 * Here x = [t - t_n]  /h and xi_i = [t_n - t_(n-i)] / h.
 *
 * The array tq is set to test quantities used in the convergence
 * test, the error test, and the selection of h at a new order.
 */

static void cpSetBDF(CPodeMem cp_mem)
{
  realtype alpha0, alpha0_hat, xi_inv, xistar_inv, hsum;
  int i,j;
  
  l[0] = l[1] = ONE;
  p[0] = p[1] = ONE;
  xi_inv = xistar_inv = ONE;
  for (i=2; i <= q; i++) {
    l[i] = ZERO;
    p[i] = ZERO;
  }
  alpha0 = alpha0_hat = -ONE;
  hsum = h;
  if (q > 1) {
    for (j=2; j < q; j++) {
      hsum += tau[j-1];
      xi_inv = h / hsum;
      alpha0 -= ONE / j;
      /* The l[i] and p[i] are coefficients of product(1 to j) (1 + x/xi_i) */
      for(i=j; i >= 1; i--) {
        l[i] += l[i-1]*xi_inv;
        p[i] += p[i-1]*xi_inv;
      }
    }
    
    /* j = q */
    alpha0 -= ONE / q;
    xistar_inv = -l[1] - alpha0;
    hsum += tau[q-1];
    xi_inv = h / hsum;
    alpha0_hat = -l[1] - xi_inv;
    for (i=q; i >= 1; i--) {
      l[i] += l[i-1]*xistar_inv;
      p[i] += p[i-1]*xi_inv;
    }
  }

  cpSetTqBDF(cp_mem, hsum, alpha0, alpha0_hat, xi_inv, xistar_inv);
}

/*
 * cpSetTqBDF
 *
 * This routine sets the test quantity array tq in the case
 * lmm_type == CP_BDF.
 */

static void cpSetTqBDF(CPodeMem cp_mem, realtype hsum, realtype alpha0,
                       realtype alpha0_hat, realtype xi_inv, realtype xistar_inv)
{
  realtype A1, A2, A3, A4, A5, A6;
  realtype C, CPrime, CPrimePrime;
  
  A1 = ONE - alpha0_hat + alpha0;
  A2 = ONE + q * A1;
  tq[2] = ABS(alpha0 * (A2 / A1));
  tq[5] = ABS((A2) / (l[q] * xi_inv/xistar_inv));
  if (qwait == 1) {
    C = xistar_inv / l[q];
    A3 = alpha0 + ONE / q;
    A4 = alpha0_hat + xi_inv;
    CPrime = A3 / (ONE - A4 + A3);
    tq[1] = ABS(CPrime / C);
    hsum += tau[q];
    xi_inv = h / hsum;
    A5 = alpha0 - (ONE / (q+1));
    A6 = alpha0_hat - xi_inv;
    CPrimePrime = A2 / (ONE - A6 + A5);
    tq[3] = ABS(CPrimePrime * xi_inv * (q+2) * A5);
  }
  tq[4] = nlscoef * tq[2];
}

/* 
 * -----------------------------------------------------------------
 * Error test
 * -----------------------------------------------------------------
 */

/*
 * cpDoErrorTest
 *
 * This routine performs the local error test. It is called to test
 * both the sates and the quadrature variables, with different values
 * of acor_norm and dsmPtr.
 *
 * The weighted local error norm dsm is loaded into *dsmPtr, and 
 * we test whether dsm <= 1.
 *
 * If the test passes, cpDoErrorTest returns CP_SUCCESS. 
 *
 * If the test fails, we undo the step just taken (call cpRestore) and 
 *
 *   - if maxnef error test failures have occurred or if ABS(h) = hmin,
 *     we return CP_ERR_FAILURE.
 *
 *   - if more than MXNEF1 error test failures have occurred, an order
 *     reduction is forced. If already at order 1, restart by reloading 
 *     zn from scratch.
 *
 *   - otherwise, set return PREDICT_AGAIN, telling cpStep to retry the
 *     step (with nflag = PREV_ERR_FAIL).
 *
 * Note that we lump all error test failures together (i.e. we do not
 * distinguish between failures due to states or to quadratures).
 */

static int cpDoErrorTest(CPodeMem cp_mem, realtype saved_t, realtype acor_norm,
                         int *nefPtr, realtype *dsmPtr)
{
  realtype dsm;

  /* If est. local error norm dsm passes test, return CP_SUCCESS */  
  dsm = acor_norm / tq[2];
  *dsmPtr = dsm; 

#ifdef CPODES_DEBUG
  printf("     dsm = %lg\n",dsm);
#endif

  if (dsm <= ONE) return(CP_SUCCESS);
  
  /* Test failed: increment counters and restore zn array */
  (*nefPtr)++;
  netf++;
  cpRestore(cp_mem, saved_t);

  /* At maxnef failures or |h| = hmin, return CP_ERR_FAILURE */
  if ((ABS(h) <= hmin*ONEPSM) || (*nefPtr == maxnef)) return(CP_ERR_FAILURE);

  /* Set etamax = 1 to prevent step size increase at end of this step */
  etamax = ONE;

  /* Set h ratio eta from dsm, rescale, and return for retry of step */
  if (*nefPtr <= MXNEF1) {
    eta = ONE / (RPowerR(BIAS2*dsm,ONE/L) + ADDON);
    eta = MAX(ETAMIN, MAX(eta, hmin / ABS(h)));
    if (*nefPtr >= SMALL_NEF) eta = MIN(eta, ETAMXF);
    cpRescale(cp_mem);
    return(PREDICT_AGAIN);
  }
  
  /* After MXNEF1 failures, force an order reduction and retry step */
  if (q > 1) {
    eta = MAX(ETAMIN, hmin / ABS(h));
    cpAdjustOrder(cp_mem,-1);
    L = q;
    q--;
    qwait = L;
    cpRescale(cp_mem);
    return(PREDICT_AGAIN);
  }

  /* If already at order 1, restart: reload zn from scratch */
  eta = MAX(ETAMIN, hmin / ABS(h));
  h *= eta;
  next_h = h;
  hscale = h;
  qwait = LONG_WAIT;
  nscon = 0;
  N_VScale(eta, zn[1], zn[1]);
  if (quadr) N_VScale(eta, znQ[1], znQ[1]);

  return(PREDICT_AGAIN);
}

/* 
 * -----------------------------------------------------------------
 * Succesful step completion functions
 * -----------------------------------------------------------------
 */

/*
 * cpCompleteStep
 *
 * This routine performs various update operations when the solution
 * has passed the local error test. 
 * We increment the step counter nst, record the values hu and qu,
 * and update the tau array.
 * The tau[i] are the last q values of h, with tau[1] the most recent.
 * The counter qwait is decremented, and if qwait == 1 (and q < qmax)
 * we save acor and tq[5] for a possible order increase.
 */

static void cpCompleteStep(CPodeMem cp_mem)
{
  int i;
  
  nst++;
  nscon++;
  hu = h;
  qu = q;

  for (i=q; i >= 2; i--)  tau[i] = tau[i-1];
  if ((q==1) && (nst > 1)) tau[2] = tau[1];
  tau[1] = h;

  /* If necessary, store Delta_n in zn[qmax] to be used in order increase
   *
   * This actually will be Delta_{n-1} in the ELTE at q+1 since it happens at
   * the next to last step of order q before a possible one at order q+1
   */

  qwait--;
  if ((qwait == 1) && (q != qmax)) {
    N_VScale(ONE, acor, zn[qmax]);
    if (quadr) N_VScale(ONE, acorQ, znQ[qmax]);
    saved_tq5 = tq[5];
    indx_acor = qmax;
  }
}

/*
 * CPprepareNextStep
 *
 * This routine handles the setting of stepsize and order for the
 * next step -- hprime and qprime.  Along with hprime, it sets the
 * ratio eta = hprime/h.  It also updates other state variables 
 * related to a change of step size or order. 
 */

static void cpPrepareNextStep(CPodeMem cp_mem, realtype dsm)
{
  /* If etamax = 1, defer step size or order changes */
  if (etamax == ONE) {
    qwait = MAX(qwait, 2);
    qprime = q;
    hprime = h;
    eta = ONE;
    return;
  }

  /* etaq is the ratio of new to old h at the current order */  
  etaq = ONE /(RPowerR(BIAS2*dsm,ONE/L) + ADDON);
  
  /* If no order change, adjust eta and acor in cpSetEta and return */
  if (qwait != 0) {
    eta = etaq;
    qprime = q;
    cpSetEta(cp_mem);
    return;
  }
  
  /* If qwait = 0, consider an order change.   etaqm1 and etaqp1 are 
     the ratios of new to old h at orders q-1 and q+1, respectively.
     cpChooseEta selects the largest; cpSetEta adjusts eta and acor */
  qwait = 2;
  etaqm1 = cpComputeEtaqm1(cp_mem);
  etaqp1 = cpComputeEtaqp1(cp_mem);  
  cpChooseEta(cp_mem); 
  cpSetEta(cp_mem);
}

/*
 * CPsetEta
 *
 * This routine adjusts the value of eta according to the various
 * heuristic limits and the optional input hmax.
 */

static void cpSetEta(CPodeMem cp_mem)
{

  /* If eta below the threshhold THRESH, reject a change of step size */
  if (eta < THRESH) {
    eta = ONE;
    hprime = h;
  } else {
    /* Limit eta by etamax and hmax, then set hprime */
    eta = MIN(eta, etamax);
    eta /= MAX(ONE, ABS(h)*hmax_inv*eta);
    hprime = h * eta;
    if (qprime < q) nscon = 0;
  }

}

/*
 * cpComputeEtaqm1
 *
 * This routine computes and returns the value of etaqm1 for a
 * possible decrease in order by 1.
 */

static realtype cpComputeEtaqm1(CPodeMem cp_mem)
{
  realtype ddn;
  
  etaqm1 = ZERO;
  if (q > 1) {
    ddn = N_VWrmsNorm(zn[q], ewt) / tq[1];
    if ( quadr && errconQ) ddn = cpQuadUpdateNorm(cp_mem, ddn, znQ[q], ewtQ);
    etaqm1 = ONE/(RPowerR(BIAS1*ddn, ONE/q) + ADDON);
  }
  return(etaqm1);
}

/*
 * cpComputeEtaqp1
 *
 * This routine computes and returns the value of etaqp1 for a
 * possible increase in order by 1.
 */

static realtype cpComputeEtaqp1(CPodeMem cp_mem)
{
  realtype dup, cquot;
  
  etaqp1 = ZERO;
  if (q != qmax) {
    cquot = (tq[5] / saved_tq5) * RPowerI(h/tau[2], L);
    N_VLinearSum(-cquot, zn[qmax], ONE, acor, tempv);
    dup = N_VWrmsNorm(tempv, ewt);
    if ( quadr && errconQ ) {
      N_VLinearSum(-cquot, znQ[qmax], ONE, acorQ, tempvQ);
      dup = cpQuadUpdateNorm(cp_mem, dup, tempvQ, ewtQ);
    }
    dup = dup / tq[3];
    etaqp1 = ONE / (RPowerR(BIAS3*dup, ONE/(L+1)) + ADDON);
  }
  return(etaqp1);
}

/*
 * cpChooseEta
 * Given etaqm1, etaq, etaqp1 (the values of eta for qprime =
 * q - 1, q, or q + 1, respectively), this routine chooses the 
 * maximum eta value, sets eta to that value, and sets qprime to the
 * corresponding value of q.  If there is a tie, the preference
 * order is to (1) keep the same order, then (2) decrease the order,
 * and finally (3) increase the order.  If the maximum eta value
 * is below the threshhold THRESH, the order is kept unchanged and
 * eta is set to 1.
 */

static void cpChooseEta(CPodeMem cp_mem)
{
  realtype etam;
  
  etam = MAX(etaqm1, MAX(etaq, etaqp1));
  
  if (etam < THRESH) {
    eta = ONE;
    qprime = q;
    return;
  }

  if (etam == etaq) {

    eta = etaq;
    qprime = q;

  } else if (etam == etaqm1) {

    eta = etaqm1;
    qprime = q - 1;

  } else {

    eta = etaqp1;
    qprime = q + 1;

    if (lmm_type == CP_BDF) {

      /* 
       * Store Delta_n in zn[qmax] to be used in order increase 
       *
       * This happens at the last step of order q before an increase
       * to order q+1, so it represents Delta_n in the ELTE at q+1
       */

      N_VScale(ONE, acor, zn[qmax]);
      if (quadr && errconQ) N_VScale(ONE, acorQ, znQ[qmax]);

    }

  }

}

/*
 * cpHandleFailure
 *
 * This routine prints error messages for all cases of failure by
 * cpHin and cpStep. It returns to CPode the value that CPode is 
 * to return to the user.
 */

static int cpHandleFailure(CPodeMem cp_mem, int flag)
{

  /* Depending on flag, print error message and return error flag */
  switch (flag) {
  case CP_ERR_FAILURE: 
    cpProcessError(cp_mem, CP_ERR_FAILURE, "CPODES", "CPode", MSGCP_ERR_FAILS, tn, h);
    break;
  case CP_CONV_FAILURE:
    cpProcessError(cp_mem, CP_CONV_FAILURE, "CPODES", "CPode", MSGCP_CONV_FAILS, tn, h);
    break;
  case CP_LSETUP_FAIL:
    cpProcessError(cp_mem, CP_LSETUP_FAIL, "CPODES", "CPode", MSGCP_SETUP_FAILED, tn);
    break;
  case CP_LSOLVE_FAIL:
    cpProcessError(cp_mem, CP_LSOLVE_FAIL, "CPODES", "CPode", MSGCP_SOLVE_FAILED, tn);
    break;
  case CP_ODEFUNC_FAIL:
    cpProcessError(cp_mem, CP_ODEFUNC_FAIL, "CPODES", "CPode", MSGCP_ODEFUNC_FAILED, tn);
    break;
  case CP_UNREC_ODEFUNC_ERR:
    cpProcessError(cp_mem, CP_UNREC_ODEFUNC_ERR, "CPODES", "CPode", MSGCP_ODEFUNC_UNREC, tn);
    break;
  case CP_REPTD_ODEFUNC_ERR:
    cpProcessError(cp_mem, CP_REPTD_ODEFUNC_ERR, "CPODES", "CPode", MSGCP_ODEFUNC_REPTD, tn);
    break;
  case CP_RTFUNC_FAIL:    
    cpProcessError(cp_mem, CP_RTFUNC_FAIL, "CPODES", "CPode", MSGCP_RTFUNC_FAILED, tn);
    break;
  case CP_TOO_CLOSE:
    cpProcessError(cp_mem, CP_TOO_CLOSE, "CPODES", "CPode", MSGCP_TOO_CLOSE);
    break;
  case CP_PROJ_FAILURE:
    cpProcessError(cp_mem, CP_PROJ_FAILURE, "CPODES", "CPode", MSGCP_PROJ_FAILS, tn, h);
    break;
  case CP_CNSTRFUNC_FAIL:
    cpProcessError(cp_mem, CP_CNSTRFUNC_FAIL, "CPODES", "CPode", MSGCP_CNSTRFUNC_FAILED, tn);
    break;
  case CP_REPTD_CNSTRFUNC_ERR:
    cpProcessError(cp_mem, CP_REPTD_CNSTRFUNC_ERR, "CPODES", "CPode", MSGCP_CNSTRFUNC_REPTD, tn);
    break;
  case CP_PROJFUNC_FAIL:
    cpProcessError(cp_mem, CP_PROJFUNC_FAIL, "CPODES", "CPode", MSGCP_PROJFUNC_FAILED, tn);
    break;
  case CP_REPTD_PROJFUNC_ERR:
    cpProcessError(cp_mem, CP_REPTD_PROJFUNC_ERR, "CPODES", "CPode", MSGCP_PROJFUNC_REPTD, tn);
    break;
  case CP_PLSETUP_FAIL:
    cpProcessError(cp_mem, CP_PLSETUP_FAIL, "CPODES", "CPode", MSGCP_PLSETUP_FAILED, tn);
    break;
  case CP_PLSOLVE_FAIL:
    cpProcessError(cp_mem, CP_PLSOLVE_FAIL, "CPODES", "CPode", MSGCP_PLSOLVE_FAILED, tn);
    break;
  default:
    return(CP_SUCCESS);   
  }

  return(flag);

}

/* 
 * -----------------------------------------------------------------
 * BDF stability limit detection functions
 * -----------------------------------------------------------------
 */

/*
 * cpBDFStab
 *
 * This routine handles the BDF Stability Limit Detection Algorithm
 * STALD.  It is called if lmm_type = CP_BDF and the SLDET option is on.
 * If the order is 3 or more, the required norm data is saved.
 * If a decision to reduce order has not already been made, and
 * enough data has been saved, cpSLdet is called.  If it signals
 * a stability limit violation, the order is reduced, and the step
 * size is reset accordingly.
 */

void cpBDFStab(CPodeMem cp_mem)
{
  int i,k, ldflag, factorial;
  realtype sq, sqm1, sqm2;
      
  /* If order is 3 or greater, then save scaled derivative data,
     push old data down in i, then add current values to top.    */

  if (q >= 3) {
    for (k = 1; k <= 3; k++)
      { for (i = 5; i >= 2; i--) ssdat[i][k] = ssdat[i-1][k]; }
    factorial = 1;
    for (i = 1; i <= q-1; i++) factorial *= i;
    sq = factorial*q*(q+1)*acnrm/tq[5];
    sqm1 = factorial*q*N_VWrmsNorm(zn[q], ewt);
    sqm2 = factorial*N_VWrmsNorm(zn[q-1], ewt);
    ssdat[1][1] = sqm2*sqm2;
    ssdat[1][2] = sqm1*sqm1;
    ssdat[1][3] = sq*sq;
  }  

  if (qprime >= q) {

    /* If order is 3 or greater, and enough ssdat has been saved,
       nscon >= q+5, then call stability limit detection routine.  */

    if ( (q >= 3) && (nscon >= q+5) ) {
      ldflag = cpSLdet(cp_mem);
      if (ldflag > 3) {
        /* A stability limit violation is indicated by
           a return flag of 4, 5, or 6.
           Reduce new order.                     */
        qprime = q-1;
        eta = etaqm1; 
        eta = MIN(eta,etamax);
        eta = eta/MAX(ONE,ABS(h)*hmax_inv*eta);
        hprime = h*eta;
        nor = nor + 1;
      }
    }
  }
  else {
    /* Otherwise, let order increase happen, and 
       reset stability limit counter, nscon.     */
    nscon = 0;
  }
}

/*
 * cpSLdet
 *
 * This routine detects stability limitation using stored scaled 
 * derivatives data. cpSLdet returns the magnitude of the
 * dominate characteristic root, rr. The presents of a stability
 * limit is indicated by rr > "something a little less then 1.0",  
 * and a positive kflag. This routine should only be called if
 * order is greater than or equal to 3, and data has been collected
 * for 5 time steps. 
 * 
 * Returned values:
 *    kflag = 1 -> Found stable characteristic root, normal matrix case
 *    kflag = 2 -> Found stable characteristic root, quartic solution
 *    kflag = 3 -> Found stable characteristic root, quartic solution,
 *                 with Newton correction
 *    kflag = 4 -> Found stability violation, normal matrix case
 *    kflag = 5 -> Found stability violation, quartic solution
 *    kflag = 6 -> Found stability violation, quartic solution,
 *                 with Newton correction
 *
 *    kflag < 0 -> No stability limitation, 
 *                 or could not compute limitation.
 *
 *    kflag = -1 -> Min/max ratio of ssdat too small.
 *    kflag = -2 -> For normal matrix case, vmax > vrrt2*vrrt2
 *    kflag = -3 -> For normal matrix case, The three ratios
 *                  are inconsistent.
 *    kflag = -4 -> Small coefficient prevents elimination of quartics.  
 *    kflag = -5 -> R value from quartics not consistent.
 *    kflag = -6 -> No corrected root passes test on qk values
 *    kflag = -7 -> Trouble solving for sigsq.
 *    kflag = -8 -> Trouble solving for B, or R via B.
 *    kflag = -9 -> R via sigsq[k] disagrees with R from data.
 */

static int cpSLdet(CPodeMem cp_mem)
{
  int i, k, j, it, kmin, kflag = 0;
  realtype rat[5][4], rav[4], qkr[4], sigsq[4], smax[4], ssmax[4];
  realtype drr[4], rrc[4],sqmx[4], qjk[4][4], vrat[5], qc[6][4], qco[6][4];
  realtype rr, rrcut, vrrtol, vrrt2, sqtol, rrtol;
  realtype smink, smaxk, sumrat, sumrsq, vmin, vmax, drrmax, adrr;
  realtype tem, sqmax, saqk, qp, s, sqmaxk, saqj, sqmin;
  realtype rsa, rsb, rsc, rsd, rd1a, rd1b, rd1c;
  realtype rd2a, rd2b, rd3a, cest1, corr1; 
  realtype ratp, ratm, qfac1, qfac2, bb, rrb;

  /* The following are cutoffs and tolerances used by this routine */

  rrcut  = RCONST(0.98);
  vrrtol = RCONST(1.0e-4);
  vrrt2  = RCONST(5.0e-4);
  sqtol  = RCONST(1.0e-3);
  rrtol  = RCONST(1.0e-2);
  
  rr = ZERO;
  
  /*  Index k corresponds to the degree of the interpolating polynomial. */
  /*      k = 1 -> q-1          */
  /*      k = 2 -> q            */
  /*      k = 3 -> q+1          */
  
  /*  Index i is a backward-in-time index, i = 1 -> current time, */
  /*      i = 2 -> previous step, etc    */
  
  /* get maxima, minima, and variances, and form quartic coefficients  */
  
  for (k=1; k<=3; k++) {
    smink = ssdat[1][k];
    smaxk = ZERO;
    
    for (i=1; i<=5; i++) {
      smink = MIN(smink,ssdat[i][k]);
      smaxk = MAX(smaxk,ssdat[i][k]);
    }
    
    if (smink < TINY*smaxk) {
      kflag = -1;  
      return(kflag);
    }
    smax[k] = smaxk;
    ssmax[k] = smaxk*smaxk;
    
    sumrat = ZERO;
    sumrsq = ZERO;
    for (i=1; i<=4; i++) {
      rat[i][k] = ssdat[i][k]/ssdat[i+1][k];
      sumrat = sumrat + rat[i][k];
      sumrsq = sumrsq + rat[i][k]*rat[i][k];
    } 
    rav[k] = PT25 * sumrat;
    vrat[k] = ABS(PT25*sumrsq - rav[k]*rav[k]);
    
    qc[5][k] = ssdat[1][k]*ssdat[3][k] - ssdat[2][k]*ssdat[2][k];
    qc[4][k] = ssdat[2][k]*ssdat[3][k] - ssdat[1][k]*ssdat[4][k];
    qc[3][k] = ZERO;
    qc[2][k] = ssdat[2][k]*ssdat[5][k] - ssdat[3][k]*ssdat[4][k];
    qc[1][k] = ssdat[4][k]*ssdat[4][k] - ssdat[3][k]*ssdat[5][k];
    
    for (i=1; i<=5; i++) {
      qco[i][k] = qc[i][k];
    }
  }                            /* End of k loop */
  
  /* Isolate normal or nearly-normal matrix case. Three quartic will
     have common or nearly-common roots in this case. 
     Return a kflag = 1 if this procedure works. If three root 
     differ more than vrrt2, return error kflag = -3.    */
  
  vmin = MIN(vrat[1],MIN(vrat[2],vrat[3]));
  vmax = MAX(vrat[1],MAX(vrat[2],vrat[3]));
  
  if(vmin < vrrtol*vrrtol) {
    if (vmax > vrrt2*vrrt2) {
      kflag = -2;  
      return(kflag);
    } else {
      rr = (rav[1] + rav[2] + rav[3])/THREE;
      
      drrmax = ZERO;
      for(k = 1;k<=3;k++) {
        adrr = ABS(rav[k] - rr);
        drrmax = MAX(drrmax, adrr);
      }
      if (drrmax > vrrt2) {
        kflag = -3;    
      }
      
      kflag = 1;

      /*  can compute charactistic root, drop to next section   */
      
    }
  } else {

    /* use the quartics to get rr. */
    
    if (ABS(qco[1][1]) < TINY*ssmax[1]) {
      kflag = -4;    
      return(kflag);
    }
    
    tem = qco[1][2]/qco[1][1];
    for(i=2; i<=5; i++) {
      qco[i][2] = qco[i][2] - tem*qco[i][1];
    }

    qco[1][2] = ZERO;
    tem = qco[1][3]/qco[1][1];
    for(i=2; i<=5; i++) {
      qco[i][3] = qco[i][3] - tem*qco[i][1];
    }
    qco[1][3] = ZERO;
    
    if (ABS(qco[2][2]) < TINY*ssmax[2]) {
      kflag = -4;    
      return(kflag);
    }
    
    tem = qco[2][3]/qco[2][2];
    for(i=3; i<=5; i++) {
      qco[i][3] = qco[i][3] - tem*qco[i][2];
    }
    
    if (ABS(qco[4][3]) < TINY*ssmax[3]) {
      kflag = -4;    
      return(kflag);
    }
    
    rr = -qco[5][3]/qco[4][3];
    
    if (rr < TINY || rr > HUNDRED) {
      kflag = -5;   
      return(kflag);
    }
    
    for(k=1; k<=3; k++) {
      qkr[k] = qc[5][k] + rr*(qc[4][k] + rr*rr*(qc[2][k] + rr*qc[1][k]));
    }  
    
    sqmax = ZERO;
    for(k=1; k<=3; k++) {
      saqk = ABS(qkr[k])/ssmax[k];
      if (saqk > sqmax) sqmax = saqk;
    } 
    
    if (sqmax < sqtol) {
      kflag = 2;
      
      /*  can compute charactistic root, drop to "given rr,etc"   */
      
    } else {

      /* do Newton corrections to improve rr.  */
      
      for(it=1; it<=3; it++) {
        for(k=1; k<=3; k++) {
          qp = qc[4][k] + rr*rr*(THREE*qc[2][k] + rr*FOUR*qc[1][k]);
          drr[k] = ZERO;
          if (ABS(qp) > TINY*ssmax[k]) drr[k] = -qkr[k]/qp;
          rrc[k] = rr + drr[k];
        } 
        
        for(k=1; k<=3; k++) {
          s = rrc[k];
          sqmaxk = ZERO;
          for(j=1; j<=3; j++) {
            qjk[j][k] = qc[5][j] + s*(qc[4][j] + 
                                      s*s*(qc[2][j] + s*qc[1][j]));
            saqj = ABS(qjk[j][k])/ssmax[j];
            if (saqj > sqmaxk) sqmaxk = saqj;
          } 
          sqmx[k] = sqmaxk;
        } 

        sqmin = sqmx[1]; kmin = 1;
        for(k=2; k<=3; k++) {
          if (sqmx[k] < sqmin) {
            kmin = k;
            sqmin = sqmx[k];
          }
        } 
        rr = rrc[kmin];
        
        if (sqmin < sqtol) {
          kflag = 3;
          /*  can compute charactistic root   */
          /*  break out of Newton correction loop and drop to "given rr,etc" */ 
          break;
        } else {
          for(j=1; j<=3; j++) {
            qkr[j] = qjk[j][kmin];
          }
        }     
      }          /*  end of Newton correction loop  */ 
      
      if (sqmin > sqtol) {
        kflag = -6;
        return(kflag);
      }
    }     /*  end of if (sqmax < sqtol) else   */
  }      /*  end of if(vmin < vrrtol*vrrtol) else, quartics to get rr. */
  
  /* given rr, find sigsq[k] and verify rr.  */
  /* All positive kflag drop to this section  */
  
  for(k=1; k<=3; k++) {
    rsa = ssdat[1][k];
    rsb = ssdat[2][k]*rr;
    rsc = ssdat[3][k]*rr*rr;
    rsd = ssdat[4][k]*rr*rr*rr;
    rd1a = rsa - rsb;
    rd1b = rsb - rsc;
    rd1c = rsc - rsd;
    rd2a = rd1a - rd1b;
    rd2b = rd1b - rd1c;
    rd3a = rd2a - rd2b;
    
    if (ABS(rd1b) < TINY*smax[k]) {
      kflag = -7;
      return(kflag);
    }
    
    cest1 = -rd3a/rd1b;
    if (cest1 < TINY || cest1 > FOUR) {
      kflag = -7;
      return(kflag);
    }
    corr1 = (rd2b/cest1)/(rr*rr);
    sigsq[k] = ssdat[3][k] + corr1;
  }
  
  if (sigsq[2] < TINY) {
    kflag = -8;
    return(kflag);
  }
  
  ratp = sigsq[3]/sigsq[2];
  ratm = sigsq[1]/sigsq[2];
  qfac1 = PT25 * (q*q - ONE);
  qfac2 = TWO/(q - ONE);
  bb = ratp*ratm - ONE - qfac1*ratp;
  tem = ONE - qfac2*bb;
  
  if (ABS(tem) < TINY) {
    kflag = -8;
    return(kflag);
  }
  
  rrb = ONE/tem;
  
  if (ABS(rrb - rr) > rrtol) {
    kflag = -9;
    return(kflag);
  }
  
  /* Check to see if rr is above cutoff rrcut  */
  if (rr > rrcut) {
    if (kflag == 1) kflag = 4;
    if (kflag == 2) kflag = 5;
    if (kflag == 3) kflag = 6;
  }
  
  /* All positive kflag returned at this point  */
  
  return(kflag);
  
}

/* 
 * -----------------------------------------------------------------
 * Internal error weight evaluation functions
 * -----------------------------------------------------------------
 */

/*
 * cpEwtSet
 *
 * This routine is responsible for setting the error weight vector ewt,
 * according to tol_type, as follows:
 *
 * (1) ewt[i] = 1 / (reltol * ABS(ycur[i]) + *abstol), i=0,...,neq-1
 *     if tol_type = CP_SS
 * (2) ewt[i] = 1 / (reltol * ABS(ycur[i]) + abstol[i]), i=0,...,neq-1
 *     if tol_type = CP_SV
 *
 * cpEwtSet returns 0 if ewt is successfully set as above to a
 * positive vector and -1 otherwise. In the latter case, ewt is
 * considered undefined.
 *
 * All the real work is done in the routines cpEwtSetSS, cpEwtSetSV.
 */

int cpEwtSet(N_Vector ycur, N_Vector weight, void *edata)
{
  CPodeMem cp_mem;
  int flag = 0;

  /* edata points to cp_mem here */
  cp_mem = (CPodeMem) edata;

  switch(tol_type) {
  case CP_SS: 
    flag = cpEwtSetSS(cp_mem, ycur, weight);
    break;
  case CP_SV: 
    flag = cpEwtSetSV(cp_mem, ycur, weight);
    break;
  }
  
  return(flag);
}

/*
 * cpEwtSetSS
 *
 * This routine sets ewt as decribed above in the case tol_type = CP_SS.
 * It tests for non-positive components before inverting. cpEwtSetSS
 * returns 0 if ewt is successfully set to a positive vector
 * and -1 otherwise. In the latter case, ewt is considered undefined.
 */

static int cpEwtSetSS(CPodeMem cp_mem, N_Vector ycur, N_Vector weight)
{
  N_VAbs(ycur, tempv);
  N_VScale(reltol, tempv, tempv);
  N_VAddConst(tempv, Sabstol, tempv);
  if (N_VMin(tempv) <= ZERO) return(-1);
  N_VInv(tempv, weight);
  return(0);
}

/*
 * cpEwtSetSV
 *
 * This routine sets ewt as decribed above in the case tol_type = CP_SV.
 * It tests for non-positive components before inverting. cpEwtSetSV
 * returns 0 if ewt is successfully set to a positive vector
 * and -1 otherwise. In the latter case, ewt is considered undefined.
 */

static int cpEwtSetSV(CPodeMem cp_mem, N_Vector ycur, N_Vector weight)
{
  N_VAbs(ycur, tempv);
  N_VLinearSum(reltol, tempv, ONE, Vabstol, tempv);
  if (N_VMin(tempv) <= ZERO) return(-1);
  N_VInv(tempv, weight);
  return(0);
}

/*
 * cpQuadEwtSet
 *
 */

static int cpQuadEwtSet(CPodeMem cp_mem, N_Vector qcur, N_Vector weightQ)
{
  int flag=0;

  switch (tol_typeQ) {
  case CP_SS: 
    flag = cpQuadEwtSetSS(cp_mem, qcur, weightQ);
    break;
  case CP_SV: 
    flag = cpQuadEwtSetSV(cp_mem, qcur, weightQ);
    break;
  }

  return(flag);

}

/*
 * cpQuadEwtSetSS
 *
 */

static int cpQuadEwtSetSS(CPodeMem cp_mem, N_Vector qcur, N_Vector weightQ)
{
  N_VAbs(qcur, tempvQ);
  N_VScale(reltolQ, tempvQ, tempvQ);
  N_VAddConst(tempvQ, SabstolQ, tempvQ);
  if (N_VMin(tempvQ) <= ZERO) return(-1);
  N_VInv(tempvQ, weightQ);

  return(0);
}

/*
 * cpQuadEwtSetSV
 *
 */

static int cpQuadEwtSetSV(CPodeMem cp_mem, N_Vector qcur, N_Vector weightQ)
{
  N_VAbs(qcur, tempvQ);
  N_VLinearSum(reltolQ, tempvQ, ONE, VabstolQ, tempvQ);
  if (N_VMin(tempvQ) <= ZERO) return(-1);
  N_VInv(tempvQ, weightQ);

  return(0);
}

/* 
 * -----------------------------------------------------------------
 * Updated WRMS norms
 * -----------------------------------------------------------------
 */

/*
 * cpQuadUpdateNorm
 *
 * Updates the norm old_nrm to account for all quadratures.
 */

static realtype cpQuadUpdateNorm(CPodeMem cp_mem, realtype old_nrm,
                                 N_Vector xQ, N_Vector wQ)
{
  realtype qnrm;

  qnrm = N_VWrmsNorm(xQ, wQ);
  if (old_nrm > qnrm) return(old_nrm);
  else                return(qnrm);
}

/*
 * cpQuadUpdateDsm
 *
 * This routine updates the local error norm dsm with quadrature
 * related information. Used only if quadratures are computed
 * with FULL error control.
 *
 * Returns the maximum over the wheighted local error norms.
 */

static realtype cpQuadUpdateDsm(CPodeMem cp_mem, realtype old_dsm, 
                                realtype dsmQ)
{
  if ( old_dsm > dsmQ ) return(old_dsm);
  else                  return(dsmQ);
}

/* 
 * -----------------------------------------------------------------
 * Error reporting functions
 * -----------------------------------------------------------------
 */

/* 
 * cpProcessError is a high level error handling function
 * - if cp_mem==NULL it prints the error message to stderr
 * - otherwise, it sets-up and calls the error hadling function 
 *   pointed to by cp_ehfun
 */

#define ehfun    (cp_mem->cp_ehfun)
#define eh_data  (cp_mem->cp_eh_data)

void cpProcessError(CPodeMem cp_mem, 
                    int error_code, const char *module, const char *fname, 
                    const char *msgfmt, ...)
{
  va_list ap;
  char msg[256];

  /* Initialize the argument pointer variable 
     (msgfmt is the last required argument to cpProcessError) */

  va_start(ap, msgfmt);

  if (cp_mem == NULL) {    /* We write to stderr */

#ifndef NO_FPRINTF_OUTPUT
    fprintf(stderr, "\n[%s ERROR]  %s\n  ", module, fname);
    fprintf(stderr, msgfmt, 0); // extra arg avoids gcc warnings
    fprintf(stderr, "\n\n");
#endif

  } else {                 /* We can call ehfun */

    /* Compose the message */

    vsprintf(msg, msgfmt, ap);

    /* Call ehfun */

    ehfun(error_code, module, fname, msg, eh_data);

  }

  /* Finalize argument processing */
  
  va_end(ap);

  return;

}

/* cpErrHandler is the default error handling function.
   It sends the error message to the stream pointed to by cp_errfp */

#define errfp    (cp_mem->cp_errfp)

static void cpErrHandler(int error_code, const char *module,
                         const char *function, char *msg, void *data)
{
  CPodeMem cp_mem;
  char err_type[10];

  /* data points to cp_mem here */

  cp_mem = (CPodeMem) data;

  if (error_code == CP_WARNING)
    sprintf(err_type,"WARNING");
  else
    sprintf(err_type,"ERROR");

#ifndef NO_FPRINTF_OUTPUT
  if (errfp!=NULL) {
    fprintf(errfp,"\n[%s %s]  %s\n",module,err_type,function);
    fprintf(errfp,"  %s\n\n",msg);
  }
#endif

  return;
}
