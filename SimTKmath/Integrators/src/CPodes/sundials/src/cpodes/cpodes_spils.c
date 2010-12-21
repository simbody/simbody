/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2006/11/24 19:09:18 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for the CPSPILS linear solvers.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "cpodes_private.h"
#include "cpodes_spils_impl.h"

/* Private constants */

#define ZERO   RCONST(0.0)
#define PT25   RCONST(0.25)
#define ONE    RCONST(1.0)

/* Algorithmic constants */

#define MAX_ITERS  3  /* max. number of attempts to recover in DQ J*v */

/* Readability Replacements */

#define ode_type  (cp_mem->cp_ode_type)
#define lrw1      (cp_mem->cp_lrw1)
#define liw1      (cp_mem->cp_liw1)
#define tq        (cp_mem->cp_tq)
#define tn        (cp_mem->cp_tn)
#define h         (cp_mem->cp_h)
#define gamma     (cp_mem->cp_gamma)
#define fe        (cp_mem->cp_fe)
#define fi        (cp_mem->cp_fi)
#define f_data    (cp_mem->cp_f_data)
#define ewt       (cp_mem->cp_ewt)
#define lmem      (cp_mem->cp_lmem)

#define ils_type  (cpspils_mem->s_type)
#define sqrtN     (cpspils_mem->s_sqrtN)   
#define ytemp     (cpspils_mem->s_ytemp)
#define yptemp    (cpspils_mem->s_yptemp)
#define x         (cpspils_mem->s_x)
#define ycur      (cpspils_mem->s_ycur)
#define ypcur     (cpspils_mem->s_ypcur)
#define fcur      (cpspils_mem->s_fcur)
#define delta     (cpspils_mem->s_delta)
#define npe       (cpspils_mem->s_npe)
#define nli       (cpspils_mem->s_nli)
#define nps       (cpspils_mem->s_nps)
#define ncfl      (cpspils_mem->s_ncfl)
#define njtimes   (cpspils_mem->s_njtimes)
#define nfes      (cpspils_mem->s_nfes)
#define last_flag (cpspils_mem->s_last_flag)

/*
 * -----------------------------------------------------------------
 * OPTIONAL INPUT and OUTPUT
 * -----------------------------------------------------------------
 */


/*
 * -----------------------------------------------------------------
 * CPSpilsSetPrecType
 * -----------------------------------------------------------------
 */

int CPSpilsSetPrecType(void *cpode_mem, int pretype)
{
  CPodeMem cp_mem;
  CPSpilsMem cpspils_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPSPILS_MEM_NULL, "CPSPILS", "CPSpilsSetPrecType", MSGS_CPMEM_NULL);
    return(CPSPILS_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(cp_mem, CPSPILS_LMEM_NULL, "CPSPILS", "CPSpilsSetPrecType", MSGS_LMEM_NULL);
    return(CPSPILS_LMEM_NULL);
  }
  cpspils_mem = (CPSpilsMem) lmem;

  /* Check for legal pretype */ 
  if ((pretype != PREC_NONE) && (pretype != PREC_LEFT) &&
      (pretype != PREC_RIGHT) && (pretype != PREC_BOTH)) {
    cpProcessError(cp_mem, CPSPILS_ILL_INPUT, "CPSPILS", "CPSpilsSetPrecType", MSGS_BAD_PRETYPE);
    return(CPSPILS_ILL_INPUT);
  }

  cpspils_mem->s_pretype = pretype;

  return(CPSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CPSpilsSetGSType
 * -----------------------------------------------------------------
 */

int CPSpilsSetGSType(void *cpode_mem, int gstype)
{
  CPodeMem cp_mem;
  CPSpilsMem cpspils_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPSPILS_MEM_NULL, "CPSPILS", "CPSpilsSetGSType", MSGS_CPMEM_NULL);
    return(CPSPILS_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(cp_mem, CPSPILS_LMEM_NULL, "CPSPILS", "CPSpilsSetGSType", MSGS_LMEM_NULL);
    return(CPSPILS_LMEM_NULL);
  }
  cpspils_mem = (CPSpilsMem) lmem;

  if (ils_type != SPILS_SPGMR) {
    cpProcessError(cp_mem, CPSPILS_ILL_INPUT, "CPSPILS", "CPSpilsSetGSType", MSGS_BAD_LSTYPE);
    return(CPSPILS_ILL_INPUT);
  }

  /* Check for legal gstype */
  if ((gstype != MODIFIED_GS) && (gstype != CLASSICAL_GS)) {
    cpProcessError(cp_mem, CPSPILS_ILL_INPUT, "CPSPILS", "CPSpilsSetGSType", MSGS_BAD_GSTYPE);
    return(CPSPILS_ILL_INPUT);
  }

  cpspils_mem->s_gstype = gstype;

  return(CPSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : CPSpilsSetMaxl
 * -----------------------------------------------------------------
 */

int CPSpilsSetMaxl(void *cpode_mem, int maxl)
{
  CPodeMem cp_mem;
  CPSpilsMem cpspils_mem;
  int mxl;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPSPILS_MEM_NULL, "CPSPILS", "CPSpilsSetMaxl", MSGS_CPMEM_NULL);
    return(CPSPILS_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(NULL, CPSPILS_LMEM_NULL, "CPSPILS", "CPSpilsSetMaxl", MSGS_LMEM_NULL);
    return(CPSPILS_LMEM_NULL);
  }
  cpspils_mem = (CPSpilsMem) lmem;

  if (ils_type == SPILS_SPGMR) {
    cpProcessError(cp_mem, CPSPILS_ILL_INPUT, "CPSPILS", "CPSpilsSetMaxl", MSGS_BAD_LSTYPE);
    return(CPSPILS_ILL_INPUT);
  }

  mxl = (maxl <= 0) ? CPSPILS_MAXL : maxl;
  cpspils_mem->s_maxl = mxl;

  return(CPSPILS_SUCCESS);
}


/*
 * -----------------------------------------------------------------
 * CPSpilsSetDelt
 * -----------------------------------------------------------------
 */

int CPSpilsSetDelt(void *cpode_mem, realtype delt)
{
  CPodeMem cp_mem;
  CPSpilsMem cpspils_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPSPILS_MEM_NULL, "CPSPILS", "CPSpilsSetDelt", MSGS_CPMEM_NULL);
    return(CPSPILS_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(cp_mem, CPSPILS_LMEM_NULL, "CPSPILS", "CPSpilsSetDelt", MSGS_LMEM_NULL);
    return(CPSPILS_LMEM_NULL);
  }
  cpspils_mem = (CPSpilsMem) lmem;

  /* Check for legal delt */
  if(delt < ZERO) {
    cpProcessError(cp_mem, CPSPILS_ILL_INPUT, "CPSPILS", "CPSpilsSetDelt", MSGS_BAD_DELT);
    return(CPSPILS_ILL_INPUT);
  }

  cpspils_mem->s_delt = (delt == ZERO) ? CPSPILS_DELT : delt;

  return(CPSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CPSpilsSetPrecSetupFn
 * -----------------------------------------------------------------
 */

int CPSpilsSetPreconditioner(void *cpode_mem, void *pset, void *psolve, void *P_data)
{
  CPodeMem cp_mem;
  CPSpilsMem cpspils_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPSPILS_MEM_NULL, "CPSPILS", "CPSpilsSetPreconditioner", MSGS_CPMEM_NULL);
    return(CPSPILS_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(cp_mem, CPSPILS_LMEM_NULL, "CPSPILS", "CPSpilsSetPreconditioner", MSGS_LMEM_NULL);
    return(CPSPILS_LMEM_NULL);
  }
  cpspils_mem = (CPSpilsMem) lmem;

  if (ode_type == CP_EXPL) {
    cpspils_mem->s_psetE = (CPSpilsPrecSetupExplFn) pset;
    cpspils_mem->s_pslvE = (CPSpilsPrecSolveExplFn) psolve;
  } else {
    cpspils_mem->s_psetI = (CPSpilsPrecSetupImplFn) pset;
    cpspils_mem->s_pslvI = (CPSpilsPrecSolveImplFn) psolve;
  }
  cpspils_mem->s_P_data = P_data;

  return(CPSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CPSpilsSetJacTimesVecFn
 * -----------------------------------------------------------------
 */

int CPSpilsSetJacTimesVecFn(void *cpode_mem, void *jtimes, void *jac_data)
{
  CPodeMem cp_mem;
  CPSpilsMem cpspils_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPSPILS_MEM_NULL, "CPSPILS", "CPSpilsSetJacTimesVecFn", MSGS_CPMEM_NULL);
    return(CPSPILS_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(cp_mem, CPSPILS_LMEM_NULL, "CPSPILS", "CPSpilsSetJacTimesVecFn", MSGS_LMEM_NULL);
    return(CPSPILS_LMEM_NULL);
  }
  cpspils_mem = (CPSpilsMem) lmem;

  if (ode_type == CP_EXPL) cpspils_mem->s_jtvE = (CPSpilsJacTimesVecExplFn) jtimes;
  else                     cpspils_mem->s_jtvI = (CPSpilsJacTimesVecImplFn) jtimes;
  cpspils_mem->s_j_data = jac_data;

  return(CPSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CPSpilsGetWorkSpace
 * -----------------------------------------------------------------
 */

int CPSpilsGetWorkSpace(void *cpode_mem, long int *lenrwLS, long int *leniwLS)
{
  CPodeMem cp_mem;
  CPSpilsMem cpspils_mem;
  int maxl;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPSPILS_MEM_NULL, "CPSPILS", "CPSpilsGetWorkSpace", MSGS_CPMEM_NULL);
    return(CPSPILS_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(cp_mem, CPSPILS_LMEM_NULL, "CPSPILS", "CPSpilsGetWorkSpace", MSGS_LMEM_NULL);
    return(CPSPILS_LMEM_NULL);
  }
  cpspils_mem = (CPSpilsMem) lmem;

  
  switch(ils_type) {
  case SPILS_SPGMR:
    maxl = cpspils_mem->s_maxl;
    *lenrwLS = lrw1*(maxl + 5) + maxl*(maxl + 4) + 1;
    *leniwLS = liw1*(maxl + 5);
    break;
  case SPILS_SPBCG:
    *lenrwLS = lrw1 * 9;
    *leniwLS = liw1 * 9;
    break;
  case SPILS_SPTFQMR:
    *lenrwLS = lrw1*11;
    *leniwLS = liw1*11;
    break;
  }


  return(CPSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CPSpilsGetNumPrecEvals
 * -----------------------------------------------------------------
 */

int CPSpilsGetNumPrecEvals(void *cpode_mem, long int *npevals)
{
  CPodeMem cp_mem;
  CPSpilsMem cpspils_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPSPILS_MEM_NULL, "CPSPILS", "CPSpilsGetNumPrecEvals", MSGS_CPMEM_NULL);
    return(CPSPILS_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(cp_mem, CPSPILS_LMEM_NULL, "CPSPILS", "CPSpilsGetNumPrecEvals", MSGS_LMEM_NULL);
    return(CPSPILS_LMEM_NULL);
  }
  cpspils_mem = (CPSpilsMem) lmem;

  *npevals = npe;

  return(CPSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CPSpilsGetNumPrecSolves
 * -----------------------------------------------------------------
 */

int CPSpilsGetNumPrecSolves(void *cpode_mem, long int *npsolves)
{
  CPodeMem cp_mem;
  CPSpilsMem cpspils_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPSPILS_MEM_NULL, "CPSPILS", "CPSpilsGetNumPrecSolves", MSGS_CPMEM_NULL);
    return(CPSPILS_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(cp_mem, CPSPILS_LMEM_NULL, "CPSPILS", "CPSpilsGetNumPrecSolves", MSGS_LMEM_NULL);
    return(CPSPILS_LMEM_NULL);
  }
  cpspils_mem = (CPSpilsMem) lmem;

  *npsolves = nps;

  return(CPSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CPSpilsGetNumLinIters
 * -----------------------------------------------------------------
 */

int CPSpilsGetNumLinIters(void *cpode_mem, long int *nliters)
{
  CPodeMem cp_mem;
  CPSpilsMem cpspils_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPSPILS_MEM_NULL, "CPSPILS", "CPSpilsGetNumLinIters", MSGS_CPMEM_NULL);
    return(CPSPILS_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(cp_mem, CPSPILS_LMEM_NULL, "CPSPILS", "CPSpilsGetNumLinIters", MSGS_LMEM_NULL);
    return(CPSPILS_LMEM_NULL);
  }
  cpspils_mem = (CPSpilsMem) lmem;

  *nliters = nli;

  return(CPSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CPSpilsGetNumConvFails
 * -----------------------------------------------------------------
 */

int CPSpilsGetNumConvFails(void *cpode_mem, long int *nlcfails)
{
  CPodeMem cp_mem;
  CPSpilsMem cpspils_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPSPILS_MEM_NULL, "CPSPILS", "CPSpilsGetNumConvFails", MSGS_CPMEM_NULL);
    return(CPSPILS_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(cp_mem, CPSPILS_LMEM_NULL, "CPSPILS", "CPSpilsGetNumConvFails", MSGS_LMEM_NULL);
    return(CPSPILS_LMEM_NULL);
  }
  cpspils_mem = (CPSpilsMem) lmem;

  *nlcfails = ncfl;

  return(CPSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CPSpilsGetNumJtimesEvals
 * -----------------------------------------------------------------
 */

int CPSpilsGetNumJtimesEvals(void *cpode_mem, long int *njvevals)
{
  CPodeMem cp_mem;
  CPSpilsMem cpspils_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPSPILS_MEM_NULL, "CPSPILS", "CPSpilsGetNumJtimesEvals", MSGS_CPMEM_NULL);
    return(CPSPILS_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(cp_mem, CPSPILS_LMEM_NULL, "CPSPILS", "CPSpilsGetNumJtimesEvals", MSGS_LMEM_NULL);
    return(CPSPILS_LMEM_NULL);
  }
  cpspils_mem = (CPSpilsMem) lmem;

  *njvevals = njtimes;

  return(CPSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CPSpilsGetNumFctEvals
 * -----------------------------------------------------------------
 */

int CPSpilsGetNumFctEvals(void *cpode_mem, long int *nfevalsLS)
{
  CPodeMem cp_mem;
  CPSpilsMem cpspils_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPSPILS_MEM_NULL, "CPSPILS", "CPSpilsGetNumRhsEvals", MSGS_CPMEM_NULL);
    return(CPSPILS_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(cp_mem, CPSPILS_LMEM_NULL, "CPSPILS", "CPSpilsGetNumRhsEvals", MSGS_LMEM_NULL);
    return(CPSPILS_LMEM_NULL);
  }
  cpspils_mem = (CPSpilsMem) lmem;

  *nfevalsLS = nfes;

  return(CPSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CPSpilsGetLastFlag
 * -----------------------------------------------------------------
 */

int CPSpilsGetLastFlag(void *cpode_mem, int *flag)
{
  CPodeMem cp_mem;
  CPSpilsMem cpspils_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPSPILS_MEM_NULL, "CPSPILS", "CPSpilsGetLastFlag", MSGS_CPMEM_NULL);
    return(CPSPILS_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(cp_mem, CPSPILS_LMEM_NULL, "CPSPILS", "CPSpilsGetLastFlag", MSGS_LMEM_NULL);
    return(CPSPILS_LMEM_NULL);
  }
  cpspils_mem = (CPSpilsMem) lmem;

  *flag = last_flag;

  return(CPSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CPSpilsGetReturnFlagName
 * -----------------------------------------------------------------
 */

char *CPSpilsGetReturnFlagName(int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case CPSPILS_SUCCESS:
    sprintf(name,"CPSPILS_SUCCESS");
    break; 
  case CPSPILS_MEM_NULL:
    sprintf(name,"CPSPILS_MEM_NULL");
    break;
  case CPSPILS_LMEM_NULL:
    sprintf(name,"CPSPILS_LMEM_NULL");
    break;
  case CPSPILS_ILL_INPUT:
    sprintf(name,"CPSPILS_ILL_INPUT");
    break;
  case CPSPILS_MEM_FAIL:
    sprintf(name,"CPSPILS_MEM_FAIL");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}

/*
 * -----------------------------------------------------------------
 * CPSPILS private functions
 * -----------------------------------------------------------------
 */


/* Additional readability Replacements */

#define pretype (cpspils_mem->s_pretype)
#define delt    (cpspils_mem->s_delt)
#define maxl    (cpspils_mem->s_maxl)
#define pslvE   (cpspils_mem->s_pslvE)
#define pslvI   (cpspils_mem->s_pslvI)
#define P_data  (cpspils_mem->s_P_data)
#define jtvE    (cpspils_mem->s_jtvE)
#define jtvI    (cpspils_mem->s_jtvI)
#define j_data  (cpspils_mem->s_j_data)

/*
 * -----------------------------------------------------------------
 * cpSpilsAtimes
 * -----------------------------------------------------------------
 * This routine generates the matrix-vector product z = Mv, where
 * M = I - gamma*J. The product J*v is obtained by calling the jtimes 
 * routine. It is then scaled by -gamma and added to v to obtain M*v.
 * The return value is the same as the value returned by jtimes --
 * 0 if successful, nonzero otherwise.
 * -----------------------------------------------------------------
 */

int cpSpilsAtimes(void *cpode_mem, N_Vector v, N_Vector z)
{
  CPodeMem   cp_mem;
  CPSpilsMem cpspils_mem;
  int jtflag;

  cp_mem = (CPodeMem) cpode_mem;
  cpspils_mem = (CPSpilsMem) lmem;

  if (ode_type == CP_EXPL) {
    jtflag = jtvE(tn, ycur, fcur, v, z, j_data, ytemp);
    njtimes++;
    if (jtflag != 0) return(jtflag);
    N_VLinearSum(ONE, v, -gamma, z, z);
  } else {
    jtflag = jtvI(tn, gamma, ycur, ypcur, fcur, v, z, j_data, ytemp, yptemp);
    njtimes++;
    if (jtflag != 0) return(jtflag);
  }

  return(0);
}

/*
 * -----------------------------------------------------------------
 * cpSpilsPSolve
 * -----------------------------------------------------------------
 * This routine interfaces between the generic SpgmrSolve routine and
 * the user's psolve routine.  It passes to psolve all required state 
 * information from cpode_mem.  Its return value is the same as that
 * returned by psolve. Note that the generic SPGMR solver guarantees
 * that CPSpilsPSolve will not be called in the case in which
 * preconditioning is not done. This is the only case in which the
 * user's psolve routine is allowed to be NULL.
 * -----------------------------------------------------------------
 */

int cpSpilsPSolve(void *cpode_mem, N_Vector r, N_Vector z, int lr)
{
  CPodeMem   cp_mem;
  CPSpilsMem cpspils_mem;
  int retval;

  cp_mem = (CPodeMem) cpode_mem;
  cpspils_mem = (CPSpilsMem)lmem;

  /* This call is counted in nps within the CPSp***Solve routine */
  if (ode_type == CP_EXPL) {
    retval = pslvE(tn, ycur, fcur, r, z, gamma, delta, lr, P_data, ytemp);
  } else {
    retval = pslvI(tn, ycur, ypcur, fcur, r, z, gamma, delta, P_data, ytemp);
  }

  return(retval);     
}

/*
 * -----------------------------------------------------------------
 * cpSpilsDQjtvExpl
 * -----------------------------------------------------------------
 * This routine generates a difference quotient approximation to
 * the Jacobian times vector for explicit ODE: Jv = f_y(t,y) * v.
 * The approximation is Jv = vnrm[f(y + v/vnrm) - f(y)], where 
 * vnrm = (WRMS norm of v) is input, i.e. WRMS norm of v/vnrm is 1.
 * -----------------------------------------------------------------
 */

int cpSpilsDQjtvExpl(realtype t, N_Vector y, N_Vector fy, 
                     N_Vector v, N_Vector Jv, void *jac_data, 
                     N_Vector tmp)
{
  CPodeMem cp_mem;
  CPSpilsMem cpspils_mem;
  realtype sig, siginv;
  int iter, retval;

  /* jac_data is cpode_mem */
  cp_mem = (CPodeMem) jac_data;
  cpspils_mem = (CPSpilsMem) lmem;

  /* Initialize perturbation to 1/||v|| */
  sig = ONE/N_VWrmsNorm(v, ewt);

  for (iter=0; iter<MAX_ITERS; iter++) {

    /* Set tmp = y + sig*v */
    N_VLinearSum(sig, v, ONE, y, tmp);

    /* Set Jv = f(tn, y+sig*v) */
    retval = fe(t, tmp, Jv, f_data); 
    nfes++;
    if (retval == 0) break;
    if (retval < 0)  return(-1);

    sig *= PT25;
  }

  if (retval > 0) return(+1);

  /* Replace Jv by (Jv - fy)/sig */
  siginv = ONE/sig;
  N_VLinearSum(siginv, Jv, -siginv, fy, Jv);

  return(0);
}

/*
 * -----------------------------------------------------------------
 * cpSpilsDQjtvImpl
 * -----------------------------------------------------------------
 */
int cpSpilsDQjtvImpl(realtype t, realtype gm, 
                     N_Vector y, N_Vector yp, N_Vector r,
                     N_Vector v, N_Vector Jv, void *jac_data,
                     N_Vector tmp1, N_Vector tmp2)
{
  CPodeMem cp_mem;
  CPSpilsMem cpspils_mem;
  N_Vector y_tmp, yp_tmp;
  realtype sig, siginv;
  int iter, retval;

  realtype dqincfac;

  /* jac_data is cp_mem */
  cp_mem = (CPodeMem) jac_data;
  cpspils_mem = (CPSpilsMem) lmem;

  
  dqincfac = ONE;


  switch(ils_type) {
  case SPILS_SPGMR:
    sig = sqrtN*dqincfac;
    break;
  case SPILS_SPBCG:
    sig = dqincfac/N_VWrmsNorm(v, ewt);
    break;
  case SPILS_SPTFQMR:
    sig = dqincfac/N_VWrmsNorm(v, ewt);
    break;
  }

  /* Rename tmp1 and tmp2 for readibility */
  y_tmp  = tmp1;
  yp_tmp = tmp2;

  for (iter=0; iter<MAX_ITERS; iter++) {

    /* Set y_tmp = y + gm*sig*v, yp_tmp = yp + sig*v. */
    N_VLinearSum(gm*sig, v, ONE, y, y_tmp);
    N_VLinearSum(sig, v, ONE, yp, yp_tmp);
    
    /* Call res for Jv = F(t, y_tmp, yp_tmp), and return if it failed. */
    retval = fi(t, y_tmp, yp_tmp, Jv, f_data);
    nfes++;
    if (retval == 0) break;
    if (retval < 0)  return(-1);

    sig *= PT25;
  }

  if (retval > 0) return(+1);

  /* Set Jv to [Jv - r]/sig and return. */
  siginv = ONE/sig;
  N_VLinearSum(siginv, Jv, -siginv, r, Jv);

  return(0);

}
