/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2006/11/29 00:05:08 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Common implementation header file for the scaled, preconditioned
 * linear solver modules.
 * -----------------------------------------------------------------
 */

#ifndef _CPSPILS_IMPL_H
#define _CPSPILS_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <cpodes/cpodes_spils.h>

/*
 * =================================================================
 *   C P S P I L S    I N T E R N A L    C O N S T A N T S
 * =================================================================
 */

/* 
 * -----------------------------------------------------------------
 * Types of iterative linear solvers 
 * -----------------------------------------------------------------
 */

#define SPILS_SPGMR   1
#define SPILS_SPBCG   2
#define SPILS_SPTFQMR 3

/*
 * -----------------------------------------------------------------
 * CPSPILS solver constants
 * -----------------------------------------------------------------
 * CPSPILS_MAXL   : default value for the maximum Krylov
 *                  dimension
 *
 * CPSPILS_MSBPRE : maximum number of steps between
 *                  preconditioner evaluations
 *
 * CPSPILS_DGMAX  : maximum change in gamma between
 *                  preconditioner evaluations
 *
 * CPSPILS_DELT   : default value for factor by which the
 *                  tolerance on the nonlinear iteration is
 *                  multiplied to get a tolerance on the linear
 *                  iteration
 * -----------------------------------------------------------------
 */
  
#define CPSPILS_MAXL   5
#define CPSPILS_MSBPRE 50
#define CPSPILS_DGMAX  RCONST(0.2)
#define CPSPILS_DELT   RCONST(0.05)

/*
 * -----------------------------------------------------------------
 * Types : CPSpilsMemRec, CPSpilsMem
 * -----------------------------------------------------------------
 * The type CPSpilsMem is pointer to a CPSpilsMemRec.
 * -----------------------------------------------------------------
 */

typedef struct {

  int s_type;           /* type of scaled preconditioned iterative LS   */

  int  s_pretype;       /* type of preconditioning                      */
  int  s_gstype;        /* type of Gram-Schmidt orthogonalization       */
  realtype s_sqrtN;     /* sqrt(N)                                      */
  realtype s_delt;      /* delt = user specified or DELT_DEFAULT        */
  realtype s_deltar;    /* deltar = delt * tq4                          */
  realtype s_delta;     /* delta = deltar * sqrtN                       */
  int  s_maxl;          /* maxl = maximum dimension of the Krylov space */

  long int s_nstlpre;   /* value of nst at the last pset call           */
  long int s_npe;       /* npe = total number of pset calls             */
  long int s_nli;       /* nli = total number of linear iterations      */
  long int s_nps;       /* nps = total number of psolve calls           */
  long int s_ncfl;      /* ncfl = total number of convergence failures  */
  long int s_njtimes;   /* njtimes = total number of calls to jtimes    */
  long int s_nfes;      /* no. of calls to f due to DQ Jacobian approx. */

  N_Vector s_ytemp;     /* temp vector passed to jtv and pslv           */
  N_Vector s_yptemp;    /* temp vector passed to jtv and pslv           */
  N_Vector s_x;         /* temp vector used by CPSpilsSolve             */
  N_Vector s_ycur;      /* CPODES current y vector in Newton Iteration  */
  N_Vector s_ypcur;     /* CPODES current y' vector in Newton Iteration */
  N_Vector s_fcur;      /* fcur = f(tn, ycur)                           */

  CPSpilsPrecSetupExplFn s_psetE; /* preconditioner setup (CP_EXPL case) */
  CPSpilsPrecSetupImplFn s_psetI; /* preconditioner setup (CP_IMPL case) */

  CPSpilsPrecSolveExplFn s_pslvE; /* preconditioner solve (CP_EXPL case) */
  CPSpilsPrecSolveImplFn s_pslvI; /* preconditioner solve (CP_IMPL case) */

  CPSpilsJacTimesVecExplFn s_jtvE; /* Jac. times vec. (CP_EXPL case)     */
  CPSpilsJacTimesVecImplFn s_jtvI; /* Jac. times vec. (CP_EXPL case)     */

  void *s_P_data;       /* data passed to psolve and pset                */
  void *s_j_data;       /* data passed to jtv                            */

  void* s_spils_mem;    /* memory used by the generic solver             */

  int s_last_flag;      /* last error flag returned by any function      */

} CPSpilsMemRec, *CPSpilsMem;

/*
 * -----------------------------------------------------------------
 * Prototypes of internal functions
 * -----------------------------------------------------------------
 */

/* Atimes and PSolve routines called by generic solver */
int cpSpilsAtimes(void *cp_mem, N_Vector v, N_Vector z);
int cpSpilsPSolve(void *cp_mem, N_Vector r, N_Vector z, int lr);

/* Difference quotient approximations for Jac times vector */
int cpSpilsDQjtvExpl(realtype t, N_Vector y, N_Vector fy, 
                     N_Vector v, N_Vector Jv, void *jac_data, 
                     N_Vector tmp); 
int cpSpilsDQjtvImpl(realtype t, realtype gm, 
                     N_Vector y, N_Vector yp, N_Vector r,
                     N_Vector v, N_Vector Jv, void *jac_data,
                     N_Vector tmp1, N_Vector tmp2);

/*
 * -----------------------------------------------------------------
 * Error Messages
 * -----------------------------------------------------------------
 */

#define MSGS_CPMEM_NULL  "Integrator memory is NULL."
#define MSGS_MEM_FAIL    "A memory request failed."
#define MSGS_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGS_BAD_LSTYPE  "Incompatible linear solver type."
#define MSGS_BAD_PRETYPE "Illegal value for pretype. Legal values are PREC_NONE, PREC_LEFT, PREC_RIGHT, and PREC_BOTH."
#define MSGS_PSOLVE_REQ  "pretype != PREC_NONE, but PSOLVE = NULL is illegal."
#define MSGS_LMEM_NULL   "Linear solver memory is NULL."
#define MSGS_BAD_GSTYPE  "Illegal value for gstype. Legal values are MODIFIED_GS and CLASSICAL_GS."
#define MSGS_BAD_DELT    "delt < 0 illegal."

#define MSGS_PSET_FAILED "The preconditioner setup routine failed in an unrecoverable manner."
#define MSGS_PSOLVE_FAILED "The preconditioner solve routine failed in an unrecoverable manner."
#define MSGS_JTIMES_FAILED "The Jacobian x vector routine failed in an unrecoverable manner."


#ifdef __cplusplus
}
#endif

#endif
