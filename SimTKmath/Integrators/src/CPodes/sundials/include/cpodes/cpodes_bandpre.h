/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2006/11/29 00:05:05 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the header file for the CPBANDPRE module, which
 * provides a banded difference quotient Jacobian-based
 * preconditioner and solver routines for use with CPSPGMR,
 * CPSPBCG, or CPSPTFQMR.
 *
 * Summary:
 * These routines provide a band matrix preconditioner based on
 * difference quotients of the ODE right-hand side function f
 * (for explicit-form ODEs) or on the residual function F (for
 * implicit-form ODEs).
 * The user supplies parameters
 *   mu = upper half-bandwidth (number of super-diagonals)
 *   ml = lower half-bandwidth (number of sub-diagonals)
 * The routines generate a band matrix of bandwidth ml + mu + 1
 * and use this to form a preconditioner for use with one of the
 * CSPILS iterative linear solvers. Although this matrix is intended
 * to approximate the Jacobian df/dy (respectively dF/dy + gamma*dF/dy,
 * it may be a very crude approximation. The true Jacobian need not
 * be banded, or its true bandwidth may be larger than ml + mu + 1,
 * as long as the banded approximation generated here is sufficiently
 * accurate to speed convergence as a preconditioner.
 *
 * Usage:
 *   The following is a summary of the usage of this module.
 *   Details of the calls to CPodeCreate, CPodeMalloc, CPSp*,
 *   and CPode are available in the User Guide.
 *   To use these routines, the sequence of calls in the user
 *   main program should be as follows:
 *
 *   #include <cpodes/cpodes_bandpre.h>
 *   #include <nvector_serial.h>
 *   ...
 *   void *bp_data;
 *   ...
 *   Set y0
 *   ...
 *   cpode_mem = CPodeCreate(...);
 *   ier = CPodeMalloc(...);
 *   ...
 *   bp_data = CPBandPrecAlloc(cpode_mem, N, mu, ml);
 *   ...
 *   flag = CPBPSptfqmr(cpode_mem, pretype, maxl, bp_data);
 *     -or-
 *   flag = CPBPSpgmr(cpode_mem, pretype, maxl, bp_data);
 *     -or-
 *   flag = CPBPSpbcg(cpode_mem, pretype, maxl, bp_data);
 *   ...
 *   flag = CPode(...);
 *   ...
 *   CPBandPrecFree(&bp_data);
 *   ...
 *   Free y0
 *   ...
 *   CPodeFree(cpode_mem);
 *
 * Notes:
 * (1) Include this file for the CPBandPrecData type definition.
 * (2) In the CPBandPrecAlloc call, the arguments N is the
 *     problem dimension.
 * (3) In the CPBPSp* call, the user is free to specify
 *     the input pretype and the optional input maxl. The last
 *     argument must be the pointer returned by CPBandPrecAlloc.
 * -----------------------------------------------------------------
 */

#ifndef _CPBANDPRE_H
#define _CPBANDPRE_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <sundials/sundials_nvector.h>

/* CPBANDPRE return values */

#define CPBANDPRE_SUCCESS        0
#define CPBANDPRE_PDATA_NULL   -11
#define CPBANDPRE_FUNC_UNRECVR -12

/*
 * -----------------------------------------------------------------
 * Function : CPBandPrecAlloc
 * -----------------------------------------------------------------
 * CPBandPrecAlloc allocates and initializes a CPBandPrecData
 * structure to be passed to CPSp* (and subsequently used
 * by CPBandPrecSetup and CPBandPrecSolve).
 *
 * The parameters of CPBandPrecAlloc are as follows:
 *
 * cpode_mem is the pointer to CPODES memory returned by CPodeCreate.
 *
 * N is the problem size.
 *
 * mu is the upper half bandwidth.
 *
 * ml is the lower half bandwidth.
 *
 * CPBandPrecAlloc returns the storage pointer of type
 * CPBandPrecData, or NULL if the request for storage cannot be
 * satisfied.
 *
 * NOTE: The band preconditioner assumes a serial implementation
 *       of the NVECTOR package. Therefore, CPBandPrecAlloc will
 *       first test for a compatible N_Vector internal
 *       representation by checking for required functions.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void *CPBandPrecAlloc(void *cpode_mem, int N, int mu, int ml);

/*
 * -----------------------------------------------------------------
 * Function : CPBPSptfqmr
 * -----------------------------------------------------------------
 * CPBPSptfqmr links the CPBANDPPRE preconditioner to the CPSPTFQMR
 * linear solver. It performs the following actions:
 *  1) Calls the CPSPTFQMR specification routine and attaches the
 *     CPSPTFQMR linear solver to the integrator memory;
 *  2) Sets the preconditioner data structure for CPSPTFQMR
 *  3) Sets the preconditioner setup routine for CPSPTFQMR
 *  4) Sets the preconditioner solve routine for CPSPTFQMR
 *
 * Its first 3 arguments are the same as for CPSptfqmr (see
 * cpsptfqmr.h). The last argument is the pointer to the CPBANDPPRE
 * memory block returned by CPBandPrecAlloc. Note that the user need
 * not call CPSptfqmr.
 *
 * Possible return values are:
 *    CPSPILS_SUCCESS      if successful
 *    CPSPILS_MEM_NULL     if the CPODES memory was NULL
 *    CPSPILS_LMEM_NULL    if the CPSPILS memory was NULL
 *    CPSPILS_MEM_FAIL     if there was a memory allocation failure
 *    CPSPILS_ILL_INPUT    if a required vector operation is missing
 *    CPBANDPRE_PDATA_NULL if the bp_data was NULL
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPBPSptfqmr(void *cpode_mem, int pretype, int maxl, void *p_data);

/*
 * -----------------------------------------------------------------
 * Function : CPBPSpbcg
 * -----------------------------------------------------------------
 * CPBPSpbcg links the CPBANDPPRE preconditioner to the CPSPBCG
 * linear solver. It performs the following actions:
 *  1) Calls the CPSPBCG specification routine and attaches the
 *     CPSPBCG linear solver to the integrator memory;
 *  2) Sets the preconditioner data structure for CPSPBCG
 *  3) Sets the preconditioner setup routine for CPSPBCG
 *  4) Sets the preconditioner solve routine for CPSPBCG
 *
 * Its first 3 arguments are the same as for CPSpbcg (see
 * cpspbcg.h). The last argument is the pointer to the CPBANDPPRE
 * memory block returned by CPBandPrecAlloc. Note that the user need
 * not call CPSpbcg.
 *
 * Possible return values are:
 *    CPSPILS_SUCCESS       if successful
 *    CPSPILS_MEM_NULL      if the CPODEs memory was NULL
 *    CPSPILS_LMEM_NULL     if the CPSPILS memory was NULL
 *    CPSPILS_MEM_FAIL      if there was a memory allocation failure
 *    CPSPILS_ILL_INPUT     if a required vector operation is missing
 *    CPBANDPRE_PDATA_NULL  if the bp_data was NULL
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPBPSpbcg(void *cpode_mem, int pretype, int maxl, void *p_data);

/*
 * -----------------------------------------------------------------
 * Function : CPBPSpgmr
 * -----------------------------------------------------------------
 * CPBPSpgmr links the CPBANDPPRE preconditioner to the CPSPGMR
 * linear solver. It performs the following actions:
 *  1) Calls the CPSPGMR specification routine and attaches the
 *     CPSPGMR linear solver to the integrator memory;
 *  2) Sets the preconditioner data structure for CPSPGMR
 *  3) Sets the preconditioner setup routine for CPSPGMR
 *  4) Sets the preconditioner solve routine for CPSPGMR
 *
 * Its first 3 arguments are the same as for CPSpgmr (see
 * cpspgmr.h). The last argument is the pointer to the CPBANDPPRE
 * memory block returned by CPBandPrecAlloc. Note that the user need
 * not call CPSpgmr.
 *
 * Possible return values are:
 *    CPSPILS_SUCCESS       if successful
 *    CPSPILS_MEM_NULL      if the CPODES memory was NULL
 *    CPSPILS_LMEM_NULL     if the CPSPILS memory was NULL
 *    CPSPILS_MEM_FAIL      if there was a memory allocation failure
 *    CPSPILS_ILL_INPUT     if a required vector operation is missing
 *    CPBANDPRE_PDATA_NULL  if the bp_data was NULL
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPBPSpgmr(void *cpode_mem, int pretype, int maxl, void *p_data);

/*
 * -----------------------------------------------------------------
 * Function : CPBandPrecFree
 * -----------------------------------------------------------------
 * CPBandPrecFree frees the memory allocated by CPBandPrecAlloc
 * in the argument bp_data.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void CPBandPrecFree(void **bp_data);

/*
 * -----------------------------------------------------------------
 * Optional output functions : CPBandPrecGet*
 * -----------------------------------------------------------------
 * CPBandPrecGetWorkSpace returns the real and integer work space used
 *                        by CPBANDPRE.
 * CPBandPrecGetNumFctEvals returns the number of calls made from
 *                          CPBANDPRE to the user's ODE function.
 *
 * The return value of CPBandPrecGet* is one of:
 *    CPBANDPRE_SUCCESS    if successful
 *    CPBANDPRE_PDATA_NULL if the bp_data memory was NULL
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPBandPrecGetWorkSpace(void *bp_data, long int *lenrwLS, long int *leniwLS);
SUNDIALS_EXPORT int CPBandPrecGetNumFctEvals(void *bp_data, long int *nfevalsBP);

/*
 * -----------------------------------------------------------------
 * The following function returns the name of the constant
 * associated with a CPBANDPRE return flag
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT char *CPBandPrecGetReturnFlagName(int flag);

#ifdef __cplusplus
}
#endif

#endif
