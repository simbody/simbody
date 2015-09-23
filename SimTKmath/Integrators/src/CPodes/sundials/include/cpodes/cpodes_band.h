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
 * This is the header file for the CPODES band linear solver, CPBAND.
 * -----------------------------------------------------------------
 */

#ifndef _CPBAND_H
#define _CPBAND_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <cpodes/cpodes_direct.h>
#include <sundials/sundials_band.h>

/*
 * -----------------------------------------------------------------
 * Function : CPBand
 * -----------------------------------------------------------------
 * A call to the CPBand function links the main CPODES integrator
 * with the CPBAND linear solver.
 *
 * cpode_mem is the pointer to the integrator memory returned by
 *           CPodeCreate.
 *
 * N is the size of the ODE system.
 *
 * mupper is the upper bandwidth of the band Jacobian
 *        approximation.
 *
 * mlower is the lower bandwidth of the band Jacobian
 *        approximation.
 *
 * The return value of CPBand is one of:
 *    CPDIRECT_SUCCESS   if successful
 *    CPDIRECT_MEM_NULL  if the CPODES memory was NULL
 *    CPDIRECT_MEM_FAIL  if there was a memory allocation failure
 *    CPDIRECT_ILL_INPUT if a required vector operation is missing
 *                       or if a bandwidth has an illegal value.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPBand(void *cpode_mem, int N, int mupper, int mlower);

#ifdef __cplusplus
}
#endif

#endif
