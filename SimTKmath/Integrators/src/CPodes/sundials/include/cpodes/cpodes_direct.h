/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2006/11/29 00:05:06 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Common header file for the direct linear solvers in CPODES.
 * -----------------------------------------------------------------
 */

#ifndef _CPDIRECT_H
#define _CPDIRECT_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <sundials/sundials_direct.h>
#include <sundials/sundials_nvector.h>

/*
 * =================================================================
 *              C P D I R E C T     C O N S T A N T S
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * CPDIRECT input constants
 * -----------------------------------------------------------------
 * fact_type: is the type of constraint Jacobian factorization used
 *   for the projection. For independent constraints (i.e. Jacobian
 *   with full row rank) any of the following three options can be
 *   used (although their computational cost varies, depending on 
 *   the number of states and constraints)
 *      CPDIRECT_LU:  use LU decomposition of G^T.
 *      CPDIRECT_QR:  use QR decomposition of G^T
 *      CPDIRECT_SC:  use Schur complement
 * If it is known (or suspected) that some constraints are redundant,
 * the following option should be used:
 *      CPDIRECT_QRP: use QR with column pivoting on G^T.  
 * -----------------------------------------------------------------
 */

/* fact_type */

#define CPDIRECT_LU   1
#define CPDIRECT_QR   2
#define CPDIRECT_SC   3
#define CPDIRECT_QRP  4

/* 
 * -----------------------------------------------------------------
 * CPDIRECT return values 
 * -----------------------------------------------------------------
 */

#define CPDIRECT_SUCCESS           0
#define CPDIRECT_MEM_NULL         -1
#define CPDIRECT_LMEM_NULL        -2
#define CPDIRECT_ILL_INPUT        -3
#define CPDIRECT_MEM_FAIL         -4

/* Additional last_flag values */

#define CPDIRECT_JACFUNC_UNRECVR  -5
#define CPDIRECT_JACFUNC_RECVR    -6

/*
 * =================================================================
 *              F U N C T I O N   T Y P E S
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Types: CPDlsDenseJacExplFn and CPDlsDenseJacImplFn
 * -----------------------------------------------------------------
 *
 * If the ODE is given in explicit form, a dense Jacobian 
 * approximation function Jac must be of type CPDlsDenseJacFn. 
 * Its parameters are:
 *
 * N   is the problem size.
 *
 * Jac is the dense matrix (of type DlsMat) that will be loaded
 *     by a CPDlsDenseJacFn with an approximation to the Jacobian 
 *     matrix J = (df_i/dy_j) at the point (t,y). 
 *
 * t   is the current value of the independent variable.
 *
 * y   is the current value of the dependent variable vector,
 *     namely the predicted value of y(t).
 *
 * fy  is the vector f(t,y).
 *
 * jac_data is a pointer to user data - the same as the jac_data
 *     parameter passed to CPDlsSetJacFn.
 *
 * tmp1, tmp2, and tmp3 are pointers to memory allocated for
 * vectors of length N which can be used by a CPDlsDenseJacFn
 * as temporary storage or work space.
 *
 * A CPDlsDenseJacFn should return 0 if successful, a positive 
 * value if a recoverable error occurred, and a negative value if 
 * an unrecoverable error occurred.
 *
 * -----------------------------------------------------------------
 *
 * If the ODE is given in implicit form, a dense Jacobian 
 * approximation function djac must be of type CPDlsDenseJacImplFn. 
 * Its parameters are:                     
 *                                                                
 * N   is the problem size, and length of all vector arguments.   
 *                                                                
 * t   is the current value of the independent variable t.        
 *                                                                
 * y   is the current value of the dependent variable vector,     
 *     namely the predicted value of y(t).                     
 *                                                                
 * yp  is the current value of the derivative vector y',          
 *     namely the predicted value of y'(t).                    
 *                                                                
 * f   is the residual vector F(tt,yy,yp).                     
 *                                                                
 * gm  is the scalar in the system Jacobian, proportional to 
 *     the step size h.
 *                                                                
 * jac_data is a pointer to user Jacobian data - the same as the    
 *     jdata parameter passed to CPDenseSetJacFn.                     
 *                                                                
 * Jac is the dense matrix (of type DenseMat) to be loaded by  
 *     an CPDenseJacImplFn routine with an approximation to the   
 *     system Jacobian matrix                                  
 *            J = dF/dy' + gamma*dF/dy                            
 *     at the given point (t,y,y'), where the ODE system is    
 *     given by F(t,y,y') = 0.  Jac is preset to zero, so only 
 *     the nonzero elements need to be loaded.  See note below.
 *                                                                
 * tmp1, tmp2, tmp3 are pointers to memory allocated for          
 *     N_Vectors which can be used by an CPDenseJacImplFn routine 
 *     as temporary storage or work space.                     
 *                                                                
 * A CPDenseJacImplFn should return                                
 *     0 if successful,                                           
 *     a positive int if a recoverable error occurred, or         
 *     a negative int if a nonrecoverable error occurred.         
 * In the case of a recoverable error return, the integrator will 
 * attempt to recover by reducing the stepsize (which changes cj).
 *
 * -----------------------------------------------------------------
 *
 * NOTE: The following are two efficient ways to load a dense Jac:         
 * (1) (with macros - no explicit data structure references)      
 *     for (j=0; j < Neq; j++) {                                  
 *       col_j = DENSE_COL(Jac,j);                                 
 *       for (i=0; i < Neq; i++) {                                
 *         generate J_ij = the (i,j)th Jacobian element           
 *         col_j[i] = J_ij;                                       
 *       }                                                        
 *     }                                                          
 * (2) (without macros - explicit data structure references)      
 *     for (j=0; j < Neq; j++) {                                  
 *       col_j = (Jac->data)[j];                                   
 *       for (i=0; i < Neq; i++) {                                
 *         generate J_ij = the (i,j)th Jacobian element           
 *         col_j[i] = J_ij;                                       
 *       }                                                        
 *     }                                                          
 * A third way, using the DENSE_ELEM(A,i,j) macro, is much less   
 * efficient in general.  It is only appropriate for use in small 
 * problems in which efficiency of access is NOT a major concern. 
 *                                                                
 * NOTE: If the user's Jacobian routine needs other quantities,   
 *     they are accessible as follows: hcur (the current stepsize)
 *     and ewt (the error weight vector) are accessible through   
 *     CPodeGetCurrentStep and CPodeGetErrWeights, respectively 
 *     (see cvode.h). The unit roundoff is available as 
 *     UNIT_ROUNDOFF defined in sundials_types.h.
 *
 * -----------------------------------------------------------------
 */
  
  
typedef int (*CPDlsDenseJacExplFn)(int N, realtype t,
                   N_Vector y, N_Vector fy, 
                   DlsMat Jac, void *jac_data,
                   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
typedef int (*CPDlsDenseJacImplFn)(int N, realtype t, realtype gm,
                   N_Vector y, N_Vector yp, N_Vector r, 
                   DlsMat Jac, void *jac_data,
                   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/*
 * -----------------------------------------------------------------
 * Types: CPDlsBandJacExplFn and CPDlsBandJacImplFn
 * -----------------------------------------------------------------
 *
 * If the ODE is given in explicit form, a band Jacobian 
 * approximation function Jac must be of type CPDlsBandJacExplFn.
 * Its parameters are:
 *
 * N is the length of all vector arguments.
 *
 * mupper is the upper half-bandwidth of the approximate banded
 * Jacobian. This parameter is the same as the mupper parameter
 * passed by the user to the linear solver initialization function.
 *
 * mlower is the lower half-bandwidth of the approximate banded
 * Jacobian. This parameter is the same as the mlower parameter
 * passed by the user to the linear solver initialization function.
 *
 * t is the current value of the independent variable.
 *
 * y is the current value of the dependent variable vector,
 *      namely the predicted value of y(t).
 *
 * fy is the vector f(t,y).
 *
 * Jac is the band matrix (of type DlsMat) that will be loaded
 * by a CPDlsBandJacFn with an approximation to the Jacobian matrix
 * Jac = (df_i/dy_j) at the point (t,y).
 * jac_data is a pointer to user data - the same as the jac_data
 *          parameter passed to CPDlsSetJacFn.
 *
 * A CPDlsBandJacExplFn should return 0 if successful, a positive
 * value if a recoverable error occurred, and a negative value if
 * an unrecoverable error occurred.
 *
 * -----------------------------------------------------------------
 *
 * If the ODE is given in explicit form, a band Jacobian 
 * approximation function Jac must be of type CPDlsBandJacImplFn.
 * Its parameters are:
 *
 * N is the problem size, and length of all vector arguments.   
 *                                                                
 * mupper is the upper bandwidth of the banded Jacobian matrix.   
 *                                                                
 * mlower is the lower bandwidth of the banded Jacobian matrix.   
 *                                                                
 * t is the current value of the independent variable t.        
 *                                                                
 * y is the current value of the dependent variable vector,     
 *    namely the predicted value of y(t).                     
 *                                                                
 * yp is the current value of the derivative vector y',          
 *    namely the predicted value of y'(t).                    
 *                                                                
 * r is the residual vector F(tt,yy,yp).                     
 *                                                                
 * gm  is the scalar in the system Jacobian, proportional to 
 *     the step size h.
 *                                                                
 * jac_data  is a pointer to user Jacobian data - the same as the    
 *    jdata parameter passed to CPBandSetJacFn.                      
 *                                                                
 * J is the band matrix (of type BandMat) to be loaded by    
 *     with an approximation to the system Jacobian matrix
 *            J = dF/dy' + gamma*dF/dy 
 *     at the given point (t,y,y'), where the ODE system is    
 *     given by F(t,y,y') = 0.  Jac is preset to zero, so only 
 *     the nonzero elements need to be loaded.  See note below.
 *                                                                
 * tmp1, tmp2, tmp3 are pointers to memory allocated for          
 *     N_Vectors which can be used by an CPBandJacImplFn routine  
 *     as temporary storage or work space.                     
 *
 * A CPDlsBandJacImplFn should return 0 if successful, a positive
 * value if a recoverable error occurred, and a negative value if
 * an unrecoverable error occurred.
 *
 * -----------------------------------------------------------------
 *
 * NOTE: Three efficient ways to load J are:
 *
 * (1) (with macros - no explicit data structure references)
 *    for (j=0; j < n; j++) {
 *       col_j = BAND_COL(Jac,j);
 *       for (i=j-mupper; i <= j+mlower; i++) {
 *         generate J_ij = the (i,j)th Jacobian element
 *         BAND_COL_ELEM(col_j,i,j) = J_ij;
 *       }
 *     }
 *
 * (2) (with BAND_COL macro, but without BAND_COL_ELEM macro)
 *    for (j=0; j < n; j++) {
 *       col_j = BAND_COL(Jac,j);
 *       for (k=-mupper; k <= mlower; k++) {
 *         generate J_ij = the (i,j)th Jacobian element, i=j+k
 *         col_j[k] = J_ij;
 *       }
 *     }
 *
 * (3) (without macros - explicit data structure references)
 *     offset = Jac->smu;
 *     for (j=0; j < n; j++) {
 *       col_j = ((Jac->data)[j])+offset;
 *       for (k=-mupper; k <= mlower; k++) {
 *         generate J_ij = the (i,j)th Jacobian element, i=j+k
 *         col_j[k] = J_ij;
 *       }
 *     }
 * Caution: Jac->smu is generally NOT the same as mupper.
 *
 * The BAND_ELEM(A,i,j) macro is appropriate for use in small
 * problems in which efficiency of access is NOT a major concern.
 *
 * NOTE: If the user's Jacobian routine needs other quantities,
 *     they are accessible as follows: hcur (the current stepsize)
 *     and ewt (the error weight vector) are accessible through
 *     CPodeGetCurrentStep and CPodeGetErrWeights, respectively
 *     (see cvode.h). The unit roundoff is available as
 *     UNIT_ROUNDOFF defined in sundials_types.h
 *
 * tmp1, tmp2, and tmp3 are pointers to memory allocated for
 * vectors of length N which can be used by a CPDlsBandJacFn
 * as temporary storage or work space.
 *
 * -----------------------------------------------------------------
 */

typedef int (*CPDlsBandJacExplFn)(int N, int mupper, int mlower,
                  realtype t, N_Vector y, N_Vector fy, 
                  DlsMat Jac, void *jac_data,
                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
typedef int (*CPDlsBandJacImplFn)(int N, int mupper, int mlower,
                  realtype t, realtype gm, 
                  N_Vector y, N_Vector yp, N_Vector r,
                  DlsMat Jac, void *jac_data,
                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/*
 * -----------------------------------------------------------------
 * Type: CPDlsDenseProjJacFn
 * -----------------------------------------------------------------
 *
 *
 * -----------------------------------------------------------------
 */

typedef int (*CPDlsDenseProjJacFn)(int Nc, int Ny, 
                   realtype t, N_Vector y, N_Vector cy,
                   DlsMat Jac, void *jac_data,
                   N_Vector tmp1, N_Vector tmp2); 



/*
 * =================================================================
 *            E X P O R T E D    F U N C T I O N S 
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Optional inputs to a CPDIRECT linear solver 
 * for implicit integration
 * -----------------------------------------------------------------
 *
 * CPDlsSetJacFn specifies the Jacobian approximation routine to be
 * used. When using dense Jacobians, a user-supplied jac routine must 
 * be of type CPDlsDenseJacFn. When using banded Jacobians, a 
 * user-supplied jac routine must be of type CPDlsBandJacFn.
 * By default, a difference quotient approximation, supplied with this 
 * solver is used.
 * CPDlsSetJacFn also specifies a pointer to user data which is 
 * passed to the user's jac routine every time it is called.
 *
 * The return value of CPDlsSetJacFn is one of:
 *    CPDIRECT_SUCCESS   if successful
 *    CPDIRECT_MEM_NULL  if the CPODES memory was NULL
 *    CPDIRECT_LMEM_NULL if the linear solver memory was NULL
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPDlsSetJacFn(void *cvode_mem, void *jac, void *jac_data);

/*
 * -----------------------------------------------------------------
 * Optional outputs from a CPDIRECT linear solver
 * for implicit integration
 * -----------------------------------------------------------------
 *
 * CPDlsGetWorkSpace   returns the real and integer workspace used
 *                     by the direct linear solver.
 * CPDlsGetNumJacEvals returns the number of calls made to the
 *                     Jacobian evaluation routine jac.
 * CPDlsGetNumFctEvals returns the number of calls to the user
 *                     f routine due to finite difference Jacobian
 *                     evaluation.
 * CPDlsGetLastFlag    returns the last error flag set by any of
 *                     the CPDIRECT interface functions.
 *
 * The return value of CPDlsGet* is one of:
 *    CPDIRECT_SUCCESS   if successful
 *    CPDIRECT_MEM_NULL  if the CPODES memory was NULL
 *    CPDIRECT_LMEM_NULL if the linear solver memory was NULL
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPDlsGetWorkSpace(void *cvode_mem, long int *lenrwLS, long int *leniwLS);
SUNDIALS_EXPORT int CPDlsGetNumJacEvals(void *cvode_mem, long int *njevals);
SUNDIALS_EXPORT int CPDlsGetNumFctEvals(void *cvode_mem, long int *nfevalsLS);
SUNDIALS_EXPORT int CPDlsGetLastFlag(void *cvode_mem, int *flag);

/*
 * -----------------------------------------------------------------
 * The following function returns the name of the constant 
 * associated with a CPDIRECT return flag
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT char *CPDlsGetReturnFlagName(int flag);

/*
 * -----------------------------------------------------------------
 * Optional I/O functions for a CPDIRECT linear solver
 * for projection
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPDlsProjSetJacFn(void *cpode_mem, void *jacP, void *jacP_data);

SUNDIALS_EXPORT int CPDlsProjGetNumJacEvals(void *cpode_mem, long int *njPevals);
SUNDIALS_EXPORT int CPDlsProjGetNumFctEvals(void *cpode_mem, long int *ncevalsLS);

#ifdef __cplusplus
}
#endif

#endif
