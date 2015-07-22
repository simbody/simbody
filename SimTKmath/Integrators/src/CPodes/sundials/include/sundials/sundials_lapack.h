/*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2006/11/29 00:05:07 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the header file for a generic package of direct matrix
 * operations for use with BLAS/LAPACK.
 * -----------------------------------------------------------------
 */

#ifndef _SUNDIALS_LAPACK_H
#define _SUNDIALS_LAPACK_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <sundials/sundials_direct.h>

/*
 * ==================================================================
 * Exported functions
 * ==================================================================
 */

/*
 * -----------------------------------------------------------------
 * Functions: LapackDenseZero and LapackBandZero
 * -----------------------------------------------------------------
 * These functions set all the elements of the M-by-N matrix A to 0.0.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void LapackDenseZero(DlsMat A);
SUNDIALS_EXPORT void LapackBandZero(DlsMat A);

/*
 * -----------------------------------------------------------------
 * Functions: LapackDenseAddI and LapackBandAddI
 * -----------------------------------------------------------------
 * These functions add 1.0 to the main diagonal (A_ii, i=1,2,...,N-1)
 * of the M-by-N matrix A (M>= N) and stores the result back in A.
 * They are typically used with square matrices.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void LapackDenseAddI(DlsMat A);
SUNDIALS_EXPORT void LapackBandAddI(DlsMat A);

/*
 * ==================================================================
 * Blas and Lapack functions
 * ==================================================================
 */

#if defined(F77_FUNC)

#define dcopy_f77       F77_FUNC(dcopy, DCOPY)
#define dscal_f77       F77_FUNC(dscal, DSCAL)
#define dgemv_f77       F77_FUNC(dgemv, DGEMV)
#define dtrsv_f77       F77_FUNC(dtrsv, DTRSV)
#define dsyrk_f77       F77_FUNC(dsyrk, DSKYR)

#define dgbtrf_f77      F77_FUNC(dgbtrf, DGBTRF)
#define dgbtrs_f77      F77_FUNC(dgbtrs, DGBTRS)
#define dgetrf_f77      F77_FUNC(dgetrf, DGETRF)
#define dgetrs_f77      F77_FUNC(dgetrs, DGETRS)
#define dgeqp3_f77      F77_FUNC(dgeqp3, DGEQP3)
#define dgeqrf_f77      F77_FUNC(dgeqrf, DGEQRF)
#define dormqr_f77      F77_FUNC(dormqr, DORMQR)
#define dpotrf_f77      F77_FUNC(dpotrf, DPOTRF)
#define dpotrs_f77      F77_FUNC(dpotrs, DPOTRS)

#else

#define dcopy_f77       dcopy_
#define dscal_f77       dscal_
#define dgemv_f77       dgemv_
#define dtrsv_f77       dtrsv_
#define dsyrk_f77       dsyrk_

#define dgbtrf_f77      dgbtrf_
#define dgbtrs_f77      dgbtrs_
#define dgeqp3_f77      dgeqp3_
#define dgeqrf_f77      dgeqrf_
#define dgetrf_f77      dgetrf_
#define dgetrs_f77      dgetrs_
#define dormqr_f77      dormqr_
#define dpotrf_f77      dpotrf_
#define dpotrs_f77      dpotrs_

#endif

/* Level-1 BLAS */

extern void dcopy_f77(int *n, const double *x, const int *inc_x, double *y, const int *inc_y);
extern void dscal_f77(int *n, const double *alpha, double *x, const int *inc_x);

/* Level-2 BLAS */

extern void dgemv_f77(const char *trans, int *m, int *n, const double *alpha, const double *a,
              int *lda, const double *x, int *inc_x, const double *beta, double *y, int *inc_y,
              int len_trans);

extern void dtrsv_f77(const char *uplo, const char *trans, const char *diag, const int *n,
              const double *a, const int *lda, double *x, const int *inc_x,
              int len_uplo, int len_trans, int len_diag);

/* Level-3 BLAS */

extern void dsyrk_f77(const char *uplo, const char *trans, const int *n, const int *k,
              const double *alpha, const double *a, const int *lda, const double *beta,
              const double *c, const int *ldc, int len_uplo, int len_trans);

/* LAPACK */

extern void dgbtrf_f77(const int *m, const int *n, const int *kl, const int *ku,
               double *ab, int *ldab, int *ipiv, int *info);

extern void dgbtrs_f77(const char *trans, const int *n, const int *kl, const int *ku, const int *nrhs,
               double *ab, const int *ldab, int *ipiv, double *b, const int *ldb,
               int *info, int len_trans);


extern void dgeqp3_f77(const int *m, const int *n, double *a, const int *lda, int *jpvt, double *tau,
               double *work, const int *lwork, int *info);

extern void dgeqrf_f77(const int *m, const int *n, double *a, const int *lda, double *tau, double *work,
               const int *lwork, int *info);

extern void dgetrf_f77(const int *m, const int *n, double *a, int *lda, int *ipiv, int *info);

extern void dgetrs_f77(const char *trans, const int *n, const int *nrhs, double *a, const int *lda,
               int *ipiv, double *b, const int *ldb, int *info, int len_trans);


extern void dormqr_f77(const char *side, const char *trans, const int *m, const int *n, const int *k,
               double *a, const int *lda, double *tau, double *c, const int *ldc,
               double *work, const int *lwork, int *info, int len_side, int len_trans);

extern void dpotrf_f77(const char *uplo, const int *n, double *a, int *lda, int *info, int len_uplo);

extern void dpotrs_f77(const char *uplo, const int *n, const int *nrhs, double *a, const int *lda,
               double *b, const int *ldb, int * info, int len_uplo);


#ifdef __cplusplus
}
#endif

#endif
