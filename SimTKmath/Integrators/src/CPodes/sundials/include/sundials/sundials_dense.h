/*
 * -----------------------------------------------------------------
 * $Revision: 1.6 $
 * $Date: 2006/11/29 00:05:07 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the header file for a generic package of DENSE matrix
 * operations, based on the DlsMat type defined in sundials_direct.h.
 *
 * There are two sets of dense solver routines listed in
 * this file: one set uses type DlsMat defined below and the
 * other set uses the type realtype ** for dense matrix arguments.
 * Routines that work with the type DlsMat begin with "Dense".
 * Routines that work with realtype** begin with "dense". 
 * -----------------------------------------------------------------
 */

#ifndef _SUNDIALS_DENSE_H
#define _SUNDIALS_DENSE_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <sundials/sundials_direct.h>

/*
 * -----------------------------------------------------------------
 * Functions: DenseGETRF and DenseGETRS
 * -----------------------------------------------------------------
 * DenseGETRF performs the LU factorization of the M by N dense
 * matrix A. This is done using standard Gaussian elimination
 * with partial (row) pivoting. Note that this applies only
 * to matrices with M >= N and full column rank.
 *
 * A successful LU factorization leaves the matrix A and the
 * pivot array p with the following information:
 *
 * (1) p[k] contains the row number of the pivot element chosen
 *     at the beginning of elimination step k, k=0, 1, ..., N-1.
 *
 * (2) If the unique LU factorization of A is given by PA = LU,
 *     where P is a permutation matrix, L is a lower trapezoidal
 *     matrix with all 1's on the diagonal, and U is an upper
 *     triangular matrix, then the upper triangular part of A
 *     (including its diagonal) contains U and the strictly lower
 *     trapezoidal part of A contains the multipliers, I-L.
 *
 * For square matrices (M=N), L is unit lower triangular.
 *
 * DenseGETRF returns 0 if successful. Otherwise it encountered
 * a zero diagonal element during the factorization. In this case
 * it returns the column index (numbered from one) at which
 * it encountered the zero.
 *
 * DenseGETRS solves the N-dimensional system A x = b using
 * the LU factorization in A and the pivot information in p
 * computed in DenseGETRF. The solution x is returned in b. This
 * routine cannot fail if the corresponding call to DenseGETRF
 * did not fail.
 * DenseGETRS does NOT check for a square matrix!
 *
 * -----------------------------------------------------------------
 * DenseGETRF and DenseGETRS are simply wrappers around denseGETRF
 * and denseGETRS, respectively, which perform all the work by
 * directly accessing the data in the DlsMat A (i.e. the field cols)
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int DenseGETRF(DlsMat A, int *p);
SUNDIALS_EXPORT void DenseGETRS(DlsMat A, int *p, realtype *b);

SUNDIALS_EXPORT int denseGETRF(realtype **a, int m, int n, int *p);
SUNDIALS_EXPORT void denseGETRS(realtype **a, int n, int *p, realtype *b);

/*
 * -----------------------------------------------------------------
 * Functions : DensePOTRF and DensePOTRS
 * -----------------------------------------------------------------
 * DensePOTRF computes the Cholesky factorization of a real symmetric
 * positive definite matrix A.
 * -----------------------------------------------------------------
 * DensePOTRS solves a system of linear equations A*X = B with a 
 * symmetric positive definite matrix A using the Cholesky factorization
 * A = L*L**T computed by DensePOTRF.
 *
 * -----------------------------------------------------------------
 * DensePOTRF and DensePOTRS are simply wrappers around densePOTRF
 * and densePOTRS, respectively, which perform all the work by
 * directly accessing the data in the DlsMat A (i.e. the field cols)
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int DensePOTRF(DlsMat A);
SUNDIALS_EXPORT void DensePOTRS(DlsMat A, realtype *b);

SUNDIALS_EXPORT int densePOTRF(realtype **a, int m);
SUNDIALS_EXPORT void densePOTRS(realtype **a, int m, realtype *b);

/*
 * -----------------------------------------------------------------
 * Functions : DenseGEQRF and DenseORMQR
 * -----------------------------------------------------------------
 * DenseGEQRF computes a QR factorization of a real M-by-N matrix A:
 * A = Q * R (with M>= N).
 * 
 * DenseGEQRF requires a temporary work vector wrk of length M.
 * -----------------------------------------------------------------
 * DenseORMQR computes the product w = Q * v where Q is a real 
 * orthogonal matrix defined as the product of k elementary reflectors
 *
 *        Q = H(1) H(2) . . . H(k)
 *
 * as returned by DenseGEQRF. Q is an M-by-N matrix, v is a vector
 * of length N and w is a vector of length M (with M>=N).
 *
 * DenseORMQR requires a temporary work vector wrk of length M.
 *
 * -----------------------------------------------------------------
 * DenseGEQRF and DenseORMQR are simply wrappers around denseGEQRF
 * and denseORMQR, respectively, which perform all the work by
 * directly accessing the data in the DlsMat A (i.e. the field cols)
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int DenseGEQRF(DlsMat A, realtype *beta, realtype *wrk);
SUNDIALS_EXPORT int DenseORMQR(DlsMat A, realtype *beta, realtype *vn, realtype *vm, 
                   realtype *wrk);

SUNDIALS_EXPORT int denseGEQRF(realtype **a, int m, int n, realtype *beta, realtype *v);
SUNDIALS_EXPORT int denseORMQR(realtype **a, int m, int n, realtype *beta,
                   realtype *v, realtype *w, realtype *wrk);

/*
 * -----------------------------------------------------------------
 * Function : DenseZero
 * -----------------------------------------------------------------
 * DenseZero sets all the elements of the M-by-N matrix A to 0.0.
 *
 * DenseZero is a wrapper around denseZero which accesses the data of
 * the DlsMat A (i.e. the field cols)
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void DenseZero(DlsMat A);
SUNDIALS_EXPORT void denseZero(realtype **a, int m, int n);

/*
 * -----------------------------------------------------------------
 * Function : DenseCopy
 * -----------------------------------------------------------------
 * DenseCopy copies the contents of the M-by-N matrix A into the
 * M-by-N matrix B.
 * 
 * DenseCopy is a wrapper around denseCopy which accesses the data
 * in the DlsMat A and B (i.e. the fields cols)
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void DenseCopy(DlsMat A, DlsMat B);
SUNDIALS_EXPORT void denseCopy(realtype **a, realtype **b, int m, int n);

/*
 * -----------------------------------------------------------------
 * Function: DenseScale
 * -----------------------------------------------------------------
 * DenseScale scales the elements of the M-by-N matrix A by the
 * constant c and stores the result back in A.
 *
 * DenseScale is a wrapper around denseScale which performs the actual
 * scaling by accessing the data in the DlsMat A (i.e. the field
 * cols).
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void DenseScale(realtype c, DlsMat A);
SUNDIALS_EXPORT void denseScale(realtype c, realtype **a, int m, int n);

/*
 * -----------------------------------------------------------------
 * Function : DenseAddI
 * -----------------------------------------------------------------
 * DenseAddI adds 1.0 to the main diagonal (A_ii, i=1,2,...,N-1) of
 * the M-by-N matrix A (M>= N) and stores the result back in A.
 * DenseAddI is typically used with square matrices.
 * DenseAddI does not check for M >= N and therefore a segmentation
 * fault will occur if M < N!
 *
 * DenseAddI is a wrapper around denseAddI which performs the actual
 * work by accessing the data in the DlsMat A (i.e. the field cols)
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void DenseAddI(DlsMat A);
SUNDIALS_EXPORT void denseAddI(realtype **a, int n);

#ifdef __cplusplus
}
#endif

#endif
