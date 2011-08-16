// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpBlas.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpoptConfig.h"
#include "IpBlas.hpp"

// Prototypes for the BLAS routines
extern "C"
{
  /** BLAS Fortran function DDOT */
  double F77_FUNC(ddot,DDOT)(ipfint *n, const double *x, ipfint *incX,
                             const double *y, ipfint *incY);
  /** BLAS Fortran function DNRM2 */
  double F77_FUNC(dnrm2,DNRM2)(ipfint *n, const double *x, ipfint *incX);
  /** BLAS Fortran function DASUM */
  double F77_FUNC(dasum,DASUM)(ipfint *n, const double *x, ipfint *incX);
  /** BLAS Fortran function IDAMAX */
  ipfint F77_FUNC(idamax,IDAMAX)(ipfint *n, const double *x, ipfint *incX);
  /** BLAS Fortran subroutine DCOPY */
  void F77_FUNC(dcopy,DCOPY)(ipfint *n, const double *x, ipfint *incX,
                             double *y, ipfint *incY);
  /** BLAS Fortran subroutine DAXPY */
  void F77_FUNC(daxpy,DAXPY)(ipfint *n, const double *alpha, const double *x,
                             ipfint *incX, double *y, ipfint *incY);
  /** BLAS Fortran subroutine DSCAL */
  void F77_FUNC(dscal,DSCAL)(ipfint *n, const double *alpha, const double *x,
                             ipfint *incX);

  void F77_FUNC(dgemv,DGEMV)(char* trans, ipfint *m, ipfint *n,
                             const double *alpha, const double *a, ipfint *lda,
                             const double *x, ipfint *incX, const double *beta,
                             double *y, ipfint *incY, int trans_len);
  void F77_FUNC(dsymv,DSYMV)(char* uplo, ipfint *n,
                             const double *alpha, const double *a, ipfint *lda,
                             const double *x, ipfint *incX, const double *beta,
                             double *y, ipfint *incY, int uplo_len);
  void F77_FUNC(dgemm,DGEMM)(char* transa, char* transb,
                             ipfint *m, ipfint *n, ipfint *k,
                             const double *alpha, const double *a, ipfint *lda,
                             const double *b, ipfint *ldb, const double *beta,
                             double *c, ipfint *ldc,
                             int transa_len, int transb_len);
  void F77_FUNC(dsyrk,DSYRK)(char* uplo, char* trans, ipfint *n, ipfint *k,
                             const double *alpha, const double *a, ipfint *lda,
                             const double *beta, double *c, ipfint *ldc,
                             int uplo_len, int trans_len);
  void F77_FUNC(dtrsm,DTRSM)(char* side, char* uplo, char* transa, char* diag,
                             ipfint *m, ipfint *n,
                             const double *alpha, const double *a, ipfint *lda,
                             const double *b, ipfint *ldb,
                             int side_len, int uplo_len,
                             int transa_len, int diag_len);
}

namespace Ipopt
{
#ifndef HAVE_CBLAS
  /* Interface to FORTRAN routine DDOT. */
  Number IpBlasDdot(Index size, const Number *x, Index incX, const Number *y,
                    Index incY)
  {
    ipfint n=size, INCX=incX, INCY=incY;

    return F77_FUNC(ddot,DDOT)(&n, x, &INCX, y, &INCY);
  }

  /* Interface to FORTRAN routine DNRM2. */
  Number IpBlasDnrm2(Index size, const Number *x, Index incX)
  {
    ipfint n=size, INCX=incX;

    return F77_FUNC(dnrm2,DNRM2)(&n, x, &INCX);
  }

  /* Interface to FORTRAN routine DASUM. */
  Number IpBlasDasum(Index size, const Number *x, Index incX)
  {
    ipfint n=size, INCX=incX;

    return F77_FUNC(dasum,DASUM)(&n, x, &INCX);
  }

  /* Interface to FORTRAN routine DASUM. */
  Index IpBlasIdamax(Index size, const Number *x, Index incX)
  {
    ipfint n=size, INCX=incX;

    return (Index) F77_FUNC(idamax,IDAMAX)(&n, x, &INCX);
  }

  /* Interface to FORTRAN routine DCOPY. */
  void IpBlasDcopy(Index size, const Number *x, Index incX, Number *y, Index incY)
  {
    ipfint N=size, INCX=incX, INCY=incY;

    F77_FUNC(dcopy,DCOPY)(&N, x, &INCX, y, &INCY);
  }

  /* Interface to FORTRAN routine DAXPY. */
  void IpBlasDaxpy(Index size, Number alpha, const Number *x, Index incX, Number *y,
                   Index incY)
  {
    ipfint N=size, INCX=incX, INCY=incY;

    F77_FUNC(daxpy,DAXPY)(&N, &alpha, x, &INCX, y, &INCY);
  }

  /* Interface to FORTRAN routine DSCAL. */
  void IpBlasDscal(Index size, Number alpha, Number *x, Index incX)
  {
    ipfint N=size, INCX=incX;

    F77_FUNC(dscal,DSCAL)(&N, &alpha, x, &INCX);
  }

  void IpBlasDgemv(bool trans, Index nRows, Index nCols, Number alpha,
                   const Number* A, Index ldA, const Number* x,
                   Index incX, Number beta, Number* y, Index incY)
  {
    ipfint M=nCols, N=nRows, LDA=ldA, INCX=incX, INCY=incY;

    char TRANS;
    if (trans) {
      TRANS = 'T';
    }
    else {
      TRANS = 'N';
    }

    F77_FUNC(dgemv,DGEMV)(&TRANS, &M, &N, &alpha, A, &LDA, x,
                          &INCX, &beta, y, &INCY, 1);
  }

  void IpBlasDsymv(Index n, Number alpha, const Number* A, Index ldA,
                   const Number* x, Index incX, Number beta, Number* y,
                   Index incY)
  {
    ipfint N=n, LDA=ldA, INCX=incX, INCY=incY;

    char UPLO='L';

    F77_FUNC(dsymv,DSYMV)(&UPLO, &N, &alpha, A, &LDA, x,
                          &INCX, &beta, y, &INCY, 1);
  }

  void IpBlasDgemm(bool transa, bool transb, Index m, Index n, Index k,
                   Number alpha, const Number* A, Index ldA, const Number* B,
                   Index ldB, Number beta, Number* C, Index ldC)
  {
    ipfint M=m, N=n, K=k, LDA=ldA, LDB=ldB, LDC=ldC;

    char TRANSA;
    if (transa) {
      TRANSA = 'T';
    }
    else {
      TRANSA = 'N';
    }
    char TRANSB;
    if (transb) {
      TRANSB = 'T';
    }
    else {
      TRANSB = 'N';
    }

    F77_FUNC(dgemm,DGEMM)(&TRANSA, &TRANSB, &M, &N, &K, &alpha, A, &LDA,
                          B, &LDB, &beta, C, &LDC, 1, 1);
  }

  void IpBlasDsyrk(bool trans, Index ndim, Index nrank,
                   Number alpha, const Number* A, Index ldA,
                   Number beta, Number* C, Index ldC)
  {
    ipfint N=ndim, K=nrank, LDA=ldA, LDC=ldC;

    char UPLO='L';
    char TRANS;
    if (trans) {
      TRANS = 'T';
    }
    else {
      TRANS = 'N';
    }

    F77_FUNC(dsyrk,DSYRK)(&UPLO, &TRANS, &N, &K, &alpha, A, &LDA,
                          &beta, C, &LDC, 1, 1);
  }

  void IpBlasDtrsm(bool trans, Index ndim, Index nrhs, Number alpha,
                   const Number* A, Index ldA, Number* B, Index ldB)
  {
    ipfint M=ndim, N=nrhs, LDA=ldA, LDB=ldB;

    char SIDE = 'L';
    char UPLO = 'L';
    char TRANSA;
    if (trans) {
      TRANSA = 'T';
    }
    else {
      TRANSA = 'N';
    }
    char DIAG = 'N';

    F77_FUNC(dtrsm,DTRSM)(&SIDE, &UPLO, &TRANSA, &DIAG, &M, &N,
                          &alpha, A, &LDA, B, &LDB, 1, 1, 1, 1);
  }

#else
  /* Interface to CBLAS routine DDOT. */
  Number IpBlasDdot(Index size, const Number *x, Index incX, const Number *y,
                    Index incY)
  {
    Not Implemented Yet!
  }

  /* Interface to CBLAS routine DNRM2. */
  Number IpBlasDnrm2(Index size, const Number *x, Index incX)
  {
    Not Implemented Yet!
  }

  /* Interface to CBLAS routine DASUM. */
  Number IpBlasDasum(Index size, const Number *x, Index incX)
  {
    Not Implemented Yet!
  }

  /* Interface to CBLAS routine DASUM. */
  Index IpBlasIdamax(Index size, const Number *x, Index incX)
  {
    Not Implemented Yet!
  }

  /* Interface to CBLAS routine DCOPY. */
  void IpBlasDcopy(Index size, const Number *x, Index incX, Number *y, Index incY)
  {
    Not Implemented Yet!
  }

  /* Interface to CBLAS routine DAXPY. */
  void IpBlasDaxpy(Index size, Number alpha, const Number *x, Index incX, Number *y,
                   Index incY)
  {
    Not Implemented Yet!
  }

  /* Interface to CBLAS routine DSCAL. */
  void IpBlasDscal(Index size, Number alpha, Number *x, Index incX)
  {
    Not Implemented Yet!
  }

#endif

} // namespace Ipopt
