// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpLapack.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Andreas Waechter              IBM    2005-12-25

#include "IpoptConfig.h"
#include "IpLapack.hpp"

// Prototypes for the LAPACK routines
extern "C"
{
  /** LAPACK Fortran subroutine DPOTRS. */
  void F77_FUNC(dpotrs,DPOTRS)(char *uplo, ipfint *n,
                               ipfint *nrhs, const double *A, ipfint *ldA,
                               double *B, ipfint *ldB, ipfint *info,
                               int uplo_len);
  /** LAPACK Fortran subroutine DPOTRF. */
  void F77_FUNC(dpotrf,DPOTRF)(char *uplo, ipfint *n,
                               double *A, ipfint *ldA,
                               ipfint *info, int uplo_len);

  /** LAPACK Fortran subroutine DSYEV */
  void F77_FUNC(dsyev,DSYEV)(char *jobz, char *uplo, ipfint *n,
                             double *A, ipfint *ldA, double *W,
                             double *WORK, ipfint *LWORK, ipfint *info,
                             int jobz_len, int uplo_len);
}

namespace Ipopt
{
  /* Interface to FORTRAN routine DPOTRS. */
  void IpLapackDpotrs(Index ndim, Index nrhs, const Number *a, Index lda,
                      Number *b, Index ldb)
  {
#ifdef COIN_HAS_LAPACK
    ipfint N=ndim, NRHS=nrhs, LDA=lda, LDB=ldb, INFO;
    char uplo = 'L';

    F77_FUNC(dpotrs,DPOTRS)(&uplo, &N, &NRHS, a, &LDA, b, &LDB, &INFO, 1);
    DBG_ASSERT(INFO==0);
#else

    std::string msg = "Ipopt has been compiled without LAPACK routine DPOTRS, but options are chosen that require this dependency.  Abort.";
    THROW_EXCEPTION(LAPACK_NOT_INCLUDED, msg);
#endif

  }

  void IpLapackDpotrf(Index ndim, Number *a, Index lda, Index& info)
  {
#ifdef COIN_HAS_LAPACK
    ipfint N=ndim, LDA=lda, INFO;

    char UPLO = 'L';

    F77_FUNC(dpotrf,DPOTRF)(&UPLO, &N, a, &LDA, &INFO, 1);

    info = INFO;
#else

    std::string msg = "Ipopt has been compiled without LAPACK routine DPOTRF, but options are chosen that require this dependency.  Abort.";
    THROW_EXCEPTION(LAPACK_NOT_INCLUDED, msg);
#endif

  }

  void IpLapackDsyev(bool compute_eigenvectors, Index ndim, Number *a,
                     Index lda, Number *w, Index& info)
  {
#ifdef COIN_HAS_LAPACK
    ipfint N=ndim, LDA=lda, INFO;

    char JOBZ;
    if (compute_eigenvectors) {
      JOBZ = 'V';
    }
    else {
      JOBZ = 'N';
    }
    char UPLO = 'L';

    // First we find out how large LWORK should be
    ipfint LWORK = -1;
    double WORK_PROBE;
    F77_FUNC(dsyev,DSYEV)(&JOBZ, &UPLO, &N, a, &LDA, w,
                          &WORK_PROBE, &LWORK, &INFO, 1, 1);
    DBG_ASSERT(INFO==0);

    LWORK = (ipfint) WORK_PROBE;
    DBG_ASSERT(LWORK>0);

    double* WORK = new double[LWORK];
    for (Index i=0; i<LWORK; i++) {
      WORK[i] = i;
    }
    F77_FUNC(dsyev,DSYEV)(&JOBZ, &UPLO, &N, a, &LDA, w,
                          WORK, &LWORK, &INFO, 1, 1);

    DBG_ASSERT(INFO>=0);
    info = INFO;

    delete [] WORK;
#else

    std::string msg = "Ipopt has been compiled without LAPACK routine DSYEV, but options are chosen that require this dependency.  Abort.";
    THROW_EXCEPTION(LAPACK_NOT_INCLUDED, msg);
#endif

  }
}
