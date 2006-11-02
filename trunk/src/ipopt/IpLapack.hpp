// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpLapack.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Andreas Waechter              IBM    2005-12-25

#ifndef __IPLAPACK_HPP__
#define __IPLAPACK_HPP__

#include "IpUtils.hpp"
#include "IpException.hpp"

namespace Ipopt
{
  DECLARE_STD_EXCEPTION(LAPACK_NOT_INCLUDED);

  /** Wrapper for LAPACK subroutine DPOTRS.  Solving a linear system
   *  given a Cholesky factorization.  We assume that the Cholesky
   *  factor is lower traiangular. */
  void IpLapackDpotrs(Index ndim, Index nrhs, const Number *a, Index lda,
                      Number *b, Index ldb);

  /** Wrapper for LAPACK subroutine DPOTRF.  Compute Cholesky
   *  factorization (lower triangular factor).  info is the return
   *  value from the LAPACK routine. */
  void IpLapackDpotrf(Index ndim, Number *a, Index lda, Index& info);

  /** Wrapper for LAPACK subroutine DSYEV.  Compute the Eigenvalue
   *  decomposition for a given matrix.  If compute_eigenvectors is
   *  true, a will contain the eigenvectors in its columns on
   *  return.  */
  void IpLapackDsyev(bool compute_eigenvectors, Index ndim, Number *a,
                     Index lda, Number *w, Index& info);

} // namespace Ipopt

#endif
