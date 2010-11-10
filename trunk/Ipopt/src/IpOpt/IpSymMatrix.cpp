// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpSymMatrix.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpSymMatrix.hpp"

namespace Ipopt
{

  SymMatrix::SymMatrix(const SymMatrixSpace* owner_space)
      :
      Matrix(owner_space),
      owner_space_(owner_space)
  {}

  SymMatrix::~SymMatrix()
  {}

  SmartPtr<const SymMatrixSpace> SymMatrix::OwnerSymMatrixSpace() const
  {
    return owner_space_;
  }

  void SymMatrix::TransMultVectorImpl(Number alpha, const Vector& x, Number beta,
                                      Vector& y) const
  {
    // Since this matrix is symetric, this is the same operation as
    // MultVector
    MultVector(alpha, x, beta, y);
  }

  SymMatrixSpace::SymMatrixSpace(Index dim)
      :
      MatrixSpace(dim,dim)
  {}

  SymMatrixSpace::~SymMatrixSpace()
  {}


  Matrix* SymMatrixSpace::MakeNew() const
  {
    return MakeNewSymMatrix();
  }

} // namespace Ipopt
