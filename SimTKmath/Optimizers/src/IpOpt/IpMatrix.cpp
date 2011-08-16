// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpMatrix.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpMatrix.hpp"

namespace Ipopt
{
  void Matrix::MultVector(Number alpha, const Vector& x, Number beta,
                          Vector& y) const
  {
    MultVectorImpl(alpha, x, beta, y);
  }

  void Matrix::TransMultVector(Number alpha, const Vector& x, Number beta,
                               Vector& y) const
  {
    TransMultVectorImpl(alpha, x, beta, y);
  }

  void Matrix::AddMSinvZ(Number alpha, const Vector& S, const Vector& Z,
                         Vector& X) const
  {
    AddMSinvZImpl(alpha, S, Z, X);
  }

  void Matrix::SinvBlrmZMTdBr(Number alpha, const Vector& S,
                              const Vector& R, const Vector& Z,
                              const Vector& D, Vector& X) const
  {
    SinvBlrmZMTdBrImpl(alpha, S, R, Z, D, X);
  }

  // Prototype for specialize methods (can and should be overloaded)
  void Matrix::AddMSinvZImpl(Number alpha, const Vector& S, const Vector& Z,
                             Vector& X) const
  {
    SmartPtr<Vector> tmp = S.MakeNew();
    tmp->AddVectorQuotient(1., Z, S, 0.);
    MultVector(alpha, *tmp, 1., X);
  }

  void Matrix::SinvBlrmZMTdBrImpl(Number alpha, const Vector& S,
                                  const Vector& R, const Vector& Z,
                                  const Vector& D, Vector& X) const
  {
    TransMultVector(alpha, D, 0., X);
    X.ElementWiseMultiply(Z);
    X.Axpy(1., R);
    X.ElementWiseDivide(S);
  }

  bool Matrix::HasValidNumbers() const
  {
    if (valid_cache_tag_ != GetTag()) {
      cached_valid_ = HasValidNumbersImpl();
      valid_cache_tag_ = GetTag();
    }
    return cached_valid_;
  }

  void Matrix::Print(SmartPtr<const Journalist> jnlst,
                     EJournalLevel level,
                     EJournalCategory category,
                     const std::string& name,
                     Index indent,
                     const std::string& prefix) const
  {
    if (IsValid(jnlst) && jnlst->ProduceOutput(level, category)) {
      PrintImpl(*jnlst, level, category, name, indent, prefix);
    }
  }

  void Matrix::Print(const Journalist& jnlst,
                     EJournalLevel level,
                     EJournalCategory category,
                     const std::string& name,
                     Index indent,
                     const std::string& prefix) const
  {
    if (jnlst.ProduceOutput(level, category)) {
      PrintImpl(jnlst, level, category, name, indent, prefix);
    }
  }

  MatrixSpace::MatrixSpace(Index nRows, Index nCols)
      :
      nRows_(nRows),
      nCols_(nCols)
  {}

  MatrixSpace::~MatrixSpace()
  {}

  bool MatrixSpace::IsMatrixFromSpace(const Matrix& matrix) const
  {
    return (matrix.OwnerSpace() == this);
  }

} // namespace Ipopt
