// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpZeroMatrix.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpZeroMatrix.hpp"

namespace Ipopt
{

  ZeroMatrix::ZeroMatrix(const MatrixSpace* owner_space)
      :
      Matrix(owner_space)
  {}

  ZeroMatrix::~ZeroMatrix()
  {}

  void ZeroMatrix::MultVectorImpl(Number alpha, const Vector &x,
                                  Number beta, Vector &y) const
  {
    //  A few sanity checks
    DBG_ASSERT(NCols()==x.Dim());
    DBG_ASSERT(NRows()==y.Dim());

    // Take care of the y part of the addition
    if( beta!=0.0 ) {
      y.Scal(beta);
    }
    else {
      y.Set(0.0);  // In case y hasn't been initialized yet
    }
  }

  void ZeroMatrix::TransMultVectorImpl(Number alpha, const Vector &x,
                                       Number beta, Vector &y) const
  {
    //  A few sanity checks
    DBG_ASSERT(NCols()==y.Dim());
    DBG_ASSERT(NRows()==x.Dim());

    // Take care of the y part of the addition
    if( beta!=0.0 ) {
      y.Scal(beta);
    }
    else {
      y.Set(0.0);  // In case y hasn't been initialized yet
    }
  }

  void ZeroMatrix::PrintImpl(const Journalist& jnlst,
                             EJournalLevel level,
                             EJournalCategory category,
                             const std::string& name,
                             Index indent,
                             const std::string& prefix) const
  {
    jnlst.Printf(level, category, "\n");
    jnlst.PrintfIndented(level, category, indent,
                         "%sZeroMatrix \"%s\" with %d row and %d columns components:\n",
                         prefix.c_str(), name.c_str(), NRows(), NCols());
  }
} // namespace Ipopt
