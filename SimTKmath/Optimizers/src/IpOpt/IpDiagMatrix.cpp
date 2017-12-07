// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpDiagMatrix.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpDiagMatrix.hpp"

namespace SimTKIpopt
{

  DiagMatrix::DiagMatrix(const SymMatrixSpace* owner_space)
      :
      SymMatrix(owner_space)
  {}

  DiagMatrix::~DiagMatrix()
  {}

  void DiagMatrix::MultVectorImpl(Number alpha, const Vector &x,
                                  Number beta, Vector &y) const
  {
    //  A few sanity checks
    DBG_ASSERT(Dim()==x.Dim());
    DBG_ASSERT(Dim()==y.Dim());
    DBG_ASSERT(IsValid(diag_));

    // Take care of the y part of the addition
    if( beta!=0.0 ) {
      y.Scal(beta);
    }
    else {
      y.Set(0.0);  // In case y hasn't been initialized yet
    }

    SmartPtr<Vector> tmp_vec = y.MakeNew();
    tmp_vec->Copy(x);
    tmp_vec->ElementWiseMultiply(*diag_);
    y.Axpy(alpha, *tmp_vec);
  }

  bool DiagMatrix::HasValidNumbersImpl() const
  {
    DBG_ASSERT(IsValid(diag_));
    return diag_->HasValidNumbers();
  }

  void DiagMatrix::PrintImpl(const Journalist& jnlst,
                             EJournalLevel level,
                             EJournalCategory category,
                             const std::string& name,
                             Index indent,
                             const std::string& prefix) const
  {
    jnlst.Printf(level, category, "\n");
    jnlst.PrintfIndented(level, category, indent,
                         "%sDiagMatrix \"%s\" with %d rows and columns, and with diagonal elements:\n",
                         prefix.c_str(), name.c_str(), Dim());
    if (IsValid(diag_)) {
      diag_->Print(&jnlst, level, category, name, indent+1, prefix);
    }
    else {
      jnlst.PrintfIndented(level, category, indent,
                           "%sDiagonal elements not set!\n", prefix.c_str());
    }
  }
} // namespace Ipopt
