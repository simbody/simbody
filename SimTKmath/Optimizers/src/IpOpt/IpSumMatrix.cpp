// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpSumMatrix.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpoptConfig.h"
#include "IpSumMatrix.hpp"

#ifdef HAVE_CSTDIO
# include <cstdio>
#else
# ifdef HAVE_STDIO_H
#  include <stdio.h>
# else
#  error "don't have header file for stdio"
# endif
#endif


// Keeps MS VC++ 8 quiet about sprintf, strcpy, etc.
#ifdef _MSC_VER
#pragma warning(disable:4996)
#endif

namespace SimTKIpopt
{

  SumMatrix::SumMatrix(const SumMatrixSpace* owner_space)
      :
      Matrix(owner_space),
      factors_(owner_space->NTerms(), 1.0),
      matrices_(owner_space->NTerms()),
      owner_space_(owner_space)
  {}

  SumMatrix::~SumMatrix()
  {}

  void SumMatrix::SetTerm(Index iterm, Number factor,
                          const Matrix& matrix)
  {
    DBG_ASSERT(iterm<owner_space_->NTerms());
    factors_[iterm] = factor;
    matrices_[iterm] = &matrix;
  }

  void SumMatrix::GetTerm(Index iterm, Number& factor, SmartPtr<const Matrix>& matrix) const
  {
    DBG_ASSERT(iterm<owner_space_->NTerms());
    factor = factors_[iterm];
    matrix = matrices_[iterm];
  }

  Index SumMatrix::NTerms() const
  {
    return owner_space_->NTerms();
  }

  void SumMatrix::MultVectorImpl(Number alpha, const Vector &x,
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

    for (Index iterm=0; iterm<NTerms(); iterm++) {
      DBG_ASSERT(IsValid(matrices_[iterm]));
      matrices_[iterm]->MultVector(alpha*factors_[iterm], x,
                                   1.0, y);
    }
  }

  void SumMatrix::TransMultVectorImpl(Number alpha, const Vector& x,
                                      Number beta, Vector& y) const
  {
    //  A few sanity checks
    DBG_ASSERT(NRows()==x.Dim());
    DBG_ASSERT(NCols()==y.Dim());

    // Take care of the y part of the addition
    if( beta!=0.0 ) {
      y.Scal(beta);
    }
    else {
      y.Set(0.0);  // In case y hasn't been initialized yet
    }

    for (Index iterm=0; iterm<NTerms(); iterm++) {
      DBG_ASSERT(IsValid(matrices_[iterm]));
      matrices_[iterm]->TransMultVector(alpha*factors_[iterm], x,
                                        1.0, y);
    }
  }

  bool SumMatrix::HasValidNumbersImpl() const
  {
    for (Index iterm=0; iterm<NTerms(); iterm++) {
      DBG_ASSERT(IsValid(matrices_[iterm]));
      if (!matrices_[iterm]->HasValidNumbers()) {
        return false;
      }
    }
    return true;
  }

  void SumMatrix::PrintImpl(const Journalist& jnlst,
                            EJournalLevel level,
                            EJournalCategory category,
                            const std::string& name,
                            Index indent,
                            const std::string& prefix) const
  {
    jnlst.Printf(level, category, "\n");
    jnlst.PrintfIndented(level, category, indent,
                         "%sSumMatrix \"%s\" of dimension %d x %d with %d terms:\n",
                         prefix.c_str(), name.c_str(), NRows(), NCols(), NTerms());
    for (Index iterm=0; iterm<NTerms(); iterm++) {
      jnlst.PrintfIndented(level, category, indent,
                           "%sTerm %d with factor %23.16e and the following matrix:\n",
                           prefix.c_str(), iterm, factors_[iterm]);
      char buffer[256];
      sprintf(buffer, "Term: %d", iterm);
      std::string name = buffer;
      matrices_[iterm]->Print(&jnlst, level, category, name, indent+1, prefix);
    }
  }

  void SumMatrixSpace::SetTermSpace(Index term_idx, const MatrixSpace& mat_space)
  {
    while(term_idx >= (Index)term_spaces_.size()) {
      term_spaces_.push_back(NULL);
    }
    term_spaces_[term_idx] = &mat_space;
  }

  SmartPtr<const MatrixSpace> SumMatrixSpace::GetTermSpace(Index term_idx) const
  {
    if (term_idx >= 0 && term_idx < (Index)term_spaces_.size()) {
      return term_spaces_[term_idx];
    }
    return NULL;
  }

  SumMatrix* SumMatrixSpace::MakeNewSumMatrix() const
  {
    DBG_ASSERT(nterms_ == (Index)term_spaces_.size());
    return new SumMatrix(this);
  }

  Matrix* SumMatrixSpace::MakeNew() const
  {
    return MakeNewSumMatrix();
  }
} // namespace Ipopt
