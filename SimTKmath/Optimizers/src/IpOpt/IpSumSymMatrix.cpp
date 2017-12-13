// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpSumSymMatrix.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpoptConfig.h"
#include "IpSumSymMatrix.hpp"

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

  SumSymMatrix::SumSymMatrix(const SumSymMatrixSpace* owner_space)
      :
      SymMatrix(owner_space),
      factors_(owner_space->NTerms(), 1.0),
      matrices_(owner_space->NTerms()),
      owner_space_(owner_space)
  {}

  SumSymMatrix::~SumSymMatrix()
  {}

  void SumSymMatrix::SetTerm(Index iterm, Number factor,
                             const SymMatrix& matrix)
  {
    DBG_ASSERT(iterm<owner_space_->NTerms());
    factors_[iterm] = factor;
    matrices_[iterm] = &matrix;
  }

  void SumSymMatrix::GetTerm(Index iterm, Number& factor, SmartPtr<const SymMatrix>& matrix) const
  {
    DBG_ASSERT(iterm<owner_space_->NTerms());
    factor = factors_[iterm];
    matrix = matrices_[iterm];
  }


  Index SumSymMatrix::NTerms() const
  {
    return owner_space_->NTerms();
  }

  void SumSymMatrix::MultVectorImpl(Number alpha, const Vector &x,
                                    Number beta, Vector &y) const
  {
    //  A few sanity checks
    DBG_ASSERT(Dim()==x.Dim());
    DBG_ASSERT(Dim()==y.Dim());

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

  bool SumSymMatrix::HasValidNumbersImpl() const
  {
    for (Index iterm=0; iterm<NTerms(); iterm++) {
      DBG_ASSERT(IsValid(matrices_[iterm]));
      if (!matrices_[iterm]->HasValidNumbers()) {
        return false;
      }
    }
    return true;
  }

  void SumSymMatrix::PrintImpl(const Journalist& jnlst,
                               EJournalLevel level,
                               EJournalCategory category,
                               const std::string& name,
                               Index indent,
                               const std::string& prefix) const
  {
    jnlst.Printf(level, category, "\n");
    jnlst.PrintfIndented(level, category, indent,
                         "%sSumSymMatrix \"%s\" of dimension %d with %d terms:\n",
                         prefix.c_str(), name.c_str(), Dim(), NTerms());
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

  void SumSymMatrixSpace::SetTermSpace
  (Index term_idx, const SymMatrixSpace& space)
  {
    while(term_idx >= (Index)term_spaces_.size()) {
      term_spaces_.push_back(NULL);
    }
    term_spaces_[term_idx] = &space;
  }

  SmartPtr<const SymMatrixSpace> SumSymMatrixSpace::GetTermSpace(Index term_idx) const
  {
    if (term_idx >= 0 && term_idx < (Index)term_spaces_.size()) {
      return term_spaces_[term_idx];
    }
    return NULL;
  }

  SumSymMatrix* SumSymMatrixSpace::MakeNewSumSymMatrix() const
  {
    DBG_ASSERT(nterms_ == (Index)term_spaces_.size());
    return new SumSymMatrix(this);
  }

  SymMatrix* SumSymMatrixSpace::MakeNewSymMatrix() const
  {
    return MakeNewSumSymMatrix();
  }
} // namespace Ipopt
