// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpLowRankUpdateSymMatrix.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Andreas Waechter                IBM    2005-12-25

#include "IpLowRankUpdateSymMatrix.hpp"

namespace SimTKIpopt
{

#ifdef IP_DEBUG
  static const Index dbg_verbosity = 0;
#endif

  LowRankUpdateSymMatrix::LowRankUpdateSymMatrix(const LowRankUpdateSymMatrixSpace* owner_space)
      :
      SymMatrix(owner_space),
      owner_space_(owner_space)
  {}

  LowRankUpdateSymMatrix::~LowRankUpdateSymMatrix()
  {}

  void LowRankUpdateSymMatrix::MultVectorImpl(Number alpha, const Vector &x,
      Number beta, Vector &y) const
  {
    DBG_START_METH("LowRankUpdateSymMatrix::MultVectorImpl",
                   dbg_verbosity);
    //  A few sanity checks
    DBG_ASSERT(Dim()==x.Dim());
    DBG_ASSERT(Dim()==y.Dim());
    DBG_ASSERT(IsValid(D_));

    SmartPtr<const Matrix> P_LR =  P_LowRank();
    if (IsNull(P_LR)) {

      // Diagonal part
      if( beta!=0.0 ) {
        SmartPtr<Vector> tmp_vec = x.MakeNewCopy();
        tmp_vec->ElementWiseMultiply(*D_);
        y.AddOneVector(alpha, *tmp_vec, beta);
      }
      else {
        y.AddOneVector(alpha, x, 0.);
        y.ElementWiseMultiply(*D_);
      }
      DBG_PRINT_VECTOR(2, "y = D*alpha*x", y);

      if (IsValid(V_)) {
        // Positive update
        V_->LRMultVector(alpha, x, 1., y);
      }

      if (IsValid(U_)) {
        // Negative update
        U_->LRMultVector(-alpha, x, 1., y);
      }
    }
    else {
      if (ReducedDiag()) {
        // Get everything into the smaller space
        SmartPtr<const VectorSpace> LR_vec_space = LowRankVectorSpace();
        SmartPtr<Vector> small_x = LR_vec_space->MakeNew();
        P_LR->TransMultVector(1., x, 0., *small_x);

        // Diagonal part in small space
        SmartPtr<Vector> small_y = LR_vec_space->MakeNew();
        small_y->Copy(*small_x);
        small_y->ElementWiseMultiply(*D_);

        if (IsValid(V_)) {
          // Positive update
          V_->LRMultVector(1., *small_x, 1., *small_y);
        }

        if (IsValid(U_)) {
          // Negative update
          U_->LRMultVector(-1., *small_x, 1., *small_y);
        }

        // Get the result back into the large space
        P_LR->MultVector(alpha, *small_y, beta, y);
      }
      else {
        // Diagonal part
        SmartPtr<Vector> tmp = x.MakeNewCopy();
        tmp->ElementWiseMultiply(*D_);
        y.AddOneVector(alpha, *tmp, beta);

        // Get x into the smaller space
        SmartPtr<const VectorSpace> LR_vec_space = LowRankVectorSpace();
        SmartPtr<Vector> small_x = LR_vec_space->MakeNew();
        P_LR->TransMultVector(1., x, 0., *small_x);

        SmartPtr<Vector> small_y = LR_vec_space->MakeNew();
        if (IsValid(V_)) {
          // Positive update
          V_->LRMultVector(1., *small_x, 0., *small_y);
        }
        else {
          small_y->Set(0.);
        }

        if (IsValid(U_)) {
          // Negative update
          U_->LRMultVector(-1., *small_x, 1., *small_y);
        }

        // Get the result back into the large space
        P_LR->MultVector(alpha, *small_y, 1., y);
      }
    }
  }

  bool LowRankUpdateSymMatrix::HasValidNumbersImpl() const
  {
    DBG_ASSERT(IsValid(D_));
    if (!D_->HasValidNumbers()) {
      return false;
    }
    if (IsValid(V_)) {
      if (!V_->HasValidNumbers()) {
        return false;
      }
    }
    if (IsValid(U_)) {
      if (!U_->HasValidNumbers()) {
        return false;
      }
    }
    return true;
  }

  void LowRankUpdateSymMatrix::PrintImpl(const Journalist& jnlst,
                                         EJournalLevel level,
                                         EJournalCategory category,
                                         const std::string& name,
                                         Index indent,
                                         const std::string& prefix) const
  {
    jnlst.Printf(level, category, "\n");
    jnlst.PrintfIndented(level, category, indent,
                         "%sLowRankUpdateSymMatrix \"%s\" with %d rows and columns:\n",
                         prefix.c_str(), name.c_str(), Dim());

    if (ReducedDiag()) {
      jnlst.PrintfIndented(level, category, indent+1,
                           "%sThis matrix has reduced diagonal.\n", prefix.c_str());
    }
    else {
      jnlst.PrintfIndented(level, category, indent+1,
                           "%sThis matrix has full diagonal.\n", prefix.c_str());
    }
    jnlst.PrintfIndented(level, category, indent+1,
                         "%sDiagonal matrix:\n", prefix.c_str());
    if (IsValid(D_)) {
      D_->Print(&jnlst, level, category, name+"-D", indent+1, prefix);
    }
    else {
      jnlst.PrintfIndented(level, category, indent,
                           "%sDiagonal matrix not set!\n", prefix.c_str());
    }

    jnlst.PrintfIndented(level, category, indent+1,
                         "%sMultiVectorMatrix V for positive update:\n", prefix.c_str());
    if (IsValid(V_)) {
      V_->Print(&jnlst, level, category, name+"-V", indent+1, prefix);
    }
    else {
      jnlst.PrintfIndented(level, category, indent,
                           "%sV matrix not set!\n", prefix.c_str());
    }

    jnlst.PrintfIndented(level, category, indent+1,
                         "%sMultiVectorMatrix U for positive update:\n", prefix.c_str());
    if (IsValid(U_)) {
      U_->Print(&jnlst, level, category, name+"-U", indent+1, prefix);
    }
    else {
      jnlst.PrintfIndented(level, category, indent,
                           "%sU matrix not set!\n", prefix.c_str());
    }
  }
} // namespace Ipopt
