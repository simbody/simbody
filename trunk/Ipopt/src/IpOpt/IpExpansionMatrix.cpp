// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpExpansionMatrix.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpExpansionMatrix.hpp"
#include "IpDenseVector.hpp"

namespace Ipopt
{

#ifdef IP_DEBUG
  static const Index dbg_verbosity = 0;
#endif

  ExpansionMatrix::ExpansionMatrix(const ExpansionMatrixSpace* owner_space)
      :
      Matrix(owner_space),
      owner_space_(owner_space)
  {}

  ExpansionMatrix::~ExpansionMatrix()
  {}

  void ExpansionMatrix::MultVectorImpl(Number alpha, const Vector &x,
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

    // See if we can understand the data
    const DenseVector* dense_x = dynamic_cast<const DenseVector*>(&x);
    DBG_ASSERT(dense_x); /* ToDo: Implement others */
    DenseVector* dense_y = dynamic_cast<DenseVector*>(&y);
    DBG_ASSERT(dense_y); /* ToDo: Implement others */

    const Index* exp_pos = ExpandedPosIndices();

    if (dense_x && dense_y) {
      Number* yvals=dense_y->Values();
      if (dense_x->IsHomogeneous()) {
        Number val = alpha * dense_x->Scalar();
        for(Index i=0; i<NCols(); i++) {
          yvals[exp_pos[i]] += val;
        }
      }
      else {
        const Number* xvals=dense_x->Values();
        for(Index i=0; i<NCols(); i++) {
          yvals[exp_pos[i]] += alpha * xvals[i];
        }
      }
    }
  }

  void ExpansionMatrix::TransMultVectorImpl(Number alpha, const Vector &x,
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

    // See if we can understand the data
    const DenseVector* dense_x = dynamic_cast<const DenseVector*>(&x);
    DBG_ASSERT(dense_x); /* ToDo: Implement others */
    DenseVector* dense_y = dynamic_cast<DenseVector*>(&y);
    DBG_ASSERT(dense_y); /* ToDo: Implement others */

    const Index* exp_pos = ExpandedPosIndices();

    if (dense_x && dense_y) {
      Number* yvals=dense_y->Values();
      if (dense_x->IsHomogeneous()) {
        Number val = alpha * dense_x->Scalar();
        for(Index i=0; i<NCols(); i++) {
          yvals[i] += val;
        }
      }
      else {
        const Number* xvals=dense_x->Values();
        for(Index i=0; i<NCols(); i++) {
          yvals[i] += alpha * xvals[exp_pos[i]];
        }
      }
    }
  }

  // Specialized method (overloaded from IpMatrix)
  void ExpansionMatrix::AddMSinvZImpl(Number alpha, const Vector& S,
                                      const Vector& Z, Vector& X) const
  {
    DBG_ASSERT(NCols()==S.Dim());
    DBG_ASSERT(NCols()==Z.Dim());
    DBG_ASSERT(NRows()==X.Dim());

    const DenseVector* dense_S = dynamic_cast<const DenseVector*>(&S);
    DBG_ASSERT(dense_S);
    const DenseVector* dense_Z = dynamic_cast<const DenseVector*>(&Z);
    DBG_ASSERT(dense_Z);
    DenseVector* dense_X = dynamic_cast<DenseVector*>(&X);
    DBG_ASSERT(dense_X);

    // if vector S is homogeneous type, call the default implementation
    // ToDo: find out how often the default implementation is called and
    // if we should implement specialized method for the homogenenous
    // case
    if (dense_S->IsHomogeneous()) {
      DBG_ASSERT(false && "dense_S is homogeneous - implement specialized option?");
      Matrix::AddMSinvZImpl(alpha, S, Z, X);
      return;
    }

    const Index* exp_pos = ExpandedPosIndices();
    const Number* vals_S = dense_S->Values();
    Number* vals_X = dense_X->Values();

    if (dense_Z->IsHomogeneous()) {
      Number val = alpha*dense_Z->Scalar();
      for(Index i=0; i<NCols(); i++) {
        vals_X[exp_pos[i]] += val/vals_S[i];
      }
    }
    else {
      const Number* vals_Z = dense_Z->Values();
      if (alpha==1.) {
        for(Index i=0; i<NCols(); i++) {
          vals_X[exp_pos[i]] += vals_Z[i]/vals_S[i];
        }
      }
      else if (alpha==-1.) {
        for(Index i=0; i<NCols(); i++) {
          vals_X[exp_pos[i]] -= vals_Z[i]/vals_S[i];
        }
      }
      else {
        for(Index i=0; i<NCols(); i++) {
          vals_X[exp_pos[i]] += alpha*vals_Z[i]/vals_S[i];
        }
      }
    }
  }

  void ExpansionMatrix::SinvBlrmZMTdBrImpl(Number alpha, const Vector& S,
      const Vector& R, const Vector& Z,
      const Vector& D, Vector& X) const
  {
    DBG_START_METH("ExpansionMatrix::SinvBlrmZMTdBrImpl",
                   dbg_verbosity);

    DBG_ASSERT(NCols()==S.Dim());
    DBG_ASSERT(NCols()==R.Dim());
    DBG_ASSERT(NCols()==Z.Dim());
    DBG_ASSERT(NRows()==D.Dim());
    DBG_ASSERT(NCols()==X.Dim());

    const DenseVector* dense_S = dynamic_cast<const DenseVector*>(&S);
    DBG_ASSERT(dense_S);
    const DenseVector* dense_R = dynamic_cast<const DenseVector*>(&R);
    DBG_ASSERT(dense_R);
    const DenseVector* dense_Z = dynamic_cast<const DenseVector*>(&Z);
    DBG_ASSERT(dense_Z);
    const DenseVector* dense_D = dynamic_cast<const DenseVector*>(&D);
    DBG_ASSERT(dense_D);
    DenseVector* dense_X = dynamic_cast<DenseVector*>(&X);
    DBG_ASSERT(dense_X);

    // if the vectors S or D are of the homogeneous type, revert to the
    // default implementation
    // ToDo: find out how often the default implementation is called and
    // if we should implement specialized method for the homogenenous
    // case

    if (dense_S->IsHomogeneous() ||
        dense_D->IsHomogeneous()) {
      DBG_ASSERT(false && "dense_S or dense_D is homogeneous - implement specialized option?");
      Matrix::SinvBlrmZMTdBrImpl(alpha, S, R, Z, D, X);
      return;
    }

    const Index* exp_pos = ExpandedPosIndices();
    const Number* vals_S = dense_S->Values();
    const Number* vals_D = dense_D->Values();
    Number* vals_X = dense_X->Values();

    if (dense_R->IsHomogeneous()) {
      Number scalar_R = dense_R->Scalar();
      if (dense_Z->IsHomogeneous()) {
        Number val = alpha*dense_Z->Scalar();
        if (val == 0.) {
          for (Index i=0; i<NCols(); i++) {
            // ToDo could treat val == 0 extra
            vals_X[i] = (scalar_R)/vals_S[i];
          }
        }
        else {
          for (Index i=0; i<NCols(); i++) {
            // ToDo could treat val == 0 extra
            vals_X[i] = (scalar_R + val*vals_D[exp_pos[i]])/vals_S[i];
          }
        }
      }
      else {
        const Number* vals_Z = dense_Z->Values();
        if (alpha==1.) {
          for (Index i=0; i<NCols(); i++) {
            vals_X[i] = (scalar_R + vals_Z[i]*vals_D[exp_pos[i]])/vals_S[i];
          }
        }
        else if (alpha==-1.) {
          for (Index i=0; i<NCols(); i++) {
            vals_X[i] = (scalar_R - vals_Z[i]*vals_D[exp_pos[i]])/vals_S[i];
          }
        }
        else {
          for (Index i=0; i<NCols(); i++) {
            vals_X[i] = (scalar_R + alpha*vals_Z[i]*vals_D[exp_pos[i]])/vals_S[i];
          }
        }
      }
    }
    else {
      const Number* vals_R = dense_R->Values();
      if (dense_Z->IsHomogeneous()) {
        Number val = alpha*dense_Z->Scalar();
        for (Index i=0; i<NCols(); i++) {
          vals_X[i] = (vals_R[i] + val*vals_D[exp_pos[i]])/vals_S[i];
        }
      }
      else {
        const Number* vals_Z = dense_Z->Values();
        if (alpha==1.) {
          for (Index i=0; i<NCols(); i++) {
            vals_X[i] = (vals_R[i] + vals_Z[i]*vals_D[exp_pos[i]])/vals_S[i];
          }
        }
        else if (alpha==-1.) {
          for (Index i=0; i<NCols(); i++) {
            vals_X[i] = (vals_R[i] - vals_Z[i]*vals_D[exp_pos[i]])/vals_S[i];
          }
        }
        else {
          for (Index i=0; i<NCols(); i++) {
            vals_X[i] = (vals_R[i] + alpha*vals_Z[i]*vals_D[exp_pos[i]])/vals_S[i];
          }
        }
      }
    }
  }

  void ExpansionMatrix::PrintImpl(const Journalist& jnlst,
                                  EJournalLevel level,
                                  EJournalCategory category,
                                  const std::string& name,
                                  Index indent,
                                  const std::string& prefix) const
  {
    jnlst.Printf(level, category, "\n");
    jnlst.PrintfIndented(level, category, indent,
                         "%sExpansionMatrix \"%s\" with %d nonzero elements:\n",
                         prefix.c_str(), name.c_str(), NCols());

    const Index* exp_pos = ExpandedPosIndices();

    for (Index i=0; i<NCols(); i++) {
      jnlst.PrintfIndented(level, category, indent,
                           "%s%s[%5d,%5d]=%23.16e  (%d)\n",
                           prefix.c_str(), name.c_str(), exp_pos[i]+1,
                           i+1, 1., i);
    }
  }

  ExpansionMatrixSpace::ExpansionMatrixSpace(Index NLargeVec,
      Index NSmallVec,
      const Index *ExpPos,
      const int offset /*= 0*/)
      :
      MatrixSpace(NLargeVec, NSmallVec),
      expanded_pos_(NULL),
      compressed_pos_(NULL)
  {
    if(NCols()>0) {
      expanded_pos_  = new Index[NCols()];
    }
    if(NRows()>0) {
      compressed_pos_ = new Index[NRows()];
    }
    for (Index j=0; j<NRows(); j++) {
      compressed_pos_[j] = -1;
    }
    for(Index i=0; i<NCols(); i++) {
      //ToDo decide for offset
      DBG_ASSERT(ExpPos[i]-offset<NRows() && ExpPos[i]-offset>=0);
      expanded_pos_[i]=ExpPos[i]-offset;
      compressed_pos_[ExpPos[i]-offset] = i;
    }
  }

} // namespace Ipopt
