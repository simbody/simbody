// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpDenseSymMatrix.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Andreas Waechter             IBM    2005-12-25

#include "IpDenseSymMatrix.hpp"
#include "IpDenseVector.hpp"
#include "IpDenseGenMatrix.hpp"
#include "IpBlas.hpp"

namespace SimTKIpopt
{

#ifdef IP_DEBUG
  static const Index dbg_verbosity = 0;
#endif

  DenseSymMatrix::DenseSymMatrix(const DenseSymMatrixSpace* owner_space)
      :
      SymMatrix(owner_space),
      owner_space_(owner_space),
      values_(new Number[NCols()*NRows()]),
      initialized_(false)
  {}

  DenseSymMatrix::~DenseSymMatrix()
  {
    delete [] values_;
  }

  void DenseSymMatrix::MultVectorImpl(Number alpha, const Vector &x,
                                      Number beta, Vector &y) const
  {
    //  A few sanity checks
    DBG_ASSERT(NCols()==x.Dim());
    DBG_ASSERT(NRows()==y.Dim());
    DBG_ASSERT(initialized_);

    // See if we can understand the data
    const DenseVector* dense_x = dynamic_cast<const DenseVector*>(&x);
    DBG_ASSERT(dense_x); /* ToDo: Implement others */
    DenseVector* dense_y = dynamic_cast<DenseVector*>(&y);
    DBG_ASSERT(dense_y); /* ToDo: Implement others */

    IpBlasDsymv(Dim(), alpha, values_, NRows(),
                dense_x->Values(), 1, beta, dense_y->Values(), 1);
  }

  void DenseSymMatrix::FillIdentity(Number factor /*=1.*/)
  {
    const Index dim = Dim();
    for (Index j=0; j<dim; j++) {
      values_[j + j*dim] = factor;
      for (Index i=j+1; i<dim; i++) {
        values_[i + j*dim] = 0.;
      }
    }
    ObjectChanged();
    initialized_ = true;
  }

  void DenseSymMatrix::AddMatrix(Number alpha, const DenseSymMatrix& A,
                                 Number beta)
  {
    DBG_ASSERT(beta==0. || initialized_);
    DBG_ASSERT(Dim()==A.Dim());

    if (alpha==0.)
      return;

    const Number* Avalues = A.Values();
    const Index dim = Dim();
    if (beta==0.) {
      for (Index j=0; j<dim; j++) {
        for (Index i=j; i<dim; i++) {
          values_[i+j*dim] = alpha*Avalues[i+j*dim];
        }
      }
    }
    else if (beta==1.) {
      for (Index j=0; j<dim; j++) {
        for (Index i=j; i<dim; i++) {
          values_[i+j*dim] += alpha*Avalues[i+j*dim];
        }
      }
    }
    else {
      for (Index j=0; j<dim; j++) {
        for (Index i=j; i<dim; i++) {
          values_[i+j*dim] = alpha*Avalues[i+j*dim] + beta*values_[i+j*dim];
        }
      }
    }
    ObjectChanged();
    initialized_ = true;
  }

  void DenseSymMatrix::HighRankUpdateTranspose(Number alpha,
      const MultiVectorMatrix& V1,
      const MultiVectorMatrix& V2,
      Number beta)
  {
    DBG_ASSERT(Dim()==V1.NCols());
    DBG_ASSERT(Dim()==V2.NCols());
    DBG_ASSERT(beta==0. || initialized_);

    const Index dim = Dim();
    if (beta==0.) {
      for (Index j=0; j<dim; j++) {
        for (Index i=j; i<dim; i++) {
          values_[i+j*dim] = alpha*V1.GetVector(i)->Dot(*V2.GetVector(j));
        }
      }
    }
    else {
      for (Index j=0; j<dim; j++) {
        for (Index i=j; i<dim; i++) {
          values_[i+j*dim] = alpha*V1.GetVector(i)->Dot(*V2.GetVector(j))
                             + beta*values_[i+j*dim];
        }
      }
    }
    initialized_ = true;
    ObjectChanged();
  }

  void DenseSymMatrix::HighRankUpdate(bool trans, Number alpha,
                                      const DenseGenMatrix& V,
                                      Number beta)
  {
    DBG_ASSERT((!trans && Dim()==V.NRows()) || (trans && Dim()==V.NCols()));
    DBG_ASSERT(beta==0. || initialized_);

    Index nrank;
    if (trans) {
      nrank = V.NRows();
    }
    else {
      nrank = V.NCols();
    }

    IpBlasDsyrk(trans, Dim(), nrank, alpha, V.Values(), V.NRows(),
                beta, values_, NRows());

    initialized_ = true;
    ObjectChanged();
  }

  void DenseSymMatrix::SpecialAddForLMSR1(const DenseVector& D,
                                          const DenseGenMatrix& L)
  {
    const Index dim = Dim();
    DBG_ASSERT(initialized_);
    DBG_ASSERT(dim==D.Dim());
    DBG_ASSERT(dim==L.NRows());
    DBG_ASSERT(dim==L.NCols());

    // First add the diagonal matrix
    const Number* Dvalues = D.Values();
    for (Index i=0; i<dim; i++) {
      values_[i+i*dim] += Dvalues[i];
    }

    // Now add the strictly-lower triagular matrix L and its transpose
    const Number* Lvalues = L.Values();
    for (Index j=0; j<dim; j++) {
      for (Index i=j+1; i<dim; i++) {
        values_[i+j*dim] += Lvalues[i+j*dim];
      }
    }
    ObjectChanged();
  }

  bool DenseSymMatrix::HasValidNumbersImpl() const
  {
    DBG_ASSERT(initialized_);
    Number sum = 0.;
    const Index dim = Dim();
    for (Index j=0; j<dim; j++) {
      sum += values_[j + j*dim];
      for (Index i=j+1; i<dim; i++) {
        sum += values_[i + j*dim];
      }
    }
    return IsFiniteNumber(sum);
  }

  void DenseSymMatrix::PrintImpl(const Journalist& jnlst,
                                 EJournalLevel level,
                                 EJournalCategory category,
                                 const std::string& name,
                                 Index indent,
                                 const std::string& prefix) const
  {
    jnlst.Printf(level, category, "\n");
    jnlst.PrintfIndented(level, category, indent,
                         "%sDenseSymMatrix \"%s\" of dimension %d (only lower triangular part printed):\n",
                         prefix.c_str(), name.c_str(), Dim());

    if (initialized_) {
      for (Index j=0; j<NCols(); j++) {
        for (Index i=j; i<NRows(); i++) {
          jnlst.PrintfIndented(level, category, indent,
                               "%s%s[%5d,%5d]=%23.16e\n",
                               prefix.c_str(), name.c_str(), i, j, values_[i+NRows()*j]);
        }
      }
    }
    else {
      jnlst.PrintfIndented(level, category, indent,
                           "The matrix has not yet been initialized!\n");
    }
  }

  DenseSymMatrixSpace::DenseSymMatrixSpace(Index nDim)
      :
      SymMatrixSpace(nDim)
  {}

} // namespace Ipopt
