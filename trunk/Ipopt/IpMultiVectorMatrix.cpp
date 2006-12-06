// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpMultiVectorMatrix.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Andreas Waechter             IBM    2005-12-24

#include "IpMultiVectorMatrix.hpp"
#include "IpDenseVector.hpp"
#include "IpDenseGenMatrix.hpp"

namespace Ipopt
{

#ifdef IP_DEBUG
  static const Index dbg_verbosity = 0;
#endif

  MultiVectorMatrix::MultiVectorMatrix(const MultiVectorMatrixSpace* owner_space)
      :
      Matrix(owner_space),
      owner_space_(owner_space),
      const_vecs_(owner_space->NCols()),
      non_const_vecs_(owner_space->NCols())
  {}

  void MultiVectorMatrix::SetVector(Index i, const Vector& vec)
  {
    DBG_ASSERT(i<NCols());
    non_const_vecs_[i] = NULL;
    const_vecs_[i] = &vec;
    ObjectChanged();
  }

  void MultiVectorMatrix::SetVectorNonConst(Index i, Vector& vec)
  {
    DBG_ASSERT(i<NCols());
    const_vecs_[i] = NULL;
    non_const_vecs_[i] = &vec;
    ObjectChanged();
  }

  void MultiVectorMatrix::MultVectorImpl(Number alpha, const Vector &x,
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
    DBG_ASSERT(dense_x);

    // We simply add all the Vectors one after the other
    if (dense_x->IsHomogeneous()) {
      Number val = dense_x->Scalar();
      for (Index i=0; i<NCols(); i++) {
        y.AddOneVector(alpha*val, *ConstVec(i), 1.);
      }
    }
    else {
      const Number* values = dense_x->Values();
      for (Index i=0; i<NCols(); i++) {
        y.AddOneVector(alpha*values[i], *ConstVec(i), 1.);
      }
    }
  }

  void MultiVectorMatrix::TransMultVectorImpl(Number alpha, const Vector &x,
      Number beta, Vector &y) const
  {
    //  A few sanity checks
    DBG_ASSERT(NCols()==y.Dim());
    DBG_ASSERT(NRows()==x.Dim());

    // See if we can understand the data
    DenseVector* dense_y = dynamic_cast<DenseVector*>(&y);
    DBG_ASSERT(dense_y);

    // Use the individual dot products to get the matrix (transpose)
    // vector product
    Number *yvals=dense_y->Values();
    if( beta!=0.0 ) {
      for (Index i=0; i<NCols(); i++) {
        yvals[i] = alpha*ConstVec(i)->Dot(x) + beta*yvals[i];
      }
    }
    else {
      for (Index i=0; i<NCols(); i++) {
        yvals[i] = alpha*ConstVec(i)->Dot(x);
      }
    }
  }

  void MultiVectorMatrix::LRMultVector(Number alpha, const Vector &x,
                                       Number beta, Vector &y) const
  {
    DBG_START_METH("MultiVectorMatrix::LRMultVector(",
                   dbg_verbosity);

    DBG_ASSERT(NRows()==x.Dim());
    DBG_ASSERT(NRows()==y.Dim());

    DBG_PRINT((1, "alpha = %e beta = %e\n", alpha, beta));
    DBG_PRINT_VECTOR(2, "x", x);

    if( beta!=0.0 ) {
      y.Scal(beta);
    }
    else {
      y.Set(0.0);
    }

    DBG_PRINT_VECTOR(2, "beta*y", y);
    for (Index i=0; i<NCols(); i++) {
      DBG_PRINT_VECTOR(2, "ConstVec(i)", *ConstVec(i));
      y.AddOneVector(alpha*ConstVec(i)->Dot(x), *ConstVec(i), 1.);
      DBG_PRINT_VECTOR(2, "y mid", y);
    }
  }

  void MultiVectorMatrix::FillWithNewVectors()
  {
    SmartPtr<const VectorSpace> vec_space = owner_space_->ColVectorSpace();
    for (Index i=0; i<NCols(); i++) {
      non_const_vecs_[i] = vec_space->MakeNew();
      const_vecs_[i] = NULL;
    }
    ObjectChanged();
  }

  void MultiVectorMatrix::ScaleRows(const Vector& scal_vec)
  {
    // Santiy checks
    DBG_ASSERT(scal_vec.Dim() == NRows());

    for (Index i=0; i<NCols(); i++) {
      Vec(i)->ElementWiseMultiply(scal_vec);
    }
    ObjectChanged();
  }

  void MultiVectorMatrix::ScaleColumns(const Vector& scal_vec)
  {
    // Santiy checks
    DBG_ASSERT(scal_vec.Dim() == NCols());

    // See if we can understand the data
    const DenseVector* dense_scal_vec =
      dynamic_cast<const DenseVector*>(&scal_vec);
    DBG_ASSERT(dense_scal_vec);

    if (dense_scal_vec->IsHomogeneous()) {
      Number val = dense_scal_vec->Scalar();
      for (Index i=0; i<NCols(); i++) {
        Vec(i)->Scal(val);
      }
    }
    else {
      const Number* values = dense_scal_vec->Values();
      for (Index i=0; i<NCols(); i++) {
        Vec(i)->Scal(values[i]);
      }
    }
    ObjectChanged();
  }

  void
  MultiVectorMatrix::AddOneMultiVectorMatrix(Number a,
      const MultiVectorMatrix& mv1,
      Number c)
  {
    DBG_ASSERT(NRows()==mv1.NRows());
    DBG_ASSERT(NCols()==mv1.NCols());

    if (c==0.) {
      FillWithNewVectors();
    }

    for (Index i=0; i<NCols(); i++) {
      Vec(i)->AddOneVector(a, *mv1.GetVector(i), c);
    }
    ObjectChanged();
  }

  void
  MultiVectorMatrix::AddRightMultMatrix(Number a,
                                        const MultiVectorMatrix& U,
                                        const Matrix& C,
                                        Number b)
  {
    DBG_ASSERT(NRows()==U.NRows());
    DBG_ASSERT(U.NCols()==C.NRows());
    DBG_ASSERT(NCols()==C.NCols());

    if (b==0.) {
      FillWithNewVectors();
    }

    // ToDo: For now, we simply use MatrixVector multiplications, but
    // we might be more efficient (at least in the non-parallel case)
    // if we used Level 3 Blas
    SmartPtr<const DenseVectorSpace> mydspace = new DenseVectorSpace(C.NRows());
    SmartPtr<DenseVector> mydvec = mydspace->MakeNewDenseVector();

    const DenseGenMatrix* dgm_C = dynamic_cast<const DenseGenMatrix*>(&C);
    DBG_ASSERT(dgm_C);
    for (Index i=0; i<NCols(); i++) {
      const Number* CValues = dgm_C->Values();
      Number* myvalues = mydvec->Values();
      for (Index j=0; j<U.NCols(); j++) {
        myvalues[j] = CValues[i*C.NRows() + j];
      }
      U.MultVector(a, *mydvec, b, *Vec(i));
    }
    ObjectChanged();
  }

  bool MultiVectorMatrix::HasValidNumbersImpl() const
  {
    for (Index i=0; i<NCols(); i++) {
      if (!ConstVec(i)->HasValidNumbers()) {
        return false;
      }
    }
    return true;
  }

  void MultiVectorMatrix::PrintImpl(const Journalist& jnlst,
                                    EJournalLevel level,
                                    EJournalCategory category,
                                    const std::string& name,
                                    Index indent,
                                    const std::string& prefix) const
  {
    jnlst.Printf(level, category, "\n");
    jnlst.PrintfIndented(level, category, indent,
                         "%sMultiVectorMatrix \"%s\" with %d columns:\n",
                         prefix.c_str(), name.c_str(), NCols());

    for (Index i=0; i<NCols(); i++) {
      if (ConstVec(i)) {
        DBG_ASSERT(name.size()<200);
        char buffer[256];
        sprintf(buffer, "%s[%2d]", name.c_str(), i);
        std::string term_name = buffer;
        ConstVec(i)->Print(&jnlst, level, category, term_name,
                           indent+1, prefix);
      }
      else {
        jnlst.PrintfIndented(level, category, indent,
                             "%sVector in column %d is not yet set!\n",
                             prefix.c_str(), i);
      }
    }
  }

  MultiVectorMatrixSpace::MultiVectorMatrixSpace(Index ncols,
      const VectorSpace& vec_space)
      :
      MatrixSpace(vec_space.Dim(), ncols),
      vec_space_(&vec_space)
  {}

} // namespace Ipopt
