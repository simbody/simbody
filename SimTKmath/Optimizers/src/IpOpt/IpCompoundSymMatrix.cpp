// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpCompoundSymMatrix.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpoptConfig.h"
#include "IpCompoundSymMatrix.hpp"
#include "IpCompoundVector.hpp"

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

  CompoundSymMatrix::CompoundSymMatrix(const CompoundSymMatrixSpace* owner_space)
      :
      SymMatrix(owner_space),
      owner_space_(owner_space),
      matrices_valid_(false)
  {
    for (Index irow=0; irow<NComps_Dim(); irow++) {
      std::vector< SmartPtr<Matrix> > row(irow+1);
      std::vector< SmartPtr<const Matrix> > const_row(irow+1);
      comps_.push_back(row);
      const_comps_.push_back(const_row);
    }
  }

  CompoundSymMatrix::~CompoundSymMatrix()
  {}

  void CompoundSymMatrix::SetComp(Index irow, Index jcol,
                                  const Matrix& matrix)
  {
    DBG_ASSERT(!matrices_valid_);
    DBG_ASSERT(irow<NComps_Dim());
    DBG_ASSERT(jcol<=irow);
    // Matrices on the diagonal must be symmetric
    DBG_ASSERT(irow!=jcol || dynamic_cast<const SymMatrix*>(&matrix));
    DBG_ASSERT(owner_space_->GetCompSpace(irow, jcol)->IsMatrixFromSpace(matrix));

    comps_[irow][jcol] = NULL;
    const_comps_[irow][jcol] = &matrix;
    ObjectChanged();
  }

  void CompoundSymMatrix::SetCompNonConst(Index irow, Index jcol,
                                          Matrix& matrix)
  {
    DBG_ASSERT(!matrices_valid_);
    DBG_ASSERT(irow < NComps_Dim());
    DBG_ASSERT(jcol <= irow);
    // Matrices on the diagonal must be symmetric
    DBG_ASSERT( irow != jcol || dynamic_cast<SymMatrix*>(&matrix));
    DBG_ASSERT(owner_space_->GetCompSpace(irow, jcol)->IsMatrixFromSpace(matrix));

    const_comps_[irow][jcol] = NULL;
    comps_[irow][jcol] = &matrix;
    ObjectChanged();
  }

  Index CompoundSymMatrix::NComps_Dim() const
  {
    return owner_space_->NComps_Dim();
  }

  void CompoundSymMatrix::MultVectorImpl(Number alpha, const Vector &x,
                                         Number beta, Vector &y) const
  {
    if (!matrices_valid_) {
      matrices_valid_ = MatricesValid();
    }
    DBG_ASSERT(matrices_valid_);

    // The vectors are assumed to be compound Vectors as well
    const CompoundVector* comp_x = dynamic_cast<const CompoundVector*>(&x);
    CompoundVector* comp_y = dynamic_cast<CompoundVector*>(&y);

    //  A few sanity checks
    if (comp_x) {
      DBG_ASSERT(NComps_Dim()==comp_x->NComps());
    }
    else {
      DBG_ASSERT(NComps_Dim() == 1);
    }
    if (comp_y) {
      DBG_ASSERT(NComps_Dim()==comp_y->NComps());
    }
    else {
      DBG_ASSERT(NComps_Dim() == 1);
    }

    // Take care of the y part of the addition
    if( beta!=0.0 ) {
      y.Scal(beta);
    }
    else {
      y.Set(0.0);  // In case y hasn't been initialized yet
    }

    for (Index irow=0; irow<NComps_Dim(); irow++) {
      SmartPtr<Vector> y_i;
      if (comp_y) {
        y_i = comp_y->GetCompNonConst(irow);
      }
      else {
        y_i = &y;
      }
      DBG_ASSERT(IsValid(y_i));

      for (Index jcol=0; jcol<=irow; jcol++) {
        SmartPtr<const Vector> x_j;
        if (comp_x) {
          x_j = comp_x->GetComp(irow);
        }
        else {
          x_j = &x;
        }
        DBG_ASSERT(IsValid(x_j));

        if (ConstComp(irow,jcol)) {
          ConstComp(irow,jcol)->MultVector(alpha, *comp_x->GetComp(jcol),
                                           1., *comp_y->GetCompNonConst(irow));
        }
      }

      for (Index jcol = irow+1; jcol < NComps_Dim(); jcol++) {
        if (ConstComp(jcol,irow)) {
          ConstComp(jcol,irow)->TransMultVector(alpha, *comp_x->GetComp(jcol),
                                                1., *comp_y->GetCompNonConst(irow));
        }
      }
    }
  }

  bool CompoundSymMatrix::HasValidNumbersImpl() const
  {
    if (!matrices_valid_) {
      matrices_valid_ = MatricesValid();
    }
    DBG_ASSERT(matrices_valid_);

    for (Index irow=0; irow<NComps_Dim(); irow++) {
      for (Index jcol=0; jcol<=irow; jcol++) {
        if (ConstComp(irow,jcol)) {
          if (!ConstComp(irow,jcol)->HasValidNumbers()) {
            return false;
          }
        }
      }
    }
    return true;
  }

  void CompoundSymMatrix::PrintImpl(const Journalist& jnlst,
                                    EJournalLevel level,
                                    EJournalCategory category,
                                    const std::string& name,
                                    Index indent,
                                    const std::string& prefix) const
  {
    jnlst.Printf(level, category, "\n");
    jnlst.PrintfIndented(level, category, indent,
                         "%sCompoundSymMatrix \"%s\" with %d rows and columns components:\n",
                         prefix.c_str(), name.c_str(), NComps_Dim());
    for (Index irow = 0; irow < NComps_Dim(); irow++ ) {
      for (Index jcol = 0; jcol <= irow; jcol++ ) {
        jnlst.PrintfIndented(level, category, indent,
                             "%sComponent for row %d and column %d:\n",
                             prefix.c_str(), irow, jcol);
        if (ConstComp(irow,jcol)) {
          DBG_ASSERT(name.size()<200);
          char buffer[256];
          sprintf(buffer, "%s[%d][%d]", name.c_str(), irow, jcol);
          std::string term_name = buffer;
          ConstComp(irow,jcol)->Print(&jnlst, level, category, term_name,
                                      indent+1, prefix);
        }
        else {
          jnlst.PrintfIndented(level, category, indent,
                               "%sThis component has not been set.\n",
                               prefix.c_str());
        }
      }
    }
  }

  bool CompoundSymMatrix::MatricesValid() const
  {
    // Check to make sure we have matrices everywhere the space has matrices
    // We already check that the matrix agrees with the block space
    // in the SetComp methods
    bool retValue = true;
    for (Index i=0; i<NComps_Dim(); i++) {
      for (Index j=0; j<=i; j++) {
        if ( (!ConstComp(i, j) && IsValid(owner_space_->GetCompSpace(i,j)))
             || (ConstComp(i, j) && IsNull(owner_space_->GetCompSpace(i,j))) ) {
          retValue = false;
          break;
        }
      }
    }

    return retValue;
  }

  CompoundSymMatrixSpace::CompoundSymMatrixSpace(Index ncomp_spaces, Index total_dim)
      :
      SymMatrixSpace(total_dim),
      ncomp_spaces_(ncomp_spaces),
      block_dim_(ncomp_spaces, -1),
      dimensions_set_(false)
  {
    for (Index irow=0; irow<ncomp_spaces_; irow++) {
      std::vector<SmartPtr<const MatrixSpace> > row(irow+1);
      std::vector< bool > allocate_row(irow+1, false);
      comp_spaces_.push_back(row);
      allocate_block_.push_back(allocate_row);
    }
  }

  void CompoundSymMatrixSpace::SetBlockDim(Index irow_jcol, Index dim)
  {
    DBG_ASSERT(!dimensions_set_ && "for now, if dimensions are set, they cannot be changed");
    DBG_ASSERT(block_dim_[irow_jcol] == -1 && "This dimension has already been set - sanity check");
    DBG_ASSERT(irow_jcol < ncomp_spaces_);
    block_dim_[irow_jcol] = dim;
  }

  Index CompoundSymMatrixSpace::GetBlockDim(Index irow_jcol) const
  {
    DBG_ASSERT(dimensions_set_ && "Cannot get block dimensions before all dimensions are set.");
    DBG_ASSERT(irow_jcol < ncomp_spaces_);
    return block_dim_[irow_jcol];
  }

  void CompoundSymMatrixSpace::SetCompSpace(Index irow, Index jcol,
      const MatrixSpace& mat_space,
      bool auto_allocate /*=false*/)
  {
    if (!dimensions_set_) {
      dimensions_set_ = DimensionsSet();
    }
    DBG_ASSERT(dimensions_set_);
    DBG_ASSERT(irow<ncomp_spaces_);
    DBG_ASSERT(jcol<=irow);
    DBG_ASSERT(IsNull(comp_spaces_[irow][jcol]));
    DBG_ASSERT(irow!=jcol || dynamic_cast<const SymMatrixSpace*> (&mat_space));
    DBG_ASSERT(block_dim_[jcol] != -1 && block_dim_[jcol] == mat_space.NCols());
    DBG_ASSERT(block_dim_[irow] != -1 && block_dim_[irow] == mat_space.NRows());

    comp_spaces_[irow][jcol] = &mat_space;
    allocate_block_[irow][jcol] = auto_allocate;
  }

  CompoundSymMatrix* CompoundSymMatrixSpace::MakeNewCompoundSymMatrix() const
  {
    if (!dimensions_set_) {
      dimensions_set_ = DimensionsSet();
    }
    DBG_ASSERT(dimensions_set_);

    CompoundSymMatrix* mat = new CompoundSymMatrix(this);
    for(Index i=0; i<NComps_Dim(); i++) {
      for (Index j=0; j<=i; j++) {
        if (allocate_block_[i][j]) {
          mat->SetCompNonConst(i, j, *GetCompSpace(i, j)->MakeNew());
        }
      }
    }

    return mat;
  }

  bool CompoundSymMatrixSpace::DimensionsSet() const
  {
    Index total_dim= 0;
    bool valid = true;
    for (Index i=0; i<ncomp_spaces_; i++) {
      if (block_dim_[i] == -1) {
        valid = false;
        break;
      }
      total_dim += block_dim_[i];
    }

    if (valid) {
      DBG_ASSERT(total_dim == Dim());
    }

    return valid;
  }

} // namespace Ipopt
