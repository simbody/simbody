// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpLowRankUpdateSymMatrix.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Andreas Waechter            IBM    2005-12-25

#ifndef __IPLOWRANKUPDATESYMMATRIX_HPP__
#define __IPLOWRANKUPDATESYMMATRIX_HPP__

#include "IpUtils.hpp"
#include "IpSymMatrix.hpp"
#include "IpMultiVectorMatrix.hpp"

namespace Ipopt
{

  /* forward declarations */
  class LowRankUpdateSymMatrixSpace;

  /** Class for symmetric matrices, represented as low-rank updates.
   *  The matrix M is represented as M = P_LR(D + V V^T - U U^T)P_LR^T
   *  (if reduced_diag is true), or M = D + P_LR(V V^T - U U^T)P_LR^T
   *  (if reduced_diag is false).  D is a diagonal matrix, and V and U
   *  are MultiVectorMatrices, and P_LR is an ExpansionMatrix.  The
   *  vectors in the low-rank update (before expansion) live in the
   *  LowRankVectorSpace.  If P_LR is NULL, P_LR is assumed to be the
   *  identity matrix.  If V or U is NULL, it is assume to be a matrix
   *  of zero columns. */
  class LowRankUpdateSymMatrix : public SymMatrix
  {
  public:

    /**@name Constructors / Destructors */
    //@{

    /** Constructor, given the corresponding matrix space. */
    LowRankUpdateSymMatrix(const LowRankUpdateSymMatrixSpace* owner_space);

    /** Destructor */
    ~LowRankUpdateSymMatrix();
    //@}

    /** Method for setting the diagonal elements (as a Vector). */
    void SetDiag(const Vector& D)
    {
      D_ = &D;
      ObjectChanged();
    }

    /** Method for getting the diagonal elements. */
    SmartPtr<const Vector> GetDiag() const
    {
      return D_;
    }

    /** Method for setting the positive low-rank update part. */
    void SetV(const MultiVectorMatrix& V)
    {
      V_ = &V;
      ObjectChanged();
    }

    /** Method for getting the positive low-rank update part. */
    SmartPtr<const MultiVectorMatrix> GetV() const
    {
      return V_;
    }

    /** Method for setting the negative low-rank update part. */
    void SetU(const MultiVectorMatrix& U)
    {
      U_ = &U;
      ObjectChanged();
    }

    /** Method for getting the negative low-rank update part. */
    SmartPtr<const MultiVectorMatrix> GetU() const
    {
      return U_;
    }

    /** Return the expansion matrix to lift the low-rank update to the
     *  higher-dimensional space. */
    SmartPtr<const Matrix> P_LowRank() const;

    /** Return the vector space in with the low-rank update vectors
     *  live. */
    SmartPtr<const VectorSpace> LowRankVectorSpace() const;

    /** Flag indicating whether the diagonal term lives in the smaller
     *  space (from P_LowRank) or in the full space. */
    bool ReducedDiag() const;

  protected:
    /**@name Methods overloaded from matrix */
    //@{
    virtual void MultVectorImpl(Number alpha, const Vector& x,
                                Number beta, Vector& y) const;

    /** Method for determining if all stored numbers are valid (i.e.,
     *  no Inf or Nan). */
    virtual bool HasValidNumbersImpl() const;

    virtual void PrintImpl(const Journalist& jnlst,
                           EJournalLevel level,
                           EJournalCategory category,
                           const std::string& name,
                           Index indent,
                           const std::string& prefix) const;
    //@}

  private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Default Constructor */
    LowRankUpdateSymMatrix();

    /** Copy Constructor */
    LowRankUpdateSymMatrix(const LowRankUpdateSymMatrix&);

    /** Overloaded Equals Operator */
    void operator=(const LowRankUpdateSymMatrix&);
    //@}

    /** corresponding matrix space */
    SmartPtr<const LowRankUpdateSymMatrixSpace> owner_space_;

    /** Vector storing the diagonal matrix D. */
    SmartPtr<const Vector> D_;

    /** Vector storing the positive low-rank update. */
    SmartPtr<const MultiVectorMatrix> V_;

    /** Vector storing the negative low-rank update. */
    SmartPtr<const MultiVectorMatrix> U_;
  };

  /** This is the matrix space for LowRankUpdateSymMatrix. */
  class LowRankUpdateSymMatrixSpace : public SymMatrixSpace
  {
  public:
    /** @name Constructors / Destructors */
    //@{
    /** Constructor, given the dimension of the matrix. */
    LowRankUpdateSymMatrixSpace(Index dim,
                                SmartPtr<const Matrix> P_LowRank,
                                SmartPtr<const VectorSpace> LowRankVectorSpace,
                                bool reduced_diag)
        :
        SymMatrixSpace(dim),
        P_LowRank_(P_LowRank),
        lowrank_vector_space_(LowRankVectorSpace),
        reduced_diag_(reduced_diag)
    {
      DBG_ASSERT(IsValid(lowrank_vector_space_));
    }

    /** Destructor */
    virtual ~LowRankUpdateSymMatrixSpace()
    {}
    //@}

    /** Overloaded MakeNew method for the SymMatrixSpace base class.
     */
    virtual SymMatrix* MakeNewSymMatrix() const
    {
      return MakeNewLowRankUpdateSymMatrix();
    }

    /** Method for creating a new matrix of this specific type. */
    LowRankUpdateSymMatrix* MakeNewLowRankUpdateSymMatrix() const
    {
      return new LowRankUpdateSymMatrix(this);
    }

    SmartPtr<const Matrix> P_LowRank() const
    {
      return P_LowRank_;
    }

    SmartPtr<const VectorSpace> LowRankVectorSpace() const
    {
      return lowrank_vector_space_;
    }

    bool ReducedDiag() const
    {
      return reduced_diag_;
    }

  private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Default Constructor */
    LowRankUpdateSymMatrixSpace();

    /** Copy Constructor */
    LowRankUpdateSymMatrixSpace(const LowRankUpdateSymMatrixSpace&);

    /** Overloaded Equals Operator */
    void operator=(const LowRankUpdateSymMatrixSpace&);
    //@}

    /** Expansion matrix to lift the low-rank approximation into a
     *  possibly higher-dimensional space.  If it is NULL, it is
     *  assume that no lift is performed. */
    SmartPtr<const Matrix> P_LowRank_;

    /** Vector space for the space in which the low-rank approximation
     *  lives. */
    SmartPtr<const VectorSpace> lowrank_vector_space_;

    /** Flag indicating whether the diagonal matrix is nonzero only in
     *  the space of V or in the full space. */
    bool reduced_diag_;
  };

  inline
  SmartPtr<const Matrix> LowRankUpdateSymMatrix::P_LowRank() const
  {
    return owner_space_->P_LowRank();
  }

  inline
  SmartPtr<const VectorSpace> LowRankUpdateSymMatrix::LowRankVectorSpace() const
  {
    return owner_space_->LowRankVectorSpace();
  }

  inline
  bool LowRankUpdateSymMatrix::ReducedDiag() const
  {
    return owner_space_->ReducedDiag();
  }

} // namespace Ipopt
#endif
