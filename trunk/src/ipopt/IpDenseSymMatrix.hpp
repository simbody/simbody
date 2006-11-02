// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpDenseSymMatrix.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Andreas Waechter             IBM    2005-12-25

#ifndef __IPDENSESYMMATRIX_HPP__
#define __IPDENSESYMMATRIX_HPP__

#include "IpUtils.hpp"
#include "IpSymMatrix.hpp"
#include "IpMultiVectorMatrix.hpp"
#include "IpDenseVector.hpp"

namespace Ipopt
{

  /** forward declarations */
  class DenseSymMatrixSpace;

  /** forward declaration so that this include file can be included
   *  from DenseGenMatrix */
  class DenseGenMatrix;

  /** Class for dense symetrix matrices.  Matrix elements are stored
   *  in one array in "Fortran" format, using BLAS "lower triangular"
   *  storage (not packed).
   */
  class DenseSymMatrix : public SymMatrix
  {
  public:

    /**@name Constructors / Destructors */
    //@{

    /** Constructor, taking the owner_space.
     */
    DenseSymMatrix(const DenseSymMatrixSpace* owner_space);

    /** Destructor */
    ~DenseSymMatrix();
    //@}

    /** Create a new DenseSymMatrix from same MatrixSpace */
    SmartPtr<DenseSymMatrix> MakeNewDenseSymMatrix() const;

    /** Retrieve the array for storing the matrix elements.  This is
     *  the non-const version, and it is assume that afterwards the
     *  calling method will set all matrix elements.  The matrix
     *  elements are stored one column after each other. */
    Number* Values()
    {
      ObjectChanged();
      initialized_ = true;
      return values_;
    }

    /** Retrieve the array that stores the matrix elements.  This is
     *  the const version, i.e., read-only.  The matrix elements are
     *  stored one column after each other. */
    const Number* Values() const
    {
      DBG_ASSERT(initialized_);
      return values_;
    }

    /** Set this matrix to be a multiple of the identity matrix. */
    void FillIdentity(Number factor=1.);

    /** Method for adding another matrix to this one.  If B is this
    *  matrix, it becomes B = alpha * A + beta * B after this call. */
    void AddMatrix(Number alpha, const DenseSymMatrix& A, Number beta);

    /** Method for adding a high-rank update to this matrix.  It
     *  computes M = alpha*op(V) op(V)^T + beta*M, where V is a
     *  DenseGenMatrix, where op(V) is V^T trans is true.  */
    void HighRankUpdate(bool trans, Number alpha, const DenseGenMatrix& V,
                        Number beta);

    /** Method for adding a high-rank update to this matrix.  It
     *  computes M = alpha*V1^T V2 + beta*M, where V1 and V2 are
     *  MultiVectorMatrices, so that V1^T V2 is symmetric.  */
    void HighRankUpdateTranspose(Number alpha,
                                 const MultiVectorMatrix& V1,
                                 const MultiVectorMatrix& V2,
                                 Number beta);

    /** Method for doing a specialized Add operation, required in the
     *  limited memory SR1 update. if M is this matrix, it computes M
     *  = M + D + L + L^T, where D is a diagonal matrix (given as a
     *  DenseVector), and L is a matrix that is assumed to be strictly
     *  lower triangular. */
    void SpecialAddForLMSR1(const DenseVector& D, const DenseGenMatrix& L);

  protected:
    /**@name Overloaded methods from Matrix base class*/
    //@{
    virtual void MultVectorImpl(Number alpha, const Vector &x, Number beta,
                                Vector &y) const;

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
    DenseSymMatrix();

    /** Copy Constructor */
    DenseSymMatrix(const DenseSymMatrix&);

    /** Overloaded Equals Operator */
    void operator=(const DenseSymMatrix&);
    //@}

    const DenseSymMatrixSpace* owner_space_;

    /** Array for storing the matrix elements (one columns after each
     *  other) */
    Number* values_;

    /** Flag indicating whether the values_ array has been initialized */
    bool initialized_;
  };

  /** This is the matrix space for DenseSymMatrix.
   */
  class DenseSymMatrixSpace : public SymMatrixSpace
  {
  public:
    /** @name Constructors / Destructors */
    //@{
    /** Constructor for matrix space for DenseSymMatrices.  Takes in
     *  dimension of the matrices.
     */
    DenseSymMatrixSpace(Index nDim);

    /** Destructor */
    ~DenseSymMatrixSpace()
    {}
    //@}

    /** Method for creating a new matrix of this specific type. */
    DenseSymMatrix* MakeNewDenseSymMatrix() const
    {
      return new DenseSymMatrix(this);
    }

    /** Overloaded MakeNew method for the MatrixSpace base class.
     */
    virtual SymMatrix* MakeNewSymMatrix() const
    {
      return MakeNewDenseSymMatrix();
    }

  };

  inline
  SmartPtr<DenseSymMatrix> DenseSymMatrix::MakeNewDenseSymMatrix() const
  {
    return owner_space_->MakeNewDenseSymMatrix();
  }

} // namespace Ipopt
#endif
