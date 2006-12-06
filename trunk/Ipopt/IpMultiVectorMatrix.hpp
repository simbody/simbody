// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpMultiVectorMatrix.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Andreas Waechter             IBM    2005-12-24

#ifndef __IPMULTIVECTORMATRIX_HPP__
#define __IPMULTIVECTORMATRIX_HPP__

#include "IpUtils.hpp"
#include "IpMatrix.hpp"

namespace Ipopt
{

  /** forward declarations */
  class MultiVectorMatrixSpace;

  /** Class for Matrices with few columns that consists of Vectors.
   *  Those matrices are for example useful in the implementation of
   *  limited memory quasi-Newton methods.
   */
  class MultiVectorMatrix : public Matrix
  {
  public:

    /**@name Constructors / Destructors */
    //@{

    /** Constructor, taking the owner_space.
     */
    MultiVectorMatrix(const MultiVectorMatrixSpace* owner_space);

    /** Destructor */
    ~MultiVectorMatrix();
    //@}

    /** Create a new MultiVectorMatrix from same MatrixSpace */
    SmartPtr<MultiVectorMatrix> MakeNewMultiVectorMatrix() const;

    /** Set a particular Vector at a given column position, replacing
     *  another vector if there has been one.  Depending on whether
     *  the Vector is const or not, it is stored in the const or
     *  non-const internal column. */
    //@{
    void SetVector(Index i, const Vector& vec);
    /* For the non-const version, keep in mind that operations that
     * change this matrix also change the Vector that has been given
     * here. */
    void SetVectorNonConst(Index i, Vector& vec);
    //@}

    /** Get a Vector in a particular column as a const Vector */
    inline SmartPtr<const Vector> GetVector(Index i) const
    {
      return ConstVec(i);
    }

    /** Get a Vector in a particular column as a non-const
     *  Vector. This is fail if the column has currently only a
     *  non-const Vector stored. */
    inline SmartPtr<Vector> GetVectorNonConst(Index i)
    {
      ObjectChanged();
      return Vec(i);
    }

    /** Method for scaling the rows of the matrix, using the
     *  ElementWiseMultiply method for each column vector. */
    void ScaleRows(const Vector& scal_vec);

    /** Method for scaling the columns of the matrix, using the Scal
     *  method for each column vector. */
    void ScaleColumns(const Vector& scal_vec);

    /** Adding another MultiVectorMatrix, using the AddOneVector
     *  methods for the individual column vectors */
    void AddOneMultiVectorMatrix(Number a, const MultiVectorMatrix& mv1,
                                 Number c);

    /** Multiplying a Matrix C (for now assumed to be a
     *  DenseGenMatrix) from the right to a MultiVectorMatrix U and
     *  adding the result to this MultiVectorMatrix V. V = a * U * C +
     *  b * V. */
    void AddRightMultMatrix(Number a, const MultiVectorMatrix& U,
                            const Matrix& C, Number b);

    /** Method for initializing all Vectors with new (uninitialized)
     *  Vectors. */
    void FillWithNewVectors();

    /** Method for adding the low-rank update matrix corresponding to
     *  this matrix to a vector. If V is this MultiVectorMatrix, the
     *  operation is y = beta*y + alpha*V*V^T*x. */
    void LRMultVector(Number alpha, const Vector &x,
                      Number beta, Vector &y) const;

    /** Vector space for the columns */
    SmartPtr<const VectorSpace> ColVectorSpace() const;

    /** Return the MultiVectorMatrixSpace */
    SmartPtr<const MultiVectorMatrixSpace> MultiVectorMatrixOwnerSpace() const;

  protected:
    /**@name Overloaded methods from Matrix base class */
    //@{
    virtual void MultVectorImpl(Number alpha, const Vector &x, Number beta,
                                Vector &y) const;

    virtual void TransMultVectorImpl(Number alpha, const Vector& x,
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
    MultiVectorMatrix();

    /** Copy Constructor */
    MultiVectorMatrix(const MultiVectorMatrix&);

    /** Overloaded Equals Operator */
    void operator=(const MultiVectorMatrix&);
    //@}

    const MultiVectorMatrixSpace* owner_space_;

    /** space for storing the const Vector's */
    std::vector<SmartPtr<const Vector> > const_vecs_;

    /** space for storing the non-const Vector's */
    std::vector<SmartPtr<Vector> > non_const_vecs_;

    /** Method for accessing the internal Vectors internally */
    //@{
    inline const Vector* ConstVec(Index i) const
    {
      DBG_ASSERT(i < NCols());
      DBG_ASSERT(IsValid(const_vecs_[i]) || IsValid(non_const_vecs_[i]));
      if (IsValid(non_const_vecs_[i])) {
        return GetRawPtr(non_const_vecs_[i]);
      }
      else {
        return GetRawPtr(const_vecs_[i]);
      }
    }

    inline Vector* Vec(Index i)
    {
      DBG_ASSERT(i < NCols());
      DBG_ASSERT(IsValid(non_const_vecs_[i]));
      return GetRawPtr(non_const_vecs_[i]);
    }
    //@}
  };

  /** This is the matrix space for MultiVectorMatrix.
   */
  class MultiVectorMatrixSpace : public MatrixSpace
  {
  public:
    /** @name Constructors / Destructors */
    //@{
    /** Constructor, given the number of columns (i.e., Vectors to be
     *  stored) and given the VectorSpace for the Vectors.
     */
    MultiVectorMatrixSpace(Index ncols,
                           const VectorSpace& vec_space);

    /** Destructor */
    ~MultiVectorMatrixSpace()
    {}
    //@}

    /** Method for creating a new matrix of this specific type. */
    MultiVectorMatrix* MakeNewMultiVectorMatrix() const
    {
      return new MultiVectorMatrix(this);
    }

    /** Overloaded MakeNew method for the MatrixSpace base class.
     */
    virtual Matrix* MakeNew() const
    {
      return MakeNewMultiVectorMatrix();
    }

    /** Accessor method for the VectorSpace for the columns */
    SmartPtr<const VectorSpace> ColVectorSpace() const
    {
      return vec_space_;
    }

  private:
    SmartPtr<const VectorSpace> vec_space_;

  };

  inline
  MultiVectorMatrix::~MultiVectorMatrix()
  {}

  inline
  SmartPtr<MultiVectorMatrix> MultiVectorMatrix::MakeNewMultiVectorMatrix() const
  {
    return owner_space_->MakeNewMultiVectorMatrix();
  }

  inline
  SmartPtr<const VectorSpace> MultiVectorMatrix::ColVectorSpace() const
  {
    return owner_space_->ColVectorSpace();
  }

  inline
  SmartPtr<const MultiVectorMatrixSpace>
  MultiVectorMatrix::MultiVectorMatrixOwnerSpace() const
  {
    return owner_space_;
  }

} // namespace Ipopt
#endif
