// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpSymMatrix.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPSYMMATRIX_HPP__
#define __IPSYMMATRIX_HPP__

#include "IpUtils.hpp"
#include "IpMatrix.hpp"

namespace Ipopt
{

  /* forward declarations */
  class SymMatrixSpace;

  /** This is the base class for all derived symmetric matrix types.
   */
  class SymMatrix : public Matrix
  {
  public:
    /** @name Constructor/Destructor */
    //@{
    /** Constructor, taking the owner_space.
     */
    SymMatrix(const SymMatrixSpace* owner_space);

    /** Destructor */
    virtual ~SymMatrix();
    //@}

    /** @name Information about the size of the matrix */
    //@{
    /** Dimension of the matrix (number of rows and columns) */
    Index Dim() const;
    //@}

    SmartPtr<const SymMatrixSpace> OwnerSymMatrixSpace() const;

  protected:
    /** @name Overloaded methods from Matrix.  Since the matrix is
     *  symmetric, it is only necessary to implement the
     *  MultVectorImpl method in a class that inherits from this base
     *  class.  If the TransMultVectorImpl is called, this base class
     *  automatically calls MultVectorImpl instead.
     */
    //@{
    virtual void TransMultVectorImpl(Number alpha, const Vector& x, Number beta,
                                     Vector& y) const;
    //@}

  private:
    /** Copy of the owner space ptr as a SymMatrixSpace instead
     *  of a MatrixSpace 
     */
    const SymMatrixSpace* owner_space_;
  };


  /** SymMatrixSpace base class, corresponding to the SymMatrix base
   *  class. */
  class SymMatrixSpace : public MatrixSpace
  {
  public:
    /** @name Constructors/Destructors */
    //@{
    /** Constructor, given the dimension (identical to the number of
     *  rows and columns).
     */
    SymMatrixSpace(Index dim);

    /** Destructor */
    virtual ~SymMatrixSpace();
    //@}

    /** Pure virtual method for creating a new matrix of this specific
     *  type. */
    virtual SymMatrix* MakeNewSymMatrix() const=0;

    /** Overloaded MakeNew method for the MatrixSpace base class.
     */
    virtual Matrix* MakeNew() const;

    /** Accessor method for the dimension of the matrices in this
     *  matrix space.
     */
    Index Dim() const;

  private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** default constructor */
    SymMatrixSpace();

    /* Copy constructor */
    SymMatrixSpace(const SymMatrixSpace&);

    /** Overloaded Equals Operator */
    SymMatrixSpace& operator=(const SymMatrixSpace&);
    //@}

  };

  /* inline methods */
  inline
  Index SymMatrixSpace::Dim() const
  {
    DBG_ASSERT(NRows() == NCols());
    return NRows();
  }

  inline
  Index SymMatrix::Dim() const
  {
    return owner_space_->Dim();
  }

} // namespace Ipopt

#endif
