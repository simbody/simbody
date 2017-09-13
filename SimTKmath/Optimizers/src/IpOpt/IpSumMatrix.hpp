// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpSumMatrix.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPSUMMATRIX_HPP__
#define __IPSUMMATRIX_HPP__

#include "IpUtils.hpp"
#include "IpMatrix.hpp"

namespace SimTKIpopt
{

  /* forward declarations */
  class SumMatrixSpace;

  /** Class for Matrices which are sum of matrices.
   *  For each term in the we store the matrix and a factor.
   */
  class SumMatrix : public Matrix
  {
  public:

    /**@name Constructors / Destructors */
    //@{
    /** Constructor, taking the owner_space.
     */
    SumMatrix(const SumMatrixSpace* owner_space);

    /** Destructor */
    virtual ~SumMatrix();
    //@}

    /** Method for setting term iterm for the sum. */
    void SetTerm(Index iterm, Number factor, const Matrix& matrix);

    /** Method for getting term iterm for the sum.  Note that counting
     *  of terms starts at 0. */
    void GetTerm(Index iterm, Number& factor, SmartPtr<const Matrix>& matrix) const;

    /** Return the number of terms */
    Index NTerms() const;

  protected:
    /**@name Methods overloaded from matrix */
    //@{
    virtual void MultVectorImpl(Number alpha, const Vector& x,
                                Number beta, Vector& y) const override;

    virtual void TransMultVectorImpl(Number alpha, const Vector& x,
                                     Number beta, Vector& y) const override;

    /** Method for determining if all stored numbers are valid (i.e.,
     *  no Inf or Nan). */
    virtual bool HasValidNumbersImpl() const override;

    virtual void PrintImpl(const Journalist& jnlst,
                           EJournalLevel level,
                           EJournalCategory category,
                           const std::string& name,
                           Index indent,
                           const std::string& prefix) const override;
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
    SumMatrix();

    /** Copy Constructor */
    SumMatrix(const SumMatrix&);

    /** Overloaded Equals Operator */
    void operator=(const SumMatrix&);
    //@}

    /** std::vector storing the factors for each term. */
    std::vector<Number> factors_;

    /** std::vector storing the matrices for each term. */
    std::vector<SmartPtr<const Matrix> > matrices_;

    /** Copy of the owner_space as a SumMatrixSpace */
    const SumMatrixSpace* owner_space_;
  };

  /** Class for matrix space for SumMatrix. */
  class SumMatrixSpace : public MatrixSpace
  {
  public:
    /** @name Constructors / Destructors */
    //@{
    /** Constructor, given the number of row and columns, as well as
     *  the number of terms in the sum.
     */
    SumMatrixSpace(Index nrows, Index ncols, Index nterms)
        :
        MatrixSpace(nrows, ncols),
        nterms_(nterms)
    {}

    /** Destructor */
    virtual ~SumMatrixSpace()
    {}
    //@}

    /** Accessor functions to get the number of terms in the sum. */
    Index NTerms() const
    {
      return nterms_;
    }

    /** Set the appropriate matrix space for each term. This must
     *  be called for each term or a runtime error will occur */
    void SetTermSpace(Index term_idx, const MatrixSpace& mat_space);

    /** Get the matrix space for a particular term */
    SmartPtr<const MatrixSpace> GetTermSpace(Index term_idx) const;

    /** Method for creating a new matrix of this specific type. */
    SumMatrix* MakeNewSumMatrix() const;

    /** Overloaded MakeNew method for the MatrixSpace base class.
     */
    virtual Matrix* MakeNew() const override;

  private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Default constructor */
    SumMatrixSpace();

    /** Copy Constructor */
    SumMatrixSpace(const SumMatrixSpace&);

    /** Overloaded Equals Operator */
    SumMatrixSpace& operator=(const SumMatrixSpace&);
    //@}

    const Index nterms_;

    std::vector< SmartPtr<const MatrixSpace> > term_spaces_;
  };

} // namespace Ipopt
#endif
