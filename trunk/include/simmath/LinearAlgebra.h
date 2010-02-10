#ifndef SimTK_LINEAR_ALGEBRA_H_
#define SimTK_LINEAR_ALGEBRA_H_
/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmath(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *

 * Portions copyright (c) 2006-2007 Stanford University and Jack Middleton.   *
 * Contributors:                                                              *
 *
 * Permission is hereby granted, free of charge, to any person obtaining      *
 * a copy of this software and associated documentation files (the            *
 * "Software"), to deal in the Software without restriction, including        *
 * without limitation the rights to use, copy, modify, merge, publish,        *
 * distribute, sublicense, and/or sell copies of the Software, and to         *
 * permit persons to whom the Software is furnished to do so, subject         *
 * to the following conditions:                                               *
 *
 * The above copyright notice and this permission notice shall be included    *
 * in all copies or substantial portions of the Software.                     *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS    *
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF                 *
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.     *
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY       *
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,       *
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE          *
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                     *
 */


/** @file
 * This is the header file that user code should include to pick up the 
 * SimTK Simmath linear algebra tools.
 */




#include <limits.h>
#include "SimTKcommon.h"
#include "SimTKmath.h"
#include "SimTKcommon/internal/BigMatrix.h"
#include "internal/common.h"

namespace SimTK {

//  default for reciprocal of the condition number
// TODO: sherm 080128 I changed this from 0.01 to a more reasonable
// value but it is still wrong because the default should depend
// on the matrix size, something like max(m,n)*eps^(7/8) where
// eps is machine precision for float or double as appropriate.
static const double DefaultRecpCondition = 1e-12;
/**
 * Abstract class for performing matrix factorizations 
 */
class SimTK_SIMMATH_EXPORT Factor {
public:

  Factor() {}
  /// creates an factorization of a matrix
  template <class ELT> Factor( Matrix_<ELT> m );
  /// solves a single right hand side using a factorization
  template <class ELT> void solve( const Vector_<ELT>& b, Vector_<ELT>& x ) const;
  /// solves multiple right hand sides using a factorization
  template <class ELT> void solve( const Matrix_<ELT>& b, Matrix_<ELT>& x ) const;
  
}; // class Factor

class FactorLURepBase;

/**
 * Class for performing LU matrix factorizations 
 */
class SimTK_SIMMATH_EXPORT FactorLU: public Factor {
    public:

    ~FactorLU();

    FactorLU();
    FactorLU( const FactorLU& c );
    FactorLU& operator=(const FactorLU& rhs);

    template <class ELT> FactorLU( const Matrix_<ELT>& m );
    /// factors a matrix
    template <class ELT> void factor( const Matrix_<ELT>& m );
    /// solves a single right hand side 
    template <class ELT> void solve( const Vector_<ELT>& b, Vector_<ELT>& x ) const;
    /// solves multiple  right hand sides 
    template <class ELT> void solve( const Matrix_<ELT>& b, Matrix_<ELT>& x ) const;

    /// returns the lower triangle of an LU factorization 
    template <class ELT> void getL( Matrix_<ELT>& l ) const;
    /// returns the upper triangle of an LU factorization 
    template <class ELT> void getU( Matrix_<ELT>& u ) const;
    /// returns the inverse of a matrix using an LU factorization
    template < class ELT > void inverse(  Matrix_<ELT>& m ) const;

    /// returns true if matrix was singular 
    bool isSingular() const;
    /// returns the first diagonal which was found to be singular
    int getSingularIndex() const;


    protected:
    class FactorLURepBase *rep;

}; // class FactorLU


class FactorQTZRepBase;
/**
 * Class to perform a QTZ (linear least squares) factorization
 */
class SimTK_SIMMATH_EXPORT FactorQTZ: public Factor {
    public:

    ~FactorQTZ();

    FactorQTZ();
    FactorQTZ( const FactorQTZ& c );
    FactorQTZ& operator=(const FactorQTZ& rhs);
    /// do QTZ factorization of a matrix
    template <typename ELT> FactorQTZ( const Matrix_<ELT>& m);
    /// do QTZ factorization of a matrix for a given reciprocal condition number
    template <typename ELT> FactorQTZ( const Matrix_<ELT>& m, double rcond );
    /// do QTZ factorization of a matrix for a given reciprocal condition number
    template <typename ELT> FactorQTZ( const Matrix_<ELT>& m, float rcond );
    /// do QTZ factorization of a matrix
    template <typename ELT> void factor( const Matrix_<ELT>& m);
    /// do QTZ factorization of a matrix for a given reciprocal condition number
    template <typename ELT> void factor( const Matrix_<ELT>& m, float rcond );
    /// do QTZ factorization of a matrix for a given reciprocal condition number
    template <typename ELT> void factor( const Matrix_<ELT>& m, double rcond );
    /// solve  for a vector x given a right hand side vector b
    template <typename ELT> void solve( const Vector_<ELT>& b, Vector_<ELT>& x ) const;
    /// solve  for an array of vectors  given multiple  right hand sides  
    template <typename ELT> void solve( const Matrix_<ELT>& b, Matrix_<ELT>& x ) const;

    template < class ELT > void inverse(  Matrix_<ELT>& m ) const;

    /// returns the rank of the matrix
    int getRank() const;
    /// returns the actual reciprocal condition number at this rank
    double getRCondEstimate() const;
//    void setRank(int rank); TBD

    protected:
    class FactorQTZRepBase *rep;
}; // class FactorQTZ
/**
 * Class to compute Eigen values and Eigen vectors of a matrix
 */
class SimTK_SIMMATH_EXPORT Eigen {
    public:

    ~Eigen();

    Eigen();
    Eigen( const Eigen& c );
    Eigen& operator=(const Eigen& rhs);

    /// create a default eigen class
    template <class ELT> Eigen( const Matrix_<ELT>& m );
    /// supply matrix which eigen values will be computed for  
    template <class ELT> void factor( const Matrix_<ELT>& m );
    /// get all the eigen values and eigen vectors of a matrix 
    template <class VAL, class VEC> void getAllEigenValuesAndVectors( Vector_<VAL>& values, Matrix_<VEC>& vectors);
    /// get all the eigen values of a matrix 
    template <class T> void getAllEigenValues( Vector_<T>& values);

    /// get a few eigen values  and eigen vectors of a symmetric matrix which are within a range of indices
    template <class VAL, class VEC> void getFewEigenValuesAndVectors( Vector_<VAL>& values, Matrix_<VEC>& vectors, int ilow, int ihi);
    /// get a few eigen vectors of a symmetric matrix which are within a range of indices
    template <class T> void getFewEigenVectors( Matrix_<T>& vectors, int ilow, int ihi );
    /// get a few eigen values of a symmetric matrix which are within a range of indices
    template <class T> void getFewEigenValues( Vector_<T>& values, int ilow, int ihi );

    /// get a few eigen values  and eigen vectors of a symmetric matrix which are within a range of eigen values
    template <class VAL, class VEC> void getFewEigenValuesAndVectors( Vector_<VAL>& values, Matrix_<VEC>& vectors, typename CNT<VAL>::TReal rlow, typename CNT<VAL>::TReal rhi);
    /// get a few eigen vectors of a symmetric matrix which are within a range of eigen values
    template <class T> void getFewEigenVectors( Matrix_<T>& vectors, typename CNT<T>::TReal rlow, typename CNT<T>::TReal rhi );
    /// get a few eigen values of a symmetric matrix which are within a range of eigen values
    template <class T> void getFewEigenValues( Vector_<T>& values, typename CNT<T>::TReal rlow, typename CNT<T>::TReal rhi );

     
    protected:
    class EigenRepBase *rep;

}; // class Eigen
/**
 * Class to compute a singular value decomposition of a matrix
 */
class SimTK_SIMMATH_EXPORT FactorSVD: public Factor {
    public:

    ~FactorSVD();
    /// default constructor  
    FactorSVD();
    /// copy  constructor  
    FactorSVD( const FactorSVD& c );
    /// copy  assign  
    FactorSVD& operator=(const FactorSVD& rhs);

    /// constructor 
    template < class ELT > FactorSVD( const Matrix_<ELT>& m );
    /// singular value decomposition of a matrix using the specified reciprocal of the condition
    /// number rcond
    template < class ELT > FactorSVD( const Matrix_<ELT>& m, float rcond );
    /// singular value decomposition of a matrix using the specified reciprocal of the condition
    /// number rcond
    template < class ELT > FactorSVD( const Matrix_<ELT>& m, double rcond );
    /// supply the matrix to do a singular value decomposition 
    template < class ELT > void factor( const Matrix_<ELT>& m );
    /// supply the matrix to do a singular value decomposition using the specified 
    /// reciprocal of the condition number rcond
    template < class ELT > void factor( const Matrix_<ELT>& m, float rcond );
    /// supply the matrix to do a singular value decomposition using the specified reciprocal of the condition
    /// reciprocal of the condition number rcond
    template < class ELT > void factor( const Matrix_<ELT>& m, double rcond );

    /// get the singular values and singular vectors of the matrix
    template < class T > void getSingularValuesAndVectors( Vector_<typename CNT<T>::TReal>& values, 
                              Matrix_<T>& leftVectors,  Matrix_<T>& rightVectors );
    /// get just the singular values of the matrix
    template < class T > void getSingularValues( Vector_<T>& values);

    /// get rank of the matrix 
    int getRank();
    /// get inverse of the matrix  using singular value decomposition (sometimes called the pseudo inverse)
    template < class ELT > void inverse(  Matrix_<ELT>& m );
    /// solve for x given a right hand side vector using the singular value decomposition
    template <class ELT> void solve( const Vector_<ELT>& b, Vector_<ELT>& x );
    /// solve for a set of x vectors  given multiple right hand side vectors 
    /// using the singular value decomposition
    template <class ELT> void solve( const Matrix_<ELT>& b, Matrix_<ELT>& x );

    protected:
    class FactorSVDRepBase *rep;

}; // class FactorSVD

} // namespace SimTK 

#endif //SimTK_LINEAR_ALGEBRA_H_
