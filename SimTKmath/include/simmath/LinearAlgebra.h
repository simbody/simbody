#ifndef SimTK_LINEAR_ALGEBRA_H_
#define SimTK_LINEAR_ALGEBRA_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
 * Authors: Jack Middleton                                                    *
 * Contributors: Michael Sherman                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

/** @file
 * This is the header file that user code should include to pick up the 
 * SimTK Simmath linear algebra tools.
 */


#include "SimTKcommon.h"
#include "simmath/internal/common.h"


namespace SimTK {

//  default for reciprocal of the condition number
// TODO: sherm 080128 I changed this from 0.01 to a more reasonable
// value but it is still wrong because the default should depend
// on the matrix size, something like max(m,n)*eps^(7/8) where
// eps is machine precision for float or double as appropriate.
static const double DefaultRecpCondition = 1e-12;

/**
 * Base class for the various matrix factorizations. 
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
