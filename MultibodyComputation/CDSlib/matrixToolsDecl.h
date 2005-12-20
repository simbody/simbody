
#ifndef __matrixToolsDecl_h__
#define __matrixToolsDecl_h__

#include <cdsVector.h>
#include <cdsList.h>
#include <cdsMatrix.h>
#include <cdsComplex.h>
#include <matrix.h>

namespace MatrixTools {

  //
  // return a column of the matrix
  //
template<class Matrix>
CDSVector<typename Matrix::ElementType>
getColumn(const Matrix&,
	  const int);

  //
  // return a row of the matrix
  //
template<class Matrix>
CDSVector<typename Matrix::ElementType>
getRow(const Matrix&,
       const int);

  //
  // set a column of the matrix
  //
template<class Matrix, class Vector>
void
setColumn(      Matrix&,
	  const int,
	  const Vector&);

  //
  // set a row of the matrix
  //
template<class Matrix, class Vector>
void
setRow(      Matrix&,
       const int,
       const Vector&);

  //
  // return matrix transpose
  //
template<class Matrix>
typename Matrix::TransposeType
transpose(const Matrix &m);

  //
  // matrix inversion
  //
template<class T> class InverseResults; //forward decl.

template<class T>
class InverseResults<FullMatrix<T> >{
public:
  CDSVector<T> work;
  int          info;
};

template<class T>
class InverseResults<SymmetricMatrix<T> >{
public:
  CDSVector<T> work;
  int          info;
};

template<class MATRIX>
MATRIX
inverse(const MATRIX& matrix,
	      InverseResults<typename MATRIX::MatrixType> ret);

template<class MATRIX> inline
MATRIX
inverse(const MATRIX& matrix)
{ return inverse(matrix, InverseResults<typename MATRIX::MatrixType>()); }



//  //
//  // return S * m * transpose(S); FIX: return type is problematic
//  //
//template<class Matrix> 
//Matrix
//orthoTransform(const Matrix& m,
//	       const Matrix& S);

template<class T> class SVDResults; //forward decl.

  //
  // perform singular value decomposition
  //

template<class Matrix>
SVDResults<typename Matrix::ElementType>
svd(const Matrix&                                  m,
    const char                                     jobu,
    const char                                     jobvt,
	SVDResults<typename Matrix::ElementType> ret);

template<class Matrix> inline
SVDResults<typename Matrix::ElementType>
svd(const Matrix&                                  m,
    const char                                     jobu='A',
    const char                                     jobvt='A')
{ return svd(m,jobu,jobvt,SVDResults<typename Matrix::ElementType>()); }

template<class T>
class SVDResults {
public:
  CDSVector<T> sigma;
  CDSMatrix<T> u;
  CDSMatrix<T> vT;
  CDSVector<T> work;
  int          info;
};

template<class T> class EigenResults; //forward decl.

  //
  // perform eigenvalue analysis
  //
template<class MATRIX>
EigenResults<typename MATRIX::MatrixType>
eigen(const MATRIX&                  m,
      const char                     jobz,
	  EigenResults<typename MATRIX::MatrixType> ret);

template<class MATRIX>
EigenResults<typename MATRIX::MatrixType>
eigen(const MATRIX&                  m,
      const char                     jobz='V')
{ return eigen(m,jobz,EigenResults<typename MATRIX::MatrixType>()); }

//  //general (nonsymmetric) matrix definition
template<class T>
class EigenResults<FullMatrix<T> >{
public:
  struct EigenPair { 
    CDS::Complex<T> value;
    CDSVector< CDS::Complex<T> > vector;
  };
  //pairs are returned in complex conjugate pairs
  CDSList< EigenPair > eigenPairs;
  CDSVector<T> work;
  int          info;
};

template<class T>
class EigenResults<SymmetricMatrix<T> > {
public:
  //currently only setup for symmetric matrices
  struct EigenPair { 
    T value;
    CDSVector<T> vector;
  };
  //pairs are returned in order of magnitude of the eigenvalue,
  // smallest to largest
public: // ??
  CDSList< EigenPair > eigenPairs;
  CDSVector<T> work;
  int          info;
};

template<class Matrix>
class EigenResults {
  //currently only setup for symmetric matrices
};

//template<class T,int SIZE, int OFFSET>
//struct EigenResults<FixedSymMatrix<T,SIZE,OFFSET> > {
//  //currently only setup for symmetric matrices
//  struct EigenPair { 
//    T value;
//    CDSVector<T> vector;
//  };
//  //pairs are returned in order of magnitude of the eigenvalue,
//  // smallest to largest
//  CDSList< EigenPair > eigenPairs;
//  CDSVector<T> work;
//  int          info;
//};
//
////template<class Matrix>
//struct EigenResults {
//  //currently only setup for symmetric matrices
//  struct EigenPair { 
//    T value;
//    CDSVector<T> vector;
//  };
//  //pairs are returned in order of magnitude of the eigenvalue,
//  // smallest to largest
//  CDSList< EigenPair > eigenPairs;
//  CDSVector<T> work;
//  int          info;
//};

  int matrixToolsTest();

};

//requirements for class Matrix:
//
// has operator()(int,int) to get, set elements
// has members cols(), rows() which return dimensions
// has member resize(int,int) to specify dimensions
//
// has typedef ElementType

#endif /* __matrixToolsDecl_h__ */
