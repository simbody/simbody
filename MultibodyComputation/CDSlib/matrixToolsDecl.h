
#ifndef __matrixToolsDecl_h__
#define __matrixToolsDecl_h__

#include "cdsVector.h"
#include "cdsList.h"
#include "cdsMatrix.h"
#include "cdsComplex.h"
#include "cdsGenericMatrix.h"

namespace MatrixTools {

  //
  // return a column of the matrix
  //
template<class MATRIX>
CDSVector<typename MATRIX::ElementType>
getColumn(const MATRIX&,
	      const int);

  //
  // return a row of the matrix
  //
template<class MATRIX>
CDSVector<typename MATRIX::ElementType>
getRow(const MATRIX&,
       const int);

  //
  // set a column of the matrix
  //
template<class MATRIX, class VECTOR>
void
setColumn(      MATRIX&,
	      const int,
	      const VECTOR&);

  //
  // set a row of the matrix
  //
template<class MATRIX, class VECTOR>
void
setRow(      MATRIX&,
       const int,
       const VECTOR&);

  //
  // return matrix transpose
  //
template<class MATRIX>
typename MATRIX::TransposeType
transpose(const MATRIX &m);

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

template<class MATRIX> MATRIX
inverse(const MATRIX& matrix,
	    InverseResults<typename MATRIX::MatrixType> ret);

template<class MATRIX> inline MATRIX
inverse(const MATRIX& matrix)
{ return inverse(matrix, InverseResults<typename MATRIX::MatrixType>()); }



//  //
//  // return S * m * transpose(S); FIX: return type is problematic
//  //
//template<class MATRIX> 
//MATRIX
//orthoTransform(const MATRIX& m,
//	       const MATRIX& S);

template<class T> class SVDResults; //forward decl.

  //
  // perform singular value decomposition
  //

template<class MATRIX>
SVDResults<typename MATRIX::ElementType>
svd(const MATRIX&                                  m,
    const char                                     jobu,
    const char                                     jobvt,
	SVDResults<typename MATRIX::ElementType> ret);

template<class MATRIX> inline
SVDResults<typename MATRIX::ElementType>
svd(const MATRIX&                                  m,
    const char                                     jobu='A',
    const char                                     jobvt='A')
{ return svd(m,jobu,jobvt,SVDResults<typename MATRIX::ElementType>()); }

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
        CDS::CDSComplex<T> value;
        CDSVector< CDS::CDSComplex<T> > vector;
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
        T            value;
        CDSVector<T> vector;
    };

    //pairs are returned in order of magnitude of the eigenvalue,
    // smallest to largest
    CDSList< EigenPair > eigenPairs;
    CDSVector<T> work;
    int          info;
};

template<class MATRIX>
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
////template<class MATRIX>
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

//requirements for class GenericMatrix:
//
// has operator()(int,int) to get, set elements
// has members cols(), rows() which return dimensions
// has member resize(int,int) to specify dimensions
//
// has typedef ElementType

#endif /* __matrixToolsDecl_h__ */
