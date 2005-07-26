/*
   CDSMatrix.h
   headers for my own Matrix template class
   9/16/95 - CDS
   1/16/96 - changed constructor - note!
   5/20/96 - changed array bounds checking - now it is implemented as
             a macro CHECKRANGE
*/

#ifndef __CDSMatrixDecl_h__
#define __CDSMatrixDecl_h__ 1

#include <cdsVector.h>
#include <assert.h>
#include "cdsExcept.h"
#include "matrix.h"

template<class T>
class CDSMatrixBase;
template<class T,int OFFSET1,int OFFSET2>
class CDSMatrix;

template<class T>
class CDSMatrixRep {
  int size1;      //# rows
  int size2;      //# columns
  T* data;
  int count;
public:
  CDSMatrixRep(const int s1,
           const int s2) : size1(s1), size2(s2) 
  { data = new T[size1*size2]; count=1; }
  CDSMatrixRep(const CDSMatrixRep* r) : size1(r->size1), size2(r->size2)
    { 
     data = new T[size1*size2]; 
     count=1; 
     for (int i=0 ; i<size1*size2 ; i++) data[i] = r->data[i];
    }
  ~CDSMatrixRep() { delete [] data; };
  friend class CDSMatrixBase<T>;
};

template<class T> 
class CDSMatrixBase {
  CDSMatrixRep<T>* rep;
protected:

  //  CDSMatrixBase();
  CDSMatrixBase(const CDSMatrixBase<T> &m);   //copy constructor
  //  template<int OFFSET1,int OFFSET2>
  //  CDSMatrixBase(const CDSMatrix<T,OFFSET1,OFFSET2> &m);   
  explicit CDSMatrixBase(const int s1=0, 
             const int s2=0);
  explicit CDSMatrixBase(const int s1, 
             const int s2,
             const T &initializer);
  explicit CDSMatrixBase(const int s1, 
             const int s2,
             const T*   a);
  template<class MAT>
  explicit CDSMatrixBase(const MAT& m);
  inline void splitRep() 
    { if (rep->count>1) {rep->count--;rep= new CDSMatrixRep<T>(rep);}}
//  template<int SIZE1,int SIZE2,int OFFSET1,int OFFSET2>
//  CDSMatrixBase(const FixedMatrix<T,SIZE1,SIZE2,OFFSET1,OFFSET2>&);
public:
  typedef T             ElementType;
  typedef FullMatrix<T> MatrixType;

  ~CDSMatrixBase ();

  CDSMatrixBase<T> &resize(const int b1,const int e1);
  int rows() const { return rep->size1; }
  int cols() const { return rep->size2; }

  T&       updData(int i, int j); 
  const T& getData(int i, int j) const;

  // access to raw pointer
  const T* pointer() const { return rep->data; }
  T* pointer() { return rep->data; }

  CDSMatrixBase<T>& operator= (const CDSMatrixBase<T>&);
  
  inline void set(const T&);
  inline void setDiag(const T&);

  template<class X>
  CDSMatrixBase<T>& scale(const X& x);
  

  //  const CDSMatrixBase<T>& operator*= (const CDSMatrixBase<T> &m);
  const CDSMatrixBase<T>& operator-= (const T &x);
  const CDSMatrixBase<T>& operator-= (const CDSMatrixBase<T> &m);
  const CDSMatrixBase<T>& operator+= (const T &x);
  const CDSMatrixBase<T>& operator+= (const CDSMatrixBase<T> &m);

  class RangeError : public CDS_NAMESPACE(out_of_range) {};
  class IOError    : public CDS_NAMESPACE(exception) {}; //FIX: should be ios_base::failure {};
};

template<class T,int OFFSET1=0,int OFFSET2=0> 
class CDSMatrix : public CDSMatrixBase<T> {
public:
  typedef CDSMatrix<T,OFFSET1,OFFSET2> TransposeType;

  CDSMatrix() : CDSMatrixBase<T>() {}
  CDSMatrix(const int s1, const int s2) : CDSMatrixBase<T>(s1,s2) {}
  CDSMatrix(const int s1, 
        const int s2,
        const T &initializer) : CDSMatrixBase<T>(s1,s2,initializer) {}

//  template<int size1,int size2,int o1,int o2>
//  CDSMatrix(const FixedMatrix<T,size1,size2,o1,o2>& fm) :
//    CDSMatrixBase<T>(fm) {}
  template<class MAT>
  CDSMatrix(const MAT& m) : CDSMatrixBase<T>(m) {}

  CDSMatrix(const CDSMatrix<T>& v) : CDSMatrixBase<T>(v) {}
  ~CDSMatrix () {}
  CDSMatrix<T>  &operator= (const CDSMatrix<T>&v) 
    { this->CDSMatrixBase<T>::operator=(v); return *this; }

  // These pass-throughs required (incorrectly) by gcc 3.4.4 (sherm).  
  T&       updData(int i, int j)       { return CDSMatrixBase<T>::updData(i,j); }
  const T& getData(int i, int j) const { return CDSMatrixBase<T>::getData(i,j); }

  int offset1() const { return OFFSET1; }
  int offset2() const { return OFFSET2; }

//  T& data(const int i,
//        const int j); 
//  const T& data(const int i,
//          const int j) const;

  //members not present in CDSMatrixBase
  T& operator() (int i1, int i2) 
    { return updData(i1-OFFSET1,i2-OFFSET2); }
  const T& operator() (int i1, int i2) const
    { return getData(i1-OFFSET1,i2-OFFSET2); }
};

//template<class T> inline
//CDSMatrixBase<T>::CDSMatrixBase() : 
//  size1(0), size2(0)
//{
// data = new T[0];
//}


template<class T> inline
T& 
CDSMatrixBase<T>::updData(const int i1,
                          const int i2)
{
 assert( i1>=0 && i1<rows() );
 assert( i2>=0 && i2<cols() );

 splitRep();

 return rep->data[i1 + i2*rows()]; //column major
} /* data */

template<class T> inline
const T&
CDSMatrixBase<T>::getData(const int i1,
                          const int i2) const
  //const version
{
 assert( i1>=0 && i1<rows() );
 assert( i2>=0 && i2<cols() );

 return rep->data[i1 + i2*rows()];
} /* const data */

template<class T,int o1, int o2>
CDSVector<T>
subCol(const CDSMatrix<T,o1,o2>& m,
       const int                 col,
       const int                 rowBeg,
       const int                 rowEnd);

template<class T>
CDSMatrix<T>
operator+(const CDSMatrix<T>& m1,
      const CDSMatrix<T>& m2);

template<class T>
CDSMatrix<T>
operator-(const CDSMatrix<T>& m1,
      const CDSMatrix<T>& m2);

template<class T>
CDSMatrix<T>
operator*(const CDSMatrix<T>& m1,
      const CDSMatrix<T>& m2);

template<class T>
CDSVector<T>
operator*(const CDSMatrixBase<T>& ,
      const CDSVectorBase<T>& v);

// moved to matrixTools
//template<class T> CDSMatrix<T> transpose(const CDSMatrix<T> &m);
//template<class T> CDSMatrix<T> inverse(const CDSMatrix<T> &m);

template<class T> 
CDSMatrix<T>
orthoTransform(const CDSMatrix<T>& m,
           const CDSMatrix<T>& S);

template<class T>
CDSMatrix<T>
blockMat22(const CDSMatrix<T>& m11,
           const CDSMatrix<T>& m12,
           const CDSMatrix<T>& m21,
           const CDSMatrix<T>& m22);

template<class T>
ostream&   
operator<<(ostream& s, const CDSMatrixBase<T>& m);

template<class T>
istream&   
operator>>(istream& s, CDSMatrixBase<T>& m);

int CDSMatrix_test();

#endif /* __CDSMatrixDecl_h__ */
