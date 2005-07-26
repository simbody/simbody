/*
   symMatrix.h
   symmetric matrix with heap memory allocation
*/

#ifndef __symMatrix_h__
#define __symMatrix_h__ 1

#include "cdsExcept.h"
#include "cdsMatrix.h"
#include "cdsVector.h"

#include "cdsIostream.h"
#include "cdsIomanip.h"
#include <cassert>

template<class T>
class SymMatrixBase;
template<class T,int OFFSET>
class SymMatrix;

template<class T>
class SymMatrixRep {
  int size;      //# rows/columns
  T* data;
  int count;
public:
  SymMatrixRep(const int s1) : size(s1)
  { data = new T[size*size]; count=1; }
  SymMatrixRep(const SymMatrixRep* r) : size(r->size)
    { 
     data = new T[size*size]; 
     count=1; 
     for (int i=0 ; i<size*(size+1)/2 ; i++) data[i] = r->data[i];
    }
  ~SymMatrixRep() { delete [] data; };
  friend class SymMatrixBase<T>;
};

template<class T> 
class SymMatrixBase {
  SymMatrixRep<T>* rep;
protected:

  //  SymMatrixBase();
  SymMatrixBase(const SymMatrixBase<T> &m);   //copy constructor
  explicit SymMatrixBase(const int size=0);
  explicit SymMatrixBase(const int size, 
             const T &initializer);
  explicit SymMatrixBase(const int size, 
             const T*   a);
  template<class MAT>
  explicit SymMatrixBase(const MAT& m);
  inline void splitRep() 
    { if (rep->count>1) {rep->count--;rep= new SymMatrixRep<T>(rep);}}
public:
  typedef T                  ElementType;
  typedef SymmetricMatrix<T> MatrixType;

  ~SymMatrixBase ();

  SymMatrixBase<T> &resize(const int b1,const int e1);
  int rows() const { return rep->size; }
  int cols() const { return rep->size; }

  T&       updData(int i, int j); 
  const T& getData(int i, int j) const;

  // access to raw pointer
  const T* pointer() const { return rep->data; }
  T* pointer() { return rep->data; }

  SymMatrixBase<T>& operator= (const SymMatrixBase<T>&);
  
  inline void set(const T&);
  inline void setDiag(const T&);

  template<class X>
  SymMatrixBase<T>& scale(const X& x);
  

  //  const SymMatrixBase<T>& operator*= (const SymMatrixBase<T> &m);
  const SymMatrixBase<T>& operator-= (const T &x);
  const SymMatrixBase<T>& operator-= (const SymMatrixBase<T> &m);
  const SymMatrixBase<T>& operator+= (const T &x);
  const SymMatrixBase<T>& operator+= (const SymMatrixBase<T> &m);

  class RangeError : public CDS_NAMESPACE(out_of_range) {};
  class IOError    : public CDS_NAMESPACE(exception) {}; //FIX: should be ios_base::failure {};
};

template<class T,int OFFSET=0> 
class SymMatrix : public SymMatrixBase<T> {
public:
  typedef SymMatrix<T,OFFSET> TransposeType;

  SymMatrix() : SymMatrixBase<T>() {}
  SymMatrix(const int s1) : SymMatrixBase<T>(s1) {}
  SymMatrix(const int s1, 
        const T &initializer) : SymMatrixBase<T>(s1,initializer) {}
  SymMatrix(const int size, 
        const T*  lowerTriRowMajorArray) :
    SymMatrixBase<T>(size,lowerTriRowMajorArray) {}

  template<class MAT>
  SymMatrix(const MAT& m) : SymMatrixBase<T>(m) {}

  SymMatrix(const SymMatrix<T>& v) : SymMatrixBase<T>(v) {}
  ~SymMatrix () {}
  SymMatrix<T>  &operator= (const SymMatrix<T>&v) 
    { this->SymMatrixBase<T>::operator=(v); return *this; }

  // silly pass-through required (erroneously) by gcc 3.4.4 (sherm).
  T&       updData(int i, int j)       { return SymMatrixBase<T>::updData(i,j); } 
  const T& getData(int i, int j) const { return SymMatrixBase<T>::getData(i,j); }

  int offset1() const { return OFFSET; }
  int offset2() const { return OFFSET; }


  //members not present in SymMatrixBase
  T& operator() (int i1, int i2) 
    { return updData(i1-OFFSET,i2-OFFSET); }
  const T& operator() (int i1, int i2) const
    { return getData(i1-OFFSET,i2-OFFSET); }
};



template<class T> inline
T& 
SymMatrixBase<T>::updData(const int i1,
                          const int j1)
{
 assert( i1>=0 && i1<rows() );
 assert( j1>=0 && j1<cols() );

 splitRep();

 int i=i1;
 int j=j1;
 if (i1<j1) { j=i1; i=j1; }

 return (rep->data[i + j*rep->size-j*(j+1)/2]); //column-major
}

template<class T> inline
const T&
SymMatrixBase<T>::getData(const int i1,
                          const int j1) const
  //const version
{
 assert( i1>=0 && i1<rows() );
 assert( j1>=0 && j1<cols() );

 int i=i1;
 int j=j1;
 if (i1<j1) { j=i1; i=j1; }

 return (rep->data[i + j*rep->size-j*(j+1)/2]); //column-major
}

template<class T>
SymMatrix<T>
operator*(const SymMatrix<T>& m1,
      const SymMatrix<T>& m2);

template<class T>
CDSVector<T>
operator*(const SymMatrixBase<T>& ,
      const CDSVectorBase<T>& v);

// moved to matrixTools
//template<class T> SymMatrix<T> transpose(const SymMatrix<T> &m);
//template<class T> SymMatrix<T> inverse(const SymMatrix<T> &m);

template<class T> 
SymMatrix<T>
orthoTransform(const SymMatrix<T>& m,
           const SymMatrix<T>& S);

template<class T>
SymMatrix<T>
blockMat22(const SymMatrix<T>& m11,
       const SymMatrix<T>& m12,
       const SymMatrix<T>& m21,
       const SymMatrix<T>& m22);

template<class T>
ostream&   
operator<<(ostream& s, const SymMatrixBase<T>& m);

template<class T>
istream&   
operator>>(istream& s, SymMatrixBase<T>& m);


//
//constructor
//

template<class T> 
SymMatrixBase<T>::SymMatrixBase(const int size) 
{
 rep = new SymMatrixRep<T>(size);
}

template<class T> 
SymMatrixBase<T>::SymMatrixBase(const int size,
                const T&  init)
{
 rep = new SymMatrixRep<T>(size);
 for (int i=0 ; i<rep->size*(rep->size+1)/2 ; i++) 
   rep->data[i] = init;
}

template<class T> 
SymMatrixBase<T>::SymMatrixBase(const int size,
                const T*  a)
  // constructor from array: not real safe
{
 rep = new SymMatrixRep<T>(size);
 for (int j=0 ; j<size  ; j++) 
   for (int i=j ; i<size  ; i++) 
     rep->data[i+size*j-j*(j+1)/2] = a[i*(i+1)/2+j];
}

//
//copy contructor
//

template<class T>
SymMatrixBase<T>::SymMatrixBase(const SymMatrixBase<T> &m)
{
 rep = m.rep;
 rep->count++;
}

template<class T>
template<class MAT>
SymMatrixBase<T>::SymMatrixBase(const MAT &m)
{
 rep = new SymMatrixRep<T>(m.rows());
 for (int i=0 ; i<rows() ; i++)
   for (int j=i ; j<cols() ; j++)
     updData(i,j) = m(m.offset1()+i,m.offset2()+j);
}

//
//destructor
//

template<class T> 
SymMatrixBase<T>::~SymMatrixBase()
{
 if (--rep->count <= 0) delete rep;
}

template<class T>
SymMatrixBase<T>& 
SymMatrixBase<T>::resize(const int rows,
             const int cols)
{
 assert( rows == cols );

 if (rows != rep->size ) {
   if (--rep->count <= 0) delete rep;
   rep = new SymMatrixRep<T>(rows);
 }
 
 //FIX: should there be a copy operation???
 return *this;
} /*resize*/

template<class T> 
SymMatrixBase<T>&
SymMatrixBase<T>::operator=(const SymMatrixBase<T> &m)
{
 if (&m == this)
   return *this;

 m.rep->count++; 
 if (--rep->count <= 0) delete rep;

 rep = m.rep;

 return *this;
}

template<class T>
inline void
SymMatrixBase<T>::set(const T &x)
{
 for (int i=0 ; i<rep->size*(rep->size+1)/2 ; i++) 
   rep->data[i] = x;
}


template<class T>
inline void
SymMatrixBase<T>::setDiag(const T &x)
{
 for (int i=0 ; i<rows() ; i++)
   updData(i,i) = x;
}


template<class T> 
template<class X>
SymMatrixBase<T>&
SymMatrixBase<T>::scale(const X& x)
{
 for (int i=0 ; i<rep->size*(rep->size+1)/2 ; i++) 
   rep->data[i] *= x;
 return *this;
}


template<class T>
const SymMatrixBase<T>&
SymMatrixBase<T>::operator+=(const SymMatrixBase<T>& m)
{
 assert( m.rep->size == rep->size );
 splitRep();
  
 for (int i=0 ; i<rep->size*(rep->size+1)/2 ; i++) 
   rep->data[i] += m.rep->data[i];
 return *this;
}

template<class T>
const SymMatrixBase<T>&
SymMatrixBase<T>::operator-=(const SymMatrixBase<T>& m)
{
 assert( m.rep->size == rep->size );
 splitRep();
  
 for (int i=0 ; i<rep->size*(rep->size+1)/2 ; i++) 
   rep->data[i] -= m.rep->data[i];
 return *this;
}

//FIX: implement
//template<class T> 
//SymMatrix<T>
//operator*(const SymMatrix<T>& m1,
//      const SymMatrix<T>& m2)
//{
// assert( m1.cols() == m2.rows() );
//
// SymMatrix<T> r(m1.rows(),m2.cols(),(T)0);
// for (int i=0 ; i<m1.rows() ; i++)
//   for (int k=0 ; k<m2.cols() ; k++)
//     for (int j=0 ; j<m1.cols() ; j++) 
//       r(i,k) += m1(i,j) * m2(j,k);
//
// return r;
//} /* operator* (matrix,matrix) */
//
//template<class T> 
//CDSVector<T>
//operator*(const SymMatrixBase<T>& m,
//      const CDSVectorBase<T>& v)
//{
// assert( m.cols() == v.size() );
// CDSVector<T> r(m.rows(),(T)0);
// 
// for (int i=0 ; i<m.rows() ; i++) 
//   for (int j=0 ; j<m.cols() ; j++) 
//     r(i) += m.data(i,j) * v.data(j);
// 
// return r;
//} /* operator* (matrix, vector) */


     

template<class T>
istream&   
operator>>(istream& s, SymMatrixBase<T>& m)
{
 FMTFLAGS_TYPE flags = s.flags();
 s.setf(ios::skipws);
 char c; 
 s>>c; if ( c!='{' ) throw SymMatrixBase<T>::IOError();
 for (int i=0 ; i<m.rows() ; i++) {
   s>>c; if ( c!='{' ) throw SymMatrixBase<T>::IOError();
   s >> m.data(i,0); // this will cause error if the matrix has a zero size!
   for (int j=1 ; j<m.cols() ; j++) {
     s >> c; if ( c!=',' ) throw SymMatrixBase<T>::IOError();
     s >> m.data(i,j);
     s>>c; if ( c!='}' ) throw SymMatrixBase<T>::IOError();
   }
   if (i<m.rows()-1) {
     s >> c; if ( c!=',' ) throw SymMatrixBase<T>::IOError();
   }
 }
 s>>c; if ( c!='}' ) throw SymMatrixBase<T>::IOError();
 s.setf(flags);
 return s;
} /* operator>> */

template<class T>
ostream&   
operator<<(ostream& s, const SymMatrixBase<T>& m)
{
 int width = s.width(); // apply width to numeric fields only
 s << "{ ";
 for (int i=0 ; i<m.rows() ; i++) {
   s << "{ ";
   s << setw(width) << m.getData(i,0); //will fail for vectors of size 0!
   for (int j=1 ; j<m.cols() ; j++) 
     s << ", " << setw(width) << m.getData(i,j) ;
   s << " }";
   if (i<m.rows()-1) s << ", ";
 }
 s << " }";
 return s;
}

int SymMatrix_test();

#endif /* __symMatrix_h__ */
