#ifndef __fixedSymMatrix_hh__
#define __fixedSymMatrix_hh__ 1


#include "sthead.h"
#include <cassert>
#include "cdsExcept.h"
#include "cdsString.h"
#include "cdsGenericMatrix.h"

#include "cdsIostream.h"
#include "cdsIomanip.h"

#include "fixedVector.h"

template<class T,int size,int OFFSET=0>
class FixedSymMatrix;
template<class T,int SIZE>
class FixedSymMatrixBase;
template<class T,int SIZE>
istream& operator>>(istream& s, 
		    FixedSymMatrixBase<T,SIZE>& x);
template<class T,int SIZE>
ostream& operator<<(ostream& s,
		    const FixedSymMatrixBase<T,SIZE>& x);

template<class T,int SIZE>
class FixedSymMatrixBase { 
protected:
  T d_[SIZE*(SIZE+1)/2]; //data sotre in lower-triangular column-major layout

public:
  typedef T                  ElementType;
  typedef SymmetricMatrix<T> MatrixType;

  FixedSymMatrixBase() {}
  explicit FixedSymMatrixBase(const T &x) 
    { for (int i=0 ; i<SIZE*(SIZE+1)/2 ; i++) d_[i] = x; }
  // constructor from a lower-triangular row-major array
  explicit FixedSymMatrixBase(const T a[SIZE*(SIZE+1)/2]) 
    { 
     for (int j=0 ; j<SIZE  ; j++) 
       for (int i=j ; i<SIZE  ; i++) 
	 d_[i+SIZE*j-j*(j+1)/2] = a[i*(i+1)/2+j];
    }
  FixedSymMatrixBase(const FixedSymMatrixBase<T,SIZE>& m) 
  {for (int i=0 ; i<SIZE*(SIZE+1)/2 ; i++) d_[i]=m.d_[i]; }

  int nrow() const { return SIZE; }
  int ncol() const { return SIZE; }

  void resize(int r,
	      int c);
  
  const T* pointer() const { return d_; }
  T* pointer() { return d_; }

  T&       updData(int i, int j); 
  const T& getData(int i, int j) const;

  FixedSymMatrixBase& operator=(const FixedSymMatrixBase&m)
    { for (int i=0 ; i<SIZE*(SIZE+1)/2 ; i++) d_[i] = m.d_[i]; return *this;}

  //template<class MATRIX>
  //FixedSymMatrixBase& operator=(const MATRIX& m)
  //{ 
  // for (int j=0 ; j<SIZE  ; j++) 
  //   for (int i=j ; i<SIZE  ; i++) 
  //     data(i,j) = m.data(i,j);
  // return *this;
  //}


  FixedSymMatrixBase& operator+=(const FixedSymMatrixBase&);
  FixedSymMatrixBase& operator-=(const FixedSymMatrixBase&);
  FixedSymMatrixBase& operator*=(const FixedSymMatrixBase<T,SIZE>&);
  
  FixedSymMatrixBase& scale(const T &);

//  template<int OFFSET> 
//  operator const FixedSymMatrix<T,SIZE,OFFSET>&() const
//    { void* p=this; return *(FixedSymMatrix<T,SIZE,OFFSET>*)p; }
//  template<int OFFSET> operator FixedSymMatrix<T,SIZE,OFFSET>&() 
//    { void* p=this; return *(FixedSymMatrix<T,SIZE,OFFSET>*)p; }

  inline void set(const T&);
  inline void setDiag(const T&);
  
  friend istream& operator>><>(istream& s, 
			       FixedSymMatrixBase<T,SIZE>& x);
  friend ostream& operator<<<>(ostream& s,
			       const FixedSymMatrixBase<T,SIZE>& x);
  class RangeError : public CDS_NAMESPACE(out_of_range) {};
  class IOError    : public CDS_NAMESPACE(exception) {}; //FIX: should be ios_base::failure {};
//  class SingularError {};
};

template<class T,int SIZE,int OFFSET>
class FixedSymMatrix : public FixedSymMatrixBase<T,SIZE> { 

public:
  typedef FixedSymMatrix<T,SIZE,OFFSET>  TransposeType;

  FixedSymMatrix() : FixedSymMatrixBase<T,SIZE>() {}
  explicit FixedSymMatrix(const T &x) : FixedSymMatrixBase<T,SIZE>(x) {} 
  // constructor from a lower-triangular row-major array
  explicit FixedSymMatrix(const T a[SIZE*(SIZE+1)/2]) : 
    FixedSymMatrixBase<T,SIZE>(a) {}
  FixedSymMatrix(const FixedSymMatrix<T,SIZE>& m) : 
    FixedSymMatrixBase<T,SIZE>(m) {}

  // silly pass-through required (erroneously) by gcc 3.4.4.
  T&       updData(int i, int j)       { return FixedSymMatrixBase<T,SIZE>::updData(i,j); }
  const T& getData(int i, int j) const { return FixedSymMatrixBase<T,SIZE>::getData(i,j); }

  int offset1() const { return OFFSET; }
  int offset2() const { return OFFSET; }

  inline T& operator()(const int,
		       const int);
  inline const T& operator()(const int,
			     const int) const;
  
};


template<class T,int SIZE>
inline T&
FixedSymMatrixBase<T,SIZE>::updData(const int i1,
				                    const int j1)
{
 assert(i1>=0 && i1<SIZE);
 assert(j1>=0 && j1<SIZE);

 int i=i1;
 int j=j1;
 if (i1<j1) { j=i1; i=j1; }

 return (d_[i + j*SIZE-j*(j+1)/2]); //column-major
}

template<class T,int SIZE>
const inline T&
FixedSymMatrixBase<T,SIZE>::getData(const int i1,
				                    const int j1) const
{
 assert(i1>=0 && i1<SIZE);
 assert(j1>=0 && j1<SIZE);

 int i=i1;
 int j=j1;
 if (i1<j1) { j=i1; i=j1; }

 return (d_[i + j*SIZE-j*(j+1)/2]); //column-major
}

template<class T,int SIZE,int OFFSET>
inline T&
FixedSymMatrix<T,SIZE,OFFSET>::operator()(const int i1,
						   const int i2)
{
 return updData(i1-OFFSET,i2-OFFSET);
}

template<class T,int SIZE,int OFFSET>
const inline T&
FixedSymMatrix<T,SIZE,OFFSET>::operator()(const int i1,
					  const int i2) const
{
 return getData(i1-OFFSET,i2-OFFSET);
}


template<class T,int SIZE> inline             //unary minus
FixedSymMatrix<T,SIZE>
operator-(const FixedSymMatrix<T,SIZE> &m)
{
 FixedSymMatrix<T,SIZE> r(m);
 r.scale(-1);
 return r;
} /* operator- */


template<class T,int SIZE> 
inline FixedSymMatrix<T,SIZE>
operator*(const FixedSymMatrix<T,SIZE>& m1,
	      const FixedSymMatrix<T,SIZE>& m2)
{
 FixedSymMatrix<T,SIZE> r((T)0);
 const int s1 = m1.nrow();
 const int s2 = m1.ncol();
 const int s3 = m2.ncol();
 for (int i=0 ; i<s1 ; i++)
   for (int k=0 ; k<s3 ; k++)
     for (int j=0 ; j<s2 ; j++) 
       r(i,k) += m1(i,j) * m2(j,k);

 return r;
} /* operator* (matrix,matrix) */


template<class T,int SIZE> 
inline FixedSymMatrixBase<T,SIZE>&
FixedSymMatrixBase<T,SIZE>::scale(const T& x)
{
 for (int i=0 ; i<SIZE*(SIZE+1)/2 ; i++)
   d_[i] *= x;
 return *this;
} /* scale */


template<class T,int SIZE> 
inline FixedSymMatrix<T,SIZE>
operator*(const T&                      x,
	  const FixedSymMatrix<T,SIZE>& m)
{
 FixedSymMatrix<T,SIZE> r(m);
 r.scale(x);
 return r;
} /* operator* (matrix,matrix) */


template<class T, int SIZE>
inline void
FixedSymMatrixBase<T,SIZE>::set(const T &x)
{
 for (int i=0 ; i<SIZE*(SIZE+1)/2 ; i++)
   d_[i] = x;
} /* set */


template<class T, int SIZE>
inline void
FixedSymMatrixBase<T,SIZE>::setDiag(const T &x)
{
 for (int i=0 ; i<SIZE ; i++)
   d_[i*SIZE - i*(i-1)/2] = x;
} /* setDiag */


template<class T,int SIZE> 
inline FixedVector<T,SIZE>
operator*(const FixedSymMatrix<T,SIZE>& m,
	  const FixedVector<T,SIZE>&    v)
{
 FixedVector<T,SIZE> r;
 r.set( (T)0 );
 
 for (int i=0 ; i<SIZE ; i++) 
   for (int j=0 ; j<SIZE ; j++) 
     r(i) += m(i,j) * v(j);
 
 return r;
} /* operator* (matrix, vector) */

template<class T,int SIZE> 
inline FixedVector<T,SIZE>
operator*(const FixedVector<T,SIZE>&    v,
	  const FixedSymMatrix<T,SIZE>& m)   // v^T * m
{
 FixedVector<T,SIZE> r;
 r.set( (T)0 );
 
 for (int j=0 ; j<SIZE ; j++) 
   for (int i=0 ; i<SIZE ; i++) 
     r(j) += v(i) * m(i,j);
 
 return r;
} /* operator* (vector, matrix) */

template<class T,int SIZE> 
inline FixedSymMatrix<T,SIZE>
operator+(const FixedSymMatrix<T,SIZE>& m1,
	  const FixedSymMatrix<T,SIZE>& m2)
{
 FixedSymMatrix<T,SIZE> ret = m1;
 ret += m2;
 return ret;
} /* operator+ (matrix, matrix) */

template<class T,int SIZE> 
inline FixedSymMatrix<T,SIZE>
operator-(const FixedSymMatrix<T,SIZE>& m1,
	  const FixedSymMatrix<T,SIZE>& m2)
{
 FixedSymMatrix<T,SIZE> ret = m1;
 ret -= m2;
 return ret;
} /* operator- (matrix, matrix) */

template<class T, int SIZE>
void
FixedSymMatrixBase<T,SIZE>::resize(int r,
				   int c)
{ 
 if ( r!=nrow() ||
      c!=ncol()   ) 
      throw CDS::exception(CDSString("FixedSymMatrix::resize: ") +
			"illegal resize operation attempted.");
} /* resize */
		
template<class T, int SIZE>
FixedSymMatrixBase<T,SIZE>&
FixedSymMatrixBase<T,SIZE>::operator+=(const FixedSymMatrixBase<T,SIZE>& m)
{
 for (int i=0 ; i<SIZE*(SIZE+1)/2 ; i++) 
   d_[i] += m.d_[i];
 return *this;
} /* operator+= */

template<class T, int SIZE>
FixedSymMatrixBase<T,SIZE>&
FixedSymMatrixBase<T,SIZE>::operator-=(const FixedSymMatrixBase<T,SIZE>& m)
{
 for (int i=0 ; i<SIZE*(SIZE+1)/2 ; i++) 
   d_[i] -= m.d_[i];
 return *this;
} /* operator-= */

//FIX: implement in efficient manner
//template<class T, int SIZE>
//FixedSymMatrixBase<T,SIZE>&
//FixedSymMatrixBase<T,SIZE>::operator*=(const FixedSymMatrixBase<T,SIZE>& m)
//{
// FixedSymMatrixBase<T,SIZE> tmp(*this);
//
// set(0.0);
// for (int i=0 ; i<S1 ; i++)
//   for (int j=0 ; j<S2 ; j++)
//     for (int k=0 ; k<S2 ; k++)
//       data(i,j) += tmp.data(i,k) * m.data(k,j);
//
// return *this;
//} /* operator*= */


     
template<class T, int SIZE, int OFFSET>
istream&   
operator>>(istream& s, FixedSymMatrix<T,SIZE,OFFSET>& x)
{
 FMTFLAGS_TYPE flags = s.flags();
 // s.setf(ios::skipws);
 char c; 
 s>>c; if ( c!='{' ) throw FixedSymMatrix<T,SIZE,OFFSET>::IOError();
 for (int i=OFFSET ; i<SIZE+OFFSET ; i++) {
   s>>c; if ( c!='{' ) throw FixedSymMatrix<T,SIZE,OFFSET>::IOError();
   s >> x(i,OFFSET); // this will cause error if the matrix has a zero size!
   for (int j=OFFSET+1 ; j<SIZE+OFFSET ; j++) {
     s >> c; if ( c!=',' ) throw FixedSymMatrix<T,SIZE,OFFSET>::IOError();
     s >> x(i,j);
     s>>c; if ( c!='}' ) throw FixedSymMatrix<T,SIZE,OFFSET>::IOError();
   }
   if (i<SIZE+OFFSET-1) {
     s >> c; if ( c!=',' ) throw FixedSymMatrix<T,SIZE,OFFSET>::IOError();
   }
 }
 s>>c; if ( c!='}' ) throw FixedSymMatrix<T,SIZE,OFFSET>::IOError();
 s.setf(flags);
 return s;
} /* operator>> */

template<class T, int SIZE>
ostream&   
operator<<(ostream& s, const FixedSymMatrixBase<T,SIZE>& x)
{
 s << "{ ";
 for (int i=0 ; i<SIZE ; i++) {
   s << "{ ";
   s << x.getData(i,0); //will fail for vectors of size 0!
   for (int j=1 ; j<SIZE ; j++) 
     s << ", " << x.getData(i,j) ;
   s << " }";
   if (i<SIZE-1) s << ", ";
 }
 s << " }";
 return s;
}

int FixedSymMatrix_test();

#endif /* __fixedMatrix_hh__ */
