#ifndef __fixedMatrixDecl_h__
#define __fixedMatrixDecl_h__ 1


#include "sthead.h"

#include "cdsExcept.h"
#include "cdsString.h"
#include "cdsGenericMatrix.h"
#include "cdsIostream.h"
#include "fixedVector.h"

// Provide template defaults for FixedMatrix class.
#include "fixedMatrixFwd.h"

#include <cassert>

template<class T,int size1,int size2>
class FixedMatrixBase;

template<class T,int size1,int size2>
istream& operator>>(istream& s, 
		    FixedMatrixBase<T,size1,size2>& x);
template<class T,int size1,int size2>
ostream& operator<<(ostream& s,
		    const FixedMatrixBase<T,size1,size2>& x);

template<class T,int size1,int size2>
class FixedMatrixBase { 
protected:
  T d_[size1*size2];

public:
  typedef T             ElementType;
  typedef FullMatrix<T> MatrixType;

  FixedMatrixBase() {}
  explicit FixedMatrixBase(const T &x) 
    { for (int i=0 ; i<size1*size2 ; i++) d_[i] = x; }
  // constructor from a row-major array
  explicit FixedMatrixBase(const T a[size1*size2]) 
    { 
     for (int i=0 ; i<size1 ; i++) 
       for (int j=0 ; j<size2 ; j++) 
	 d_[i+size1*j] = a[i*size2+j]; 
    }
  FixedMatrixBase(const FixedMatrixBase<T,size1,size2>& m) 
    {for (int i=0 ; i<size1*size2 ; i++) d_[i]=m.d_[i]; }

  int rows() const { return size1; }
  int cols() const { return size2; }

  void resize(int r,
	      int c);
  
  const T* pointer() const { return d_; }
  T* pointer() { return d_; }

  T&       data(int i, int j); 
  const T& data(int i, int j) const;

  FixedMatrixBase& operator=(const FixedMatrixBase&m)
    { for (int i=0 ; i<size1*size2 ; i++) d_[i] = m.d_[i]; return *this;}


  FixedMatrixBase& operator+=(const FixedMatrixBase&);
  FixedMatrixBase& operator-=(const FixedMatrixBase&);
  FixedMatrixBase& operator*=(const FixedMatrixBase<T,size1,size2>&);
//  const FixedMatrix& operator/=(const T &);
  
  FixedMatrixBase& scale(const T &);
//
//  void copyArray(const T* a);
//

//  template<int noffset1,int noffset2>   //change offsets
//  operator FixedMatrix<T,size1,size2,noffset1,noffset2>&()
//    { void* p=this; return *(FixedMatrix<T,size1,size2,noffset1,noffset2>*)p; }

//  template<int o1,int o2> 
//  operator const FixedMatrix<T,size1,size2,o1,o2>&() const
//    { void* p=(void*)this; return *(FixedMatrix<T,size1,size2,o1,o2>*)p; }
//  template<int o1,int o2> operator FixedMatrix<T,size1,size2,o1,o2>&() 
//    { void* p=this; return *(FixedMatrix<T,size1,size2,o1,o2>*)p; }

  //    { void* p=this; return *(FixedMatrix<T,size1,size2,o1,o2>*)p; }
//  template<int o1,int o2> operator FixedMatrix<T,size1,size2,o1,o2>() 
//    { void* p=this; return *(FixedMatrix<T,size1,size2,o1,o2>*)p; }

  inline void set(const T&);
  inline void setDiag(const T&);
  //  FixedMatrix transpose();
  
//
  friend istream& operator>><>(istream& s, 
			     FixedMatrixBase<T,size1,size2>& x);
  friend ostream& operator<<<>(ostream& s,
			     const FixedMatrixBase<T,size1,size2>& x);
  class RangeError : public CDS_NAMESPACE(out_of_range) {};
  class IOError    : public CDS_NAMESPACE(exception) {}; //FIX: should be ios_base::failure {};
//  class SingularError {};
};

template<class T,int size1,int size2,int OFFSET1,int OFFSET2>
class FixedMatrix : public FixedMatrixBase<T,size1,size2> { 

public:
  typedef FixedMatrix<T,size2,size1,OFFSET2,OFFSET1>  TransposeType;

  FixedMatrix() : FixedMatrixBase<T,size1,size2>() {}
  explicit FixedMatrix(const T &x) : FixedMatrixBase<T,size1,size2>(x) {} 
  // constructor from a row-major array
  explicit FixedMatrix(const T a[size1*size2]) : 
    FixedMatrixBase<T,size1,size2>(a) {}
  // this def. required (and works) for gcc 3.3
  //FixedMatrix(const FixedMatrix& m) : 
  //  FixedMatrixBase<T,size1,size2>(m) {}
  // this def required for gcc<3.3
  FixedMatrix(const FixedMatrixBase<T,size1,size2>& m) : 
    FixedMatrixBase<T,size1,size2>(m) {}

//  
//  FixedMatrix& operator=(const FixedMatrix&m)
//    { for (int i=0 ; i<size1*size2 ; i++) d_[i] = m.d_[i]; return *this;}
//
//  template<int noffset1,int noffset2>   //change offsets
//  operator FixedMatrix<T,size1,size2,noffset1,noffset2>&()
//    { void* p=this; return *(FixedMatrix<T,size1,size2,noffset1,noffset2>*)p; }

//  template<int o1,int o2> operator const FixedMatrix<T,size1,size2,o1,o2>&() const
//    { void* p=this; return *(FixedMatrix<T,size1,size2,o1,o2>*)p; }
  //    { void* p=this; return *(FixedMatrix<T,size1,size2,o1,o2>*)p; }
//  template<int o1,int o2> 
//  operator FixedMatrix<T,size1,size2,o1,o2>() 
//    { void* p=this; return *(FixedMatrix<T,size1,size2,o1,o2>*)p; }
//  template<int o1,int o2> operator FixedMatrix<T,size1,size2,o1,o2>() 
//    { void* p=this; return *(FixedMatrix<T,size1,size2,o1,o2>*)p; }

  int offset1() const { return OFFSET1; }
  int offset2() const { return OFFSET2; }

  T&       data(int i, int j)       { return FixedMatrixBase<T,size1,size2>::data(i,j); } 
  const T& data(int i, int j) const { return FixedMatrixBase<T,size1,size2>::data(i,j); }

  inline T& operator()(int, int);
  inline const T& operator()(int, int) const;
  
};

//template<class T,int size>
//FixedMatrix<T,size,size> 
//inverse(const FixedMatrix<T,size,size>&);

//template<class T,int size1,int size2>
//FixedMatrix<T,size2,size1> 
//transpose(const FixedMatrix<T,size1,size2>&);

template<class T,int SIZE1,int SIZE2>
inline T&
FixedMatrixBase<T,SIZE1,SIZE2>::data(const int i1,
				     const int i2)
{
 assert(i1>=0 && i1<SIZE1);
 assert(i2>=0 && i2<SIZE2);

 return (d_[i1 + SIZE1*i2]); //column-major
} /* data */

template<class T,int SIZE1,int SIZE2>
const inline T&
FixedMatrixBase<T,SIZE1,SIZE2>::data(const int i1,
				     const int i2) const
{
 assert(i1>=0 && i1<SIZE1);
 assert(i2>=0 && i2<SIZE2);

 return (d_[i1 + SIZE1*i2]); //column-major
} /* data const */

template<class T,int SIZE1,int SIZE2,int OFFSET1,int OFFSET2>
inline T&
FixedMatrix<T,SIZE1,SIZE2,OFFSET1,OFFSET2>::operator()(const int i1,
						       const int i2)
{
 return data(i1-OFFSET1,i2-OFFSET2);
} /* operator() */

template<class T,int SIZE1,int SIZE2,int OFFSET1,int OFFSET2>
const inline T&
FixedMatrix<T,SIZE1,SIZE2,OFFSET1,OFFSET2>::operator()(const int i1,
						       const int i2) const
{
 return data(i1-OFFSET1,i2-OFFSET2);
} /* operator() */

//generate a 2x2 block matrix
template<class T, int m1, int m3, int n1, int n2>
inline FixedMatrix<T,m1+m3,n1+n2>
blockMat22(const FixedMatrix<T,m1,n1>& m11,
	       const FixedMatrix<T,m1,n2>& m12,
	       const FixedMatrix<T,m3,n1>& m21,
	       const FixedMatrix<T,m3,n2>& m22)
{
 FixedMatrix<T,m1+m3,n1+n2> ret;
 for (int i=0 ; i<m1 ; i++)
   for (int j=0 ; j<n1 ; j++)
     ret(i,j) = m11(i,j);
 for (l_int i=0 ; i<m1 ; i++)
   for (int j=0 ; j<n2 ; j++)
     ret(i,j+n1) = m12(i,j);
 for (l_int i=0 ; i<m3 ; i++)
   for (int j=0 ; j<n1 ; j++)
     ret(i+m1,j) = m21(i,j);
 for (l_int i=0 ; i<m3 ; i++)
   for (int j=0 ; j<n2 ; j++)
     ret(i+m1,j+n1) = m22(i,j);
 return ret;
} /* blockMat22 */

//template<class T, int m1, int m3, int n1, int n2>
//FixedMatrix<T,m1+m3,n1+n2>
//blockMat22(const FixedMatrix<T,m1,n1>& m11,
//	     const FixedMatrix<T,m1,n2>& m12,
//	     const FixedMatrix<T,m3,n1>& m21,
//	     const FixedMatrix<T,m3,n2>& m22);

//generate a 1x2 block matrix
template<class T, int m1, int m2, int m3>
FixedMatrix<T,m1,m2+m3>
blockMat12(const FixedMatrix<T,m1,m2>&, 
	       const FixedMatrix<T,m1,m3>&);

//generate a 2x1 block matrix
template<class T, int m1, int n1, int m2>
FixedMatrix<T,m1+m2,n1>
blockMat21(const FixedMatrix<T,m1,n1>&, 
	       const FixedMatrix<T,m2,n1>&);

////build a matrix by placing row-vectors on top of one-another
//template<class T, int m1, int m2, int m3>
//FixedMatrix<T,m1,m2+m3>
//catRow(const FixedVector<T,d1>&, 
//	 const FixedVector<T,d1>&);


//template<class T,int size> inline
//FixedMatrix<T,size>
//operator-(const FixedMatrix<T,size> &v1,
//	    const FixedMatrix<T,size> &v2)
//{
// FixedMatrix<T,size> r;
// for (int i=0 ; i<size ; i++)
//   r[i] = v1[i] - v2[i];
// return r;
//} /* operator- */
//

template<class T,int s1,int s2> inline             //unary minus
FixedMatrix<T,s1,s2>
operator-(const FixedMatrix<T,s1,s2> &m)
{
 FixedMatrix<T,s1,s2> r(m);
 r.scale(-1);
 return r;
}


template<class T,int s1,int s2,int s3> 
inline FixedMatrix<T,s1,s3>
operator*(const FixedMatrix<T,s1,s2>& m1,
	  const FixedMatrix<T,s2,s3>& m2)
{
 FixedMatrix<T,s1,s3> r((T)0);
 for (int i=0 ; i<s1 ; i++)
   for (int k=0 ; k<s3 ; k++)
     for (int j=0 ; j<s2 ; j++) 
       r(i,k) += m1(i,j) * m2(j,k);

 return r;
}

//
// perform generalized orthogonal transform on m: S * m * transpose(S)
//
template<class T,int s1,int s2> 
inline FixedMatrix<T,s2,s2>
orthoTransform(const FixedMatrix<T,s1,s1>& m,
	           const FixedMatrix<T,s2,s1>& S)
{
 // return S * m * transpose(S);
// FixedMatrix<T,s2,s1> dum((T)0);
// for (int i=0 ; i<s2 ; i++) 
//   for (int k=0 ; k<s1 ; k++) 
//     for (int j=0 ; j<s1 ; j++) 
//	 dum(i,k) += S(i,j) * m(j,k);

 FixedMatrix<T,s2,s1> dum = S * m;
 FixedMatrix<T,s2,s2> ret((T)0);
 for (int i=0 ; i<s2 ; i++) 
   for (int k=0 ; k<s2 ; k++) 
     for (int j=0 ; j<s1 ; j++) 
       ret(i,k) += dum(i,j) * S(k,j);
 
// for (int i=0 ; i<s2 ; i++)            slower?
//   for (int l=0 ; l<s2 ; l++) {
//     for (int j=0 ; j<s1 ; j++) 
//	 for (int k=0 ; k<s1 ; k++) 
//	   ret(i,l) += S(i,j) * m(j,k) * S(l,k);
//   }
 return ret;
}

template<class T,int s1,int s2> 
inline FixedMatrixBase<T,s1,s2>&
FixedMatrixBase<T,s1,s2>::scale(const T& x)
{
 for (int i=0 ; i<s1*s2 ; i++)
   d_[i] *= x;
 return *this;
}


template<class T,int s1,int s2> 
inline FixedMatrix<T,s1,s2>
operator*(const T&                    x,
	      const FixedMatrix<T,s1,s2>& m)
{
 FixedMatrix<T,s1,s2> r(m);
 r.scale(x);
 return r;
}


template<class T, int s1, int s2>
inline void
FixedMatrixBase<T,s1,s2>::set(const T &x)
{
 for (int i=0 ; i<s1*s2 ; i++)
   d_[i] = x;
}


template<class T, int s1, int s2>
inline void
FixedMatrixBase<T,s1,s2>::setDiag(const T &x)
{
 for (int i=0 ; i<((s1<s2)?s1:s2) ; i++)
   d_[i + s1*i] = x;
}

//template<class T, int s1, int s2>
//inline FixedMatrix<T,s2,s1>
//transpose(const FixedMatrix<T,s1,s2> &m)
//{
// FixedMatrix<T,s2,s1> r;
// for (int i=0 ; i<s1 ; i++)
//   for (int j=0 ; j<s2 ; j++) 
//     r(j,i) = m(i,j);
// return r;
//} /* transpose */


template<class T,int s1,int s2> 
inline FixedVector<T,s1>
operator*(const FixedMatrix<T,s1,s2>& m,
	  const FixedVector<T,s2>&    v)
{
 FixedVector<T,s1> r;
 r.set( (T)0 );
 
 for (int i=0 ; i<s1 ; i++) 
   for (int j=0 ; j<s2 ; j++) 
     r(i) += m(i,j) * v(j);
 
 return r;
} /* operator* (matrix, vector) */

template<class T,int s1,int s2> 
inline FixedVector<T,s2>
operator*(const FixedVector<T,s1>&    v,
	  const FixedMatrix<T,s1,s2>& m)   // v^T * m
{
 FixedVector<T,s2> r;
 r.set( (T)0 );
 
 for (int j=0 ; j<s2 ; j++) 
   for (int i=0 ; i<s1 ; i++) 
     r(j) += v(i) * m(i,j);
 
 return r;
} /* operator* (vector, matrix) */

template<class T,int S1,int S2> 
inline FixedMatrix<T,S1,S2>
operator+(const FixedMatrix<T,S1,S2>& m1,
	  const FixedMatrix<T,S1,S2>& m2)
{
 FixedMatrix<T,S1,S2> ret = m1;
 ret += m2;
 return ret;
} /* operator+ (matrix, matrix) */

template<class T,int S1,int S2> 
inline FixedMatrix<T,S1,S2>
operator-(const FixedMatrix<T,S1,S2>& m1,
	  const FixedMatrix<T,S1,S2>& m2)
{
 FixedMatrix<T,S1,S2> ret = m1;
 ret -= m2;
 return ret;
} /* operator- (matrix, matrix) */


//template<class T,int s1> 
//inline 
//FixedVector<T,s1>
//operator+(const FixedVector<T,s1>& v1,
//	    const FixedVector<T,s1>& v2)
//{
// FixedVector<T,s1> r= v1;
// r += v2;
// return r;
//} /* operator+ (vector, vector) */


//template<class T,int size1,int size2,int size3> inline
//FixedMatrix<T,size1,size3>
//operator*(const FixedMatrix<T,size1,size2>& m1,
//	    const FixedMatrix<T,size2,size3>& m2)
//{
// FixedMatrix<T,size1,size3> r;
// for (int i=0 ; i<size1 ; i++)
//   for (int j=0 ; j<size3 ; j++)
//     for (int k=0 ; j<size2 ; j++)
//	 r(i,j) = m1(i,k) * m2(k,j);
// return r;
//} /* operator* */


int FixedMatrix_test();

#endif /* __fixedMatrixDecl_h__ */
