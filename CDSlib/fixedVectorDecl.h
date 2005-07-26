#ifndef __fixedVectorDecl_h__
#define __fixedVectorDecl_h__ 1

//
// Numeric Vector class with fixed dimension with local storage
//
//

#include <sthead.h>

#include <cdsIostream.h>

#include <assert.h>
#include <cdsExcept.h>
#include <vector.h>

// This gives the proper declaration of FixedMatrix with default
// template parameters.
#include "fixedMatrixFwd.h"

template<class T,int size,int OFFSET>
class FixedVector;


template<class T,int SIZE>
class FixedVectorBase {

protected:

  T d_[SIZE];

public:
  typedef T ElementType;

  FixedVectorBase() {}
  explicit FixedVectorBase(const T &x) 
    { for (int i=0 ; i<SIZE ; i++) d_[i] = x; }
  explicit FixedVectorBase(const T a[SIZE])
    { for (int i=0 ; i<SIZE ; i++) d_[i] = a[i]; }
  FixedVectorBase(const FixedVectorBase<T,SIZE>& v) 
    {for (int i=0 ; i<SIZE ; i++) d_[i]=v.d_[i]; }
  //  template<int OFFSET>
  //  FixedVectorBase(const CDSVector<T,OFFSET>& v);
  template<class VEC>
  FixedVectorBase(const CDS::Vector<VEC>& v) //generic constructor from vector
  { 
   assert(SIZE==v.size()); 
   for (int i=0 ; i<SIZE ; i++) d_[i]=v(v.offset()+i); 
  }
  //
  int size() const  { return SIZE;   }
  //  
  const T& operator()(const int i) const { return d_[i]; }
  //
  FixedVectorBase<T,SIZE>& operator=(const FixedVectorBase<T,SIZE>& v)
  { for (int i=0 ; i<SIZE ; i++) d_[i] = v.d_[i]; return *this; }
  template<class VEC>
  FixedVectorBase<T,SIZE>& operator=(const CDS::Vector<VEC>&);
  //
  inline const FixedVectorBase& operator+=(const FixedVectorBase&);
  const FixedVectorBase& operator-=(const FixedVectorBase&);
  const FixedVectorBase& operator/=(const T &);
  template<class T2>
  inline const FixedVectorBase& operator*=(const T2 &);
  //
  void copyArray(const T* a);
  void set(const T& x) { *this = FixedVectorBase<T,SIZE>(x); }

  // access to the raw array
  const T* getData() const { return d_; }
  T*       updData()       { return d_; }

//  template<int OFFSET> operator FixedVector<T,SIZE,OFFSET>&() 
//    { void* p=this; return *(FixedVector<T,SIZE,OFFSET>*)p; }
//  template<int OFFSET> operator FixedVector<T,SIZE,OFFSET> const&() 
//    { void* p=this; return *(FixedVector<T,SIZE,OFFSET>*)p; }

  //

  class RangeError : public CDS_NAMESPACE(out_of_range) {};
  class IOError    : public CDS_NAMESPACE(exception)    {};  //FIX: should inherit from ios_base::failure
  //
};

template<class T,int SIZE,int OFFSET=0>
class FixedVector : public FixedVectorBase<T,SIZE> {

public:
  FixedVector() : FixedVectorBase<T,SIZE>() {}
  explicit FixedVector(const T &x) : FixedVectorBase<T,SIZE>(x) {}
  //  explicit FixedVector(const T a[SIZE]) : FixedVectorBase<T,SIZE>(a) {}
  explicit FixedVector(const T a[SIZE]) : FixedVectorBase<T,SIZE>(a) {}
  FixedVector(const FixedVector<T,SIZE>& v) : 
    FixedVectorBase<T,SIZE>(v) {}
  FixedVector(const FixedVectorBase<T,SIZE>& v) : 
    FixedVectorBase<T,SIZE>(v) {}
  template<class VECTOR>
  FixedVector(const CDS::Vector<VECTOR>& v) :
    FixedVectorBase<T,SIZE>(v) {}

  // access to the raw array
  const T* getData() const { return FixedVectorBase<T,SIZE>::getData(); }
  T*       updData()       { return FixedVectorBase<T,SIZE>::updData(); }

  // conversion to generic vector
  CDS::Vector<FixedVector<T,SIZE,OFFSET> > vector()
  { return CDS::Vector<FixedVector<T,SIZE,OFFSET> >(*this); }
  const CDS::Vector<FixedVector<T,SIZE,OFFSET> > vector() const
  { return CDS::Vector<FixedVector<T,SIZE,OFFSET> >(*this); }

  int offset() const { return OFFSET; }
  //
  T& operator()(const int&);
  T& operator[](const int &index) { return (*this)(index); }
  const T& operator()(const int&) const;
  const T& operator[](const int &index) const { return (*this)(index); }
  //
  FixedVector<T,SIZE,OFFSET>& operator=(const FixedVector<T,SIZE,OFFSET>& v)
  { FixedVectorBase<T,SIZE>::operator=(v); return *this; }

  template<int msize1,int msize2,int moff1, int moff2>
  static FixedVector subCol(const FixedMatrix<T,msize1,msize2,moff1,moff2>& m,
                  int col,
                const int o1,
                const int o2);

};

template<class T,int SIZE> 
ostream& operator<<(ostream&,const FixedVectorBase<T,SIZE> &);


template<class T,int SIZE,int OFFSET> inline
T& 
FixedVector<T,SIZE,OFFSET>::operator()(const int &index)
{
 assert( index>=OFFSET );
 assert( index<SIZE+OFFSET );

 return (updData()[index-OFFSET]);
}

template<class T,int SIZE,int OFFSET> inline
const T& 
FixedVector<T,SIZE,OFFSET>::operator()(const int &index) const
{
 assert( index>=OFFSET );
 assert( index<SIZE+OFFSET );

 return (getData()[index-OFFSET]);
}


template<class T,int SIZE1,int SIZE2>
FixedVector<T,SIZE1+SIZE2>
blockVec(const FixedVector<T,SIZE1>& v1,
     const FixedVector<T,SIZE2>& v2);


template<class T,int SIZE> inline
bool
operator==(const FixedVectorBase<T,SIZE> &v1,
       const FixedVectorBase<T,SIZE> &v2)
{
  for (int i=0 ; i<SIZE ; i++)
   if ( !( v1(i) == v2(i) ) )
 return 0;
 return 1;
} /* operator== */

template<class T,int SIZE> inline
bool
operator!=(const FixedVectorBase<T,SIZE> &v1,
       const FixedVectorBase<T,SIZE> &v2)
{
 return !(v1==v2);
} /* operator!= */

template<class T,int SIZE> inline
FixedVector<T,SIZE>
operator-(const FixedVector<T,SIZE> &v1,
      const FixedVector<T,SIZE> &v2)
{
 FixedVector<T,SIZE> r;
 for (int i=0 ; i<SIZE ; i++)
   r(i) = v1(i) - v2(i);
 return r;
} /* operator- */

template<class T,int SIZE> inline 
FixedVector<T,SIZE>
operator-(const FixedVectorBase<T,SIZE>& v)  //unary -
{ 
 FixedVectorBase<T,SIZE> ret = v; //; ret = v; 
 ret *= -1;
 return ret;
} /* operator- (unary) */

template<class T,int SIZE> inline
FixedVector<T,SIZE>
operator*(const T&                   x,
      const FixedVector<T,SIZE>& v)
{
 FixedVector<T,SIZE> r;
 for (int i=0 ; i<SIZE ; i++)
   r(i) = x * v(i);
 return r;
}

template<class T,int SIZE> inline
FixedVector<T,SIZE>
operator*(const FixedVector<T,SIZE>& v,
      const T&                   x)
      
{
 return x * v;
}

template<class T,int SIZE> inline
FixedVector<T,SIZE>
operator/(const FixedVector<T,SIZE>& v,
      const T&                   x)
{
 FixedVector<T,SIZE> r;
 for (int i=0 ; i<SIZE ; i++)
   r(i) = v(i) / x;;
 return r;
} /* operator/ */

template<class T,int SIZE> inline
T
abs2(const FixedVector<T,SIZE>& v)
{
 //onst FixedVector<T,SIZE>& vr = v;
 T ret = v(0)*v(0);
 for (int i=1 ; i<SIZE ; i++)
   ret += v(i)*v(i);
 return ret;
} /* abs2 */ 

template<class T,int SIZE> inline
T
dot(const FixedVector<T,SIZE>& v1,
    const FixedVector<T,SIZE>& v2) 
{ T ret=0; 
 for (int i=0;i<SIZE;i++) ret += v1(i)*v2(i);
 return ret;
} /* dot */

template<class T,int SIZE> inline
T
operator*(const FixedVector<T,SIZE>& v1,
      const FixedVector<T,SIZE>& v2)
{ return dot(v1,v2); }

template<class T,int SIZE> inline
FixedVector<T,SIZE>
operator+(const FixedVector<T,SIZE>& v1,
      const FixedVector<T,SIZE>& v2)
{ FixedVector<T,SIZE> ret(v1); ret += v2; return ret; }

template<class T,int SIZE>
T
norm(const FixedVector<T,SIZE>& v);

template<class T, int SIZE>
template<class T2>
inline const FixedVectorBase<T,SIZE>&
FixedVectorBase<T,SIZE>::operator*=(const T2 &x)
  //
  // return v * x
  //
{
 for (int i=0 ; i<SIZE ; i++)
   d_[i] *= x;

  return *this;
} /* operator*= (T2) */

template<class T, int SIZE>
inline const FixedVectorBase<T,SIZE>&
FixedVectorBase<T,SIZE>::operator+=(const FixedVectorBase<T,SIZE>& v)
{
 for (int i=0 ; i<SIZE ; i++)
   d_[i] += v.d_[i];

 return *this;
} /* operator+= */

//template<class T, int SIZE, int OFFSET>
//template<int bsize, int boffset>
//inline FixedVector<T,SIZE,OFFSET>
//FixedVector<T,SIZE,OFFSET>::subVec(const FixedVector<T,bsize,boffset>& v1,
//                     const int o1,
//                     const int o2)
//  //
//  // takes CDSVector as argument because the indices specifying range depend
//  // on the offset.
//  //
//{ 
// FixedVector<T,SIZE,OFFSET> ret;
// for (int i=o1 ; i<=o2 ; i++)
//   ret(i-o1+OFFSET) = v1(i);
// return ret;
//} /* subVec */

int FixedVector_test();

#endif /* __fixedVectorDecl_h__ */
