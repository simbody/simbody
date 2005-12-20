/*
 *  CDSVector.hh
 *  headers for my own GenericVector template class
 *
 *  This is for big arrays which can dynamically change size. For small 
 *  arrays of fixed size, use FixedVector.
 *
 *  10/1/99 - removed CDSArray class- only really need one
 *  10/8/99 - changed offset to a template parameter
 *  10/28/99- first shot at reference counting
 */

#ifndef __CDSVectorDecl_h__
#define __CDSVectorDecl_h__ 1

//#define CDSVECTOR_DEBUG_ALLOC 1

#include <cdsAlloc.h>
#include <cdsMath.h>
#include "cdsGenericVector.h"

template<class T, class ALLOC>             //forward declarations
class CDSVectorBase;
template<class T, int offset, class ALLOC>
class CDSVector; 

#include <cdsIostream.h>

#include <cassert>
#include <cstddef>  //for size_t

template<class T, class ALLOC>
class CDSVectorRep {
  int size;
  T* data;
  int count;
#ifdef CDSVECTOR_DEBUG_ALLOC
  static T forDiags;
#endif /* CDSVECTOR_DEBUG_ALLOC */
public:
  CDSVectorRep(const int s) : 
    size(s), data((T*)ALLOC::alloc(size*sizeof(T))), count(1)
  { for (int i=0 ; i<size ; i++) new(data+i) T(); }

  CDSVectorRep(const CDSVectorRep* r) :
    size(r->size), data((T*)ALLOC::alloc(size*sizeof(T))), count(1)
  { for (int i=0 ; i<size ; i++) new(data+i) T(r->data[i]); }

  ~CDSVectorRep() 
  { for (int i=0 ; i<size ; i++) data[i].~T(); ALLOC::free(data); }

  friend class CDSVectorBase<T,ALLOC>;
#ifdef CDSVECTOR_DEBUG_ALLOC
  void* operator new(size_t);
  void  operator delete(void*);
  void* operator new(size_t,void* location) { return location; }
  void  operator delete(void*, void*) { };
  void* operator new[](size_t);
  void  operator delete[](void*);
#endif /* CDSVECTOR_DEBUG_ALLOC */
};

template<class T, class ALLOC=CDS::DefaultAlloc>
class CDSVectorBase {

protected:
  CDSVectorRep<T,ALLOC>* rep;

  CDSVectorBase();
  CDSVectorBase(const int size);
  CDSVectorBase(const int size,const T& i); //initialize w/ i
  inline void splitRep() 
    { if (rep->count>1) {
      rep->count--;rep= new CDSVectorRep<T,ALLOC>(rep);}}
public:
  CDSVectorBase(const CDSVectorBase<T,ALLOC> &v); // this should be private:
                                            // public for blockVec: friendship?
  template<class VEC>
  explicit CDSVectorBase(const VEC& v); //generic constructor from vector
  ~CDSVectorBase();
  
  CDSVectorBase<T,ALLOC> &resize(const int size); 
  int size()   const {return rep->size;}

  T&       updData(int i)       {return rep->data[i];}
  const T& getData(int i) const {return rep->data[i];}

  const T* pointer() const { return rep->data; }
  T* pointer() { return rep->data; }

  //copy array into an Array
  inline CDSVectorBase<T,ALLOC> &copy(      T const* array,
				const int      size);

  //void     AllocateMem(); //??
  CDSVectorBase<T,ALLOC>  &operator=(const CDSVectorBase<T,ALLOC>&);
  void set(const T&);

  //
  // convert to CDSVector
//  template<int offset> 
//  operator CDSVector<T,offset,ALLOC>&() 
//    { void* p=this; return *(CDSVector<T,offset,ALLOC>*)p; }

//  template<int offset> operator CDSVector<T,offset>() 
//    { void* p=this; return *(CDSVector<T,offset>*)p; }
//  template<int offset> operator const CDSVector<T,offset>() const
//    { const void* p=this; return *(CDSVector<T,offset>*)p; }

  const T& operator()(const int i) const { return rep->data[i]; }

  //
  // math operations
  //  const FixedVectorBase& operator*=(const T &);
  template <class T1> 
  CDSVectorBase<T,ALLOC> &operator*=(const T1 &x);
  //
  CDSVectorBase<T,ALLOC>& operator+=(const CDSVectorBase<T,ALLOC> &v);
  CDSVectorBase<T,ALLOC>& operator-=(const CDSVectorBase<T,ALLOC> &v);

#ifdef CDSVECTOR_DEBUG_ALLOC 
  static T forDiags;
  //#endif /* CDSVECTOR_DEBUG_ALLOC */
  void* operator new(size_t);
  void* operator new(size_t,void* location) { return location; }
  void  operator delete(void*,void*) { }
  void  operator delete(void*);
  void* operator new[](size_t);
  void  operator delete[](void*);
#endif /* CDSVECTOR_DEBUG_ALLOC */
};

template<class T, int offset_=0, class ALLOC=CDS::DefaultAlloc>
class CDSVector : public CDSVectorBase<T,ALLOC> {
public:
  typedef T ElementType;

  CDSVector() : CDSVectorBase<T,ALLOC>() {}
  explicit CDSVector(const int size) : CDSVectorBase<T,ALLOC>(size) {}
  CDSVector(const int size,const T& i) : CDSVectorBase<T,ALLOC>(size,i) {}
  CDSVector(const CDSVectorBase<T,ALLOC> &v) : 
    CDSVectorBase<T,ALLOC>(v) {}
  CDSVector(const CDSVector<T,offset_,ALLOC> &v) : 
    CDSVectorBase<T,ALLOC>(v) {}
  ~CDSVector () {}
  CDSVector<T,offset_,ALLOC> &operator= (const CDSVectorBase<T,ALLOC>&v) 
    { ((CDSVectorBase<T,ALLOC>&) *this)=v; return *this; }

  template<class VEC>  //generic constructor from vector
  explicit CDSVector(const VEC& v) : CDSVectorBase<T,ALLOC>(v) {} 
//  CDSVector<T>  &operator= (const CDSVector<T>&v) 
//    { this->CDSVectorBase<T>::operator=(v); return *this; }

  // These pass-throughs to the base class were needed to get this to go through gcc 3.4.4 (sherm).
  void     splitRep()           {CDSVectorBase<T,ALLOC>::splitRep();}
  int      size()         const {return CDSVectorBase<T,ALLOC>::size();}
  T&       updData(int i)       {return CDSVectorBase<T,ALLOC>::updData(i);}
  const T& getData(int i) const {return CDSVectorBase<T,ALLOC>::getData(i);}

  //assignment from generic GenericVector
  template<class VEC>
  CDSVector<T,offset_,ALLOC> &operator=(const CDS::GenericVector<VEC>& v);

  //convert to generic vector
  CDS::GenericVector<CDSVector<T,offset_,ALLOC> > vector()
  { return CDS::GenericVector<CDSVector<T,offset_,ALLOC> >(*this); }

  //members not present in CDSVectorBase
  int offset() const {return offset_;}
  inline T  &operator() (const int &index);
  inline T  &operator[] (const int &index) { return (*this)(index); }
  inline const T  &operator() (const int &index) const;
  inline const T  &operator[] (const int &index) const 
  { return (*this)(index); }

//  // SHOULDN'T BE NECESSARY
//  // convert to CDSVector with different offset
//  template<int offset> operator CDSVectorBase<T>() 
//    { void* p=this; return *(CDSVector<T,offset>*)p; }
//  template<int offset> operator CDSVector<T,offset> const&() 
//    { void* p=this; return *(CDSVector<T,offset>*)p; }
};

class stream;
template<class T,class ALLOC> 
ostream& operator<<(ostream&,const CDSVectorBase<T,ALLOC> &);

template<class T,class ALLOC>  inline 
CDSVectorBase<T,ALLOC>::CDSVectorBase() 
{
 rep = new CDSVectorRep<T,ALLOC>(0);
}

template<class T,class ALLOC>  inline 
CDSVectorBase<T,ALLOC>& 
CDSVectorBase<T,ALLOC>::copy(      T const*   array,
		       const int   	size)
{
 splitRep();
 resize(size);
 for (int i=0 ; i<size ; i++)
   rep->data[i] = array[i];
 return *this;
} /* copy */


template<class T,int offset_,class ALLOC> inline
T& 
CDSVector<T,offset_,ALLOC>::operator()(const int &index)
{
 assert( index>=offset_ && index<offset_+size() );

 splitRep();
 
 return updData( index-offset_ );
}

template<class T,int offset_,class ALLOC> inline
const T& 
CDSVector<T,offset_,ALLOC>::operator()(const int &index) const
{
 assert( index>=offset_ && index<offset_+size() );

 return getData( index-offset_ );
}

template<class T,class ALLOC>
template<class T1>
CDSVectorBase<T,ALLOC>& 
CDSVectorBase<T,ALLOC>::operator*=(const T1 &x)
{
 splitRep();
 T *p = rep->data;

 for (int i=0 ; i<size() ; i++)
   *p++ *= x;
// for (int i=0 ; i<size() ; i++)
//   p[i] *= x;
 return *this;
} /* *= */

//template<class T, int moff1, int moff2>
//class CDSMatrix;
//
//template<class T, int moff1, int moff2>
//CDSVector<T>
//subCol(const CDSMatrix<T,moff1,moff2>& m,
//	       int col,
//	 const int o1,
//	 const int o2);

template<class T,class ALLOC>
CDSVector<T,0,ALLOC>
blockVec(const CDSVectorBase<T,ALLOC>& v1,
	 const CDSVectorBase<T,ALLOC>& v2);

//template<class T,int offset>
//CDSVector<T>
//subVec(const CDSVector<T,offset>& v1,
//	 const int o1,
//	 const int o2)
//  //
//  // takes CDSVector as argument because the indices specifying range depend
//  // on the offset.
//  //
//{ 
// CDSVector<T> ret(o2-o1+1);
// for (int i=o1 ; i<=o2 ; i++)
//   ret(i-o1) = v1(i);
// return ret;
//} /* blockVec */

template<class T,class ALLOC>
CDSVectorBase<T,ALLOC>
operator+(const CDSVectorBase<T,ALLOC>& v1,
	  const CDSVectorBase<T,ALLOC>& v2)
{
 return CDSVectorBase<T,ALLOC>(v1) += v2;
}

template<class T,class ALLOC>
CDSVectorBase<T,ALLOC>
operator-(const CDSVectorBase<T,ALLOC>& v1,
	  const CDSVectorBase<T,ALLOC>& v2)
{
 return CDSVectorBase<T,ALLOC>(v1) -= v2;
}


template<class T,class ALLOC>
CDSVectorBase<T,ALLOC>
operator-(const CDSVectorBase<T,ALLOC>& v)  //unary -
{
 return CDSVectorBase<T,ALLOC>(v) *= -1;
} /* operator- (unary) */

//template<class T,class ALLOC>
//T
//operator*(const CDSVectorBase<T>& v1,
//	    const CDSVectorBase<T>& v2)
//{ return dot(v1,v2); }

template<class T,class ALLOC>
CDSVector<T,0,ALLOC>
operator*(const T&               s,
	  const CDSVectorBase<T,ALLOC>& v)
{
 return CDSVectorBase<T,ALLOC>(v) *= s;
}

//template<class T> inline
//double
//dot(const T& v1,
//    const T& v2) 
//{ return v1*v2;
//} /* dot */

template<class T,class ALLOC> inline
double
dot(const CDSVectorBase<T,ALLOC>& v1,
    const CDSVectorBase<T,ALLOC>& v2) 
{ double ret=0; 
 for (int i=0;i<v1.size();i++) ret += v1(i) * v2(i);
 return ret;
} /* dot */

template<class T,class ALLOC> inline
double
abs2(const CDSVectorBase<T,ALLOC>& v) { return dot(v,v); }

#endif /* __CDSVectorDecl_h__ */
