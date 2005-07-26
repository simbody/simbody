#ifndef __CDSVectorDefs_h__
#define __CDSVectorDefs_h__ 1

//
//constructor
//

template<class T,class ALLOC> 
CDSVectorBase<T,ALLOC>::CDSVectorBase(const int size) 
{
 rep = new CDSVectorRep<T,ALLOC>(size);
}

template<class T,class ALLOC> 
CDSVectorBase<T,ALLOC>::CDSVectorBase(const int size,
				      const T&  init)
{
 rep = new CDSVectorRep<T,ALLOC>(size);
 for (int i=0 ; i<size ; i++)
   rep->data[i] = init;
}

//
//copy contructor
//

template<class T,class ALLOC>
CDSVectorBase<T,ALLOC>::CDSVectorBase(const CDSVectorBase<T,ALLOC> &v)
{
 rep = v.rep;
 rep->count++;
}

//
//generic contructor
//

template<class T,class ALLOC>
template<class VEC>
CDSVectorBase<T,ALLOC>::CDSVectorBase(const VEC &v)
{
 rep = new CDSVectorRep<T,ALLOC>(v.size());
 for (int i=0 ; i<v.size() ; i++)
   rep->data[i] = v(v.offset()+i);
}

//
//destructor
//

template<class T,class ALLOC> 
CDSVectorBase<T,ALLOC>::~CDSVectorBase()
{
 if (--rep->count <= 0) delete rep;
}

template<class T,class ALLOC>
CDSVectorBase<T,ALLOC>& 
CDSVectorBase<T,ALLOC>::resize(const int size)
{
 if (size == rep->size) 
   return *this;

 CDSVectorRep<T,ALLOC>* oldRep = rep;

 rep = new CDSVectorRep<T,ALLOC>(size);
 for (int i=0 ; i<CDSMath::min(size,oldRep->size) ; i++)
   rep->data[i] = oldRep->data[i];

 if (--oldRep->count <= 0) delete oldRep;
 
 return *this;
} /*resize*/

template<class T,class ALLOC> 
CDSVectorBase<T,ALLOC>& 
CDSVectorBase<T,ALLOC>::operator=(const CDSVectorBase<T,ALLOC> &x)
{
 if (&x == this)
   return *this;

 x.rep->count++; 
 if (--rep->count <= 0) delete rep;

 rep = x.rep;

 return *this;
}

template<class T,class ALLOC> 
void
CDSVectorBase<T,ALLOC>::set(const T &x)
{
 splitRep();
 for (int i=0 ; i<rep->size ; i++)
   rep->data[i] = x;
} /* set */


//
//operator+= (vector)
//

template<class T,class ALLOC>
CDSVectorBase<T,ALLOC> &CDSVectorBase<T,ALLOC>::operator+=(const CDSVectorBase<T,ALLOC> &v)
{
 assert( v.rep->size == rep->size );

 splitRep();
  
 for (int i=0 ; i<rep->size ; i++)
   rep->data[i] += v.rep->data[i];
 return *this;
}

template<class T,class ALLOC>
CDSVectorBase<T,ALLOC> &CDSVectorBase<T,ALLOC>::operator-=(const CDSVectorBase<T,ALLOC> &v)
{
 assert( v.rep->size == rep->size );

 splitRep();
  
 for (int i=0 ; i<rep->size ; i++)
   rep->data[i] -= v.rep->data[i];
 return *this;
}

template<class T,int OFFSET,class ALLOC> 
template<class VEC>
CDSVector<T,OFFSET,ALLOC>& 
CDSVector<T,OFFSET,ALLOC>::operator=(const CDS::Vector<VEC>& v)
{
 resize( v.size() );

 for (int i=0 ; i<v.size() ; i++)
   (*this)[i] = v(i);

 return *this;
} /* assignment from generic vector */


template<class T,class ALLOC>
CDSVector<T,0,ALLOC>
blockVec(const CDSVectorBase<T,ALLOC>& v1,
	 const CDSVectorBase<T,ALLOC>& v2)
{ 
 CDSVector<T,0,ALLOC> ret( v1.size() + v2.size() );
// const CDSVector<T,ALLOC>& v1r = v1;
// const CDSVector<T,ALLOC>& v2r = v2;
 for (int i=0 ; i<v1.size() ; i++)
   ret(i) = v1(i);
 for (l_int i=0 ; i<v2.size() ; i++)
   ret(v1.size()+i) = v2(i);
 return ret;
} /* blockVec */

template<class T,class ALLOC>
ostream&   
operator<<(ostream& s, const CDSVectorBase<T,ALLOC>& x)
{
 s << "{ ";
 if (x.size()) {
   s << x(0);
   for (int i=1 ; i<x.size() ; i++)
     s << ", " << x(i) ;
 }
 s << " }";
 return s;
}

#ifdef CDSVECTOR_DEBUG_ALLOC 
#include <typeinfo.h>

template<class T,class ALLOC> T CDSVectorRep<T,ALLOC>::forDiags;
template<class T,class ALLOC> T CDSVectorBase<T,ALLOC>::forDiags;

template<class T,class ALLOC>
void*
CDSVectorBase<T,ALLOC>::operator new(size_t s)
{
 cout << "CDSVectorBase<" << typeid(forDiags).name() 
      << ">: operator new allocating " << s << " bytes\n";
 return CDSAlloc::new(s);
} /* operator new */

template<class T,class ALLOC>
void
CDSVectorBase<T,ALLOC>::operator delete(void* p)
{
 cout << "CDSVectorBase<" << typeid(forDiags).name() 
      << ">: operator delete freeing " << p << "\n";
 CDSAlloc::delete p;
} /* operator delete */

template<class T,class ALLOC>
void*
CDSVectorBase<T,ALLOC>::operator new[](size_t s)
{
 cout << "CDSVectorBase<" << typeid(forDiags).name() 
      << ">: operator new[] allocating " << s << " bytes\n";
 return CDSAlloc::operator new[](s);
} /* operator new[] */

template<class T,class ALLOC>
void
CDSVectorBase<T,ALLOC>::operator delete[](void* p)
{
 cout << "CDSVectorBase<" << typeid(forDiags).name() 
      << ">: operator delete[] freeing " << p << "\n";
 CDSAlloc::delete [] p;
} /* operator delete[] */

template<class T,class ALLOC>
void*
CDSVectorRep<T,ALLOC>::operator new(size_t s)
{
 void* p = CDSAlloc::operator new(s);
 cout << "CDSVectorBase<" << typeid(forDiags).name() 
      << ">: operator new allocating " << s << " bytes: " << p << '\n';
 return p;
} /* operator new */

template<class T,class ALLOC>
void
CDSVectorRep<T,ALLOC>::operator delete(void* p)
{
 cout << "CDSVectorBase<" << typeid(forDiags).name() 
      << ">: operator delete freeing " << p << "\n";
 CDSAlloc::delete p;
} /* operator delete */

template<class T,class ALLOC>
void*
CDSVectorRep<T,ALLOC>::operator new[](size_t s)
{
 cout << "CDSVectorBase<" << typeid(forDiags).name() 
      << ">: operator new[] allocating " << s << " bytes\n";
 return CDSAlloc::operator new[](s);
} /* operator new[] */

template<class T,class ALLOC>
void
CDSVectorRep<T,ALLOC>::operator delete[](void* p)
{
 cout << "CDSVectorBase<" << typeid(forDiags).name() 
      << ">: operator delete[] freeing " << p << "\n";
 CDSAlloc::delete [] p;
} /* operator delete[] */

#else /* not debug */

//template<class T, class ALLOC>
//void*
//CDSVectorBase<T,ALLOC>::operator new(size_t s)
//{
// return ALLOC::alloc(s);
//} /* operator new */
//
//template<class T,class ALLOC>
//void
//CDSVectorBase<T,ALLOC>::operator delete(void* p)
//{
// ALLOC::free(p);
//} /* operator delete */
//
//template<class T,class ALLOC>
//void*
//CDSVectorBase<T,ALLOC>::operator new[](size_t s)
//{
// return ALLOC::alloc(s);
//} /* operator new[] */
//
//template<class T,class ALLOC>
//void
//CDSVectorBase<T,ALLOC>::operator delete[](void* p)
//{
// ALLOC::free(p);
//} /* operator delete[] */
//
//template<class T,class ALLOC>
//void*
//CDSVectorRep<T,ALLOC>::operator new(size_t s)
//{
// return ALLOC::alloc(s);
//} /* operator new */
//
//template<class T,class ALLOC>
//void
//CDSVectorRep<T,ALLOC>::operator delete(void* p)
//{
// ALLOC::free(p);
//} /* operator delete */
//
//template<class T,class ALLOC>
//void*
//CDSVectorRep<T,ALLOC>::operator new[](size_t s)
//{
// return ALLOC::alloc(s);
//} /* operator new[] */
//
//template<class T,class ALLOC>
//void
//CDSVectorRep<T,ALLOC>::operator delete[](void* p)
//{
// ALLOC::free(p);
//} /* operator delete[] */

#endif /* CDSVECTOR_DEBUG_ALLOC */

int CDSVector_test();

#endif /* __CDSVectorDefs_h__ */
