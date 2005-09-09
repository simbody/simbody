/*
   CDSList.h
   headers for my own List template class
   9/1/95 - CDS
   see below for docs
*/

#ifndef __cdslist_hh__
#define __cdslist_hh__ 1

#include <sthead.h> 
#include <cdsIostream.h> 
#include <cdsMath.h>
#include <vector.h>

#include <cassert>

template<class T>
class CDSList; 

template<class T>
class CDSListRep {
  CDSListRep();
  CDSListRep(const CDSListRep&);
  void operator=(const CDSListRep&);
public:
  int size;  //number of elements used
  int asize;  //numer of elements allocated

  T*  data;
  int count;
  CDSListRep(const int s,
	     const int a);
  CDSListRep(const int a);  //doesn't call T's constructor
  CDSListRep(const CDSListRep* r);
  ~CDSListRep();

};


template<class T> 
class CDSList {
  CDSListRep<T>* rep;
protected:
  inline void splitRep() 
  { if (rep->count>1) {rep->count--;rep= new CDSListRep<T>(rep);}}
public:


  CDSList (const int s,             //requires default constructor
	   const int a=10) : 
    rep(new CDSListRep<T>(s,a)) {}
  CDSList () :                      //copy constructor required eventually
    rep(new CDSListRep<T>(10)) {}
  
  ~CDSList ();
  
  CDSList (const  CDSList<T>&); //copy constructor

  CDSList &resize(const int s);
  void reserve(const int size);

  CDSList& operator= (const  CDSList<T>&);
  template<class VECTOR>
  CDSList& operator= (const CDS::Vector<VECTOR>&);

  // convert to generic vector
  CDS::Vector<CDSList<T> > vector() 
     { return CDS::Vector<CDSList<T> >(*this); }
  const CDS::Vector<CDSList<T> > vector() const
     { return CDS::Vector<CDSList<T> >(*this); }

  const T  &operator[] (const int &index) const;
        T  &operator[] (const int &index);

  const T* pointer() const { return rep->data; }

  // these so CDSList can be used in a generic vector
  typedef T ElementType;
  int offset() const { return 0; }
  const T& operator()(const int &index) const { return (*this)[index]; }
        T& operator()(const int &index)       { return (*this)[index]; }


  void append(const T member);
  const T& prepend(const T& member);  //slower currently
  //  T& append(T& member);
  bool contains(const T& member) const;
  void remove(const int i);
  int  getIndex(const T& member) const;
  //void setBlockSize(const int);

  T pop();  //remove and return the last item

  int  size() const {return rep->size;}

  void reverse();
  // if for a,b return 1 if a>b 
  //                   0 if a==b
  //                  -1 if a<b
  static int stdComparer(const T& a,
			 const T& b) {return ( a>b?1 : (a==b?0:-1) );}
  template<class Comparer>
  void sort(Comparer comparer);
  template<class Comparer>
  void sort_si(Comparer comparer);
  template<class Comparer>
  void sort_shell(Comparer comparer);
  template<class Comparer>
  void sort_heap(Comparer comparer);
  //  template<class T2> 
//  friend ostream& operator<<(ostream& os, 
//			       const CDSList<T,DefaultSize>& v);

  const CDSListRep<T>* getRep() const { return rep; }

#ifdef DEBUG_ALLOC
  void* operator new(size_t t,void* p) { return p; }
  void* operator new(size_t);
  void  operator delete(void*);
  void* operator new[](size_t);
  void  operator delete[](void*);
#endif /* DEBUG_ALLOC */
};
  
template <class T> 
ostream& 
operator<<(ostream& os, const CDSList<T>& v);

//constructors

template<class T>
CDSListRep<T>::CDSListRep(const int s,
			  const int a) :
  size(s), asize(a), count(1)
{
 if ( size>asize )
   asize += size;

 //extra 1 for roundoff if sizeof(T)<sizeof(void*)
 data = (T*)new void*[ 1+asize*sizeof(T)/sizeof(void*) ]; 

 T* p=(T*)data;
 for (int i=0 ; i<size ; i++,p++)
   new((void*)p) T();
} /* CDSListRep::CDSListRep */

template<class T>
CDSListRep<T>::CDSListRep(const int a) :
  size(0), asize(a), count(1)
{
 data = (T*)new void*[1+asize*sizeof(T)/sizeof(void*)];

} /* CDSListRep::CDSListRep */

template<class T>
CDSListRep<T>::CDSListRep(const CDSListRep* r) : 
  size(r->size), asize(r->asize), count(1)
{ 
 data = (T*)new void*[ 1+asize*sizeof(T)/sizeof(void*) ]; 

 for (int i=0 ; i<size ; i++)
   new(data+i) T( r->data[i] );
} /* CDSListRep::CDSListRep */

template<class T>
CDSListRep<T>::~CDSListRep()
{
 assert(count==0);
 for (int i=0 ; i<size ; i++)
   data[i].~T();

 delete [] (void**)data;
} /* CDSListRep::~CDSListRep */


template<class T>
CDSList<T>::~CDSList()
{
 assert( rep->count>0 );
 if (--rep->count <= 0) delete rep;
} /* CDSListRep::~CDSListRep */

template<class T>
inline T&
CDSList<T>::operator[](const int &index)  //non-const version
{
 assert( index>=0  && index<size() );

 splitRep();
 return (rep->data[index]);
}

template<class T>
inline const T& 
CDSList<T>::operator[](const int &index) const  //const version
{
 assert( index>=0  && index<size() );

 return (rep->data[index]);
}

//
// CDSList container class:
//   contained objects must have valid copy constructor
//
// in addition
//   default contructor is required if following methods are used:
//     list constructor with arguments
//     resize
//
//   assignment operator required for
//     nonconst operator[]
//     remove
//     reverse
//     sort
//     
//   equality operator is required for
//     getIndex
//     contains
//
//


template<class T> 
CDSList<T>& 
CDSList<T>::resize(const int nsize)
  // calls default constructor
{
 if (nsize == rep->size)
   return *this;

 //case count>1: always do copy
 if (rep->count>1) {
   CDSListRep<T>* newRep = new CDSListRep<T>(nsize,rep->asize);
   for (int i=0 ; i<CDSMath::min(nsize,size()) ; i++) 
     new(newRep->data+i) T(rep->data[i]);
   //if nsize>size, call constructors
   for (l_int i=size() ; i<nsize ; i++)
     new(newRep->data+i) T();
   rep->count--;
   rep = newRep;
 } else {
   //case count==1: do realloc only if we don't have enough allocated space
   if ( nsize>rep->asize ) {
     rep->asize += nsize;
     
     T* ndata = (T*)new void*[1+rep->asize*sizeof(T)/sizeof(void*) ];

     for (int i=0 ; i<size() ; i++) {
       new(ndata+i) T( rep->data[i] );
       rep->data[i].~T();
     }
     
     delete [] (void**)rep->data;
     rep->data = ndata;
   } 
   //if nsize<size, call destructors
   for (int i=nsize ; i<rep->size ; i++)
     rep->data[i].~T();
   //if nsize>size, call constructors
   for (l_int i=rep->size ; i<nsize ; i++)
     new(rep->data+i) T();
   rep->size = nsize;
 }

 return *this;
} /*resize*/

 
template<class T> 
void
CDSList<T>::reserve(const int nasize)
  // calls copy constructor
{
 if (nasize <= rep->asize) 
   return;

 splitRep();

 int asize = rep->asize+nasize; //Consider this.

 T* ndata = (T*)new void*[1+asize*sizeof(T)/sizeof(void*) ];
 for (int i=0 ; i<size() ; i++) {
   new(ndata+i) T(rep->data[i]);
   rep->data[i].~T();
 }
 delete [] (void**)rep->data;
 rep->data = ndata;
 rep->asize = asize;
} /*reserve*/

 
template<class T> 
CDSList<T>::CDSList(const CDSList<T> &l)
  //
  // copy constructor
  //
{
 rep = l.rep;
 rep->count++;

} /* CDSList(const CDSList&) */

template<class T> 
CDSList<T>& 
CDSList<T>::operator=(const CDSList<T> &l)
  //
  // assignment of list to list
  //
{
 if (&l==this)  //self-assignment
   return *this;

 l.rep->count++;
 if (--rep->count <= 0) delete rep;
 rep = l.rep;

 return *this;
}

template<class T> 
template<class VECTOR>
CDSList<T>& 
CDSList<T>::operator=(const CDS::Vector<VECTOR>& v)
  //
  // assignment of vector to list
  //
{
 resize(v.size());

 for (int i=0 ; i<size() ; i++)
   (*this)[i] = v(i);

 return *this;
}

template<class T>
void
CDSList<T>::append(const T member) // this calls the copy constructor
  // member not passed by reference: 
  //   if member is another element of this list, a reference can be
  //   rendered invalid by a call to reserve.
{
 splitRep();

 reserve(size()+1);
 new( rep->data+size() ) T(member);
 rep->size++;
} /* append */
 

template<class T>
const T& 
CDSList<T>::prepend(const T& member) // this calls the copy constructor
{
 splitRep();

 reserve(size()+1);
 for (int i=size() ; i>0 ; i--)  {//ugly: two copies could happen
   new( rep->data+i ) T( rep->data[i-1] );
   rep->data[i-1].~T();
 }
 new( rep->data ) T(member);
 rep->size++;

 return member;
}
 

template<class T> 
int
CDSList<T>::getIndex(const T& member) const
{
 for (int i=0 ; i<size() ; i++)
   if ( rep->data[i] == member )
     return i;

 return -1;
}
 

template<class T> 
bool
CDSList<T>::contains(const T& member) const
{
 if ( getIndex(member) >= 0 )
   return 1;
 else
   return 0;
}
 

template<class T>
void 
CDSList<T>::remove(const int i)
  //
  // remove an element of the list
  //
{
 assert( i<size() && i>=0 );

 splitRep();

 for (int j=i+1 ; j<size() ; j++)
   rep->data[j-1] = rep->data[j];

 rep->data[size()-1].~T();
 rep->size--;
} /* remove */
 
template<class T>
T
CDSList<T>::pop()
  //
  // remove the last element of the list and return it
  //
{
 T ret = rep->data[size()-1];
 remove( size()-1 );
 return ret;
} /* pop */
 
template <class T> 
ostream& 
operator<<(ostream& os, const CDSList<T>& v)
{
 os << '{';
 for (int i=0 ; i<v.size() ; i++) {
   os << ' ' << v[i];
   if ( i<v.size()-1 ) os << ',';
 }
 os << " }";
 return os;
}

template<class T>
void
CDSList<T>::reverse()
{
 splitRep();
 T tmp;
 for (int i=0 ; i<size()/2 ; i++) {
   tmp = rep->data[i];
   rep->data[i] = rep->data[size()-i-1];
   rep->data[size()-i-1] = tmp;
 }
} /* reverse */

template<class T>
template<class Comparer>
void
CDSList<T>::sort(Comparer c)
{
 splitRep();

 if ( size() < 50) 
   sort_si(c);
 else if ( size() < 1000 )
   sort_shell(c);
 else //size() >=1000
   sort_heap(c);
} /* sort */

template<class T>
template<class Comparer>
void
CDSList<T>::sort_si(Comparer compare) //straight insertion
{
 for (int j=1 ; j<size() ; j++) {
   T a = rep->data[j];
   int i=j-1;
   while ( i>=0 &&
	   compare(rep->data[i],a)>0 ) {
     rep->data[i+1] = rep->data[i];
     i--;
   }
   rep->data[i+1] = a;
 }
} /* sort_si */

template<class T>
template<class Comparer>
void
CDSList<T>::sort_shell(Comparer compare) //shell's method
{
 const double ALN2I = 1.442695022;
 const double TINY  = 1e-5;

 T* data1 = rep->data-1;     //ugly: routine uses offset 1
 int lognb2 = (int)(log((double)size()) * ALN2I + TINY);
 int m=size();
 for (int nn=1 ; nn<=lognb2 ; nn++) {
   m >>= 1;
   for (int j=m+1 ; j<=size() ; j++) {
     int i=j-m;
     T t = data1[j];
     while ( i>=1 && 
	     compare(data1[i],t)>0 ) {
       data1[i+m] = data1[i];
       i -= m;
     }
     data1[i+m] = t;
   }
 }
} /* sort_shell */

template<class T>
template<class Comparer>
void
CDSList<T>::sort_heap(Comparer compare) //heap sort
{
 T* data1 = rep->data-1;     //ugly: routine uses offset 1

 int l = (size()>>1)+1;
 int ir = size();
 for (;;) {
   T rra;
   if (l>1)
     rra = data1[--l];
   else {
     rra=data1[ir];
     data1[ir]=data1[1];
     if (--ir == 1) {
       data1[1] = rra;
       break;
     }
   }
   int i=l;
   int j=l<<1;
   while (j<=ir) {
     if (j<ir && 
	 compare(data1[j],data1[j+1])<0 )
       ++j;
     if ( compare(rra,data1[j])<0 ) {
       data1[i]=data1[j];
       j += (i=j);
     } else
       j=ir+1;
   }
   data1[i] = rra;
 }
} /* sort_heap */

#ifdef DEBUG_ALLOC

template<class T>
void*
CDSList<T>::operator new(size_t s)
{
 cout << "CDSList: operator new allocating " << s << " bytes\n";
 return ::operator new(s);
} /* operator new */

template<class T>
void
CDSList<T>::operator delete(void* p)
{
 cout << "CDSList: operator delete freeing " << p << "\n";
 ::delete p;
} /* operator delete */

template<class T>
void*
CDSList<T>::operator new[](size_t s)
{
 cout << "CDSList: operator new[] allocating " << s << " bytes\n";
 return ::operator new[](s);
} /* operator new[] */

template<class T>
void
CDSList<T>::operator delete[](void* p)
{
 cout << "CDSList: operator delete[] freeing " << p << "\n";
 ::delete [] p;
} /* operator delete[] */

#endif /* DEBUG_ALLOC */

int CDSList_test();

#endif /*__cdslist_hh__*/
