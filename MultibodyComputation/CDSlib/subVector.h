
#ifndef __subVector_hh__
#define __subVector_hh__

#include <assert.h>
#include "cdsGenericVector.h"

#include <iostream>
#include <iomanip>
using namespace std;

//
// subVector class (should work with CDSVector and FixedVectors)
//
// note: to use with const vectors, const-qualify template argument V
//

template<class V>
class SubVector {
  V&  v;
  int offset_;
  int size_;
  SubVector(); //default constructor inaccessible
public:
  typedef typename V::ElementType ElementType;

  SubVector(V&  v,
        int offset,
        int size) : 
    v(v), offset_(offset), size_(size) 
  {
   OSF_SUBVECTOR_HACK;
   assert(offset>=v.offset());
   assert(offset+size<=v.offset()+v.size());
  }
  //int size()   const { cout << "size: " << size_ << endl; return size_;}
  int size()   const { return size_;}
  int offset() const { return 0;}    //indexing into v is zero-offset

  // assignment from generic vector
  template<class VEC>
  void operator=(const CDS::GenericVector<VEC>& vright) {
   assert(vright.size() == size());
   for (int i=0 ; i<vright.size() ; i++)
     v(offset_+i) = vright(vright.offset()+i);
  }

  // conversion to generic vector
  CDS::GenericVector<SubVector<V> > vector()
  { return CDS::GenericVector<SubVector<V> >(*this); }

  const CDS::GenericVector<SubVector<V> > vector() const
  { return CDS::GenericVector<SubVector<V> >(*this); }

  //template<class VEC>
  //SubVector<V>&
  //operator=(const CDSVector<VEC>& vright) {
  // assert(vright.size() == size());
  // for (int i=0 ; i<vright.size() ; i++)
  //   v(offset_+i) = vright(vright.offset()+i);
  // return *this;
  //}
  SubVector<V>&
  operator=(const SubVector<V>& vright) {
   assert(vright.size() == size());
   for (int i=0 ; i<vright.size() ; i++)
     v(offset_+i) = vright(vright.offset()+i);
   return *this;
  }
  // zero-offset index operator
  const ElementType& operator()(int i) const {       
   assert(i>=0 && i<size_);
   return v(offset_+i);
  }
  template<class VEC>           //operator+=
  SubVector<V>&
  operator+=(const CDS::GenericVector<VEC>& vr) {
   assert(vr.size() == size());
   for (int i=0 ; i<vr.size() ; i++)
     v(offset_+i) += vr(vr.offset()+i);
   return *this;
  }
  template<class V2>           //operator-=
  SubVector<V>&
  operator-=(const V2& vr) {
   assert(vr.size() == size());
   for (int i=0 ; i<vr.size() ; i++)
     v(offset_+i) -= vr(vr.offset()+i);
   return *this;
  }
   
};

template<class V>
ostream& operator<<(      ostream&      s,
            const SubVector<V>& v)
{
 int width = s.width(); // apply width to numeric fields only
 s << setw(0) << "{ ";
 if ( v.size() ) {
   s << setw(width) << v(0); 
   for (int i=1 ; i<v.size() ; i++)
     s << ", " << setw(width) << v(i) ;
 }
 s << " }";
 return s;
} /* operator<< */

template<class V,class ATOMTYPE=double>
class ConstSubVector {
  const V&  v;
  int offset_;
  int size_;
public:
  ConstSubVector(const V&  v,
               int offset,
               int size) : 
    v(v), offset_(offset), size_(size) 
  {
   assert(offset>=v.offset());
   assert(offset+size<=v.offset()+v.size());
  }
  int size()   const { return size_;}
  int offset() const { return 0;}    //indexing into v is zero-offset

  // assignment from a vector
//  template<class V2>
//  void operator=(const V2& vright) {
//   assert(vright.size() == size());
//   for (int i=0 ; i<vright.size() ; i++)
//     v(offset_+i) = vright(vright.offset()+i);
//  }
  const ATOMTYPE& operator()(int i) const {       
   assert(i>=0 && i<size_);
   return v(offset_+i);
  }
};

template<class V,class ATOMTYPE>
ostream& operator<<(ostream&                          s,
                    const ConstSubVector<V,ATOMTYPE>& v)
{
 int width = s.width(); // apply width to numeric fields only
 s << setw(0) << "{ ";
 if ( v.size() ) {
     s << setw(width) << v(0); 
   for (int i=1 ; i<v.size() ; i++)
       s << ", " << setw(width) << v(i) ;
 }
 s << " }";
 return s;
}

template<class VectorType, class V>
VectorType
operator+(const SubVector<V>& v1,
      const VectorType&   v2)
{
 assert( v1.size() == v2.size() );

 VectorType ret = v2;
 for (int i=0 ; i<ret.size() ; i++)
   ret(i) += v1(i);
 return ret;
}


#endif /* __subVector_hh__ */
