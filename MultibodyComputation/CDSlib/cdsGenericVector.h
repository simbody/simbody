
#ifndef __vector_hh__
#define __vector_hh__

#include <sthead.h>

CDS_BEGIN_NAMESPACE

//
// Generic Vector template class
//
// wraps other container classes to allow definition of generic functions.
//
//

template<class VECTOR>
struct GenericVector {
  VECTOR& v;
  explicit GenericVector(VECTOR& v) : v(v) {}
  explicit GenericVector(const VECTOR& v) : v( (VECTOR&)v ) {}

  //type of data in the vector
  typedef typename VECTOR::ElementType ElementType;

  int size() const { return v.size(); }
  int offset() const { return v.offset(); }

  //access to vector elements
  typename VECTOR::ElementType& operator()(int i) { return v(i); }
  const typename VECTOR::ElementType& operator()(int i) const { return v(i); }
  
} ;

CDS_END_NAMESPACE

#endif /* __vector_hh__ */
