#ifndef SimTK_SimTKCOMMON_ARRAY_H_
#define SimTK_SimTKCOMMON_ARRAY_H_

/* Copyright (c) 2005-6 Stanford University and Michael Sherman.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "SimTKcommon/internal/common.h"

#include <cstddef>
#include <cassert>
#include <ostream>


namespace SimTKimpl {
	
class ArrayHelperImpl;

/**
 * Non-templatized helper class for ArrayBase<T>.
 */

class SimTK_SimTKCOMMON_EXPORT ArrayHelper {
public:		
	// no default constructor
	explicit ArrayHelper(const TypeManipulatorT&, ptrdiff_t n=0);
	ArrayHelper(const TypeManipulatorT&, ptrdiff_t n, const void* ival, bool repeat);
	ArrayHelper(const ArrayHelper&);

    // Sub-array constructors, read only and writable
    ArrayHelper(const ArrayHelper&, ptrdiff_t offset, ptrdiff_t length);
    ArrayHelper(ArrayHelper&,       ptrdiff_t offset, ptrdiff_t length);    
    
	~ArrayHelper();
	ArrayHelper& operator=(const ArrayHelper&);

    void      reverse();                 // reverse the order of all elements
	
    // use signed integers to avoid signed/unsigned mismatch misery
	ptrdiff_t capacity() const;
	ptrdiff_t size() const;	
		
	const void* operator[](ptrdiff_t i) const; 
	void* operator[](ptrdiff_t i);
	
	void push_back(const void*);
	void pop_back();
	void clear();
	void resize(ptrdiff_t n, const void* x=0);
    void reserve(ptrdiff_t n);
	
private:
	ArrayHelperImpl* impl;
    
    ArrayHelperImpl&        updImpl()       {assert(impl); return *impl;}
    const ArrayHelperImpl&  getImpl() const {assert(impl); return *impl;}
};

} //namespace SimTKimpl

namespace SimTK {
template <class T> class Array;
template <class T> class ArrayView;

template <class T> class ArrayBase {
public:
    ArrayBase() : ah(type()) { }
    explicit ArrayBase(ptrdiff_t n) : ah(type(),n) { }
    ArrayBase(ptrdiff_t n, const T& ival)  : ah(type(),n,&ival,true) { }
    ArrayBase(ptrdiff_t n, const T* ivals) : ah(type(),n,ivals,false) { }
        
    // Sub-array constructors, read only and writable, and their operator equivalents
    ArrayBase(const ArrayBase& a, ptrdiff_t offset, ptrdiff_t length) : ah(a.ah,offset,length) { }
    ArrayBase(ArrayBase&       a, ptrdiff_t offset, ptrdiff_t length) : ah(a.ah,offset,length) { }

    // Default copy, assignment, destructor made explicit here for convenient debugging.
    // These are "deep" copies.
    ArrayBase(const ArrayBase& src) : ah(src.ah) { }
    ~ArrayBase() { }
    ArrayBase& operator=(const ArrayBase& src) {
        if (this != &src)
            ah = src.ah;
        return *this;
    }

    // Assign the same value to every element.
    ArrayBase& operator=(const T& x) {
        for (ptrdiff_t i=0; i < size(); ++i)
            (*this)[i] = x;
        return *this;
    }

    typedef T*       iterator;
    typedef const T* const_iterator;
    typedef T&       reference;
    typedef const T& const_reference;

    reference       front()       {assert(size()); return (*this)[0];}
    const_reference front() const {assert(size()); return (*this)[0];}
    reference       back()        {assert(size()); return (*this)[size()-1];}
    const_reference back()  const {assert(size()); return (*this)[size()-1];}

    iterator       begin()       {return &front();}
    const_iterator begin() const {return &front();}
    iterator       end()         {return &back() + 1;}
    const_iterator end()   const {return &back() + 1;}

    void reverse() {ah.reverse();}

    ptrdiff_t capacity() const { return ah.capacity(); }
    ptrdiff_t size()     const { return ah.size(); }

    const T& operator[](ptrdiff_t i) const { return getAs(ah[i]); }
    T&       operator[](ptrdiff_t i)       { return updAs(ah[i]); }
    
    // Append to end of array.
    ArrayBase& operator+=(const T& x) { push_back(x); return *this; }
        
    void push_back(const T& x)           { ah.push_back(&x); }
    void pop_back()                      { ah.pop_back(); }
    void clear()                         { ah.clear(); }
    void resize(ptrdiff_t n)             { ah.resize(n); }
    void resize(ptrdiff_t n, const T& x) { ah.resize(n,&x); }
    void reserve(ptrdiff_t n)            { ah.reserve(n); }
    
private:
    SimTKimpl::ArrayHelper  ah;
    
    static const SimTKimpl::TypeManipulatorT& type()  
        { return SimTKimpl::MakeTypeManipulator<T>::getTypeManipulatorT(); }
    static const T& getAs(const void* p) { return SimTKimpl::MakeTypeManipulator<T>::getAs(p); }
    static T&       updAs(void* p)       { return SimTKimpl::MakeTypeManipulator<T>::updAs(p); }    
};  

/*
 * This class is a duplicate of Array and can be cast to an Array with no
 * harm. However, this provides an alternate type of temporary object
 * so that we can force copy construction when appropriate. Assume a is
 * an Array. Then we want identical behavior for all of these:
 *      Array b(a); Array b=a; {Array b; b=a; } // b is a COPY of a
 *  and Array b(a(2,3)); Array b=a(2,3)); { Array b; b=a(2,3); }
 *          // b is a 3-element COPY of a[2],a[3],a[4].
 * This won't happen without special handling because the compiler is
 * permitted to skip the copy constructor if a(2,3) were to return an Array,
 * in which case Array b(a(2,3)) MIGHT end up a view rather than a copy.
 * The same applies for functions returning Array when used to construct
 * another Array -- the compiler is free to avoid the temporary by constructing
 * directly into the destination, meaning you could get an unwanted view.
 */
   
template <class T> class ArrayView : public ArrayBase<T> {
    typedef ArrayBase<T> Base;
public:   
    // no default constructor
    ArrayView(const ArrayView& v) : Base(v, 0, v.size()) { } // an identical view
    explicit ArrayView(const Base& l) : Base(l, 0, l.size()) { }  // a view of the whole thing
    ~ArrayView() { } // convenient for debugging
       
    // Sub-array constructors, read only and writable, and their operator equivalents
    ArrayView(const ArrayView& a, ptrdiff_t offset, ptrdiff_t length) : Base(a,offset,length) { }
    ArrayView(ArrayView&       a, ptrdiff_t offset, ptrdiff_t length) : Base(a,offset,length) { }
    ArrayView operator()(ptrdiff_t offset, ptrdiff_t length) const { return ArrayView(*this, offset, length); }
    ArrayView operator()(ptrdiff_t offset, ptrdiff_t length)       { return ArrayView(*this, offset, length); } 
    
    // Shallow assignment only (ArrayBase understands)
    ArrayView& operator=(const ArrayView& v) { Base::operator=(v); return *this; }
    ArrayView& operator=(const Base& b) { Base::operator=(b); return *this; }
    
    // Conversions
    operator const Array<T>&() const { return *reinterpret_cast<const Array<T>*>(this); }
    operator Array<T>&()             { return *reinterpret_cast<Array<T>*>(this); }         
private:
    // NO DATA MEMBERS ALLOWED
    ArrayView() { assert(false); } // default construction not allowed     
};

/**
 * Container class like std::vector<T> but with hidden implementation.
 * This container provides functionality like STL's vector<T> but with
 * an opaque implementation suitable for passing across a binary interface.
 *
 * Like std::vector<T>, Array<T> stores items in consecutive memory locations
 * (subject to packing rules).
 */
template <class T> class Array : public ArrayBase<T> {
    typedef ArrayBase<T> Base;
public:
	Array() : Base() { }
	explicit Array(ptrdiff_t n) : Base(n) { }
      
    // Default copy, assignment, destructor made explicit here for convenient debugging.
    // These are "deep" copies.
    Array(const Array& src) : Base(src) { }
    ~Array() { }
    Array& operator=(const Array& src) { Base::operator=(src); return *this; }
    // Assign to every element.
    Array& operator=(const T& x)       { Base::operator=(x); return *this; }

    // These allocate new space and fill it with copies of the supplied object(s)
    Array(ptrdiff_t n, const T& ival)  : Base(n,ival) { }
    Array(ptrdiff_t n, const T* ivals) : Base(n,ivals) { }
        
    // Sub-array constructors, read only and writable, and their operator equivalents
    Array(const Array& a, ptrdiff_t offset, ptrdiff_t length) : Base(a,offset,length) { }
    Array(Array&       a, ptrdiff_t offset, ptrdiff_t length) : Base(a,offset,length) { }
    
    const ArrayView<T> operator()(ptrdiff_t offset, ptrdiff_t length) const 
        { return ArrayView<T>(*this, offset, length); }
    ArrayView<T>       operator()(ptrdiff_t offset, ptrdiff_t length)       
        { return ArrayView<T>(*this, offset, length); }	
    
    // Append an element to the end of the Array. TODO: is this confusing?
    Array& operator+=(const T& x) { Base::operator+=(x); return *this; }
        
    // Conversions
    operator const ArrayView<T>&() const 
        { return *reinterpret_cast<const ArrayView<T>*>(this); }
    operator ArrayView<T>&()             
        { return *reinterpret_cast<ArrayView<T>*>(this); }
private:
    // NO DATA MEMBERS ALLOWED   
};	

/// If the type T supports an "==" operator, you can instantiate
/// this method to find the first element of an Array<T> which matches the
/// supplied test element. The index is returned if found, otherwise -1.
template <class T> inline ptrdiff_t
findFirstOf(const Array<T>& a, const T& test) {
    for (ptrdiff_t i=0; i<a.size(); ++i)
        if (a[i] == test) return i;
    return -1;
} 
template <class T> inline ptrdiff_t
findFirstOf(const ArrayView<T>& a, const T& test) 
  { return findFirstOf((const Array<T>&)a,test); }

/// If the type T supports an "==" operator, you can instantiate
/// this method to find the indices of all elements of an Array<T> at which 
/// the elements match the supplied test element. The returned index
/// list will have length 0 if there are no matches.

template <class T> inline Array<ptrdiff_t>
findAllOf(const Array<T>& a, const T& test) {
    Array<ptrdiff_t> matches;
    for (ptrdiff_t i=0; i<a.size(); ++i)
        if (a[i] == test) matches.push_back(i);
    return matches;
} 

template <class T> inline Array<ptrdiff_t>
findAllOf(const ArrayView<T>& a, const T& test) { 
    return findAllOf((const Array<T>&)a,test); 
}

/// If the type T supports an "==" operator, you can instantiate
/// this method to find out if an Array<T> contains a particular
/// test element.
template <class T> inline bool
contains(const Array<T>& a, const T& test) {
    return findFirstOf(a,test) != -1;
}
template <class T> inline bool
contains(const ArrayView<T>& a, const T& test) {
    return findFirstOf(a,test) != -1;
}

template <class T> inline std::ostream&
operator<<(std::ostream& s, const Array<T>& a) {
    for (ptrdiff_t i=0; i<a.size(); ++i)
        s << a[i] << std::endl;
    return s;
} 
template <class T> inline std::ostream&
operator<<(std::ostream& s, const ArrayView<T>& a) {
    return s << (const Array<T>&)a;
}
} // namespace SimTK
  

#endif // SimTK_SimTKCOMMON_ARRAY_H_
