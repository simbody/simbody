#ifndef SimTK_SimTKCOMMON_OLDARRAY_H_
#define SimTK_SimTKCOMMON_OLDARRAY_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-7 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

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
	explicit ArrayHelper(const TypeManipulatorT&, int n=0);
	ArrayHelper(const TypeManipulatorT&, int n, const void* ival, bool repeat);
	ArrayHelper(const ArrayHelper&);

    // Sub-array constructors, read only and writable
    ArrayHelper(const ArrayHelper&, int offset, int length);
    ArrayHelper(ArrayHelper&,       int offset, int length);    
    
	~ArrayHelper();
	ArrayHelper& operator=(const ArrayHelper&);

    void reverse(); // reverse the order of all elements
	
    // use signed integers to avoid signed/unsigned mismatch misery
	int capacity() const;
	int size() const;	
		
	const void* operator[](int i) const; 
	void* operator[](int i);
	
	void push_back(const void*);
	void pop_back();
	void clear();
	void resize(int n, const void* x=0);
    void reserve(int n);
	
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
    explicit ArrayBase(int n) : ah(type(),n) { }
    ArrayBase(int n, const T& ival)  : ah(type(),n,&ival,true) { }
    ArrayBase(int n, const T* ivals) : ah(type(),n,ivals,false) { }
        
    // Sub-array constructors, read only and writable, and their operator equivalents
    ArrayBase(const ArrayBase& a, int offset, int length) : ah(a.ah,offset,length) { }
    ArrayBase(ArrayBase&       a, int offset, int length) : ah(a.ah,offset,length) { }

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
        for (int i=0; i < size(); ++i)
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

    // If the Array is empty, begin() and end() both return null so 
    // begin()==end() will return true.
    iterator       begin()       {return size() ? &front() : 0;}
    const_iterator begin() const {return size() ? &front() : 0;}
    iterator       end()         {return size() ? &back() + 1 : 0;}
    const_iterator end()   const {return size() ? &back() + 1 : 0;}

    void reverse() {ah.reverse();}

    int capacity() const { return ah.capacity(); }
    int size()     const { return ah.size(); }

    const T& operator[](int i) const {return getAs(ah[i]);}
    T&       operator[](int i)       {return updAs(ah[i]);}
    
    // Append to end of array.
    ArrayBase& operator+=(const T& x) { push_back(x); return *this; }
        
    void push_back(const T& x)           {ah.push_back(&x);}
    void pop_back()                      {ah.pop_back();}
    void clear()                         {ah.clear();}
    void resize(int n)                   {ah.resize(n);}
    void resize(int n, const T& x)       {ah.resize(n,&x);}
    void reserve(int n)                  {ah.reserve(n);}
    
private:
    SimTKimpl::ArrayHelper  ah;
    
    static const SimTKimpl::TypeManipulatorT& type() {
        return SimTKimpl::MakeTypeManipulator<T>::getTypeManipulatorT();
    }
    static const T& getAs(const void* p) {return SimTKimpl::MakeTypeManipulator<T>::getAs(p);}
    static T&       updAs(void* p)       {return SimTKimpl::MakeTypeManipulator<T>::updAs(p);}    
};  

/**
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
    ArrayView(const ArrayView& a, int offset, int length) : Base(a,offset,length) { }
    ArrayView(ArrayView&       a, int offset, int length) : Base(a,offset,length) { }
    ArrayView operator()(int offset, int length) const { return ArrayView(*this, offset, length); }
    ArrayView operator()(int offset, int length)       { return ArrayView(*this, offset, length); } 
    
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
 * (subject to ordinary C++ packing rules).
 */
template <class T> class Array : public ArrayBase<T> {
    typedef ArrayBase<T> Base;
public:
	Array() : Base() { }
	explicit Array(int n) : Base(n) { }
      
    // Default copy, assignment, destructor made explicit here for convenient debugging.
    // These are "deep" copies.
    Array(const Array& src) : Base(src) { }
    ~Array() { }
    Array& operator=(const Array& src) {Base::operator=(src); return *this;}
    // Assign to every element.
    Array& operator=(const T& x)       {Base::operator=(x); return *this;}

    // These allocate new space and fill it with copies of the supplied object(s)
    Array(int n, const T& ival)  : Base(n,ival) { }
    Array(int n, const T* ivals) : Base(n,ivals) { }
        
    // Sub-array constructors, read only and writable, and their operator equivalents
    Array(const Array& a, int offset, int length) : Base(a,offset,length) { }
    Array(Array&       a, int offset, int length) : Base(a,offset,length) { }
    
    const ArrayView<T> operator()(int offset, int length) const 
      { return ArrayView<T>(*this, offset, length); }
    ArrayView<T>       operator()(int offset, int length)       
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
template <class T> inline int
findFirstOf(const Array<T>& a, const T& test) {
    for (int i=0; i<a.size(); ++i)
        if (a[i] == test) return i;
    return -1;
} 
template <class T> inline int
findFirstOf(const ArrayView<T>& a, const T& test) 
  { return findFirstOf((const Array<T>&)a,test); }

/// If the type T supports an "==" operator, you can instantiate
/// this method to find the indices of all elements of an Array<T> at which 
/// the elements match the supplied test element. The returned index
/// list will have length 0 if there are no matches.

template <class T> inline Array<int>
findAllOf(const Array<T>& a, const T& test) {
    Array<int> matches;
    for (int i=0; i<a.size(); ++i)
        if (a[i] == test) matches.push_back(i);
    return matches;
} 

template <class T> inline Array<int>
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
    for (int i=0; i<a.size(); ++i)
        s << a[i] << std::endl;
    return s;
} 
template <class T> inline std::ostream&
operator<<(std::ostream& s, const ArrayView<T>& a) {
    return s << (const Array<T>&)a;
}

} // namespace SimTK
  
#endif // SimTK_SimTKCOMMON_OLDARRAY_H_
