#ifndef SimTK_SimTKCOMMON_LIST_H_
#define SimTK_SimTKCOMMON_LIST_H_

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
#include "SimTKcommon/internal/Array.h"
#include "SimTKcommon/internal/Concretize.h"

#include <complex>

namespace SimTK {
     
template <class T> class List;
template <class T> class ListView; 

template <class T> class ListBase : public ArrayBase< Concretize<T> > {
    typedef Concretize<T> Element;
    typedef ArrayBase<Element> Base;
public:
    ListBase() { }
    ListBase(size_t n) : Base(n) { }

    // default copy, assignment, destructor made explicit here for convenient debugging
    ListBase(const ListBase& lb) : Base(lb) { }
    ~ListBase() { }
    ListBase& operator=(const ListBase& lb) { Base::operator=(lb); return *this; }
    
    // These make copies of the supplied objects    
    ListBase(size_t n, const T& initVal) : Base(n,Element(&initVal)) { }
    ListBase(size_t n, const T* initVal) : Base(n,Element(initVal))  { }
    ListBase(size_t n, const T* const* initVals) : Base(n) 
        { for (size_t i=0; i<n; ++i) Base::operator[](i) = *initVals[i]; }

    // This steals all the space and zeroes out the passed-in array of pointers.        
    ListBase(size_t n, T** pp) : Base(n) 
        { for (size_t i=0;i<n;++i) Base::operator[](i).replace(&pp[i]); }
    
    // Sub-ListBase constructors, read only and writable.
    // These share space with the original ListBase.
    ListBase(const ListBase& l, size_t offset, size_t length) : Base(l,offset,length) { }
    ListBase(ListBase&       l, size_t offset, size_t length) : Base(l,offset,length) { }
    
    Element& operator[](size_t i)       { return Base::operator[](i); }
    const T& operator[](size_t i) const { return Base::operator[](i).getRef(); }    
        
    void push_back(const T& x)   { Base::push_back(Element(&x)); } // copies x
    void push_back(T* p)         { Base::push_back(Element(p));  } // steals p
    void push_back(T*& pr)       { Base::push_back(Element(pr));  } // steals p and clears it
    ListBase& operator+=(const T& x) { push_back(x); }
    ListBase& operator+=(T*  p)      { push_back(p); }
    ListBase& operator+=(T*& pr)     { push_back(pr); }    
        
    bool isEmpty(size_t i) const { return Base::operator[](i).isEmpty(); }
};

// Partial specialization for pointers -- these are just Arrays of pointers
template <class T> class ListBase<T*>: public ArrayBase<T*> {
    typedef T* Element;
    typedef ArrayBase<Element> Base;
public:
    ListBase() { }
    ListBase(size_t n) : Base(n) { }
    ListBase(size_t n, Element initVal) : Base(n,initVal) { }     
    ListBase(size_t n, const Element* initVals) : Base(n,initVals) { }
    ListBase(const ListBase& l, size_t offset, size_t length) : Base(l,offset,length) { }
    ListBase(ListBase&       l, size_t offset, size_t length) : Base(l,offset,length) { }
    
    // default copy, assignment, destructor made explicit here for convenient debugging
    ListBase(const ListBase& lb) : Base(lb) { }
    ~ListBase() { }
    ListBase& operator=(const ListBase& lb) { Base::operator=(lb); return *this; }
 
    // inherit most operators and methods from Base   
    
    ListBase& operator+=(Element x) { Base::operator+=(x); return *this; }
    bool isEmpty(size_t) const { return false; }
};

// ListBase template specializations for known concrete objects
#define SimTK_LIST_SPECIALIZE(T)                                    \
template <> class ListBase< T > : public ArrayBase< T > {           \
    typedef ArrayBase< T > Base;                                    \
public:                                                             \
    ListBase() { }                                                  \
    ListBase(size_t n) : Base(n) { }                                \
    ListBase(size_t n, const T& initVal) : Base(n,initVal) { }      \
    ListBase(size_t n, const T* initVals): Base(n,initVals) { }     \
    ListBase(const ListBase& l, size_t offset, size_t length) : Base(l,offset,length) { }   \
    ListBase(ListBase&       l, size_t offset, size_t length) : Base(l,offset,length) { }   \
    /* default copy, assignment, destructor; inherit most operators & methods */            \
    ListBase& operator+=(const T& x) { Base::operator+=(x); return *this; } \
    bool isEmpty(size_t) const { return false; }                            \
};

// Built-in types
SimTK_LIST_SPECIALIZE(bool);            SimTK_LIST_SPECIALIZE(signed char); 
SimTK_LIST_SPECIALIZE(char);            SimTK_LIST_SPECIALIZE(unsigned char);
SimTK_LIST_SPECIALIZE(short);           SimTK_LIST_SPECIALIZE(int); 
SimTK_LIST_SPECIALIZE(long); 
SimTK_LIST_SPECIALIZE(unsigned short);  SimTK_LIST_SPECIALIZE(unsigned int); 
SimTK_LIST_SPECIALIZE(unsigned long);
SimTK_LIST_SPECIALIZE(float);           SimTK_LIST_SPECIALIZE(double); 
SimTK_LIST_SPECIALIZE(long double);
SimTK_LIST_SPECIALIZE(std::complex<float>);
SimTK_LIST_SPECIALIZE(std::complex<double>); 
SimTK_LIST_SPECIALIZE(std::complex<long double>);


/**
 * ListView is a 'dummy' class which is completely interchangeable with List. 
 * We use this to give tight control over deep and shallow copying, following 
 * these rules: 
 *      (1) any copy construction or initialization into a List does a deep 
 *          copy of the source, and
 *      (2) any copy construction or initialization into a ListView does a
 *          shallow copy.
 * However, List& objects are interchangeable with ListView& objects. 
 * So passing a ListView object to a List& function parameter is allowed.
 */    
template <class T> class ListView : public ListBase<T> {
    typedef ListBase<T> Base;
public:
    // no default constructor
    ListView(const ListView& v) : Base(v, 0, v.size()) { } // an identical view
    explicit ListView(const Base& l) : Base(l, 0, l.size()) { }  // a view of the whole thing
    ~ListView() { } // convenient for debugging
       
    // Sub-array constructors, read only and writable, and their operator equivalents
    ListView(const Base& a, size_t offset, size_t length) : Base(a,offset,length) { }
    ListView(Base&       a, size_t offset, size_t length) : Base(a,offset,length) { }
    ListView operator()(size_t offset, size_t length) const { return ListView(*this, offset, length); }
    ListView operator()(size_t offset, size_t length)       { return ListView(*this, offset, length); } 
        
    // Shallow assignment only (ListBase understands)
    ListView& operator=(const ListView& v) { Base::operator=(v); return *this; }
    ListView& operator=(const Base& b) { Base::operator=(b); return *this; }
    
    // Conversions
    operator const List<T>&() const { return *reinterpret_cast<const List<T>*>(this); }
    operator List<T>&()             { return *reinterpret_cast<List<T>*>(this); }
private:
    // NO DATA MEMBERS ALLOWED
    ListView() { assert(false); } // default construction not allowed
};

/**
 * @brief Container class with hidden implementation.
 * This container makes randomly accessible lists of any object, abstract or
 * concrete but as a consequence can make no guarantees about adjacent List
 * items being adjacent in memory.
 *
 * Because the implementation is opaque, List<T> is an acceptable data type
 * for passing across a binary interface, while STL's vector<T> is not.
 */
template <class T> class List : public ListBase<T> {
    typedef ListBase<T> Base;
public:
	List() : Base() { }
	explicit List(size_t n) : Base(n) { }
    ~List() { } // convenient for debugging
    
    // default copy, assignment, destructor (these will be "deep" copies)
    
    // These allocate new space and fill it with copies of the supplied objects    
	List(size_t n, const T& initVal) : Base(n,initVal) { }
    List(size_t n, const T* initVal) : Base(n,initVal)  { }
    List(size_t n, const T* const* initVals) : Base(n,initVals) { }

    // This allocates a new List but steals the initializing objects rather
    // than copying them. For safety it zeroes out the passed-in array of pointers.        
    List(size_t n, T** pp) : Base(n,pp) { }
    
    // Sub-list constructors, read only and writable, and their operator equivalents.
    // These share space with the original List.
    List(const List& l, size_t offset, size_t length) : Base(l,offset,length) { }
    List(List&       l, size_t offset, size_t length) : Base(l,offset,length) { }
    
    const ListView<T> operator()(size_t offset, size_t length) const 
        { return ListView<T>(*this, offset, length); }
    ListView<T>       operator()(size_t offset, size_t length)       
        { return ListView<T>(*this, offset, length); }    
        
    List& operator+=(const T& x) { Base::operator+=(x); return *this; }
    List& operator+=(T*  p)      { Base::operator+=(p); return *this; }
    List& operator+=(T*& x)      { Base::operator+=(x); return *this; }   
        
    // Conversions
    operator const ListView<T>&() const 
        { return *reinterpret_cast<const ListView<T>*>(this); }
    operator ListView<T>&()             
        { return *reinterpret_cast<ListView<T>*>(this); }
private:
    // NO DATA MEMBERS ALLOWED
};

/// If the type T supports an "==" operator, you can instantiate
/// this method to find the first element of an List<T> which matches the
/// supplied test element. The index is returned if found, otherwise -1.
template <class T> inline ptrdiff_t
findFirstOf(const List<T>& l, const T& test) {
    for (ptrdiff_t i=0; i<l.size(); ++i)
        if (l[i] == test) return i;
    return -1;
} 

template <class T> inline ptrdiff_t
findFirstOf(const ListView<T>& lv, const T& test) {
    return findFirstOf((const List<T>&)lv, test);
}

/// If the type T supports an "==" operator, you can instantiate
/// this method to find the indices of all elements of an List<T> at which 
/// the elements match the supplied test element. The returned index
/// list will have length 0 if there are no matches.

template <class T> inline Array<ptrdiff_t>
findAllOf(const List<T>& l, const T& test) {
    Array<ptrdiff_t> matches;
    for (ptrdiff_t i=0; i<l.size(); ++i)
        if (l[i] == test) matches.push_back(i);
    return matches;
} 

template <class T> inline Array<ptrdiff_t>
findAllOf(const ListView<T>& lv, const T& test) {
    return findAllOf((const List<T>&)lv, test);
}

/// If the type T supports an "==" operator, you can instantiate
/// this method to find out if an List<T> contains a particular
/// test element.
template <class T> inline bool
contains(const List<T>& l, const T& test) {
    return findFirstOf(l,test) != -1;
}

template <class T> inline bool
contains(const ListView<T>& lv, const T& test) {
    return contains((const List<T>&)lv, test);
}

} // namespace SimTK


#endif // SimTK_SimTKCOMMON_LIST_H_
