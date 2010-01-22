#ifndef SimTK_SimTKCOMMON_ARRAY_H_
#define SimTK_SimTKCOMMON_ARRAY_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
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

/**
 * This file defines the SimTK::Array_<T> class and related support classes.
 */

#include "SimTKcommon/internal/common.h"

#include <algorithm>
#include <iterator>
#include <ostream>

namespace SimTK {

// We want the index_type and size_type for ordinary integral types to
// be either both signed or both unsigned so that an index value can
// be compared against a container's size() method without a warning.
// Also, there has to be a signed difference_type that can hold the
// difference between any two valid indices. This means we can't use
// the full range of unsigned types.

template <class X> struct IndexTraits {
    typedef X                           index_type;
    typedef typename X::size_type       size_type;
    typedef typename X::difference_type difference_type;
    // We require that max_index()+1 fit in size_type and that
    // -max_index() and max_index() fit in difference type.
    static const size_type max_size = X::max_size;
    static const char* index_name() {return X::index_name();}
};

// If max_size is m, then indices range from 0..m-1, so index differences
// range from 1-m to m-1. If the signed difference type has the same number 
// of bits as the unsigned index, we have to limit m so that m-1 is 
// representable in the signed difference type.
template <> struct IndexTraits<unsigned> {
    typedef unsigned    index_type;
    typedef unsigned    size_type;
    typedef int         difference_type;
    static const size_type max_size = 0x7fffffffU;
    static const char* index_name() {return "unsigned";}
};

template <> struct IndexTraits<int> {
    typedef int             index_type;
    typedef int             size_type;
    typedef int             difference_type;
    static const size_type  max_size = 0x7fffffff;
    static const char* index_name() {return "int";}
};

template <> struct IndexTraits<unsigned short> {
    typedef unsigned short  index_type;
    typedef unsigned short  size_type;
    typedef short           difference_type;
    static const size_type  max_size = 0x7fffU;
    static const char* index_name() {return "unsigned short";}
};

template <> struct IndexTraits<short> {
    typedef short     index_type;
    typedef short     size_type;
    typedef short     difference_type;
    static const size_type max_size = 0x7fff;
    static const char* index_name() {return "short";}
}; 

template <> struct IndexTraits<unsigned long long> {
    typedef unsigned long long  index_type;
    typedef unsigned long long  size_type;
    typedef long long           difference_type;
    static const size_type max_size = 0x7fffffffffffffffULL;
    static const char* index_name() {return "unsigned long long";}
};

template <> struct IndexTraits<long long> {
    typedef long long   index_type;
    typedef long long   size_type;
    typedef long long   difference_type;
    static const size_type max_size = 0x7fffffffffffffffLL;
    static const char* index_name() {return "long long";}
};

// A container using unsigned char as an index should use unsigned char
// as its size, meaning the max size is 255 and the max index must be
// 254. Then the difference type must hold -254..254 which takes a short.
template <> struct IndexTraits<unsigned char> {
    typedef unsigned char index_type;
    typedef unsigned char size_type;
    typedef short         difference_type;
    static const size_type max_size = 255;
    static const char* index_name() {return "unsigned char";}
};

// A container using signed char as an index should used signed char as
// its size also, so the max size is 127 meaning the max index is 126
// and the difference range is -126..126 which fits in a signed char.
template <> struct IndexTraits<signed char> {
    typedef signed char index_type;
    typedef signed char size_type;
    typedef signed char difference_type;
    static const size_type max_size = 127;
    static const char* index_name() {return "signed char";}
};

// We won't use the top bit of a char index so sizes are 0 to 127
// and index differences -126..126 which fits in a signed char.
template <> struct IndexTraits<char> {
    typedef char        index_type;
    typedef char        size_type;
    typedef signed char difference_type;
    static const size_type max_size = 127;
    static const char* index_name() {return "char";}
};

/**
 * Array_<T> is a plug-compatible replacement for std::vector<T> with
 * a number of advantages:
 *  - Most important, it is safe to pass through an API to a binary library 
 *    without worrying about compiler version or Release/Debug compatibility
 *    issues. It is guaranteed to have exactly the same memory layout in
 *    all releases, across compiler versions, and between Release and Debug
 *    builds. This allows us to use Array_<T> in the SimTK API where use
 *    of std::vector<T> would be desirable but problematic.
 *  - You can specify an optional index type which can be used to provide
 *    type-safe indexing (i.e. the Array can only be indexed by indices of
 *    a particular type, like MobilizedBodyIndex). 
 *  - It is extremely fast in Release builds with zero overhead, inline
 *    methods (but with few safety checks, too!). Microsoft in its wisdom
 *    decided that STL contains should still range check in Release builds
 *    for safety, but that makes them too slow for some purposes (and also
 *    corrupts the promise of generic programming but that's another story).
 *  - It has some dangerous extensions that permit the expert user
 *    to construct objects directly into the Array without having to copy them,
 *    a big win for complicated objects and even bigger for those that don't
 *    have copy constructors!
 *  - It has a constant-time "fast erase" operation you can use if you
 *    don't mind the Array being reordered after the erase.
 *  - It supports all standard types, methods, and iterators of 
 *    std::vector<T> so it works smoothly with the STL containers and
 *    algorithms.
 *  - It is convertible to and from std::vector<T>, although that requires
 *    copying of the elements.
 *
 * Basically this is just an std::vector without the nanny-state coddling.
 */
template <class T, class X=int, int MX=IndexTraits<X>::max_size> 
class Array_ {
public:
    typedef T           value_type;
    typedef T*          pointer;
    typedef const T*    const_pointer;
    typedef T&          reference;
    typedef const T&    const_reference;
    typedef T*          iterator;
    typedef const T*    const_iterator;

	typedef std::reverse_iterator<iterator>         reverse_iterator;
	typedef std::reverse_iterator<const_iterator>   const_reverse_iterator;

    typedef X                                       index_type;
    typedef IndexTraits<X>                          index_traits;
    typedef typename index_traits::size_type        size_type;
    typedef typename index_traits::difference_type  difference_type;

    /// Return the maximum allowable size for this container.
    size_type max_size() const {return MX;}
    const char* index_name() const {return index_traits::index_name();}

    //TODO
    // insert, assign
    // check standard
    // more raw methods

    /// Default constructor allocates no heap space and is very fast.
    Array_() : data(0), nUsed(0), nAllocated(0) {}

    /// Construct an Array containing n default-constructed elements.
    /// T's default constructor is called exactly n times. If n is zero
    /// no space will be allocated.
    explicit Array_(size_type n) {
        SimTK_SIZECHECK(n, max_size(), "Array_<T>::ctor(n)");
        allocateNoConstruct(n);
        defaultConstruct(data, data+n);
        nUsed = n;
    }

    /// Construct an Array containing n elements each set to a copy of the
    /// given initial value. T's copy constructor will be called exactly
    /// n times. If n is zero no space will be allocated.
    Array_(size_type n, const T& initVal) {
        SimTK_SIZECHECK(n, max_size(), "Array_<T>::ctor(n,T)");
        allocateNoConstruct(n);
        copyConstruct(data, data+n, initVal);
        nUsed = n;
    }

    /// Construct an Array from a range of values identified by a begin
    /// and an end iterator. This is templatized so can be used with any
    /// source type T2 that is assignment-compatible with this Array's element
    /// type T.
    template <class T2>
    Array_(const T2* b, const T2* e) {
        SimTK_ERRCHK(b!=0 && e!=0, "Array_::ctor(b,e)",
            "One or both iterators was null.");
        SimTK_ERRCHK(b <= e, "Array_::ctor(b,e)",
            "Iterators were out of order.");
        SimTK_ERRCHK3(isSizeOK(e-b), "Array_::ctor(b,e)",
            "Source has %llu elements but this Array is limited to %llu"
            " elements (index type %s).",
            (unsigned long long)(e-b), ullMaxSize(), index_name());
        const size_type sz(e-b);
        allocateNoConstruct(sz);
        copyConstruct(data, data+sz, b);
        nUsed = sz;
    }

    /// Copy constructor allocates exactly as much memory as is
    /// in use in the source (not its capacity) and copy constructs 
    /// the elements so that T's copy constructor will be called exactly
    /// src.size() times. If the source is empty, no heap space will be
    /// allocated.
    Array_(const Array_& src) {
        nUsed = src.nUsed;
        allocateNoConstruct(nUsed);
        copyConstruct(data, data+nUsed, src.data);
    }

    /// Construct this Array as a copy of another Array with different
    /// template parameters. This will work as long as the source is not
    /// larger than will fit here, and as long as the source element type T2
    /// is assignment compatible with this Array's element type T. One of T's
    /// constructors will be called exactly src.size() times; the particular
    /// constructor is whichever one best matches T(T2).
    template <class T2, class X2, int MX2>
    Array_(const Array_<T2,X2,MX2>& src) {
        SimTK_ERRCHK_ALWAYS(isSizeOK(src.nUsed), "Array_<T>::ctor(Array_<T2>)",
            "Source Array was too big to fit in this one.");
        const size_type sz(src.nUsed);
        allocateNoConstruct(sz);
        copyConstruct(data, data+sz, src.data);
        nUsed = sz;
    }

    /// Copy assignment is defined to destruct 
    /// the current contents and then <em>copy construct</em>
    /// the source into the destination.Note that this is not equivalent to 
    /// elementwise assignment. 
    Array_& operator=(const Array_& src) {
        if (this != &src) {
            eraseAll(); // all elements destructed; space unchanged
            if (   nAllocated < src.nUsed 
                || nAllocated/2 > std::max(minAlloc(), src.nUsed)) 
            {   reallocateNoDestructOrConstruct(src.nUsed); }
            copyConstruct(data, data+src.nUsed, src.data);
            nUsed = src.nUsed;
        }
        return *this;
    }

    template <class T2, class X2, int MX2>
    Array_& operator=(const Array_<T2,X2,MX2>& src) {
        SimTK_ERRCHK_ALWAYS(isSizeOK(src.nUsed), "Array_<T>::operator=()",
            "Source Array was too big to fit in this one.");
        const size_type len(src.nUsed);

        eraseAll(); // all elements destructed; space unchanged
        if (   nAllocated < src.nUsed 
            || nAllocated/2 > std::max(minAlloc(), src.nUsed)) 
        {   reallocateNoDestructOrConstruct(src.nUsed); }
        copyConstruct(data, data+src.nUsed, src.data);
        nUsed = src.nUsed;

        return *this;
    }

    // Array destructor calls the destructor for each element
    // before freeing everything.
    ~Array_() {
        eraseAll(); // each element is destructed; nUsed=0
        deallocateNoDestruct(); // free data; nAllocated=0 
    }

    // Constant time exchange of data with another Array_<T>.
    void swap(Array_& other) {
        std::swap(nAllocated,other.nAllocated);
        std::swap(nUsed,other.nUsed);
        std::swap(data,other.data);
    }

    size_type size() const {return nUsed;}
    size_type capacity() const {return nAllocated;}
    bool empty() const {return nUsed==0;}

    const T& front() const {assert(!empty()); return data[0];}
    const T& back() const {assert(!empty()); return data[nUsed-1];}
    T& front() {assert(!empty()); return data[0];}
    T& back() {assert(!empty()); return data[nUsed-1];}
    const_iterator begin() const {return data;}
    const_iterator end() const {return &data[nUsed];} // one past end
    iterator begin() {return data;}
    iterator end() {return &data[nUsed];} // one past end

    const_reverse_iterator rbegin() const 
    {   return const_reverse_iterator(end()); }
    const_reverse_iterator rend() const 
    {   return const_reverse_iterator(begin()); }
    reverse_iterator rbegin()
    {   return reverse_iterator(end()); }
    reverse_iterator rend()
    {   return reverse_iterator(begin()); }

    const T& operator[](index_type i) const {
        SimTK_INDEXCHECK(i,nUsed,"Array_<T>::operator[]() const");
        return data[i];
    }
    T& operator[](index_type i) {
        SimTK_INDEXCHECK(i,nUsed,"Array_<T>::operator[]()");
        return data[i];
    }

    const T& at(index_type i) const {
        SimTK_INDEXCHECK_ALWAYS(i,nUsed,"Array_<T>::at() const");
        return data[i];
    }
    T& at(index_type i) {
        SimTK_INDEXCHECK_ALWAYS(i,nUsed,"Array_<T>::at()");
        return data[i];
    }

    // Erase all the elements currently in this Array. Size is zero after
    // this call, and capacity may or may not change depending on whether
    // a lot of heap space is in use.
    void clear() {
        eraseAll();
        if (nAllocated > minAlloc())
            deallocateNoDestruct();
    }

    // Erase elements in range [b,e), packing in any later elements into
    // the newly-available space. Calls T's destructor once for
    // each erased element; calls T's copy constructor once for each
    // copied element; calls T's destructor once for each element vacated
    // by copying. Length is reduced by the number of elements erased but
    // the capacity does not change.
    void erase(T* b, T* e) {
        SimTK_ERRCHK(begin() <= b && b <= e && e <= end(),
            "Array<T>::erase(b,e)", "Iterators out of range.");
        T* const first = b;
        T* const last1 = e; // last + 1

        // Destruct the elements we're erasing.
        while (b != last1)
            (b++)->~T(); // destruct

        // If there are any stragglers, compress them into the gap.
        b = first;
        while (e != end()) {
            new(b++) T(*e); // copy construct
            (e++)->~T(); // destruct
        }
        nUsed -= (last1-first);
    }

    // Erase just one element. This is equivalent to erase(p,p+1) but
    // faster. Note that that means p cannot be end(). Length 
    // is reduced by 1.
    void erase(T* p) {
        SimTK_ERRCHK(begin() <= p && p < end(),
            "Array<T>::erase(p)", "Iterator must point to a valid element.");

        p->~T(); // destruct
        while (++p != end()) {
            new(p-1) T(*p); // copy construct
            p->~T();        // destruct source
        }
        --nUsed;
    }

    // Same as erase(begin(),end()) but faster. T's destructor is called
    // once for each element. Upon return nUsed is zero but the capacity
    // is unchanged.
    void eraseAll() {
        T* p = begin();
        while (p != end())
            (p++)->~T(); // destruct
        nUsed = 0;
    }

    // Erase one element and move the last one in its place. This avoids
    // having to compress the elements so this runs in constant time:
    // the element is destructed; then if it wasn't the last element the
    // copy constructor is used to copy the last element into the vacated
    // space, and the destructor is called to clear the last element. The
    // nUsed is reduced by 1 but the capacity does not change.
    void eraseFast(T* p) {
        assert(begin() <= p && p < end());
        p->~T(); // destruct
        if (p+1 != end()) {
            new(p) T(back()); // copy last element
            back().~T(); // destruct last
        }
        --nUsed;
    }

    /// This method increases the size of the Array by 1 element at the
    /// end and initializes it by copy constructing it from the given
    /// element, just like the std::vector::push_back() method. If 
    /// capacity() > size(), that's all that will happen. If capacity()==size(), there is no room for another element so
    /// we'll allocate more space and move all the elements there by 
    /// calling T's copy constructor for each element, and calling its 
    /// destructor once for each of the original elements. We return an
    /// iterator pointing to the new element, which is also the element
    /// whose reference is returned by back() after this call.
    /// @return An iterator pointing to the newly added element.
    T* push_back(const T& elt) {
        if (nAllocated == nUsed)
            growByAtLeast(1,"Array_<T>::push_back(elt)");
        T* const p = data + nUsed++;
        copyConstruct(p, elt);
        return p;
    }

    /// This dangerous method increases the Array's size by 1 element at
    /// the end but doesn't perform any construction so the memory is
    /// filled with garbage. You must immediately construct into this space,
    /// using code like:
    /// <pre>   new(a.raw_push_back()) MyConstructor(...args...); </pre>
    /// This is a substantial performance improvement when the element type
    /// is something complicated since the constructor is called once and
    /// not copied.
    /// @return An iterator pointing at the unconstructed element.
    T* raw_push_back() {
        if (nAllocated == nUsed)
            growByAtLeast(1,"Array_<T>::raw_push_back()");
        return data + nUsed++;
    }

    void pop_back() {
        SimTK_ERRCHK(!empty(), "Array_<T>::pop_back()", "Array was empty.");
        back().~T(); // destruct last element
        --nUsed;
    }

    void resize(size_type sz) {
        if (sz == nUsed) return;
        if (sz == 0) {clear(); return;}
        if (sz < nUsed) {
            erase(&data[sz], end());
            return;
        }
        // sz > nUsed
        reserve(sz);
        defaultConstruct(data+nUsed, data+sz);
        nUsed = sz;
    }

    void resize(size_type sz, const T& initVal) {
        if (sz == nUsed) return;
        if (sz == 0) {clear(); return;}
        if (sz < nUsed) {
            erase(&data[sz], end());
            return;
        }
        // sz > nUsed
        reserve(sz);
        copyConstruct(data+nUsed, data+sz, initVal);
        nUsed = sz;
    }

    void reserve(size_type newCapacity) {
        if (nAllocated >= newCapacity)
            return;
        assert(newCapacity > nUsed);
        T* newdata = allocN(newCapacity); // no construction yet
        copyConstructThenDestructSource(newdata, newdata+nUsed, data);
        freeN(data);
        data = newdata;
        nAllocated = newCapacity;
    }

private:
    // This method is used when we have already decided we need to make
    // room for some new elements by reallocation. We'll issue an error
    // message if this violates the max_size restriction (we can afford
    // to do that even in a Release build since we don't expect to grow
    // very often). Otherwise we'll allocate some more space and copy
    // construct the existing elements into the new space. We will generally
    // allocate more new space than requested, at least doubling the
    // capacity, unless that would exceed max_size(). There is also a
    // minimum amount of allocation we'll do if the current size is zero
    // or very small. Note that this method does not change the current size.
    void growByAtLeast(size_type n, const char* methodName) {
        SimTK_ERRCHK3_ALWAYS(isGrowthOK(n), methodName,
            "Can't grow this Array by %llu element(s) because it would"
            " then exceed its max_size of %llu (index type %s).",
            (unsigned long long)n, ullMaxSize(), index_name());

        // Be careful not to overflow size_type as you could if you 
        // double nAllocated rather than halving max_size.
        const size_type halfMax = max_size()/2;
        size_type newCapacity = 
            nAllocated <= halfMax 
                ? size_type(2*nAllocated) : max_size();
        newCapacity = std::max(newCapacity, minAlloc());
        T* newdata = allocN(newCapacity); // no construction yet
        copyConstructThenDestructSource(newdata, newdata+nUsed, data);
        freeN(data);
        data = newdata;
        nAllocated = newCapacity;
    }


    void allocateNoConstruct(size_type n) 
    {   data = allocN(n); nAllocated=n; nUsed=0; }
    void deallocateNoDestruct() 
    {   freeN(data); data=0; nAllocated=nUsed=0; }
    void reallocateNoDestructOrConstruct(size_type n)
    {   deallocateNoDestruct(); allocateNoConstruct(n); }

    // This sets the smallest allocation we'll do when growing.
    size_type minAlloc() const 
    {   return std::min(max_size(), size_type(4)); }
  
    // Allocate without construction. Returns a null pointer if asked
    // to allocate 0 elements. In Debug mode we'll fill memory with 
    // all 1's as a bug catcher.
    static T* allocN(size_type n) {
        if (n==0) return 0;
        unsigned char* newdata = new unsigned char[n * sizeof(T)];
        #ifndef NDEBUG
            unsigned char* b=newdata; 
            const unsigned char* e=newdata+(n*sizeof(T));
            while (b != e) *b++ = 0xff;
        #endif
        return reinterpret_cast<T*>(newdata);
    }

    // Free memory without calling destructors. Nothing happens if passed
    // a null pointer.
    static void freeN(T* p) {
        delete[] reinterpret_cast<char*>(p);
    }

    // default construct one element
    static void defaultConstruct(T* p) {new(p) T();}
    // default construct range [b,e)
    static void defaultConstruct(T* b, T* const e) 
    {   while (b!=e) new(b++) T(); }

    // copy construct one element from a given value
    static void copyConstruct(T* p, const T& v) {new(p) T(v);}
    // copy construct range [b,e) with repeats of a given value
    static void copyConstruct(T* b, const T* const e, const T& v)
    {   while(b!=e) new(b++) T(v); }
    // copy construct range [b,e) from sequence of source values
    static void copyConstruct(T* b, const T* const e, const T* src)
    {   while(b!=e) new(b++) T(*src++); }

    // Templatized copy construct will work if the source elements are
    // assignment compatible with the destination elements.
    template <class T2>
    static void copyConstruct(T* b, const T* const e, const T2* src)
    {   while(b!=e) new(b++) T(*src++); }

    // Copy construct range [b,e] from sequence of source values and
    // destruct the source after it is copied. It's better to alternate
    // copying and destructing than to do this in two passes since we
    // will already have touched the memory.
    static void copyConstructThenDestructSource
       (T* b, const T* const e, T* src)
    {   while(b!=e) {new(b++) T(*src); src++->~T();} }


    // destruct one element
    static void destruct(T* p) {p->~T();}
    // destruct range [b,e)
    static void destruct(T* b, const T* const e)
    {   while(b!=e) b++->~T(); }

    // Check that a source object's size will fit in the Array being
    // careful to avoid overflow and warnings in the comparison.
    template <class S> 
    bool isSizeOK(S srcSz) const
    {   return (unsigned long long)srcSz <= ullMaxSize(); }

    template <class S>
    bool isGrowthOK(S n) const
    {   return isSizeOK(ullCapacity() + (unsigned long long)n); }

    // Return the current capacity (i.e. the number of elements we 
    // can fit in the currently allocated space) in the widest possible
    // integer type so that we can check for overflow without getting
    // derailed by wraparound that can happen with small index types.
    unsigned long long ullCapacity() const
    {   return (unsigned long long)capacity(); }

    unsigned long long ullMaxSize() const
    {   return (unsigned long long)max_size(); }

private:
    T*        data;
    size_type nUsed;
    size_type nAllocated;
};

/// Output a human readable representation of an Array to an std::ostream
/// (like std::cout). The format is '{' \e elements '}' where \e elements
/// is a space-separated list of the Array's contents output by invoking
/// the "<<" operator on the elements. This function will not compile if
/// the element type does not support the "<<" operator. No newline is
/// issued before or after the output.
/// @relates Array_
template <class T, class X, int MX> inline std::ostream&
operator<<(std::ostream& o, const Array_<T,X,MX>& a) {
    o << '{';
    if (!a.empty()) {
        o << a[X(0)];
        for (X i(1); i<a.size(); ++i)
            o << ' ' << a[i];
    }
    return o << '}';
} 


} // namespace SimTK
  
#endif // SimTK_SimTKCOMMON_ARRAY_H_
