#ifndef SimTK_SimTKCOMMON_ARRAY_H_
#define SimTK_SimTKCOMMON_ARRAY_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-13 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

/** @file
 * This file defines the Array_<T,X> class and related support classes
 * including base classes ArrayViewConst_<T,X> and ArrayView_<T,X>, and
 * helper class ArrayIndexTraits<X>.
 */

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/ExceptionMacros.h"
#include "SimTKcommon/internal/Serialize.h"

#include <algorithm>
#include <iterator>
#include <vector>
#include <ostream>
#include <climits>
#include <typeinfo>

namespace SimTK {

// These are the classes defined in this header.
template <class X>                   struct ArrayIndexTraits;
template <class T, class X=unsigned> class  ArrayViewConst_;
template <class T, class X=unsigned> class  ArrayView_;
template <class T, class X=unsigned> class  Array_;

// NOTE: I have attempted to force the compilers to inline certain trivial
// methods here because I observed Visual C++ 2013 fail to inline operator[]
// in a performance-critical method (getCacheEntry() to be specific). It is
// essential that there be no overhead introduced by the Array_ classes, which
// Simbody uses extensively specifically because std::vector was too slow.
// (sherm 20140404).

//==============================================================================
//                           CLASS ArrayIndexTraits
//==============================================================================

/** This templatized type is used by the Array_<T,X> classes to obtain the
information they need to use the class X as an index class for the array. 
There must be a specialization here providing ArrayIndexTraits for each of
the built-in integral types that is suitable for use as an index. Any other
type X will qualify as an index if it defines the following members:

  - typedef size_type
  - typedef difference_type
  - static size_type max_size()
  - operator size_type() const (conversion from X to size_type)

max_size() determines the largest number of elements that a container may hold
if its index type is X. 

size_type must be an integral type large enough to hold
all the values from 0 to max_size(), including all the index values (which 
range from 0 to max_size()-1). size_type may be signed or unsigned; it is the 
type returned by size(), max_size(), capacity(), etc. and may be compared 
directly against an index of type X without producing a compiler warning. 

difference_type is a signed integral type that can hold all possible 
differences between two indices, that is, values between -(max_size()-1) and
+(max_size()-1). In most cases we use an integral type with the same number of 
bits for size_type and difference_type but when the index type is very small 
(bool, unsigned char, or unsigned short) we want to allow the full range (2, 
255, or 65535 elements, resp.) in which case we need a wider type to hold the 
differences. 

The conversion operator ensures that we can write size_type(i) for an index
i of type X. An explicit conversion member does not need to be present as long
as the conversion size_type(i) already works as it does for all the integral 
type specializations. 

The SimTK type-generating macro SimTK_DEFINE_UNIQUE_INDEX_TYPE() provides the
necessary members so that these types can be used directly as index types for
Array_ objects with no further preparation. For example, you can make an
Array_<int,MobilizedBodyIndex> that stores ints that can be indexed only via
MobilizedBodyIndex indices. 

@tparam X   A type suitable for use as an Array_ index.
@see Array_
**/
template <class X> struct ArrayIndexTraits {
    /** The signed or unsigned integral type to which an object of index type
    X can be converted without producing any compiler warnings. **/
    typedef typename X::size_type       size_type;
    /** A signed integral type large enough to hold the full range of 
    possible signed differences i-j between two indices i and j of type X. **/
    typedef typename X::difference_type difference_type;
    /** The maximum allowable size for any Array_<T,X> that uses this type X
    as its index type. **/
    static size_type max_size() {return X::max_size();}
};

/** Specialization of ArrayIndexTraits for \c unsigned (that is, \c unsigned
\c int) used as an index. **/
template <> struct ArrayIndexTraits<unsigned> {
    typedef unsigned        size_type;
    typedef int             difference_type;
    static size_type        max_size() {return (unsigned)INT_MAX;}
};

/** Specialization of ArrayIndexTraits for (signed) \c int used as an index. **/
template <> struct ArrayIndexTraits<int> {
    typedef int             size_type;
    typedef int             difference_type;
    static size_type        max_size() {return INT_MAX;}
};

/** Specialization of ArrayIndexTraits for \c unsigned \c long used as an index. 
@warning
Different 64 bit platforms have different lengths for long. In particular, 
64 bit MSVC++ has sizeof(long)==sizeof(int) while 64 bit gcc has 
sizeof(long)==sizeof(long long). We recommend that you avoid using long
and unsigned long (ever, not just here) and instead use int or long long (or 
their unsigned versions) which are unambiguously 32 or 64 bits, resp. **/
template <> struct ArrayIndexTraits<unsigned long> {
    typedef unsigned long       size_type;
    typedef long                difference_type;
    static size_type            max_size() {return (unsigned long)LONG_MAX;}
};

/** Specialization of ArrayIndexTraits for (signed) \c long used as an index. 
@warning
Different 64 bit platforms have different lengths for long. In particular, 
64 bit MSVC++ has sizeof(long)==sizeof(int) while 64 bit gcc has 
sizeof(long)==sizeof(long long). We recommend that you avoid using long
and unsigned long (ever, not just here) and instead use int or long long (or 
their unsigned versions) which are unambiguously 32 or 64 bits, resp. **/
template <> struct ArrayIndexTraits<long> {
    typedef long                size_type;
    typedef long                difference_type;
    static size_type            max_size() {return LONG_MAX;}
};

/** Specialization of ArrayIndexTraits for \c unsigned \c short used as an 
index. We don't have any bits to spare here so we want to allow the full 
65535 elements for an unsigned short indexed container. That means the index
difference range is -65534..+65534 which doesn't fit in a short so we have to 
use an int for difference_type to accommodate the whole range. **/
template <> struct ArrayIndexTraits<unsigned short> {
    typedef unsigned short      size_type;
    typedef int                 difference_type;
    static size_type            max_size() {return USHRT_MAX;}
};

/** Specialization of ArrayIndexTraits for (signed) \c short used as an 
index. In contrast to unsigned short, here the max size is 32767 so the index 
difference range -32766..+32766 still fits in a short so we don't need a wider
type for difference_type. **/
template <> struct ArrayIndexTraits<short> {
    typedef short               size_type;
    typedef short               difference_type;
    static size_type            max_size() {return SHRT_MAX;}
}; 


/** Specialization of ArrayIndexTraits for \c unsigned \c char used as
an index. Here we don't have any bits to spare and we want to use the full
max size of 255. The max index must then be 254, so the difference_type must 
hold -254..254 which takes a short. **/
template <> struct ArrayIndexTraits<unsigned char> {
    typedef unsigned char       size_type;
    typedef short               difference_type;
    static size_type            max_size() {return UCHAR_MAX;} // not CHAR_MAX
};

/** Specialization of ArrayIndexTraits for \c signed \c char used as
an index. In contrast with the unsigned char case which allows 255 elements, 
the max size here is 127 meaning the max index is 126 and the difference range
is -126..126 which still fits in a signed char so we don't need a wider type 
for difference_type. **/
template <> struct ArrayIndexTraits<signed char> {
    typedef signed char         size_type;
    typedef signed char         difference_type;
    static size_type            max_size() {return SCHAR_MAX;}
};

/** Specialization of ArrayIndexTraits for \c char used as
an index. The C++ standard does not specify whether \c char is a signed or
unsigned type; here we'll limit its max size to 127 so that we don't have
to use the non-standard high bit. That means it behaves just like the signed
char case; if you want the full range to a size of 255, use unsigned char
instead. **/
template <> struct ArrayIndexTraits<char> {
    typedef char                size_type;
    typedef signed char         difference_type;
    static size_type            max_size() {return (char)SCHAR_MAX;}
};

/** Specialization of ArrayIndexTraits for \c bool used as an index. OK, this 
seems unlikely but it works fine -- you get a container that can hold only two
elements indexed by either \c false (0) or \c true (1). You'll get warnings
if you try to index with an ordinary int. If anyone finds a legitimate use for
this, please post to the forum and let us know! **/
template <> struct ArrayIndexTraits<bool> {
    typedef unsigned char       size_type;
    typedef signed char         difference_type;
    static size_type            max_size() {return 2;}
};

/** Specialization of ArrayIndexTraits for \c unsigned \c long \c long used as
an index. This only makes sense in a 64-bit compilation. **/ 
template <> struct ArrayIndexTraits<unsigned long long> {
    typedef unsigned long long  size_type;
    typedef long long           difference_type;
    static size_type            max_size() 
                                    {return (unsigned long long)LLONG_MAX;}
};

/** Specialization of ArrayIndexTraits for \c long \c long used as
an index. This only makes sense in a 64-bit compilation. **/ 
template <> struct ArrayIndexTraits<long long> {
    typedef long long           size_type;
    typedef long long           difference_type;
    static size_type            max_size() {return LLONG_MAX;}
};

// Don't show this in Doxygen.
/** @cond **/
// This helper class decides what integral type we should use to best pack
// the index type's size_type representation. The idea is to pack the whole
// Array_ structure into 8 bytes on a 32 bit machine, 16 bytes on a 64 bit
// machine, using the largest integral type that will work, giving a layout
// like this:          |       data pointer     |
//                     |   nUsed   | nAllocated |

// The default implementation just uses the integral type itself.
template <class Integral, class is64Bit> struct ArrayIndexPackTypeHelper 
{   typedef Integral packed_size_type;};

// On 32 bit machine, pack anything smaller than a short into a short.
template<> struct ArrayIndexPackTypeHelper<bool,FalseType> 
{   typedef unsigned short packed_size_type;};
template<> struct ArrayIndexPackTypeHelper<char,FalseType> 
{   typedef unsigned short packed_size_type;};
template<> struct ArrayIndexPackTypeHelper<unsigned char,FalseType> 
{   typedef unsigned short packed_size_type;};
template<> struct ArrayIndexPackTypeHelper<signed char,FalseType> 
{   typedef short packed_size_type;};

// On 64 bit machine, pack anything smaller than an int into an int.
template<> struct ArrayIndexPackTypeHelper<bool,TrueType> 
{   typedef unsigned int packed_size_type;};
template<> struct ArrayIndexPackTypeHelper<char,TrueType> 
{   typedef unsigned int packed_size_type;};
template<> struct ArrayIndexPackTypeHelper<unsigned char,TrueType> 
{   typedef unsigned int packed_size_type;};
template<> struct ArrayIndexPackTypeHelper<signed char,TrueType> 
{   typedef int packed_size_type;};
template<> struct ArrayIndexPackTypeHelper<unsigned short,TrueType> 
{   typedef unsigned int packed_size_type;};
template<> struct ArrayIndexPackTypeHelper<short,TrueType> 
{   typedef int packed_size_type;};

template <class Integral> struct ArrayIndexPackType
{   typedef typename ArrayIndexPackTypeHelper<Integral,Is64BitPlatformType>
                        ::packed_size_type  packed_size_type;};
/** @endcond **/






//==============================================================================
//                            CLASS ArrayViewConst_
//==============================================================================
/** This Array_ helper class is the base class for ArrayView_ which is the
base class for Array_; here we provide only the minimal read-only "const"
functionality required by any Array_ object, and shallow copy semantics. The 
ability to write is added by the ArrayView_ class, and the additional ability 
to reallocate, insert, erase, etc. is added by the Array_ class. 

This class is particularly useful for recasting existing const data into a
const Array_ without copying. For example a const std::vector can be passed
to a const Array& argument by an implicit, near-zero cost conversion to an
ArrayViewConst_ which can then convert to a const Array&. 

An ArrayViewConst_ is given all the data it is going to have at the time it is 
constructed (except when it is being accessed from the derived Array_ class 
that has more capability). The contents and size of a ArrayViewConst_ cannot be 
changed after construction. In particular, the default copy assignment operator
is suppressed. The destructor simply disconnects the ArrayViewConst_ handle 
from the data it was referencing; no element destruction or heap deallocation 
occurs. 

@tparam T 
    The type of object to be stored in this container. 
@tparam X 
    The type to be used for indexing this container, with default unsigned
    (not size_t). Any integral type may be used, as well as user types that 
    satisfy the requirements discussed with class ArrayIndexTraits. 
@see Array_, ArrayView_, ArrayIndexTraits **/
template <class T, class X> class ArrayViewConst_ {
public:


//------------------------------------------------------------------------------
/** @name                        Typedefs

Types required of STL containers, plus index_type which is an extension, and
packed_size_type which is an implementation detail. **/
/*@{*/
/** The type of object stored in this container. **/
typedef T           value_type;
/** The index type (an extension). **/ 
typedef X           index_type;
/** A writable pointer to a value_type. **/
typedef T*          pointer;
/** A const pointer to a value_type. **/
typedef const T*    const_pointer;
/** A writable value_type reference. **/
typedef T&          reference;
/** A const value_type reference. **/
typedef const T&    const_reference;
/** A writable iterator for this container (same as pointer here). **/
typedef T*          iterator;
/** A const iterator for this container (same as const_pointer here). **/
typedef const T*    const_iterator;
/** A writable reverse iterator for this container. **/
typedef std::reverse_iterator<iterator>                 reverse_iterator;
/** A const reverse iterator for this container. **/
typedef std::reverse_iterator<const_iterator>           const_reverse_iterator;
/** An integral type suitable for all indices and sizes for this array. **/
typedef typename ArrayIndexTraits<X>::size_type         size_type;
/** A signed integral type that can represent the difference between any two
legitimate index values for this array. **/
typedef typename ArrayIndexTraits<X>::difference_type   difference_type;
/** The integral type we actually use internally to store size_type values. **/
typedef typename ArrayIndexPackType<size_type>::packed_size_type 
                                                        packed_size_type;
/*@}    End of typedefs **/


//------------------------------------------------------------------------------
/** @name         Construction, conversion, and destruction

Constructors here are limited to those that don't allocate new data, and can
only accept const data to reference. Copy assignment is suppressed. **/
/*@{*/

/** Default constructor allocates no heap space and is very fast. **/
ArrayViewConst_() : pData(0), nUsed(0), nAllocated(0) {}

/** Copy constructor is shallow; the constructed const array object will be
referencing the original source data. However, if the source is zero length, 
this will result in a default-constructed array view handle with a null data
pointer, even if the source had some unused data allocated.

@param[in]  src 
    The object whose data will be referenced. 
@par Complexity:
    Constant time; extremely fast. **/
ArrayViewConst_(const ArrayViewConst_& src) 
:   pData(0), nUsed(src.nUsed), nAllocated(0) {
    if (nUsed) pData = const_cast<T*>(src.pData);
} 

// Copy assignment is suppressed.

/** Construct an ArrayViewConst_<T> by referencing (sharing) a given range of 
const data [first,last1), without copying that data. This will work as long as 
the size of the source data does not exceed the array's max_size. The resulting
object is not resizeable but can be used to read elements of the original data.
This will becomes invalid if the original data is destructed or resized, but 
there is no way for the ArrayViewConst_ class to detect that.

@param[in]  first   
    A pointer to the first data element to be referenced.
@param[in]  last1   
    A pointer to the position one element past the last one in the range to be
    referenced.
@remarks
  - If the source data is empty, the resulting ArrayViewConst_ will also 
    be empty and will look as though it had been default-constructed. 
  - You can break the connection between the array handle and the data it
    was constructed from by calling disconnect().
@pre first <= last1, last1-first <= max_size()
@par Complexity:
    Dirt cheap. There will be no construction, destruction, or heap allocation
    performed.
@see disconnect() **/
ArrayViewConst_(const T* first, const T* last1) 
:   pData(0),nUsed(0),nAllocated(0) { 
    if (last1==first) return; // empty

    SimTK_ERRCHK((first&&last1)||(first==last1), 
        "ArrayViewConst_<T>(first,last1)", 
        "One of the source pointers was null (0); either both must be"
        " non-null or both must be null.");

    SimTK_ERRCHK3(this->isSizeOK(last1-first), 
        "ArrayViewConst_<T>(first,last1)",
        "The source data's size %llu is too big for this array which"
        " is limited to %llu elements by its index type %s.",
        this->ull(last1-first), ullMaxSize(), indexName());

    pData = const_cast<T*>(first); 
    nUsed = packed_size_type(last1-first); 
    // nAllocated is already zero
}

/** Construct a ArrayViewConst_<T> by referencing (sharing) the data in a const 
std::vector<T>, without copying the data; this is also an implicit conversion. 
This will work as long as the size of the vector does not exceed the array's 
max_size. The resulting array object is not resizeable but can be used to read
elements of the original std::vector. The array becomes invalid if the original
std::vector is destructed or resized, but there is no way for the array class
to detect that.

@param[in]  src
    The std::vector<T> whose data will be referenced by the constructed 
    ArrayViewConst_ handle.
@remarks
  - If the source std::vector is empty, the resulting array will also be empty 
    and will look as though it had been default-constructed. It will therefore
    not have any connection to the source vector.
  - This is quite dangerous to use since the connection between the array
    and the vector is tenuous and subject to the vector remaining untouched 
    during the lifetime of the array handle. There is no reference 
    counting; destructing the vector leaves the array referring to garbage. Be
    careful!
  - You can break the connection between the array view and the vector it was 
    constructed from by calling disconnect().
@pre src.size() <= max_size()
@par Complexity:
    Dirt cheap. There will be no construction, destruction, or heap allocation
    performed.
@see disconnect() **/
template <class A>
ArrayViewConst_(const std::vector<T,A>& src) 
:   pData(0),nUsed(0),nAllocated(0) { 
    if (src.empty()) return;

    SimTK_ERRCHK3(this->isSizeOK(src.size()),
        "ArrayViewConst_<T>::ctor(std::vector<T>)",
        "The source std::vector's size %llu is too big for this array which"
        " is limited to %llu elements by its index type %s.",
        this->ull(src.size()), ullMaxSize(), indexName());

    pData = const_cast<T*>(&src.front()); 
    nUsed = packed_size_type(src.size()); 
    // nAllocated is already zero
}
/** This is an implicit conversion to const ArrayView_<T,X>&, which is 
harmless since the const result won't permit writing on the elements. **/
operator const ArrayView_<T,X>&() const
{   return *reinterpret_cast<const ArrayView_<T,X>*>(this); }
/** This is an implicit conversion to const Array_<T,X>&, which is harmless
since the const result can't be used to write on or resize the data. **/
operator const Array_<T,X>&() const
{   return *reinterpret_cast<const Array_<T,X>*>(this); }

/** Disconnect this array handle from any data to which it refers, 
restoring it to the condition it would be in if it had just been 
default-constructed. The data pointer will simply be set to null; we'll assume
the owner will clean things up later. In either case the size() and capacity() 
will be zero after this call and data() will return null (0). **/
void disconnect() {
    SimTK_ASSERT(nAllocated==0,
        "ArrayViewConst_::deallocate(): called on an owner Array_");
    nUsed = 0;
    pData = 0;
}

/** The destructor just disconnects the array view handle from its data; see
disconnect() for more information. @see disconnect() **/
~ArrayViewConst_() {
    disconnect();
}
/*@}    End of constructors, etc. **/


//------------------------------------------------------------------------------
/** @name                       Size and capacity 

These methods examine the number of elements (size) or the amount of allocated
heap space (capacity). See the derived Array_<T,X> class for methods that can 
\e change the size or capacity. **/
/*@{*/

/** Return the current number of elements stored in this array. **/
size_type size() const {return size_type(nUsed);}
/** Return the maximum allowable size for this array. **/
size_type max_size() const 
{   return ArrayIndexTraits<X>::max_size(); }
/** Return true if there are no elements currently stored in this array. This
is equivalent to the tests begin()==end() or size()==0. **/
bool empty() const {return nUsed==0;}
/** Return the number of elements this array can currently hold without
requiring reallocation. The value returned by capacity() is always greater 
than or equal to size(), even if the data is not owned by this array in
which case we have capacity()==size() and the array is not reallocatable. **/
size_type capacity() const 
{   return size_type(nAllocated?nAllocated:nUsed); }
/** Return the amount of heap space owned by this array; this is the same
as capacity() for owner arrays but is zero for non-owners. 
@note There is no equivalent of this method for std::vector. **/
size_type allocated() const {return size_type(nAllocated);}
/** Does this array own the data to which it refers? If not, it can't be
resized, and the destructor will not free any heap space nor call any element
destructors. If the array does not refer to any data it is considered to be
an owner since it is resizeable. 
@note There is no equivalent of this method for std::vector. **/
bool isOwner() const {return nAllocated || pData==0;}
/*}*/


//------------------------------------------------------------------------------
/** @name                  Read-only element access

These methods provide read-only (const) access to individual elements that are
currently present in the array. The derived ArrayView_<T,X> class adds the
non-const versions of these methods. **/
/*@{*/

/** Select an element by its index, returning a const reference. Note that only 
a value of the array's templatized index type is allowed (default is unsigned).
This will be range-checked in a Debug build but not in Release.
@pre 0 <= \a i < size()
@par Complexity:
    Constant time. **/

SimTK_FORCE_INLINE const T& operator[](index_type i) const {
    SimTK_INDEXCHECK(size_type(i),size(),"ArrayViewConst_<T>::operator[]()");
    return pData[i];
}
/** Same as operator[] but always range-checked, even in a Release build.  
@pre 0 <= \a i < size()
@par Complexity:
    Constant time. **/
const T& at(index_type i) const {
    SimTK_INDEXCHECK_ALWAYS(size_type(i),size(),"ArrayViewConst_<T>::at()");
    return pData[i];
}
/** Same as the const form of operator[]; exists to provide a non-operator
method for element access in case that's needed. **/
SimTK_FORCE_INLINE const T& getElt(index_type i) const {
    SimTK_INDEXCHECK(size_type(i),size(),"ArrayViewConst_<T>::getElt()");
    return pData[i];
}
/** Return a const reference to the first element in this array, which must
not be empty (we'll check in a Debug build but not Release).
@pre The array is not empty.
@par Complexity:
    Constant time. **/
SimTK_FORCE_INLINE const T& front() const 
{   SimTK_ERRCHK(!empty(), "ArrayViewConst_<T>::front()", "Array was empty.");
    return pData[0]; }
/** Return a const reference to the last element in this array, which must
not be empty (we'll check in a Debug build but not Release).
@pre The array is not empty.
@par Complexity:
    Constant time. **/
SimTK_FORCE_INLINE const T& back() const 
{   SimTK_ERRCHK(!empty(), "ArrayViewConst_<T>::back()", "Array was empty.");
    return pData[nUsed-1]; }

/** Select a contiguous subarray of the elements of this array and create 
another ArrayViewConst_ that refers only to those element (without copying). 
@param[in]      index
    The index of the first element to be included in the subarray; this can
    be one past the end of the array if \a length is zero. 
@param[in]      length
    The length of the subarray to be produced.
@return
    A new ArrayViewConst_<T,X> object referencing the original data.
@note 
    If \a length==0 the returned array will be in a default-constructed,
    all-zero and null state with no connection to the original data.
@pre \a index >= 0, \a length >= 0
@pre \a index + \a length <= size()
@pre We'll validate preconditions in Debug builds but not Release.
@par Complexity:
    Dirt cheap; no element construction or destruction or heap allocation
    is required. **/ 
ArrayViewConst_ operator()(index_type index, size_type length) const {
    const size_type ix(index);
    SimTK_ERRCHK2(isSizeInRange(ix, size()), "ArrayViewConst_<T>(index,length)",
        "For this operator, we must have 0 <= index <= size(), but"
        " index==%llu and size==%llu.", this->ull(ix), ullSize());
    SimTK_ERRCHK2(isSizeInRange(length, size_type(size()-ix)), 
        "ArrayViewConst_<T>(index,length)", 
        "This operator requires 0 <= length <= size()-index, but"
        " length==%llu and size()-index==%llu.",this->ull(length),this->ull(size()-ix));

    return ArrayViewConst_(pData+ix, pData+ix+length);
}
/** Same as const form of operator()(index,length); exists to provide 
non-operator access to that functionality in case it is needed. **/
ArrayViewConst_ getSubArray(index_type index, size_type length) const
{   return (*this)(index,length); }

/*@}    End of element access. **/


//------------------------------------------------------------------------------
/** @name                   Iterators (const only)

These methods deal in iterators, which are STL generalized pointers. For this
class, iterators are just ordinary const pointers to T, and you may depend on 
that. By necessity, reverse iterators can't be just pointers; however, they 
contain an ordinary iterator (i.e. a pointer) that can be obtained by calling 
the reverse iterator's base() method. **/
/*@{*/

/** Return a const pointer to the first element of this array if any, otherwise
cend(), which may be null (0) in that case but does not have to be. This method
is from the proposed C++0x standard; there is also an overloaded begin() from
the original standard that returns a const pointer. **/
const T* cbegin() const {return pData;}
/** Return a const pointer to what would be the element just after the last one
in the array; this may be null (0) if there are no elements but doesn't have to
be. This method is from the proposed C++0x standard; there is also an 
overloaded end() from the original standard that returns a const pointer. **/
const T* cend() const {return pData + nUsed;}
/** The const version of begin() is the same as cbegin(). **/
const T* begin() const {return pData;}
/** The const version of end() is the same as cend(). **/
const T* end() const {return pData + nUsed;}

/** Return a const reverse iterator pointing to the last element in the array 
or crend() if the array is empty. **/
const_reverse_iterator crbegin() const 
{   return const_reverse_iterator(cend()); }
/** Return the past-the-end reverse iterator that tests equal to a reverse
iterator that has been incremented past the front of the array. You cannot 
dereference this iterator. **/
const_reverse_iterator crend() const 
{   return const_reverse_iterator(cbegin()); }
/** The const version of rbegin() is the same as crbegin(). **/
const_reverse_iterator rbegin() const {return crbegin();} 
/** The const version of rend() is the same as crend(). **/
const_reverse_iterator rend() const {return crend();}

/** Return a const pointer to the first element of the array, or possibly
(but not necessarily) null (0) if the array is empty.
@note
    cdata() does not appear to be in the C++0x standard although it would seem
    obvious in view of the cbegin() and cend() methods that had to be added. 
    The C++0x overloaded const data() method is also available. **/
const T* cdata() const {return pData;}
/** The const version of the data() method is identical to cdata(). **/
const T* data() const {return pData;}
/*@}    End of iterators. **/


//------------------------------------------------------------------------------
                                 protected:
//------------------------------------------------------------------------------
// The remainder of this class is for the use of the ArrayView_<T,X> and
// Array_<T,X> derived classes only and should not be documented for users to 
// see. 
                                     
// Don't let doxygen see any of this.
/** @cond **/
packed_size_type psize() const {return nUsed;}
packed_size_type pallocated() const {return nAllocated;}

// These provide direct access to the data members for our trusted friends.
void setData(const T* p)        {pData = const_cast<T*>(p);}
void setSize(size_type n)       {nUsed = packed_size_type(n);}
void incrSize()                 {++nUsed;}
void decrSize()                 {--nUsed;}
void setAllocated(size_type n)  {nAllocated = packed_size_type(n);}

// Check whether a given size is the same as the current size of this array,
// avoiding any compiler warnings due to mismatched integral types.
template <class S> 
bool isSameSize(S sz) const
{   return ull(sz) == ullSize(); }

// Check that a source object's size will fit in the array being careful to
// avoid overflow and warnings in the comparison.
template <class S> 
bool isSizeOK(S srcSz) const
{   return ull(srcSz) <= ullMaxSize(); }

// This is identical in function to std::distance() (reports how many 
// elements lie between two iterators) but avoids any slow 
// Release-build bugcatchers that Microsoft may have felt compelled to add.
// The implementation is specialized for random access iterators because
// they can measure distance very fast.
template<class Iter> static
typename std::iterator_traits<Iter>::difference_type
iterDistance(const Iter& first, const Iter& last1) {
    return iterDistanceImpl(first,last1,
                typename std::iterator_traits<Iter>::iterator_category());
}

// Generic slow implementation for non-random access iterators. This is fine
// for forward and bidirectional iterators, but it will *consume* input
// iterators so is useless for them.
template<class Iter> static
typename std::iterator_traits<Iter>::difference_type
iterDistanceImpl(const Iter& first, const Iter& last1, std::input_iterator_tag) {
    typename std::iterator_traits<Iter>::difference_type d = 0;
    for (Iter src=first; src != last1; ++src, ++d)
        ;
    return d;
}

// Fast specialization for random access iterators (including ordinary
// pointers) -- just subtract.
template<class Iter> static
typename std::iterator_traits<Iter>::difference_type
iterDistanceImpl(const Iter& first, const Iter& last1, 
                 std::random_access_iterator_tag) {
    return last1 - first;
}

// This method attempts to determine whether any elements in the iterator range
// [first,last1) overlap with the elements stored in this array. This is used 
// for error checks for operations where source is not permitted to overlap the
// destination. For random access iterators (including ordinary pointers), we 
// can answer this question definitively because we expect the data to be 
// consecutive in memory. For other kinds of iterators, we will just assume
// there is no overlap. Note that null ranges do not overlap even if the
// pair of equal iterators points within the other range -- what matters is
// the number of overlapping elements.
template<class Iter> bool
overlapsWithData(const Iter& first, const Iter& last1) {
    return overlapsWithDataImpl(first,last1,
                typename std::iterator_traits<Iter>::iterator_category());
}

// This is a partial specialization of the above where the data is given
// with ordinary pointers.
template <class T2> bool
overlapsWithData(const T2* first, const T2* last1) {
    // Find the start and end+1 of the alleged overlap region. There is
    // overlap iff end+1 > start. Note that this works if either range 
    // is [0,0) or [p,p), or if last1 is illegally less than first (we just
    // want to report no overlap in that case -- it is someone else's business
    // to complain).
    const T* obegin = std::max(cbegin(), (const T*)first);
    const T* oend1  = std::min(cend(),   (const T*)last1);

    return obegin < oend1;
}

// This is the generic implementation for any type of input iterator other than
// random access (i.e., bidirectional, forward, or input) -- assume no overlap.
template<class Iter> bool
overlapsWithDataImpl(const Iter&, const Iter&, std::input_iterator_tag) 
{   return false; }

// Here we can actually test for overlap since we have random access iterators.
// We convert them to pointers and then look for memory overlap.
template<class Iter> bool
overlapsWithDataImpl(const Iter& first, const Iter& last1, 
                     std::random_access_iterator_tag) {
    // We must check that the input iterators span a non-zero range before
    // assuming we can dereference them.
    if (last1 <= first)
        return false; // zero or malformed source range: no overlap

    // We now know we can dereference first and last1-1 (can't safely 
    // dereference last1 but we can use pointer arithmetic to point past
    // the (last-1)th element in memory). We then take the dereferenced
    // object's address to get ordinary pointers that we can use to 
    // watch for illegal overlap.
    return overlapsWithData(&*first, &*(last1-1)); // use pointer overload
}

// Cast an integral type to maximal-width unsigned long long to avoid accidental
// overflows that might otherwise occur due to wraparound that can happen 
// with small index types.
template <class S>
static unsigned long long ull(S sz)
{   return (unsigned long long)sz; }

// Return size(), capacity(), and max_size() cast to unsigned long long.
unsigned long long ullSize()     const {return ull(size());}
unsigned long long ullCapacity() const {return ull(capacity());}
unsigned long long ullMaxSize()  const {return ull(max_size());}

// Useful in error messages for explaining why something was too big.
const char* indexName() const {return NiceTypeName<X>::name();}

/** @endcond **/

private:
//------------------------------------------------------------------------------
//                               DATA MEMBERS
// These are the only data members and this layout is guaranteed not to change
// from release to release. If data is null, then nUsed==nAllocated==0.

T*                  pData;      // ptr to data referenced here, or 0 if none
packed_size_type    nUsed;      // number of elements currently present (size)
packed_size_type    nAllocated; // heap allocation; 0 if pData is not owned

ArrayViewConst_& operator=(const ArrayViewConst_& src); // suppressed
};






//==============================================================================
//                            CLASS ArrayView_
//==============================================================================
/** This Array_ helper class is the base class for Array_, extending 
ArrayViewConst_ to add the ability to modify elements, but not the ability to 
change size or reallocate. 

@tparam T 
    The type of object to be stored in this container. 
@tparam X 
    The type to be used for indexing this container, with default unsigned
    (not size_t). Any integral type may be used, as well as user types that 
    satisfy the requirements discussed with class ArrayIndexTraits. 
@see Array_, ArrayViewConst_, ArrayIndexTraits **/
template <class T, class X> class ArrayView_ : public ArrayViewConst_<T,X> {
typedef ArrayViewConst_<T,X> CBase;
public:
//------------------------------------------------------------------------------
/** @name                        Typedefs

Types required of STL containers, plus index_type which is an extension, and
packed_size_type which is an implementation detail. **/
/*{*/
typedef T                                               value_type;
typedef X                                               index_type;
typedef T*                                              pointer;
typedef const T*                                        const_pointer;
typedef T&                                              reference;
typedef const T&                                        const_reference;
typedef T*                                              iterator;
typedef const T*                                        const_iterator;
typedef std::reverse_iterator<iterator>                 reverse_iterator;
typedef std::reverse_iterator<const_iterator>           const_reverse_iterator;
typedef typename ArrayIndexTraits<X>::size_type         size_type;
typedef typename ArrayIndexTraits<X>::difference_type   difference_type;
typedef typename ArrayIndexPackType<size_type>::packed_size_type 
                                                        packed_size_type;
/*}*/


//------------------------------------------------------------------------------
/** @name       Construction, conversion, and destruction

Constructors here are limited to those that don't allocate new data, however
they can reference writable data. **/
/*@{*/

/** Default constructor allocates no heap space and is very fast. **/
ArrayView_() : CBase() {}

/** Copy constructor is shallow. **/
ArrayView_(const ArrayView_& src) : CBase(src) {}

/** Construct from a range of writable memory. **/
ArrayView_(T* first, const T* last1) : CBase(first,last1) {} 

/** Construct to reference memory owned by a writable std::vector. **/
template <class A>
ArrayView_(std::vector<T,A>& v) : CBase(v) {}

/** Implicit conversion of const ArrayView_ to const Array_& (zero cost). **/
operator const Array_<T,X>&() const 
{   return *reinterpret_cast<const Array_<T,X>*>(this); }  

/** Implicit conversion of non-const ArrayView_ to Array_& (zero cost). **/
operator Array_<T,X>&() 
{   return *reinterpret_cast<Array_<T,X>*>(this); } 

/** Forward to base class disconnect() method -- clears the handle without
doing anything to the data. */
void disconnect() {this->CBase::disconnect();}

/** The destructor just disconnects the array view handle from its data; see
ArrayViewConst_<T,X>::disconnect() for more information. **/
~ArrayView_() {this->CBase::disconnect();}
/*@}    End of construction, etc. **/


//------------------------------------------------------------------------------
/** @name                     Assignment

Assignment is permitted only if the source and destination are the same 
size. The semantics here are different than for a resizeable Array_
object: here the meaning is elementwise assignment rather than destruction
followed by copy construction. That is, if our elements are of type T, and
the source elements are of type T2, we will use the operator of T that
best matches the signature T::operator=(const T2&) to perform the assignments.
When the source also has type T, this is just T's copy assignment operator. 
We never perform any element destruction or construction here. **/
/*@{*/

/** Copy assignment; source must be the same size as this array. **/
ArrayView_& operator=(const ArrayView_& src) {
    if (&src != this)
        avAssignIteratorDispatch(src.cbegin(), src.cend(),
                                 std::random_access_iterator_tag(),
                                 "ArrayView_<T>::operator=(ArrayView_<T>)");
    return *this;
}


/** Assignment from any other array object is allowed as long as the number
of elements matches and the types are assignment compatible. **/
template <class T2, class X2>
ArrayView_& operator=(const ArrayViewConst_<T2,X2>& src) {
    if ((const void*)&src != (void*)this)
        avAssignIteratorDispatch(src.cbegin(), src.cend(),
                                 std::random_access_iterator_tag(),
                                 "ArrayView_<T>::operator=(Array_<T2>)");
    return *this;
}

// Help out dumb compilers struggling to match the template arguments and
// promote the Array_ or ArrayView_ to ArrayConstView_ at the same time.

/** Assignment from any other array object is allowed as long as the number
of elements matches and the types are assignment compatible. **/
template <class T2, class X2>
ArrayView_& operator=(const ArrayView_<T2,X2>& src)
{   return *this = static_cast<const ArrayViewConst_<T2,X2>&>(src); }
/** Assignment from any other array object is allowed as long as the number
of elements matches and the types are assignment compatible. **/
template <class T2, class X2>
ArrayView_& operator=(const Array_<T2,X2>& src)
{   return *this = static_cast<const ArrayViewConst_<T2,X2>&>(src); }

/** Assignment from any std::vector object is allowed as long as the number
of elements matches and the types are assignment compatible. **/
template <class T2, class A2>
ArrayView_& operator=(const std::vector<T2,A2>& src) {
    avAssignIteratorDispatch(src.begin(), src.end(),
                             std::random_access_iterator_tag(),
                             "ArrayView_<T>::operator=(std::vector<T2>)");
    return *this;
}

/** Fill assignment -- all elements are set to fillValue. @see fill() **/
ArrayView_& operator=(const T& fillValue) 
{   fill(fillValue); return *this; }

/** Assign the supplied fill value to each element of this array, using T's
copy assignment operator for each element. Note that this also serves to allow
fill from an object whose type T2 is different from T, as long as there is a 
constructor T(T2) that works since that can be invoked (implicitly or 
explicitly) to convert the T2 object to type T prior to the call. **/ 
ArrayView_& fill(const T& fillValue) {
    for (T* d = begin(); d != end(); ++d)
        *d = fillValue; // using T::operator=(T)
    return *this;
}

/** This is the same as fill() but has the usual std::vector signature for
compatibility; it will only work if the given number of elements is the same
as this array's (fixed) size. **/
void assign(size_type n, const T& fillValue) {
    SimTK_ERRCHK2(n == size(), "ArrayView_<T>::assign(n,value)",
        "Assignment to an ArrayView is permitted only if the source"
        " is the same size. Here n==%llu element(s) but the"
        " ArrayView has a fixed size of %llu.", 
        this->ull(n), this->ull(size()));

    fill(fillValue);
}

/** Assign to this array to make it a copy of the elements in range 
[first,last1) given by ordinary pointers, provided that the range is the same
size as the array. It is not allowed for the source range to include any of the 
elements currently in the array. The source elements can be 
of a type T2 that may be the same or different than this array's element type 
T as long as there is a T=T2 assignment operator that works. Note that although
the source arguments are pointers, those may be iterators for some container 
depending on implementation details of the container. Specifically, any 
ArrayViewConst_, ArrayView_, or Array_ iterator is an ordinary pointer.

@param[in]      first 
    A pointer to the first element to be copied.
@param[in]      last1 
    A pointer to the element one past the last element to be copied.
@pre last1-first == size()
@par Complexity:
    The T=T2 assignment operator will be called exactly size() times. **/
template <class T2>
void assign(const T2* first, const T2* last1) {
    const char* methodName = "ArrayView_<T>::assign(T2* first, T2* last1)";
    SimTK_ERRCHK((first&&last1)||(first==last1), methodName, 
        "One of the source pointers was null (0); either both must be"
        " non-null or both must be null.");
    // Valid pointers are random access iterators.
    avAssignIteratorDispatch(first, last1, std::random_access_iterator_tag(),
                             methodName);
}

/** Assign to this array to make it a copy of the elements in range 
[first,last1) given by non-pointer iterators (the pointer case is handled 
with a specialized assign() variant). It is not allowed for this range to 
include any of the elements currently in the array. The source elements can be 
of a type T2 that may be the same or different than this array's element type 
T as long as there is a T=T2 operator that works.

The source must have the same number of elements as the current (fixed) size
of this ArrayView. For input_iterators we'll be happy if we get enough elements
and won't insist that the input stream is empty after that. For forward_ and
bidirectional_iterators we'll copy the elements and complain at the end if
there are too few or too many. For random_access_iterators we'll check in
advance since we can do that fast.

@param[in]      first 
    An iterator pointing to the first element to be copied.
@param[in]      last1 
    An iterator pointing to the element one past the last element to be copied.

@remarks
This variant of assign() will not be called when the iterators are forward 
iterators from ArrayViewConst_, ArrayView_, or Array_ objects since those are 
ordinary pointers. 

@pre last1 is reachable from first
@pre distance(first,last1)==size()
@par Complexity:
    The T=T2 assignment operator will be called exactly size() times. **/

// Watch out for integral types matching this signature -- they must be
// forwarded to the assign(n, fillValue) signature instead.
template <class Iter>
void assign(const Iter& first, const Iter& last1)
{   avAssignDispatch(first,last1,typename IsIntegralType<Iter>::Result(),
                     "ArrayView_<T>::assign(Iter first, Iter last1)"); }
/*@}    End of assignment. */


//------------------------------------------------------------------------------
/** @name                     Element access

These methods provide read and write access to individual elements that are 
currently present in the array; the ArrayViewConst_<T,X> base class provides the
read-only (const) methods. **/
/*@{*/

/** Select an element by its index, returning a const reference. Note that only 
a value of the array's templatized index type is allowed (default is unsigned).
This will be range-checked in a Debug build but not in Release.
@pre 0 <= \a i < size()
@par Complexity:
    Constant time. **/
SimTK_FORCE_INLINE const T& operator[](index_type i) const 
{   return this->CBase::operator[](i); }

/** Select an element by its index, returning a writable (lvalue) reference. 
Note that only a value of the Array's templatized index type is allowed 
(default is unsigned). This will be range-checked in a Debug build but not 
in Release. 
@pre 0 <= \a i < size()
@par Complexity:
    Constant time. **/
SimTK_FORCE_INLINE T& operator[](index_type i) 
{   return const_cast<T&>(this->CBase::operator[](i)); }

/** Same as operator[] but always range-checked, even in a Release build.  
@pre 0 <= \a i < size()
@par Complexity:
    Constant time. **/
const T& at(index_type i) const {return this->CBase::at(i);}

/** Same as operator[] but always range-checked, even in a Release build.  
@pre 0 <= \a i < size()
@par Complexity:
    Constant time. **/
T& at(index_type i) {return const_cast<T&>(this->CBase::at(i));}

/** Same as the const form of operator[]; exists to provide a non-operator
method for element access in case that's needed. **/
SimTK_FORCE_INLINE const T& getElt(index_type i) const 
{   return this->CBase::getElt(i); }
/** Same as the non-const form of operator[]; exists to provide a non-operator
method for element access in case that's needed. **/
SimTK_FORCE_INLINE T& updElt(index_type i) 
{   return const_cast<T&>(this->CBase::getElt(i)); }

/** Return a const reference to the first element in this array, which must
not be empty.
@pre The array is not empty.
@par Complexity:
    Constant time. **/
SimTK_FORCE_INLINE const T& front() const {return this->CBase::front();} 

/** Return a writable reference to the first element in this array, which must
not be empty.
@pre The array is not empty.
@par Complexity:
    Constant time. **/
SimTK_FORCE_INLINE T& front() {return const_cast<T&>(this->CBase::front());}

/** Return a const reference to the last element in this array, which must
not be empty.
@pre The array is not empty.
@par Complexity:
    Constant time. **/
SimTK_FORCE_INLINE const T& back() const {return this->CBase::back();}

/** Return a writable reference to the last element in this array, which must
not be empty.
@pre The array is not empty.
@par Complexity:
    Constant time. **/
SimTK_FORCE_INLINE T& back() {return const_cast<T&>(this->CBase::back());}

/** Select a contiguous subarray of the elements of this array and create 
another ArrayView_ that refers only to those element (without copying). 
@param[in]      index
    The index of the first element to be included in the subarray; this can
    be one past the end of the array if \a length is zero. 
@param[in]      length
    The length of the subarray to be produced.
@return
    A new ArrayView_<T,X> object referencing the original data.
@note 
    If \a length==0 the returned array will be in a default-constructed,
    all-zero and null state with no connection to the original data.
@pre \a index >= 0, \a length >= 0
@pre \a index + \a length <= size()
@pre We'll validate preconditions in Debug builds but not Release.
@par Complexity:
    Dirt cheap; no element construction or destruction or heap allocation
    is required. **/ 
ArrayView_ operator()(index_type index, size_type length) {
    const size_type ix(index);
    SimTK_ERRCHK2(isSizeInRange(ix, size()), "ArrayView_<T>(index,length)",
        "For this operator, we must have 0 <= index <= size(), but"
        " index==%llu and size==%llu.", this->ull(ix), ullSize());
    SimTK_ERRCHK2(isSizeInRange(length, size_type(size()-ix)), 
        "ArrayView_<T>(index,length)", 
        "This operator requires 0 <= length <= size()-index, but"
        " length==%llu and size()-index==%llu.",this->ull(length),this->ull(size()-ix));

    return ArrayView_(data()+ix, data()+ix+length);
}
/** Same as non-const operator()(index,length); exists to provide non-operator
access to that functionality in case it is needed. **/
ArrayView_ updSubArray(index_type index, size_type length)
{   return (*this)(index,length); }
/*@}    End of element access. **/


//------------------------------------------------------------------------------
/** @name                       Iterators

These methods deal in iterators, which are STL generalized pointers. For this
class, iterators are just ordinary pointers to T, and you may depend on that.
By necessity, reverse iterators can't be just pointers; however, they contain 
an ordinary iterator (i.e. a pointer) that can be obtained by calling the 
reverse iterator's base() method. **/
/*@{*/

/** Return a const pointer to the first element of this array if any, otherwise
end(), which may be null (0) in that case but does not have to be. This method
is from the proposed C++0x standard; there is also an overloaded begin() from
the original standard that returns a const pointer. **/
SimTK_FORCE_INLINE const T* cbegin() const {return this->CBase::cbegin();}
/** The const version of begin() is the same as cbegin(). **/
SimTK_FORCE_INLINE const T* begin() const {return this->CBase::cbegin();}
/** Return a writable pointer to the first element of this array if any,
otherwise end(). If the array is empty, this \e may return null (0) but does 
not have to -- the only thing you can be sure of is that begin() == end() for 
an empty array. **/
SimTK_FORCE_INLINE T* begin() {return const_cast<T*>(this->CBase::cbegin());}

/** Return a const pointer to what would be the element just after the last one
in the array; this may be null (0) if there are no elements but doesn't have to
be. This method is from the proposed C++0x standard; there is also an 
overloaded end() from the original standard that returns a const pointer. **/
SimTK_FORCE_INLINE const T* cend() const {return this->CBase::cend();}
/** The const version of end() is the same as cend(). **/
SimTK_FORCE_INLINE const T* end() const {return this->CBase::cend();}
/** Return a writable pointer to what would be the element just after the last
one in this array. If the array is empty, this \e may return null (0) but does 
not have to -- the only thing you can be sure of is that begin()==end() for an 
empty array. **/
SimTK_FORCE_INLINE T* end() {return const_cast<T*>(this->CBase::cend());}

/** Return a const reverse iterator pointing to the last element in the array 
or crend() if the array is empty. **/
const_reverse_iterator crbegin() const 
{   return this->CBase::crbegin(); }
/** The const version of rbegin() is the same as crbegin(). **/
const_reverse_iterator rbegin() const 
{   return this->CBase::crbegin(); } 
/** Return a writable reverse iterator pointing to the last element in the
array or rend() if the array is empty. **/
reverse_iterator rbegin() {return reverse_iterator(end());}

/** Return the past-the-end reverse iterator that tests equal to a reverse
iterator that has been incremented past the front of the array. You cannot 
dereference this iterator. **/
const_reverse_iterator crend() const 
{   return this->CBase::crend(); }
/** The const version of rend() is the same as crend(). **/
const_reverse_iterator rend() const 
{   return this->CBase::crend(); }
/** Return a writable past-the-end reverse iterator that tests equal to a 
reverse iterator that has been incremented past the front of the array. You
cannot dereference this iterator. **/
reverse_iterator rend() {return reverse_iterator(begin());}

/** Return a const pointer to the first element of the array, or possibly
(but not necessarily) null (0) if the array is empty.
@note
    cdata() does not appear to be in the C++0x standard although it would seem
    obvious in view of the cbegin() and cend() methods that had to be added. 
    The C++0x overloaded const data() method is also available. **/
SimTK_FORCE_INLINE const T* cdata() const {return this->CBase::cdata();}
/** The const version of the data() method is identical to cdata().
@note This method is from the proposed C++0x std::vector. **/
SimTK_FORCE_INLINE const T* data() const {return this->CBase::cdata();}
/** Return a writable pointer to the first allocated element of the array, or
a null pointer if no space is associated with the array.
@note This method is from the proposed C++0x std::vector. **/
SimTK_FORCE_INLINE T* data() {return const_cast<T*>(this->CBase::cdata());}
/*@}    End of iterators. */


//------------------------------------------------------------------------------
/** @name                   Size and capacity 

These methods report the number of elements (size) or the amount of allocated 
heap space (capacity) or both but cannot be used to change size. **/
/*@{*/

// Note: these have to be explicitly forwarded to the base class methods
// in order to keep gcc from complaining. Note that the "this->" is 
// apparently necessary in order to permit delayed definition of templatized 
// methods. Doxygen picks up the comments from the base class.

SimTK_FORCE_INLINE size_type size()      const {return this->CBase::size();}
size_type max_size()  const {return this->CBase::max_size();}
bool      empty()     const {return this->CBase::empty();}
size_type capacity()  const {return this->CBase::capacity();}
size_type allocated() const {return this->CBase::allocated();}
bool      isOwner()   const {return this->CBase::isOwner();}
/*@}    End of size and capacity. **/


//------------------------------------------------------------------------------
                                   private:
//------------------------------------------------------------------------------
// no data members are allowed

//------------------------------------------------------------------------------
//                       ARRAY VIEW ASSIGN DISPATCH
// This is the assign() implementation for ArrayView_ when the class that 
// matched the alleged InputIterator template argument turned out to be one of 
// the integral types in which case this should match the assign(n, fillValue) 
// signature.
template <class IntegralType>
void avAssignDispatch(IntegralType n, IntegralType v, TrueType isIntegralType,
                      const char*) 
{   assign(size_type(n), value_type(v)); }

// This is the assign() implementation for ArrayView_ when the class that 
// matched the alleged InputIterator template argument is NOT an integral type 
// and may very well be an iterator. 
template <class InputIterator> 
void avAssignDispatch(const InputIterator& first, const InputIterator& last1, 
                      FalseType isIntegralType, const char* methodName) 
{   avAssignIteratorDispatch(first, last1, 
        typename std::iterator_traits<InputIterator>::iterator_category(),
        methodName); }

// This is the assign() implementation for a plain input_iterator
// (i.e., not a forward, bidirectional, or random access iterator). These
// have the unfortunate property that we can't count the elements in advance.
// Here we're going to complain if there aren't enough; but will simply stop
// when we get size() elements and not insist that the input stream reached
// the supplied last1 iterator. Semantics is elementwise assignment.
template <class InputIterator>
void avAssignIteratorDispatch(const InputIterator& first, 
                              const InputIterator& last1, 
                              std::input_iterator_tag, 
                              const char* methodName) 
{
    T* p = begin();
    InputIterator src = first;
    while (src != last1 && p != end())
        *p++ = *src++; // call T's assignment operator

    // p now points just after the last element that was copied.
    const size_type nCopied = size_type(p - begin());
    SimTK_ERRCHK2_ALWAYS(nCopied == size(), methodName,
        "The supplied input_iterator provided only %llu elements but this"
        " ArrayView has a fixed size of %llu elements.",
        this->ull(nCopied), ullSize());

    // We don't care if there are still more input elements available.
}

// This is the assign() implementation that works for forward and bidirectional
// iterators, but is not used for random_access_iterators. Here we'll count
// the elements as we copy them and complain at the end if there were too
// few or too many.
template <class ForwardIterator>
void avAssignIteratorDispatch(const ForwardIterator& first, 
                              const ForwardIterator& last1,
                              std::forward_iterator_tag, 
                              const char* methodName) 
{
    T* p = begin();
    ForwardIterator src = first;
    while (src != last1 && p != end())
        *p++ = *src++; // call T's assignment operator

    // p now points just after the last element that was copied.
    const size_type nCopied = size_type(p - begin());
    SimTK_ERRCHK2_ALWAYS(nCopied == size(), methodName,
        "The supplied forward_ or bidirectional_iterator source range provided"
        " only %llu elements but this ArrayView has a fixed size of"
        " %llu elements.", this->ull(nCopied), ullSize());

    // Make sure we ran out of source elements.
    SimTK_ERRCHK1_ALWAYS(src == last1, methodName,
        "The supplied forward_ or bidirectional_iterator source range"
        " contained too many elements; this ArrayView has a fixed size of"
        " %llu elements.", ullSize());
}

// This is the assign() implementation that works for random_access_iterators
// including ordinary pointers. Here we check the number of elements in advance
// and complain if the source and destination aren't the same size. The 
// copying loop can be done faster in this case.
template <class RandomAccessIterator>
void avAssignIteratorDispatch(const RandomAccessIterator& first, 
                              const RandomAccessIterator& last1,
                              std::random_access_iterator_tag, 
                              const char* methodName) 
{
    SimTK_ERRCHK2_ALWAYS(this->isSameSize(last1-first), methodName,
        "Assignment to an ArrayView is permitted only if the source"
        " is the same size. Here the source had %llu element(s) but the"
        " ArrayView has a fixed size of %llu.", 
        this->ull(last1-first), this->ull(size()));

    SimTK_ERRCHK_ALWAYS(!this->overlapsWithData(first,last1), methodName,
        "Source range can't overlap with the destination data.");

    T* p = begin();
    RandomAccessIterator src = first;
    while (p != end())
        *p++ = *src++; // call T's assignment operator
}


//------------------------------------------------------------------------------
// The following private methods are protected methods in the ArrayViewConst_ 
// base class, so they should not need repeating here. However, we explicitly 
// forward to the base methods to avoid gcc errors. The gcc complaint
// is due to their not depending on any template parameters; the "this->"
// apparently fixes that problem.

packed_size_type psize()      const 
{   return this->CBase::psize(); }
packed_size_type pallocated() const 
{   return this->CBase::pallocated(); }

// This just cast sizes to unsigned long long so that we can do comparisons
// without getting warnings.
unsigned long long ullSize()     const 
{   return this->CBase::ullSize(); }
unsigned long long ullCapacity() const 
{   return this->CBase::ullCapacity(); }
unsigned long long ullMaxSize()  const 
{   return this->CBase::ullMaxSize(); }
// This is the index type name and is handy for error messages to explain
// why some size was too big.
const char* indexName() const   {return this->CBase::indexName();}
};





//==============================================================================
//                               CLASS Array_
//==============================================================================
/** The SimTK::Array_<T> container class is a plug-compatible replacement for 
the C++ standard template library (STL) std::vector<T> class, but with some
important advantages in performance, and functionality, and binary 
compatibility.

@tparam T 
    The type of object to be stored in this container. 
@tparam X 
    The type to be used for indexing this container, with default unsigned
    (not size_t). Any integral type may be used, as well as user types that 
    satisfy the requirements discussed with class ArrayIndexTraits. 

@par Performance:
There are several performance and memory footprint problems with the C++ 
standard STL design in general, and with Microsoft's implementation in 
particular, that are addressed here. Microsoft in its wisdom decided that STL 
containers should still do runtime range checks in Release builds for safety, 
but that makes them too slow for use in some high-performance contexts (and 
also breaks the promise of generic programming but that's another rant). In 
practice, VC++9 std::vector runs about half speed for simple operations like 
indexing and push_back. Attempting to disable these runtime checks with 
_SECURE_SCL breaks binary compatibility. In contrast the performance of this 
Array_<T> class on any platform is indistinguishable from what you would get 
by managing your own heap-allocated arrays.

@par
Regarding memory footprint, the typical implementation of std::vector uses
three pointers: 12 bytes for 32 bit machines; 24 bytes for 64 bit machines.
Microsoft somehow manages to trump this with 20 to 24 bytes on a 32 bit
machine -- I don't know what they do on a 64 bit machine but I'm not 
optimistic! Array_ instead uses one pointer and two lengths for a total size 
as little as 8 bytes on 32 bits and 16 on 64 bits; see below for details.

@par
Some nuts and bolts:

- We promise that no heap allocation occurs when an empty Array_<T> object 
  is declared (that is, when an Array_<T> is default-constructed); in
  that case both begin() and end() are null.
- Array_<T> methods are extremely fast in Release builds with zero overhead, 
  inline, unchecked methods. The implementations of inline methods are kept
  small to ensure that they are actually inlined in practice; and generated 
  assembly code was examined to make sure.
- There are some dangerous extensions provided that permit the expert user
  to construct objects directly into the array without having to copy them,
  a big win for complicated objects and even bigger for those that don't
  have copy constructors!
- There is a constant-time eraseFast() method you can use if you don't mind the
  array being reordered after the erase. This avoids the extremely expensive
  "compress" activity required by the standard erase() method.
- The optional index-type template parameter can be used to reduce the memory
  footprint to as little as 8 bytes on a 32 bit machine (e.g., a 32 bit 
  pointer and two shorts).
- The default size_type for an Array_<T> is a 32-bit unsigned integer rather 
  than a size_t. On a 64-bit machine that keeps the overhead down substantially
  since the structure is then one 64-bit pointer and two 32-bit integers, 
  fitting tightly into a cleanly alignable 16 bytes.


@par Functionality:
For the most part Array_<T> is a plug-compatible replacement for std::vector<T>,
and everything that both classes can do is done with an identical API. However,
there are a few additions and subtractions:

- This class always uses the default new/delete allocator; there is no option
  to specify your own as there is in std::vector.
- Instead of an allocator, the second template argument X to Array_<T,X> is an 
  optional index type which can be used to provide type-safe indexing (i.e. the
  array can only be indexed by indices of a particular type, like 
  MobilizedBodyIndex). This has zero performance cost if the index is an 
  integral type or class consisting of only an integral value such as those
  produced by the SimTK_DEFINE_UNIQUE_INDEX_TYPE macro.
- You can create uninitialized slots in the array and construct directly into
  them rather than having to construct a temporary object which must then be
  copied into the array.
- You can create Array_<T> objects that reference existing data, including
  the contents of std::vectors.
- Where possible this class implements the new std::vector features proposed
  for the C++0x standard (see below).

@par Compatibility:
Included here are binary compatibility issues and compatibility with the C++
standard STL objects.

- Most important, it is safe to pass an Array_<T> through an API to a binary 
  library without worrying about compiler version or Release/Debug compatibility
  issues. For a given compiler (e.g. gcc or Microsoft cl) and word size (64 bit
  vs. 32 bit), Array_<T> has an extremely stable memory layout that is preserved 
  across compiler versions, and between Release and Debug builds. This allows us
  to use Array_<T> in the SimTK API where use of std::vector<T> would be 
  desirable but problematic.
- It supports all standard types, methods, iterators, and operators of the 
  C++98 standard std::vector and the C++0x proposed improvements other than
  those requiring rvalue references, so it works smoothly with all STL 
  containers and algorithms.
- It is convertible to and from std::vector, usually without copying the 
  elements. It is easy to provide APIs that accept either Array_<T> or 
  std::vector<T>; the std::vector's data is referenced by an Array_ handle
  that is used to convey the data across the API without binary compatibility
  problems.

@see Array_, ArrayViewConst_, ArrayIndexTraits **/
template <class T, class X> class Array_ : public ArrayView_<T,X> {
    typedef ArrayView_<T,X>      Base;
    typedef ArrayViewConst_<T,X> CBase;
public:


//------------------------------------------------------------------------------
/** @name                        Typedefs

Types required of STL containers, plus index_type which is an extension, and
packed_size_type which is an implementation detail. **/

// Doxygen picks up individual descriptions from the base class.
/*{*/
typedef T                                               value_type;
typedef X                                               index_type;
typedef T*                                              pointer;
typedef const T*                                        const_pointer;
typedef T&                                              reference;
typedef const T&                                        const_reference;
typedef T*                                              iterator;
typedef const T*                                        const_iterator;
typedef std::reverse_iterator<iterator>                 reverse_iterator;
typedef std::reverse_iterator<const_iterator>           const_reverse_iterator;
typedef typename ArrayIndexTraits<X>::size_type         size_type;
typedef typename ArrayIndexTraits<X>::difference_type   difference_type;
typedef typename ArrayIndexPackType<size_type>::packed_size_type 
                                                        packed_size_type;
/*}*/

//------------------------------------------------------------------------------
/** @name        Construction, conversion and destruction

A variety of constructors are provided for this class, including all those
required by the C++ standard for std::vector implementations, plus additional
ones providing smooth conversions between Array_<T> and std::vector<T> objects.
**/
/*{*/

/** Default constructor allocates no heap space and is very fast. **/
Array_() : Base() {}

/** Construct an array containing \a n default-constructed elements. T's default 
constructor (if any) is called exactly \a n times. If \a n is zero no heap space 
will be allocated; although in that case it is preferable to use the default 
constructor if you can since that will be somewhat faster. **/
explicit Array_(size_type n) : Base() {
    SimTK_SIZECHECK(n, max_size(), "Array_<T>::ctor(n)");
    allocateNoConstruct(n);
    defaultConstruct(data(), data()+n);
    setSize(n);
}

/** Construct an array containing \a n elements each set to a copy of the given 
initial value. T's copy constructor will be called exactly \a n times. If \a n
is zero no space will be allocated. **/
Array_(size_type n, const T& initVal) : Base() {
    SimTK_SIZECHECK(n, max_size(), "Array_<T>::ctor(n,T)");
    setSize(n);
    allocateNoConstruct(size());
    fillConstruct(begin(), cend(), initVal);
}
/** Construct an Array_<T> from a range [first,last1) of values identified by a 
pair of iterators. 
@note
The standard requires that if an integral type matches this signature, it must
behave as the Array_(size_type,value_type) constructor.
@par Complexity:
The performance of this constructor depends on the type
of iterator: 
- random_access_iterator: n=(last1-first); a single space allocation; 
  n calls to T's copy constructor. 
- forward or bidirectional iterator: must increment from first to last1 to
  determine n; otherwise same as random access.
- input iterator: can't determine n in advance; expect log n reallocations
  during construction as we "push back" one input element at a time.
**/
template <class InputIter>
Array_(const InputIter& first, const InputIter& last1) : Base() {
    ctorDispatch(first,last1,typename IsIntegralType<InputIter>::Result());
}

/** Construct an Array_<T> from a range [first,last1) of values identified by a 
pair of ordinary pointers to elements of type T2 (where T2 might be the same as
T but doesn't have to be). This is templatized so can be used with any source 
type T2 for which there is a working conversion constructor T(T2), provided
that the number of source elements does not exceed the array's max_size(). **/
template <class T2>
Array_(const T2* first, const T2* last1) : Base() {
    SimTK_ERRCHK((first&&last1)||(first==last1), "Array_<T>(first,last1)", 
        "Pointers must be non-null unless they are both null.");
    SimTK_ERRCHK3(this->isSizeOK(last1-first), "Array_<T>(first,last1)",
        "Source has %llu elements but this array is limited to %llu"
        " elements by its index type %s.",
        this->ull(last1-first), ullMaxSize(), indexName());

    setSize(size_type(last1-first));
    allocateNoConstruct(size());
    copyConstruct(begin(), cend(), first);
}

/** Construct an Array_<T> by copying from an std::vector<T2>, where T2 may
be the same type as T but doesn't have to be. This will work as long as the 
size of the vector does not exceed the array's max_size(), and provided there 
is a working T(T2) conversion constructor. **/
template <class T2>
explicit Array_(const std::vector<T2>& v) : Base() { 
    if (v.empty()) return;

    SimTK_ERRCHK3(this->isSizeOK(v.size()), "Array_<T>::ctor(std::vector<T2>)",
        "The source std::vector's size %llu is too big for this array which"
        " is limited to %llu elements by its index type %s.",
        this->ull(v.size()), ullMaxSize(), indexName());

    // Call the above constructor, making sure to use pointers into the
    // vector's data rather than the iterators begin() and end() in case
    // they are different types.
    new (this) Array_(&v.front(), (&v.back())+1);
}

/** Copy constructor allocates exactly as much memory as is in use in the 
source (not its capacity) and copy constructs the elements so that T's copy 
constructor will be called exactly src.size() times. If the source is empty, 
no heap space will be allocated. **/
Array_(const Array_& src) : Base() {
    setSize(src.size());
    allocateNoConstruct(size());
    copyConstruct(begin(), cend(), src.data());
}

/** Construct this Array_<T,X> as a copy of another Array_<T2,X2> where T2!=T
or X2!=X. This will work as long as the source is not larger than will fit 
here, and as long as the source element type T2 is assignment compatible with 
this array's element type T. One of T's constructors will be called exactly 
src.size() times; the particular constructor is whichever one best matches 
T(T2). **/
template <class T2, class X2>
Array_(const Array_<T2,X2>& src) : Base() {
    new (this) Array_(src.begin(), src.cend()); // see above
}

/** Construct an Array_<T> by referencing (sharing) a given range of data
[first,last1), without copying that data; better to use the corresponding
ArrayView_<T> constructor if you can. This is very fast but can be 
dangerous -- it is most useful for argument passing where the array handle 
will be discarded immediately after use. Note that this is available only if 
you have write access to the data because there is no way to construct 
a non-writable array. This will work as long as the size of the data does 
not exceed the array's max_size. The resulting array object is not resizeable
but can be used to read and write elements of the original data. The
array is invalid if the original data is destructed or resized, but there is
no way for the array class to detect that.

@remarks
  - If the source data is empty, the resulting array will also 
    be empty and will look just like a default-constructed array. It will
    therefore not have any connection to the source and will be an
    ordinary resizable array.
  - This is quite dangerous to use since the connection between the array and
    the data is tenuous and subject to the data remaining untouched during
    the lifetime of the array handle. There is no reference counting;
    destructing the original data would leave the array referring to garbage. 
    Be careful!
  - You can break the connection between the array and the data it was 
    constructed from by calling deallocate().

@par Complexity:
    Dirt cheap. There will be no construction, destruction, or heap allocation
    performed.
@see deallocate() **/
Array_(T* first, const T* last1, const DontCopy&) : Base(first,last1) {}

/** Construct an Array_<T> by referencing (sharing) the data in an 
std::vector<T>, without copying the data; better to use the ArrayView_<T>
constructor instead if you can. This is very fast but can be 
dangerous -- it is most useful for argument passing where the array handle 
will be discarded immediately after use. Note that this is available only if 
you have write access to the std::vector because there is no way to construct 
a non-writable array. This will work as long as the size of the vector does 
not exceed the array's max_size. The resulting array object is not resizeable
but can be used to read and write elements of the original std::vector. The
array is invalid if the original std::vector is destructed or resized.

@remarks
  - If the source std::vector is empty, the resulting array will also 
    be empty and will look just like a default-constructed array. It will
    therefore not have any connection to the source vector and will be an
    ordinary resizable array.
  - This is quite dangerous to use since the connection between the array and
    the vector is tenuous and subject to the vector remaining untouchged during
    the lifetime of the array handle. There is no reference counting;
    destructing the vector leaves the array referring to garbage. Be careful!
  - You can break the connection between the array and the vector it was 
    constructed from by calling deallocate().

@par Complexity:
    Dirt cheap. There will be no construction, destruction, or heap allocation
    performed.
@see deallocate() **/
template <class A>
Array_(std::vector<T,A>& v, const DontCopy&) : Base(v) {}

/** The destructor performs a deallocate() operation which may result in 
element destruction and freeing of heap space; see deallocate() for more
information. @see deallocate() **/
~Array_() {
    deallocate();
}

/** Empty this array of its contents, returning the array to its 
default-constructed, all-zero state. If this array is the owner of its data,
the destructor (if any) is called for each data element and the array's
allocated heap space is freed. If it is a non-owner the array handle is
cleaned out using disconnect() but the referenced data is untouched.
@note There is no equivalent to this method for std::vector.
@return A reference to the now-empty, default-constructed array, ready for
reassignment. **/
Array_& deallocate() {
    if (allocated()) { // owner with non-zero allocation
        clear(); // each element is destructed; size()=0; allocated() unchanged
        deallocateNoDestruct(); // free data(); allocated()=0
    }
    this->Base::disconnect(); // clear the handle
    return *this;
}

/*@}    End of constructors, etc. **/


//------------------------------------------------------------------------------
/** @name           Assignment methods and operators

These methods put new data values in an existing array, but the meaning of
assignment is subtly different for resizeable (owner) arrays and fixed 
(non-owner) arrays. The standard std::vector type is always an owner so the
non-owner description here is an extension applying only to Array_.

For the normal case of resizeable arrays, assignment does not have an 
elementwise definition because the source will typically have a different 
number of elements than the array's current size. So regardless of the actual 
numbers, assignment in the resizeable case is defined as it is for std::vector:
first clear the array by erasing (destructing) all the current elements in the
array, then reserve sufficient heap space to hold a copy of the source, then 
use appropriate constructors of type T (most commonly T's copy constructor 
T(T)) to initialize each element to be a copy of the corresponding source 
element. T's assignment operators are never used in this case.

For fixed arrays, the source must have the same number of elments as are 
currently in the array and the meaning is conventional elementwise assignment;
that is, an appropriate assignment operator of type T (most commonly T's copy 
assignment operator T=T) is used to change the value of each existing element. 

So there are different requirements on the value type T for owner and non-owner
assignments to type T2: for owner assignment T must have a constructor T(T2)
available; for non-owner assignment, T must have an assignment operator T=T2 
available; .

@remarks
- When reallocating the destination array, we may reuse the existing heap 
allocation if it is sufficient and not \e too big; otherwise we'll reallocate 
before copying. 
- The fill() method here has elementwise assignment semantics regardless of
whether the array is an owner or non-owner. **/
/*@{*/

/** Set this array to be \a n copies of the supplied \a fillValue. Note that 
this serves to allow fill from an object whose type T2 is different from T, as
long as there is a constructor T(T2) that works since that can be invoked 
(implicitly or explicitly) to convert the T2 object to type T prior to the
call. If this is a non-owner array then \a n must be the same as the current
size(); consider using the fill() method instead.
@param[in] n            The number of elements to be in the result.
@param[in] fillValue    The value to which to initialize each element.

@pre \a n <= max_size()
@pre for non-owner, n==size()
@par Complexity:
For a non-owner with \a n==size(), there will be exactly \a n calls to T's
copy assignment operator. For an owner, there will be size() calls to T's
destructor (if it has one), possibly a heap reallocation (but with no element
copying), followed by \a n calls to T's copy constructor. 
@see fill() **/ 
void assign(size_type n, const T& fillValue) {
    SimTK_ERRCHK3(this->isSizeOK(n), "Array_<T>::assign(n,value)",
        "Requested size %llu is too big for this array which is limited"
        " to %llu elements by its index type %s.",
        this->ull(n), ullMaxSize(), indexName());

    SimTK_ERRCHK2(isOwner() || n==size(), "Array_<T>::assign(n,value)",
        "Requested size %llu is not allowed because this is a non-owner"
        " array of fixed size %llu.", this->ull(n), this->ull(size()));

    if (!isOwner())
        this->Base::fill(fillValue);
    else {
        clear(); // all elements destructed; allocation unchanged
        reallocateIfAdvisable(n); // change size if too small or too big
        fillConstruct(data(), cdata()+n, fillValue);
        setSize(n);
    }
}

/** Assign all current elements of the array to the same \a fillValue. This is
similar to assign(size(),fillValue) but the semantics are subtly different.
Here we use repeated application of T's copy assignment operator T=fillValue,
whereas the assign() semantics are to first destruct all the existing elements,
then allocate if necessary, then use the copy constructor to initialize the
new elements. Note that you can use this to fill from a source type T2 that
is different from T as long as there exists a suitable constructor T(T2) that
can be used to create the type T \a fillValue from the original T2 source.
@note Unlike other assignment methods, the behavior of fill() is identical for
owner and non-owner arrays.

@param[in] fillValue    The value to which all existing elements are set.
@par Complexity:
Just size() calls to T's copy assignment operator. **/
void fill(const T& fillValue) {this->Base::fill(fillValue);}


/** Assign to this array to to make it a copy of the elements in range 
[first,last1) given by ordinary pointers. It is not allowed for this range to 
include any of the elements currently in the array. The source elements can be 
of a type T2 that may be the same or different than this array's element type 
T as long as there is a working constructor T(T2) (for owner arrays) or a 
working assignment operator T=T2 (for non-owner arrays). Note that although the
source arguments are pointers, those may be iterators for some container 
depending on implementation details of the container. Specifically, any 
Array_<T2>::iterator or const_iterator is an ordinary pointer.

@param[in] first    A pointer to the first source element to be copied.
@param[in] last1    A pointer to one element past the last source element.

@par Complexity:
For non-owner arrays, n=last1-first must equal the current size() in which
case there will be exactly size() calls to the T=T2 assignment operator.
For owner arrays, say the array initially has capacity c, and the 
source provides n new elements. If type T has a destructor, it will be called 
exactly size() times. Reallocation will then occur if c < n and may occur if 
c >> n to avoid leaving a lot of unused space. Then the constructor T(T2) will 
be called exactly n times. **/
template <class T2>
void assign(const T2* first, const T2* last1) {
    const char* methodName = "Array_<T>::assign(T2* first, T2* last1)";
    SimTK_ERRCHK((first&&last1)||(first==last1), methodName, 
        "Pointers must be non-null unless they are both null.");
    SimTK_ERRCHK(!this->overlapsWithData(first,last1), methodName,
        "Source range can't overlap the current array contents.");
    // Pointers are random access iterators.
    assignIteratorDispatch(first,last1,std::random_access_iterator_tag(),
                           methodName);
}


/** Assign this array from a range [first,last1) given by non-pointer 
iterators. See the assign(first,last1) method with pointer arguments for a
relevant discussion.

@remarks
  - For a non-owner array this is only allowed if we can calculate the number of
    source elements, and if that number is exactly the same as the current 
    size().
  - See Complexity discussion below for behavior for the different kinds of
    iterators that might be supplied.
  - It is not permitted for any of the source elements to overlap in memory
    with the initial contents of the array.

@param[in]      first    
    An iterator pointing to the first source element to be copied.
@param[in]      last1    
    A iterator pointing one element past the last source element.

@pre last1-first <= max_size()
@pre for non-owner array, last1-first == size()
@par Complexity:
For a non-owner array, this is only allowed if n=last1-first equals the
current size(), in which case we'll perform exactly n calls to the appropriate
assignment operator of element type T. For owner arrays, if we can determine 
how many elements n=last1-first the source contains in advance, we'll do only 
a single allocation here and call one of T's constructors exactly n times after 
just size() destructor calls needed to erase the original data. If the 
iterators are random access iterators, calculating n is a fast constant-time
operation. For forward or bidirectional iterators, we have to advance through
the iterators once to count the source elements prior to allocating space,
adding an O(n) cost. For input iterators, we can't count them in advance so
we just have to add elements as we find them using push_back() meaning we may
need to reallocate log(n) times, calling the destructor and copy constructor
each time to move the elements around. 
@see assign(T2* first, T2* last1) **/
template <class Iter>
void assign(const Iter& first, const Iter& last1) {
    assignDispatch(first,last1,typename IsIntegralType<Iter>::Result(),
                   "Array_<T>::assign(Iter first, Iter last1)");
}

/** Copy assignment operator destructs the current contents of this array and 
then makes it a copy of the source array by repeated calls to the element 
type's copy constructor. At most one reallocation of heap space occurs that 
may result in this array having a larger or smaller capacity, although of 
course it will be at least as large as the source. **/
Array_& operator=(const Array_& src) {
    if (this != &src)
        assignIteratorDispatch(src.begin(), src.end(), 
                               std::random_access_iterator_tag(),
                               "Array_<T>::operator=(Array_<T>)");
    return *this;
}

/** This is assignment from a source array whose element type T2 and/or index 
type X2 are different from this array's T and X. This will work as long as 
this array can accommodate all the elements in the source and T2 is assignment
compatible with T. See discussion for the copy assignment operator for more 
information. */
template <class T2, class X2>
Array_& operator=(const Array_<T2,X2>& src) {
    assignIteratorDispatch(src.begin(), src.end(), 
                           std::random_access_iterator_tag(),
                           "Array_<T>::operator=(Array_<T2,X2>)");
    return *this;
}


/** This is assignment from a source std::vector<T2>. This will work as long as 
this array can accommodate all the elements in the source and T2 is assignment
compatible with T. See discussion for the copy assignment operator for more 
information. */
template <class T2, class A>
Array_& operator=(const std::vector<T2,A>& src) {
    assignIteratorDispatch(src.begin(), src.end(), 
                           std::random_access_iterator_tag(),
                           "Array_<T>::operator=(std::vector)");
    return *this;
}

/** This is a specialized algorithm providing constant time exchange of data 
with another array that has identical element and index types. This is \e much 
faster than using the std::swap() algorithm on the arrays since that would
involve O(n) copying operations. This method makes no calls to any constructors
or destructors. This is allowable even for non-owner arrays; the non-owner
attribute will follow the non-owned data. **/
void swap(Array_& other) {
    T* const pTmp=data(); setData(other.data()); other.setData(pTmp);
    size_type nTmp=size(); setSize(other.size()); other.setSize(nTmp);
    nTmp=allocated(); setAllocated(other.allocated()); other.setAllocated(nTmp);
}

/** This dangerous extension allows you to supply your own already-allocated
heap space for use by this array, which then becomes the owner of the supplied
heap space. Any memory currently associated with the array is deallocated; 
see deallocate() for more information. 
@see deallocate(), shareData() **/
Array_& adoptData(T* newData, size_type dataSize, 
                  size_type dataCapacity) 
{
    SimTK_SIZECHECK(dataCapacity, max_size(), "Array_<T>::adoptData()");
    SimTK_ERRCHK2(dataSize <= dataCapacity, "Array_<T>::adoptData()", 
        "Specified data size %llu was greater than the specified data"
        " capacity of %llu.", this->ull(dataSize), this->ull(dataCapacity));
    SimTK_ERRCHK(newData || dataCapacity==0, "Array_<T>::adoptData()",
        "A null data pointer is allowed only if the size and capacity are"
        " specified as zero.");
    SimTK_ERRCHK(!this->overlapsWithData(newData, newData+dataSize), 
        "Array_<T>::adoptData()",
        "The new data can't overlap with the old data.");

    deallocate();
    setData(newData);
    setSize(dataSize);
    setAllocated(dataCapacity);
    return *this;
}
/** A variant of adoptData() that assumes the capacity is the same as the
current size. @see adoptData(data,size,capacity) **/
Array_& adoptData(T* newData, size_type dataSize) 
{   return adoptData(newData, dataSize, dataSize); }


/** This dangerous extension allows you to make this array handle refer to
someone else's data without copying it. Any memory currently associated
with the array is deallocated; see deallocate() for more information. This
method makes the array a fixed-size, non-owner array that cannot be 
reallocated, and no element destruction nor heap deallocation will occur when
the handle is subsequently destructed or deallocated.
@note
  - A null (0) pointer is allowed for the pointer as long as \a dataSize==0,
    however in that case the array handle ends up deallocated (that is, 
    indistinguishable from a default-constructed array) so is resizeable.
  - This is implemented by setting the nAllocated data member to zero while
    the nUsed data member is set to the given \a dataSize.
@see deallocate(), adoptData() **/
Array_& shareData(T* newData, size_type dataSize) {
    SimTK_SIZECHECK(dataSize, max_size(), "Array_<T>::shareData()");
    SimTK_ERRCHK(newData || dataSize==0, "Array_<T>::shareData()",
        "A null data pointer is allowed only if the size is zero.");
    SimTK_ERRCHK(!this->overlapsWithData(newData, newData+dataSize), 
        "Array_<T>::shareData()",
        "The new data can't overlap with the old data.");

    deallocate();
    setData(newData);
    setSize(dataSize);
    setAllocated(0); // indicates shared data
    return *this;
}

/** Same as shareData(data,size) but uses a pointer range [first,last1) to
identify the data to be referenced. **/
Array_& shareData(T* first, const T* last1) {
    SimTK_ERRCHK3(this->isSizeOK(last1-first), "Array_<T>::shareData(first,last1)",
        "Requested size %llu is too big for this array which is limited"
        " to %llu elements by its index type %s.",
        this->ull(last1-first), ullMaxSize(), indexName());
    return shareData(first, size_type(last1-first));
}

/*@}    End of assignment. **/


//------------------------------------------------------------------------------
/** @name                   Size and capacity 

These methods examine and alter the number of elements (size) or the amount of 
allocated heap space (capacity) or both. **/
/*@{*/

// Note: these have to be explicitly forwarded to the base class methods
// in order to keep gcc from complaining. Note that the "this->" is 
// apparently necessary in order to permit delayed definition of templatized 
// methods.

/** Return the current number of elements stored in this array. **/
SimTK_FORCE_INLINE size_type size() const {return this->CBase::size();}
/** Return the maximum allowable size for this array. **/
size_type max_size() const {return this->CBase::max_size();}
/** Return true if there are no elements currently stored in this array. This
is equivalent to the tests begin() == end() or size()==0. **/
bool empty() const {return this->CBase::empty();}
/** Return the number of elements this array can currently hold without
requiring reallocation. The value returned by capacity() is always greater 
than or equal to size(), even if the data is not owned by this array in
which case we have capacity() == size() and the array is not reallocatable. **/
size_type capacity() const {return this->CBase::capacity();}

/** Change the size of this Array, preserving all the elements that will still 
fit, and default constructing any new elements that are added. This is not
allowed for non-owner arrays unless the requested size is the same as the 
current size. **/
void resize(size_type n) {
    if (n == size()) return;

    SimTK_ERRCHK2(isOwner(), "Array_<T>::resize(n)",
        "Requested size change to %llu is not allowed because this is a"
        " non-owner array of fixed size %llu.", this->ull(n), this->ull(size()));

    if (n < size()) {
        erase(data()+n, cend());
        return;
    }
    // n > size()
    reserve(n);
    defaultConstruct(data()+size(), cdata()+n); // data() has changed
    setSize(n);
}

/** Change the size of this array, preserving all the elements that will still 
fit, and initializing any new elements that are added by repeatedly copy-
constructing from the supplied value. This is not allowed for non-owner arrays
unless the requested size is the same as the current size. **/
void resize(size_type n, const T& initVal) {
    if (n == size()) return;

    SimTK_ERRCHK2(isOwner(), "Array_<T>::resize(n,value)",
        "Requested size change to %llu is not allowed because this is a"
        " non-owner array of fixed size %llu.", this->ull(n), this->ull(size()));

    if (n < size()) {
        erase(data()+n, cend());
        return;
    }
    // n > size()
    reserve(n);
    fillConstruct(data()+size(), cdata()+n, initVal);
    setSize(n);
}

/** Ensure that this array has enough allocated capacity to hold the indicated 
number of elements. No heap reallocation will occur after this until the array
is grown beyond this capacity, meaning that adding elements will not invalidate
any iterators or element addresses until that point. This method will never 
reduce the capacity of the array. It is OK to call this on a non-owner array
as long as you are not asking for an increase in capacity. **/
void reserve(size_type n) {
    if (capacity() >= n)
        return;

    SimTK_ERRCHK2(isOwner(), "Array_<T>::reserve()",
        "Requested capacity change to %llu is not allowed because this is a"
        " non-owner array of fixed size %llu.", this->ull(n), this->ull(size()));

    T* newData = allocN(n); // no construction yet
    copyConstructThenDestructSource(newData, newData+size(), data());
    freeN(data());
    setData(newData);
    setAllocated(n);
}

/** Request that the capacity of this array be reduced to the minimum necessary
to hold the number of elements currently in use. In practice no shrinkage will 
occur if the current size is just slightly too big, unless the current size is
exactly zero in which case we guarantee to deallocate all heap space associated
with this array leaving a null data pointer and begin()==end()==0, exactly as
though the array had just been default-constructed. Otherwise you can check
capacity() afterwards to see what happened. If the capacity() is reduced by 
this method, then all the elements will have been moved to new locations so 
existing iterators and references into the array will become invalid.

@note
  - This method is from the proposed C++0x standard for std::vector, except for
    the guaranteed behavior for a zero-size container.
  - It is OK to call this on a non-owner array but it has no effect since
    capacity()==size() already in that case.

@par Complexity:
    If the capacity is reduced, there will be one call to T's copy constructor 
    and destructor (if any) for each element currently in the array. Otherwise
    this is very fast. **/
void shrink_to_fit() {
    // Allow 25% slop, but note that if size()==0 this will always reallocate
    // unless capacity is already zero.
    if (capacity() - size()/4 <= size()) // avoid overflow if size() near max
        return;
    T* newData = allocN(size());
    copyConstructThenDestructSource(newData, newData+size(), data());
    deallocateNoDestruct(); // data()=0, allocated()=0, size() unchanged
    setData(newData);
    setAllocated(size());
}

/** Return the amount of heap space owned by this array; this is the same
as capacity() for owner arrays but is zero for non-owners. 
@note There is no equivalent of this method for std::vector. **/
size_type allocated() const 
{   return this->CBase::allocated(); }
/** Does this array own the data to which it refers? If not, it can't be
resized, and the destructor will not free any heap space nor call any element
destructors. If the array does not refer to \e any data it is considered to be
an owner and it is resizeable. 
@note There is no equivalent of this method for std::vector. **/
bool isOwner() const {return this->CBase::isOwner();}
/*@}    End of size and capacity. **/


//------------------------------------------------------------------------------
/** @name                       Iterators

These methods deal in iterators, which are STL generalized pointers. For this
class, iterators are just ordinary pointers to T, and you may depend on that.
By necessity, reverse iterators can't be just pointers; however, they contain 
an ordinary iterator (i.e. a pointer) that can be obtained by calling the 
reverse iterator's base() method. **/
/*@{*/

/** Return a const pointer to the first element of this array if any, otherwise
cend(), which may be null (0) in that case but does not have to be. This method
is from the proposed C++0x standard; there is also an overloaded begin() from
the original standard that returns a const pointer. **/
SimTK_FORCE_INLINE const T* cbegin() const {return this->CBase::cbegin();}
/** The const version of begin() is the same as cbegin(). **/
SimTK_FORCE_INLINE const T* begin() const {return this->CBase::cbegin();}
/** Return a writable pointer to the first element of this array if any,
otherwise end(). If the array is empty, this \e may return null (0) but does 
not have to -- the only thing you can be sure of is that begin() == end() for 
an empty array. **/
SimTK_FORCE_INLINE T* begin() {return this->Base::begin();}

/** Return a const pointer to what would be the element just after the last one
in the array; this may be null (0) if there are no elements but doesn't have to
be. This method is from the proposed C++0x standard; there is also an 
overloaded end() from the original standard that returns a const pointer. **/
SimTK_FORCE_INLINE const T* cend() const {return this->CBase::cend();}
/** The const version of end() is the same as cend(). **/
SimTK_FORCE_INLINE const T* end() const {return this->CBase::cend();}
/** Return a writable pointer to what would be the element just after the last
one in this array. If the array is empty, this \e may return null (0) but does 
not have to -- the only thing you can be sure of is that begin()==end() for an 
empty array. **/
SimTK_FORCE_INLINE T* end() {return this->Base::end();}

/** Return a const reverse iterator pointing to the last element in the array 
or crend() if the array is empty. **/
const_reverse_iterator crbegin() const 
{   return this->CBase::crbegin(); }
/** The const version of rbegin() is the same as crbegin(). **/
const_reverse_iterator rbegin() const 
{   return this->CBase::crbegin(); } 
/** Return a writable reverse iterator pointing to the last element in the
array or rend() if the array is empty. **/
reverse_iterator rbegin() {return this->Base::rbegin();}

/** Return the past-the-end reverse iterator that tests equal to a reverse
iterator that has been incremented past the front of the array. You cannot 
dereference this iterator. **/
const_reverse_iterator crend() const 
{   return this->CBase::crend(); }
/** The const version of rend() is the same as crend(). **/
const_reverse_iterator rend() const 
{   return this->CBase::crend(); }
/** Return a writable past-the-end reverse iterator that tests equal to a 
reverse iterator that has been incremented past the front of the array. You
cannot dereference this iterator. **/
reverse_iterator rend() {return this->Base::rend();}

/** Return a const pointer to the first element of the array, or possibly
(but not necessarily) null (0) if the array is empty.
@note
    cdata() does not appear to be in the C++0x standard although it would seem
    obvious in view of the cbegin() and cend() methods that had to be added. 
    The C++0x overloaded const data() method is also available. **/
SimTK_FORCE_INLINE const T* cdata() const {return this->CBase::cdata();}
/** The const version of the data() method is identical to cdata().
@note This method is from the proposed C++0x std::vector. **/
SimTK_FORCE_INLINE const T* data() const {return this->CBase::cdata();}
/** Return a writable pointer to the first allocated element of the array, or
a null pointer if no space is associated with the array.
@note This method is from the proposed C++0x std::vector. **/
SimTK_FORCE_INLINE T* data() {return this->Base::data();}
/*@}*/

/** @name                     Element access

These methods provide read and write access to individual elements, or groups
of elements, that are currently present in the array. **/
/*@{*/

/** Select an element by its index, returning a const reference. Note that only 
a value of the Array's templatized index type is allowed (default is unsigned).
This will be range-checked in a Debug build but not in Release.
@pre 0 <= \a i < size()
@par Complexity:
    Constant time. **/
SimTK_FORCE_INLINE const T& operator[](index_type i) const 
{   return this->CBase::operator[](i); }

/** Select an element by its index, returning a writable (lvalue) reference. 
Note that only a value of the Array's templatized index type is allowed 
(default is unsigned). This will be range-checked in a Debug build but not 
in Release. 
@pre 0 <= \a i < size()
@par Complexity:
    Constant time. **/
SimTK_FORCE_INLINE T& operator[](index_type i) {return this->Base::operator[](i);}

/** Same as operator[] but always range-checked, even in a Release build.  
@pre 0 <= \a i < size()
@par Complexity:
    Constant time. **/
const T& at(index_type i) const {return this->CBase::at(i);}

/** Same as operator[] but always range-checked, even in a Release build.  
@pre 0 <= \a i < size()
@par Complexity:
    Constant time. **/
T& at(index_type i) {return const_cast<T&>(this->Base::at(i));}

/** Same as the const form of operator[]; exists to provide a non-operator
method for element access in case that's needed. **/
SimTK_FORCE_INLINE const T& getElt(index_type i) const 
{   return this->CBase::getElt(i); }
/** Same as the non-const form of operator[]; exists to provide a non-operator
method for element access in case that's needed. **/
SimTK_FORCE_INLINE T& updElt(index_type i) {return this->Base::updElt(i);}

/** Return a const reference to the first element in this array, which must
not be empty.
@pre The array is not empty.
@par Complexity:
    Constant time. **/
SimTK_FORCE_INLINE const T& front() const {return this->CBase::front();} 

/** Return a writable reference to the first element in this array, which must
not be empty.
@pre The array is not empty.
@par Complexity:
    Constant time. **/
SimTK_FORCE_INLINE T& front() {return const_cast<T&>(this->Base::front());}

/** Return a const reference to the last element in this array, which must
not be empty.
@pre The array is not empty.
@par Complexity:
    Constant time. **/
SimTK_FORCE_INLINE const T& back() const {return this->CBase::back();}

/** Return a writable reference to the last element in this array, which must
not be empty.
@pre The array is not empty.
@par Complexity:
    Constant time. **/
SimTK_FORCE_INLINE T& back() {return const_cast<T&>(this->Base::back());}

/** Select a subrange of this const array by starting index and length, and
return a ArrayViewConst_ referencing that data without copying it. **/
ArrayViewConst_<T,X> operator()(index_type index, size_type length) const
{   return CBase::operator()(index,length); }
/** Same as const form of operator()(index,length); exists to provide 
non-operator access to that functionality in case it is needed. **/
ArrayViewConst_<T,X> getSubArray(index_type index, size_type length) const
{   return CBase::getSubArray(index,length); }

/** Select a subrange of this array by starting index and length, and
return an ArrayView_ referencing that data without copying it. **/
ArrayView_<T,X> operator()(index_type index, size_type length)
{   return Base::operator()(index,length); }
/** Same as non-const operator()(index,length); exists to provide non-operator
access to that functionality in case it is needed. **/
ArrayView_<T,X> updSubArray(index_type index, size_type length)
{   return Base::updSubArray(index,length); }
/*@}    End of element access. **/


//------------------------------------------------------------------------------
/**@name                Element insertion and removal

These are methods that change the number of elements in the array by insertion
or erasure. **/
/*@{*/

/** This method increases the size of the Array by one element at the end and 
initializes that element by copy constructing it from the given value. If 
capacity() > size(), that's all that will happen. If capacity()==size(), there
is no room for another element so we'll allocate more space and move all the 
elements there. A reference to the just-inserted element can be obtained using
the back() method after the call to push_back().
@param[in]      value
    An object of type T from which the new element is copy-constructed.

@remarks
  - If you are appending a default-constructed object of type T, consider using
    the alternate non-standard but safe push_back() method rather than 
    push_back(T()). The non-standard method default-constructs the new element 
    internally. That avoids a call to the copy constructor which can be 
    expensive for some objects, and nonexistent for others.
  - If you are constructing the source object with a non-default constructor,
    and the object is expensive or impossible to default-construct and/or 
    copy-construct, consider using the non-standard and dangerous method 
    raw_push_back() which enables you to construct the new element in place. 

@par Complexity:
    Constant time if no reallocation is required; otherwise the current 
    contents of the array must be copied to new space, costing one call to T's
    copy constructor and destructor (if any) for each element currently in the
    array. Either way there is also one call to T's copy constructor to 
    construct the new element from the supplied value. **/
void push_back(const T& value) {
    if (pallocated() == psize())
        growAtEnd(1,"Array_<T>::push_back(value)");
    copyConstruct(end(), value);
    incrSize();
}

/** This is a non-standard version of push_back() that increases the size of the
array by one default-constructed element at the end. This avoids having to 
default-construct the argument to the standard push_back(value) method which 
then has to copy-construct it into the array. By carefully avoiding 
reallocation and using this form of push_back() you can use the Array_<T> class
to hold objects of type T even if T has no copy constructor, which is 
prohibited by the standard std::vector<T> definition. 

@par Complexity:
    Same as the standard push_back(value) method except without the final
    call to T's copy constructor.
@see push_back(value) 
**/
void push_back() {
    if (pallocated() == psize())
        growAtEnd(1,"Array_<T>::push_back()");
    defaultConstruct(end());
    incrSize();
}

/** This dangerous method increases the Array's size by one element at the end 
but doesn't perform any construction so the memory is filled with garbage. You 
must immediately construct into this space, using code like:
@code       
    new(a.raw_push_back()) MyConstructor(...args...);       
@endcode
This is a substantial performance improvement when the element type is something
complicated since the constructor is called once and not copied; it can also be
used for objects that have neither default nor copy constructors.
@return 
    An iterator (pointer) pointing at the unconstructed element. 
@par Complexity:
    Same as ordinary push_back().
@see push_back(value), push_back() 
**/
T* raw_push_back() {
    if (pallocated() == psize())
        growAtEnd(1,"Array_<T>::raw_push_back()");
    T* const p = end();
    incrSize();
    return p;
}

/** Remove the last element from this array, which must not be empty. The 
element is destructed, not returned. The array's size() is reduced by one. **/
void pop_back() {
    SimTK_ERRCHK(!empty(), "Array_<T>::pop_back()", "Array was empty.");
    destruct(&back());
    decrSize();
}

/** Erase elements in range [first,last1), packing in any later elements into 
the newly-available space and reducing the array's size by the number of 
elements erased. Capacity is unchanged. If the range is empty nothing happens.

@pre begin() <= \a first <= \a last1 <= end()
@param      first
    Points to the first element that will be erased.
@param      last1
    Points one element past the last element to be erased.
@return
    An iterator pointing to the new location of the element immediately
    following the erased ones, or end() if there are none. Either way, this is 
    the same memory address as the passed-in \a first argument since there can
    be no reallocation here.
@par Complexity:
    Calls T's destructor once for each erased element and calls T's copy 
    constructor and destructor once for each element that has to be moved. **/
T* erase(T* first, const T* last1) {
    SimTK_ERRCHK(begin() <= first && first <= last1 && last1 <= end(),
    "Array<T>::erase(first,last1)", "Pointers out of range or out of order.");

    const size_type nErased = size_type(last1-first);
    SimTK_ERRCHK(isOwner() || nErased==0, "Array<T>::erase(first,last1)",
        "No elements can be erased from a non-owner array.");

    if (nErased) {
        destruct(first, last1); // Destruct the elements we're erasing.
        moveElementsDown(first+nErased, nErased); // Compress followers into the gap.
        setSize(size()-nErased);
    }
    return first;
}

/** Erase just one element, moving all subsequent elements down one slot and
reducing the array's size by one. This is equivalent to erase(p,p+1) but faster;
that means \a p cannot be end() because end()+1 is not defined. Capacity is 
unchanged.

@note If you don't mind the elements being reordered, you can erase an element
in constant time using the non-standard extension eraseFast().

@param      p
    Points to the element that will be erased; \a p cannot be end().
@return
    A pointer to the element that replaced the one at \a p, or end() if \a p 
    was the last element. Either way, this is the same memory address as the 
    erased element had since there can be no reallocation here.
@pre begin() <= \a p < end()
@par Complexity:
    Calls T's destructor once for the erased element and calls T's copy 
    constructor and destructor once for each element that has to be moved. 
@see eraseFast() **/
T* erase(T* p) {
    SimTK_ERRCHK(begin() <= p && p < end(),
        "Array<T>::erase(p)", "Pointer must point to a valid element.");
    SimTK_ERRCHK(isOwner(), "Array<T>::erase(p)",
        "No elements can be erased from a non-owner array.");

    destruct(p);              // Destruct the element we're erasing.
    moveElementsDown(p+1, 1); // Compress followers into the gap.
    decrSize();
    return p;
}


/** Be careful with this non-standard extension; it erases one element and 
then moves the last one in its place which changes the element order
from what it was before (unlike the standard erase() method). This avoids
having to compress the elements so this runs in constant time:
the element is destructed; then if it wasn't the last element the
copy constructor is used to copy the last element into the vacated
space, and the destructor is called to clear the last element. The
size is reduced by 1 but the capacity does not change. 

@param      p
    Points to the element that will be erased; \a p cannot be end().
@return
    A pointer to the element that replaced the one at \a p, or end() if \a p 
    was the last element. Either way, this is the same memory address as the 
    erased element had since there can be no reallocation here.
@pre begin() <= \a p < end()
@par Complexity:
    Calls T's destructor once for the erased element and calls T's copy 
    constructor and destructor once for each element that has to be moved.
@see erase() **/
T* eraseFast(T* p) {
    SimTK_ERRCHK(begin() <= p && p < end(),
        "Array<T>::eraseFast(p)", "Pointer must point to a valid element.");
    SimTK_ERRCHK(isOwner(), "Array<T>::eraseFast(p)",
        "No elements can be erased from a non-owner array.");

    destruct(p);
    if (p+1 != end()) 
        moveOneElement(p, &back());
    decrSize();
    return p;
}

/** Erase all the elements currently in this array without changing the 
capacity; equivalent to erase(begin(),end()) but a little faster. Size is 
zero after this call. T's destructor is called exactly once for each element 
in the array.

@par Complexity:
    O(n) if T has a destructor; constant time otherwise. **/
void clear() {
    SimTK_ERRCHK(isOwner() || empty(), "Array_<T>::clear()", 
        "clear() is not allowed for a non-owner array.");
    destruct(begin(), end());
    setSize(0);
}


/** Insert \a n copies of a given value at a particular location within this 
array, moving all following elements up by \a n positions.

@param[in]      p        
    Where to insert the new elements. This must be an iterator (pointer) that 
    is valid for this array, that is, begin() <= \a p <= end().
@param[in]      n
    How many copies of the given \a value to insert. Nothing happens if 
    \a n is zero.
@param[in]      value    
    A value of the element type that is copied into the newly-created elements
    using T's copy constructor.
@return         
    A pointer to the first of the newly-created elements in the array. This 
    will be different from \a p if reallocation occurred, otherwise it is the
    same as \a p was on entry.

@pre begin() <= \a p <= end()
@par Complexity:
    If size() + \a n > capacity() then the array must be reallocated, resulting
    in size() copy constructor/destructor call pairs to move the old data to 
    the new location. Otherwise, the m=(end()-\a p) elements above the insertion 
    point must be moved up \a n positions resulting in m copy/destruct pairs.
    Then there are n additional copy constructor calls to construct the new 
    elements from the given value. 
**/
T* insert(T* p, size_type n, const T& value) {
    T* const gap = insertGapAt(p, n, "Array<T>::insert(p,n,value)");
    // Copy construct into the inserted elements and note the size change.
    fillConstruct(gap, gap+n, value);
    setSize(size()+n);
    return gap;
}

/** Insert a new element at a given location within this array, initializing 
it to a copy of a given value and moving all following elements up one 
position. This is identical to insert(\a p,1,\a value) but slightly faster; see
that method for full documentation. **/
T* insert(T* p, const T& value)  {
    T* const gap = insertGapAt(p, 1, "Array<T>::insert(p,value)");
    // Copy construct into the inserted element and note the size change.
    copyConstruct(gap, value);
    incrSize();
    return gap;
}

/** Insert elements in a range [first,last1) into this array at a given position
\a p, moving all following elements up by n=(last1-first) positions. This 
variant of insert() takes iterators which are ordinary pointers, although the
source elements do not have to be of type T as long as there is a constructor
T(T2) that works.

@param[in]      p        
    Where to insert the new elements. This must be an iterator (pointer) that 
    is valid for this array, that is, begin() <= \a p <= end().
@param[in]      first
    This is a pointer to the first element of the source to be copied.
@param[in]      last1    
    This points one element past the last element of the source to be copied.
@return         
    A pointer to the first of the newly-created elements in the array. This 
    will be different from \a p if reallocation occurred, otherwise it is the
    same as \a p was on entry.

@pre begin() <= \a p <= end()
@pre first <= last1
@pre The range [first,last1) does not include any of the current contents 
of this array.
@par Complexity:
    If capacity() < size()+n then the array must be reallocated, resulting in 
    size() calls to T's copy constructor and destructor (if any)to move the old
    data to the new location. Otherwise, the m=(end()-\a p) elements above the 
    insertion point must be moved up n positions resulting in m copy/destruct 
    pairs. Then there are n additional copy constructor calls to construct the 
    new elements from the given value. **/
template <class T2>
T* insert(T* p, const T2* first, const T2* last1) {
    const char* methodName = "Array_<T>::insert(T* p, T2* first, T2* last1)";
    SimTK_ERRCHK((first&&last1) || (first==last1), methodName, 
        "One of first or last1 was null; either both or neither must be null.");
    SimTK_ERRCHK(!this->overlapsWithData(first,last1), methodName,
        "Source range can't overlap with the current array contents.");
    // Pointers are random access iterators.
    return insertIteratorDispatch(p, first, last1,
                                  std::random_access_iterator_tag(),
                                  methodName);
}

/** Insert elements in a range [first,last1) where the range is given by
non-pointer iterators. **/
template <class Iter>
T* insert(T* p, const Iter& first, const Iter& last1) {
    return insertDispatch(p, first, last1,
                          typename IsIntegralType<Iter>::Result(),
                          "Array_<T>::insert(T* p, Iter first, Iter last1)");
}
/*@}    End of insertion and erase. **/



//------------------------------------------------------------------------------
                                  private:
//------------------------------------------------------------------------------


// This method is used when we have already decided we need to make room for 
// some new elements by reallocation, by creating a gap somewhere within the
// existing data. We'll issue an error message if this would violate the  
// max_size restriction (we can afford to do that even in a Release build since
// we don't expect to grow very often). Otherwise we'll allocate some more space
// and copy construct the existing elements into the new space, leaving a gap 
// where indicated. Note that this method does not change the current size but
// does change the capacity.
//
// The gapPos must point within the existing data with null OK if the array
// itself is null, and end() being OK any time although you should use the
// more specialized growAtEnd() method if you know that's what's happening.
//
// Don't call this with a gapSz of zero.
T* growWithGap(T* gapPos, size_type gapSz, const char* methodName) {
    assert(gapSz > 0); // <= 0 is a bug, not a user error

    // Note that gapPos may be null if begin() and end() are also.
    SimTK_ERRCHK(begin() <= gapPos && gapPos <= end(), methodName, 
        "Given insertion point is not valid for this array.");

    // Get some new space of a reasonable size.
    setAllocated(calcNewCapacityForGrowthBy(gapSz, methodName));
    T* newData   = allocN(allocated());

    // How many elements will be before the gap?
    const size_type nBefore = (size_type)(gapPos-begin());

    // Locate the gap in the new space allocation.
    T* newGap    = newData + nBefore;
    T* newGapEnd = newGap  + gapSz; // one past the last element in the gap

    // Copy elements before insertion point; destruct source as we go.
    copyConstructThenDestructSource(newData,   newGap,        data());
    // Copy/destruct elements at and after insertion pt; leave gapSz gap.
    copyConstructThenDestructSource(newGapEnd, newData+size(), gapPos);

    // Throw away the old data and switch to the new.
    freeN(data()); setData(newData);
    return newGap;
}

// Same as growWithGap(end(), n, methodName); see above.
void growAtEnd(size_type n, const char* methodName) {
    assert(n > 0); // <= 0 is a bug, not a user error
    // Get some new space of a reasonable size.
    setAllocated(calcNewCapacityForGrowthBy(n, methodName));
    T* newData   = allocN(allocated());
    // Copy all the elements; destruct source as we go.
    copyConstructThenDestructSource(newData, newData+size(), data());
    // Throw away the old data and switch to the new.
    freeN(data()); setData(newData);
}

// This method determines how much we should increase the array's capacity
// when asked to insert n elements, due to an insertion or push_back. We will
// generally allocate more new space than requested, in anticipation of
// further insertions. This has to be based on the current size so that only
// log(n) reallocations are performed to insert n elements one at a time. Our
// policy here is to at least double the capacity unless that would exceed 
// max_size(). There is also a minimum amount of allocation we'll do if the 
// current size is zero or very small. 
size_type calcNewCapacityForGrowthBy(size_type n, const char* methodName) const {
    SimTK_ERRCHK3_ALWAYS(isGrowthOK(n), methodName,
        "Can't grow this Array by %llu element(s) because it would"
        " then exceed the max_size of %llu set by its index type %s.",
        (unsigned long long)n, ullMaxSize(), indexName());

    // At this point we know that capacity()+n <= max_size().
    const size_type mustHave = capacity() + n;

    // Be careful not to overflow size_type as you could if you 
    // double capacity() rather than halving max_size().
    const size_type wantToHave = capacity() <= max_size()/2 
                                    ? 2*capacity() 
                                    : max_size();

    const size_type newCapacity = std::max(std::max(mustHave,wantToHave),
                                           minAlloc());
    return newCapacity;
}

// Insert an unconstructed gap of size n beginning at position p. The return
// value is a pointer to the first element in the gap. It will be p if no
// reallocation occurred, otherwise it will be pointing into the new data.
// On return size() will be unchanged although allocated() may be bigger.
T* insertGapAt(T* p, size_type n, const char* methodName) {
    // Note that p may be null if begin() and end() are also.
    SimTK_ERRCHK(begin() <= p && p <= end(), methodName, 
        "Given insertion point is not valid for this Array.");

    if (n==0) return p; // nothing to do

    SimTK_ERRCHK_ALWAYS(isOwner(), methodName,
        "No elements can be inserted into a non-owner array.");

    // Determine the number of elements before the insertion point and
    // the number at or afterwards (those must be moved up by one slot).
    const size_type before = (size_type)(p-begin()), after = (size_type)(end()-p);

    // Grow the container if necessary. Note that if we have to grow we
    // can create the gap at the same time we copy the old elements over
    // to the new space.
    if (capacity() >= size()+n) {
        moveElementsUp(p, n); // leave a gap at p
    } else { // need to grow
        setAllocated(calcNewCapacityForGrowthBy(n, methodName));
        T* newdata = allocN(allocated());
        // Copy the elements before the insertion point, and destroy source.
        copyConstructThenDestructSource(newdata, newdata+before, data());
        // Copy the elements at and after the insertion point, leaving a gap
        // of n elements.
        copyConstructThenDestructSource(newdata+before+n,
                                        newdata+before+n+after,
                                        p); // i.e., pData+before
        p = newdata + before; // points into newdata now
        freeN(data());
        setData(newdata);
    }

    return p;
}

//------------------------------------------------------------------------------
//                           CTOR DISPATCH
// This is the constructor implementation for when the class that matches
// the alleged InputIterator type turns out to be one of the integral types
// in which case this should be the ctor(n, initValue) constructor.
template <class IntegralType> void
ctorDispatch(IntegralType n, IntegralType v, TrueType isIntegralType) {
    new(this) Array_(size_type(n), value_type(v));
}

// This is the constructor implementation for when the class that matches
// the alleged InputIterator type is NOT an integral type and may very well
// be an iterator. In that case we split into iterators for which we can
// determine the number of elements in advance (forward, bidirectional,
// random access) and input iterators, for which we can't. Note: iterator
// types are arranged hierarchically random->bi->forward->input with each
// deriving from the one on its right, so the forward iterator tag also
// matches bi and random.
template <class InputIterator> void
ctorDispatch(const InputIterator& first, const InputIterator& last1, 
             FalseType isIntegralType) 
{   ctorIteratorDispatch(first, last1, 
        typename std::iterator_traits<InputIterator>::iterator_category()); }

// This is the slow generic ctor implementation for any iterator that can't
// make it up to forward_iterator capability (that is, an input_iterator).
// The issue here is that we can't advance the iterator to count the number
// of elements before allocating because input_iterators are consumed when
// reference so we can't go back to look. That means we may have to reallocate
// memory log n times as we "push back" these elements onto the array.
template <class InputIterator> void
ctorIteratorDispatch(const InputIterator& first, const InputIterator& last1, 
                     std::input_iterator_tag) 
{
    InputIterator src = first;
    while (src != last1) {
        // We can afford to check this always since we are probably doing I/O.
        // Throwing an exception in a constructor is tricky, though -- this
        // won't go through the Array_ destructor although it will call the
        // Base (ArrayView_) destructor. Since we have already allocated
        // some space, we must call deallocate() manually.
        if (size() == max_size()) {
            deallocate();
            SimTK_ERRCHK2_ALWAYS(!"too many elements",
                "Array_::ctor(InputIterator first, InputIterator last1)",
                "There were still source elements available when the array"
                " reached its maximum size of %llu as determined by its index"
                " type %s.", ullMaxSize(), indexName());
        }
        push_back(*src++);
    }
}

// This is the faster constructor implementation for iterator types for which
// we can calculate the number of elements in advance. This will be optimal
// for a random access iterator since we can count in constant time, but for
// forward or bidirectional we'll have to advance n times to count and then
// go back again to do the copy constructions.
template <class ForwardIterator> void
ctorIteratorDispatch(const ForwardIterator& first, const ForwardIterator& last1, 
                     std::forward_iterator_tag) 
{
    typedef typename std::iterator_traits<ForwardIterator>::difference_type
        difference_type;
    // iterDistance() is constant time for random access iterators, but 
    // O(last1-first) for forward and bidirectional since it has to increment 
    // to count how far apart they are.
    const difference_type nInput = this->iterDistance(first,last1);

    SimTK_ERRCHK(nInput >= 0, 
        "Array_(ForwardIterator first, ForwardIterator last1)", 
        "Iterators were out of order.");

    SimTK_ERRCHK3(this->isSizeOK(nInput), 
        "Array_(ForwardIterator first, ForwardIterator last1)",
        "Source has %llu elements but this array is limited to %llu"
        " elements by its index type %s.",
        this->ull(nInput), ullMaxSize(), indexName());

    const size_type n = size_type(nInput);
    setSize(n);
    allocateNoConstruct(n);
    copyConstruct(data(), data()+n, first);
}

//------------------------------------------------------------------------------
//                           INSERT DISPATCH
// This is the insert() implementation for when the class that matches
// the alleged InputIterator type turns out to be one of the integral types
// in which case this should be the insert(p, n, initValue) constructor.
template <class IntegralType> 
T* insertDispatch(T* p, IntegralType n, IntegralType v, 
                  TrueType isIntegralType, const char*) 
{   return insert(p, size_type(n), value_type(v)); }

// This is the insert() implementation for when the class that matches
// the alleged InputIterator type is NOT an integral type and may very well
// be an iterator. See ctorDispatch() above for more information.
template <class InputIterator> 
T* insertDispatch(T* p, const InputIterator& first, const InputIterator& last1, 
                  FalseType isIntegralType, const char* methodName) 
{   return insertIteratorDispatch(p, first, last1, 
        typename std::iterator_traits<InputIterator>::iterator_category(),
        methodName); }

// This is the slow generic insert implementation for any iterator that can't
// make it up to forward_iterator capability (that is, an input_iterator).
// See ctorIteratorDispatch() above for more information.
template <class InputIterator> 
T* insertIteratorDispatch(T* p, InputIterator first, InputIterator last1, 
                          std::input_iterator_tag, const char* methodName) 
{
    size_type nInserted = 0;
    while (first != last1) {
        // We can afford to check this always since we are probably doing I/O.
        SimTK_ERRCHK2_ALWAYS(size() < max_size(), methodName,
            "There were still source elements available when the array"
            " reached its maximum size of %llu as determined by its index"
            " type %s.", ullMaxSize(), indexName());
        p = insert(p, *first++);  // p may now point to reallocated memory
        ++p; ++nInserted;
    }
    // p now points just after the last inserted element; subtract the
    // number inserted to get a pointer to the first inserted element.
    return p-nInserted;
}

// This is the faster constructor implementation for iterator types for which
// we can calculate the number of elements in advance. This will be optimal
// for a random access iterator since we can count in constant time, but for
// forward or bidirectional we'll have to advance n times to count and then
// go back again to do the copy constructions.
template <class ForwardIterator>
T* insertIteratorDispatch(T* p, const ForwardIterator& first, 
                                const ForwardIterator& last1,
                                std::forward_iterator_tag,
                                const char* methodName) 
{
    typedef typename std::iterator_traits<ForwardIterator>::difference_type
        difference_type;
    // iterDistance() is constant time for random access iterators, but 
    // O(last1-first) for forward and bidirectional since it has to increment 
    // to count how far apart they are.
    const difference_type nInput = this->iterDistance(first,last1);

    SimTK_ERRCHK(nInput >= 0, methodName, "Iterators were out of order.");

    SimTK_ERRCHK3(isGrowthOK(nInput), methodName,
        "Source has %llu elements which would make this array exceed the %llu"
        " elements allowed by its index type %s.",
        this->ull(nInput), ullMaxSize(), indexName());

    const size_type n = size_type(nInput);
    p = insertGapAt(p, n, methodName);
    copyConstruct(p, p+n, first);
    setSize(size()+n);
    return p;
}

//------------------------------------------------------------------------------
//                           ASSIGN DISPATCH
// This is the assign() implementation for when the class that matches
// the alleged InputIterator type turns out to be one of the integral types
// in which case this should be the assign(n, initValue) constructor.
template <class IntegralType>
void assignDispatch(IntegralType n, IntegralType v, TrueType isIntegralType,
                    const char* methodName) 
{   assign(size_type(n), value_type(v)); }

// This is the assign() implementation for when the class that matches
// the alleged InputIterator type is NOT an integral type and may very well
// be an iterator. See ctorDispatch() above for more information.
template <class InputIterator> 
void assignDispatch(const InputIterator& first, const InputIterator& last1, 
                    FalseType isIntegralType, const char* methodName) 
{   assignIteratorDispatch(first, last1, 
        typename std::iterator_traits<InputIterator>::iterator_category(),
        methodName); }

// This is the slow generic implementation for a plain input_iterator
// (i.e., not a forward, bidirectional, or random access iterator). These
// have the unfortunate property that we can't count the elements in advance.
template <class InputIterator>
void assignIteratorDispatch(const InputIterator& first, 
                            const InputIterator& last1, 
                            std::input_iterator_tag, 
                            const char* methodName) 
{
    // TODO: should probably allow this and just blow up when the size()+1st
    // element is seen.
    SimTK_ERRCHK_ALWAYS(isOwner(), methodName,
        "Assignment to a non-owner array can only be done from a source"
        " designated with forward iterators or pointers because we"
        " must be able to verify that the source and destination sizes"
        " are the same.");

    clear(); // TODO: change space allocation here?
    InputIterator src = first;
    while (src != last1) {
        // We can afford to check this always since we are probably doing I/O.
        SimTK_ERRCHK2_ALWAYS(size() < max_size(), methodName,
            "There were still source elements available when the array"
            " reached its maximum size of %llu as determined by its index"
            " type %s.", ullMaxSize(), indexName());

        push_back(*src++);
    }
}

// This is the faster implementation that works for forward, bidirectional,
// and random access iterators including ordinary pointers. We can check here that the 
// iterators are in the right order, and that the source is not too big to
// fit in this array. Null pointer checks should be done prior to calling,
// however, since iterators in general aren't pointers.
template <class ForwardIterator>
void assignIteratorDispatch(const ForwardIterator& first, 
                            const ForwardIterator& last1,
                            std::forward_iterator_tag, 
                            const char* methodName) 
{
    typedef typename std::iterator_traits<ForwardIterator>::difference_type
        IterDiffType;
    // iterDistance() is constant time for random access iterators, but 
    // O(last1-first) for forward and bidirectional since it has to increment 
    // to count how far apart they are.
    const IterDiffType nInput = this->iterDistance(first,last1);

    SimTK_ERRCHK(nInput >= 0, methodName, "Iterators were out of order.");

    SimTK_ERRCHK3(this->isSizeOK(nInput), methodName,
        "Source has %llu elements but this Array is limited to %llu"
        " elements by its index type %s.",
        this->ull(nInput), ullMaxSize(), indexName());

    const size_type n = size_type(nInput);
    if (isOwner()) {
        // This is an owner Array; assignment is considered deallocation
        // followed by copy construction.

        clear(); // all elements destructed; allocation unchanged
        reallocateIfAdvisable(n); // change allocation if too small or too big
        copyConstruct(data(), cdata()+n, first);
        setSize(n);
    } else {
        // This is a non-owner Array. Assignment can occur only if the
        // source is the same size as the array, and the semantics are of
        // repeated assignment using T::operator=() not destruction followed
        // by copy construction.
        SimTK_ERRCHK2(n == size(), methodName,
            "Source has %llu elements which does not match the size %llu"
            " of the non-owner array it is being assigned into.",
            this->ull(n), ullSize());

        T* p = begin();
        ForwardIterator src = first;
        while (src != last1)
            *p++ = *src++; // call T's assignment operator
    }
}

// We are going to put a total of n elements into the Array (probably
// because of an assignment or resize) and we want the space allocation
// to be reasonable. That means of course that the allocation must be 
// *at least* n, but we also don't want it to be too big. Our policy
// here is that if it is currently less than twice what we need we
// won't reallocate, otherwise we'll shrink the space. When changing
// the size to zero or something very small we'll treat the Array as
// though its current size is minAlloc, meaning we won't reallocate
// if the existing space is less than twice minAlloc.
// nAllocated will be set appropriately; size() is not touched here.
// No constructors or destructors are called.
void reallocateIfAdvisable(size_type n) {
    if (allocated() < n || allocated()/2 > std::max(minAlloc(), n)) 
        reallocateNoDestructOrConstruct(n);
}


void allocateNoConstruct(size_type n) 
{   setData(allocN(n)); setAllocated(n); } // size() left unchanged
void deallocateNoDestruct() 
{   freeN(data()); setData(0); setAllocated(0); } // size() left unchanged
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

// Free memory without calling T's destructor. Nothing happens if passed
// a null pointer.
static void freeN(T* p) {
    delete[] reinterpret_cast<char*>(p);
}

// default construct one element
static void defaultConstruct(T* p) {new(p) T();}
// default construct range [b,e)
static void defaultConstruct(T* b, const T* e) 
{   while (b!=e) new(b++) T(); }

// copy construct range [b,e) with repeats of a given value
static void fillConstruct(T* b, const T* e, const T& v)
{   while(b!=e) new(b++) T(v); }

// copy construct one element from a given value
static void copyConstruct(T* p, const T& v) {new(p) T(v);}
// copy construct range [b,e) from sequence of source values
static void copyConstruct(T* b, const T* e, T* src)
{   while(b!=e) new(b++) T(*src++); }
// Templatized copy construct will work if the source elements are
// assignment compatible with the destination elements.
template <class InputIterator>
static void copyConstruct(T* b, const T* e, InputIterator src)
{   while(b!=e) new(b++) T(*src++); }

// Copy construct range [b,e] from sequence of source values and
// destruct the source after it is copied. It's better to alternate
// copying and destructing than to do this in two passes since we
// will already have touched the memory.
static void copyConstructThenDestructSource(T* b, const T* e, T* src)
{   while(b!=e) {new(b++) T(*src); src++->~T();} }

// We have an element at from that we would like to move into the currently-
// unconstructed slot at to. Both from and to are expected to be pointing to
// elements within the currently allocated space. From's slot will be left
// unconstructed.
void moveOneElement(T* to, T* from) {
    assert(data() <= to   && to   < data()+allocated());
    assert(data() <= from && from < data()+allocated());
    copyConstruct(to, *from); 
    destruct(from);
}


// Move elements from p to end() down by n>0 places to fill an unconstructed 
// gap beginning at p-n. Any leftover space at the end will be unconstructed.
void moveElementsDown(T* p, size_type n) {
    assert(n > 0);
    for (; p != end(); ++p)
        moveOneElement(p-n,p);
}

// Move elements from p to end() up by n>0 places to make an unconstructed gap
// at [p,p+n). Note that this has to be done backwards so that we don't
// write on any elements until after they've been copied.
void moveElementsUp(T* p, size_type n) {
    assert(n > 0);
    T* src = end(); // points one past source element (begin()-1 not allowed)
    while (src != p) {
        --src; // now points to source
        moveOneElement(src+n, src);;
    }
}

// destruct one element
static void destruct(T* p) {p->~T();}
// destruct range [b,e)
static void destruct(T* b, const T* e)
{   while(b!=e) b++->~T(); }

// Check that growing this array by n elements wouldn't cause it to exceed
// its allowable maximum size.
template <class S>
bool isGrowthOK(S n) const
{   return this->isSizeOK(ullCapacity() + this->ull(n)); }

// The following private methods are protected methods in the ArrayView base 
// class, so they should not need repeating here. Howevr, we explicitly 
// forward to the Base methods to avoid gcc errors. The gcc complaint
// is due to their not depending on any template parameters; the "this->"
// apparently fixes that problem.

// These provide direct access to the data members.
packed_size_type psize()      const {return this->CBase::psize();}
packed_size_type pallocated() const {return this->CBase::pallocated();}

void setData(const T* p)       {this->CBase::setData(p);}
void setSize(size_type n)      {this->CBase::setSize(n);}
void incrSize()                {this->CBase::incrSize();}
void decrSize()                {this->CBase::decrSize();}
void setAllocated(size_type n) {this->CBase::setAllocated(n);}
// This just cast sizes to unsigned long long so that we can do comparisons
// without getting warnings.
unsigned long long ullSize()     const 
{   return this->CBase::ullSize(); }
unsigned long long ullCapacity() const 
{   return this->CBase::ullCapacity(); }
unsigned long long ullMaxSize()  const 
{   return this->CBase::ullMaxSize(); }
// This is the index type name and is handy for error messages to explain
// why some size was too big.
const char* indexName() const   {return this->CBase::indexName();}
};


// This "private" static method is used to implement ArrayView's 
// fillArrayViewFromStream() and Array's readArrayFromStream() namespace-scope
// static methods, which are in turn used to implement ArrayView's and 
// Array's stream extraction operators ">>". This method has to be in the 
// header file so that we don't need to pass streams through the API, but it 
// is not intended for use by users and has no Doxygen presence, unlike 
// fillArrayFromStream() and readArrayFromStream() and (more commonly)
// the extraction operators.
template <class T, class X> static inline 
std::istream& readArrayFromStreamHelper
   (std::istream& in, bool isFixedSize, Array_<T,X>& out)
{
    // If already failed, bad, or eof, set failed bit and return without 
    // touching the Array.
    if (!in.good()) {in.setstate(std::ios::failbit); return in;}

    // If the passed-in Array isn't an owner, then we have to treat it as
    // a fixed size ArrayView regardless of the setting of the isFixedSize
    // argument.
    if (!out.isOwner())
        isFixedSize = true; // might be overriding the argument here

    // numRequired will be ignored unless isFixedSize==true.
    const typename Array_<T,X>::size_type 
        numRequired = isFixedSize ? out.size() : 0;

    if (!isFixedSize)
        out.clear(); // We're going to replace the entire contents of the Array.

    // Skip initial whitespace. If that results in eof this may be a successful
    // read of a 0-length, unbracketed Array. That is OK for either a
    // variable-length Array or a fixed-length ArrayView of length zero.
    std::ws(in); if (in.fail()) return in;
    if (in.eof()) {
        if (isFixedSize && numRequired != 0)
            in.setstate(std::ios_base::failbit); // zero elements not OK
        return in;
    }
    
    // Here the stream is good and the next character is non-white.
    assert(in.good());

    // Use this for raw i/o (peeks and gets).
    typename       std::iostream::int_type ch;
    const typename std::iostream::int_type EOFch = 
        std::iostream::traits_type::eof();

    // Now see if the sequence is bare or surrounded by (), [], or {}.
    bool lookForCloser = true;
    char openBracket, closeBracket;
    ch = in.peek(); if (in.fail()) return in;
    assert(ch != EOFch); // we already checked above

    openBracket = (char)ch;
    if      (openBracket=='(') {in.get(); closeBracket = ')';}
    else if (openBracket=='[') {in.get(); closeBracket = ']';}
    else if (openBracket=='{') {in.get(); closeBracket = '}';}
    else lookForCloser = false;

    // If lookForCloser is true, then closeBracket contains the terminating
    // delimiter, otherwise we're not going to quit until eof.

    // Eat whitespace after the opening bracket to see what's next.
    if (in.good()) std::ws(in);

    // If we're at eof now it must be because the open bracket was the
    // last non-white character in the stream, which is an error.
    if (!in.good()) {
        if (in.eof()) {
            assert(lookForCloser); // or we haven't read anything that could eof
            in.setstate(std::ios::failbit);
        }
        return in;
    }

    // istream is good and next character is non-white; ready to read first
    // value or terminator.

    // We need to figure out whether the elements are space- or comma-
    // separated and then insist on consistency.
    bool commaOK = true, commaRequired = false;
    bool terminatorSeen = false;
    X nextIndex(0);
    while (true) {
        char c;

        // Here at the top of this loop, we have already successfully read 
        // n=nextIndex values of type T. For fixed-size reads, it might be
        // the case that n==numRequired already, but we still may need to
        // look for a closing bracket before we can declare victory.
        // The stream is good() (not at eof) but it might be the case that 
        // there is nothing but white space left; we don't know yet because
        // if we have satisfied the fixed-size count and are not expecting
        // a terminator then we should quit without absorbing the trailing
        // white space.
        assert(in.good());

        // Look for closing bracket before trying to read value.
        if (lookForCloser) {
            // Eat white space to find the closing bracket.
            std::ws(in); if (!in.good()) break; // eof?
            ch = in.peek(); assert(ch != EOFch);
            if (!in.good()) break;
            c = (char)ch;
            if (c == closeBracket) {   
                in.get(); // absorb the closing bracket
                terminatorSeen = true; 
                break; 
            }
            // next char not a closing bracket; fall through
        }

        // We didn't look or didn't find a closing bracket. The istream is good 
        // but we might be looking at white space.

        // If we already got all the elements we want, break for final checks.
        if (isFixedSize && (nextIndex == numRequired))
            break; // that's a full count.

        // Look for comma before value, except the first time.
        if (commaOK && nextIndex != 0) {
            // Eat white space to find the comma.
            std::ws(in); if (!in.good()) break; // eof?
            ch = in.peek(); assert(ch != EOFch);
            if (!in.good()) break;
            c = (char)ch;
            if (c == ',') {
                in.get(); // absorb comma
                commaRequired = true; // all commas from now on
            } else { // next char not a comma
                if (commaRequired) // bad, e.g.: v1, v2, v3 v4 
                {   in.setstate(std::ios::failbit); break; }
                else commaOK = false; // saw: v1 v2 (no commas now)
            }
            if (!in.good()) break; // might be eof
        }

        // No closing bracket yet; don't have enough elements; skipped comma 
        // if any; istream is good; might be looking at white space.
        assert(in.good());

        // Now read in an element of type T.
        // The extractor T::operator>>() will ignore leading white space.
        if (!isFixedSize)
            out.push_back(); // grow by one (default consructed)
        in >> out[nextIndex]; if (in.fail()) break;
        ++nextIndex;

        if (!in.good()) break; // might be eof
    }

    // We will get here under a number of circumstances:
    //  - the fail bit is set in the istream, or
    //  - we reached eof
    //  - we saw a closing brace
    //  - we got all the elements we wanted (for a fixed-size read)
    // Note that it is possible that we consumed everything except some
    // trailing white space (meaning we're not technically at eof), but
    // for consistency with built-in operator>>()'s we won't try to absorb
    // that trailing white space.

    if (!in.fail()) {
        if (lookForCloser && !terminatorSeen)
            in.setstate(std::ios::failbit); // missing terminator

        if (isFixedSize && nextIndex != numRequired)
            in.setstate(std::ios::failbit); // wrong number of values
    }

    return in;
}



//------------------------------------------------------------------------------
//                          RELATED GLOBAL OPERATORS
//------------------------------------------------------------------------------
// These are logically part of the Array_<T,X> class but are not actually 
// class members; that is, they are in the SimTK namespace.

// Some of the serialization methods could have been member functions but 
// then an attempt to explicitly instantiate the whole Array_<T> class for
// some type T would fail if T did not support the requisite I/O operations
// even if those operations were never used. This came up when Chris Bruns was
// trying to wrap Array objects for Python, which requires explicit 
// instantiation.

/**@name             Array_<T> serialization and I/O
These methods are at namespace scope but are logically part of the Array
classes. These deal with reading and writing Arrays from and to streams,
which places an additional requirement on the element type T: the element 
must support the same operation you are trying to do on the Array as a 
whole. **/
/*@{*/

/** Specialize writeUnformatted() for Array_<E,X> to delegate to element type
E, with spaces separating the elements. 
@relates SimTK::Array_ **/
template <class T, class X> inline void
writeUnformatted(std::ostream& o, const Array_<T,X>& v) {
    for (X i(0); i < v.size(); ++i) {
        if (i != 0) o << " ";
        writeUnformatted(o, v[i]);
    }
}


/** Specialize writeFormatted() for Array_<E,X> to delegate to element type
E, with surrounding parentheses and commas separating the elements. 
@relates SimTK::Array_ **/
template <class T, class X> inline void
writeFormatted(std::ostream& o, const Array_<T,X>& v) {
    o << '(';
    for (X i(0); i < v.size(); ++i) {
        if (i != 0) o << ',';
        writeFormatted(o, v[i]);
    }
    o << ')';
}

/** Output a human readable representation of an array to an std::ostream
(like std::cout). The format is ( \e elements ) where \e elements is a 
comma-separated list of the Array's contents output by invoking the "<<" 
operator on the elements. This function will not compile if the element type 
does not support the "<<" operator. No newline is issued before
or after the output. @relates Array_ **/
template <class T, class X> inline 
std::ostream&
operator<<(std::ostream& o, 
           const ArrayViewConst_<T,X>& a) 
{
    o << '(';
    if (!a.empty()) {
        o << a.front();
        for (const T* p = a.begin()+1; p != a.end(); ++p)
            o << ',' << *p;
    }
    return o << ')';
} 


/** Specialization of readUnformatted() for variable-length Array_<T,X>; 
continues reading whitespace-separated tokens until error or eof. 
@relates Array_ **/
template <class T, class X> inline bool 
readUnformatted(std::istream& in, Array_<T,X>& v) {
    v.clear();
    T element;
    std::ws(in); // Make sure we're at eof if stream is all whitespace.
    while (!in.eof() && readUnformatted(in, element))
        v.push_back(element);
    return !in.fail(); // eof is expected and ok
}

/** Specialization of readUnformatted() for fixed-length ArrayView_<T,X>; reads
whitespace-separated tokens until the expected number have been read. If fewer 
are available we fail. 
@relates ArrayView_ **/
template <class T, class X> inline bool 
readUnformatted(std::istream& in, ArrayView_<T,X>& v) {
    for (X i(0); i < v.size(); ++i)
        if (!readUnformatted(in, v[i])) return false;
    return true;
}


/** Specialization of readFormatted() for variable-length Array_<T,X>; 
uses readArrayFromStream() to consume an appropriately-formatted array
until error, closing parenthesis or bracket, or eof.
@see readArrayFromStream() for details
@relates Array_ **/
template <class T, class X> inline bool 
readFormatted(std::istream& in, Array_<T,X>& v) {
    return !readArrayFromStream(in,v).fail();
}

/** Specialization of readFormatted() for fixed-length ArrayView_<T,X>; uses
fillArrayViewFromStream() to consume an appropriately-formatted fixed-size
array.
@see fillArrayViewFromStream() for details
@relates ArrayView_ **/
template <class T, class X> inline bool 
readFormatted(std::istream& in, ArrayView_<T,X>& v) {
    return !fillArrayViewFromStream(in,v).fail();
}

/** Read in an Array_<T> from a stream, as a sequence of space-separated or
comma-separated values optionally surrounded by parentheses (), square 
brackets [], or curly braces {}. We will continue to read elements of 
type T from the stream until we find a reason to stop, using type T's stream 
extraction operator>>() to read in each element and resizing the Array as
necessary. If the data is bracketed, we'll read until we hit the closing 
bracket. If it is not bracketed, we'll read until we hit eof() or get an error
such as the element extractor setting the stream's fail bit due to bad 
formatting. On successful return, the stream will be positioned right after 
the final read-in element or terminating bracket, and the stream's status will 
be good() or eof(). We will not consume trailing whitespace after bracketed 
elements; that means the stream might actually be empty even if we don't 
return eof(). If you want to know whether there is anything else in the 
stream, follow this call with the STL whitespace skipper std::ws() like this:
@code
    if (readArrayFromStream(in,array) && !in.eof()) 
        std::ws(in); // might take us to eof
    if (in.fail()) {...} // probably a formatting error
    else {
        // Here if the stream is good() then there is more to read; if the
        // stream got used up the status is guaranteed to be eof().
    }
@endcode
A compilation error will occur if you try to use this method on an Array_<T>
for a type T for which there is no stream extraction operator>>(). 
@note If you want to fill an owner Array_<T> with a fixed amount of data from
the stream, resize() the array to the appropriate length and then use 
fillArrayFromStream() instead. @see fillArrayFromStream()
@relates Array_ **/
template <class T, class X> static inline 
std::istream& readArrayFromStream(std::istream& in, Array_<T,X>& out)
{   return readArrayFromStreamHelper<T,X>(in, false /*variable sizez*/, out); }



/** Read in a fixed number of elements from a stream into an Array. We expect 
to read in exactly size() elements of type T, using type T's stream extraction 
operator>>(). This will stop reading when we've read size() elements, or set 
the fail bit in the stream if we run out of elements or if any element's 
extract operator sets the fail bit. On successful return, all size() elements 
will have been set, the stream will be positioned right after the final 
read-in element or terminating bracket, and the stream's status will be good()
or eof(). We will not consume trailing whitespace after reading all the 
elements; that means the stream might actually be empty even if we don't 
return eof(). If you want to know whether there is anything else in the 
stream, follow this call with std::ws() like this:
@code
    if (fillArrayFromStream(in,array))
        if (!in.eof()) std::ws(in); // might take us to eof
    if (in.fail()) {...} // deal with I/O or formatting error
    // Here if the stream is good() then there is more to read; if the
    // stream got used up the status is guaranteed to be eof().
@endcode
A compilation error will occur if you try to use this method on an Array_<T>
for a type T for which there is no stream extraction operator>>().
@note If you want to read in a variable number of elements and have the 
Array_<T> resized as needed, use readArrayFromStream() instead.
@see readArrayFromStream()
@relates Array_ **/
template <class T, class X> static inline 
std::istream& fillArrayFromStream(std::istream& in, Array_<T,X>& out)
{   return readArrayFromStreamHelper<T,X>(in, true /*fixed size*/, out); }

/** Read in a fixed number of elements from a stream into an ArrayView. See
fillArrayFromStream() for more information; this works the same way.
@see fillArrayFromStream()
@relates ArrayView_ **/
template <class T, class X> static inline 
std::istream& fillArrayViewFromStream(std::istream& in, ArrayView_<T,X>& out)
{   return readArrayFromStreamHelper<T,X>(in, true /*fixed size*/, out); }




/** Read an Array_<T> from a stream as a sequence of space- or comma-separated
values of type T, optionally delimited by parentheses, brackets, or braces.
The Array_<T> may be an owner (variable size) or a view (fixed size n). In
the case of an owner, we'll read all the elements in brackets or until eof if
there are no brackets. In the case of a view, there must be exactly n elements
in brackets, or if there are no brackets we'll consume exactly n elements and
then stop. Each element is read in with its own operator ">>" so this won't 
work if no such operator is defined for type T.
@relates Array_ **/
template <class T, class X> inline
std::istream& operator>>(std::istream& in, Array_<T,X>& out) 
{   return readArrayFromStream<T,X>(in, out); }

/** Read a (fixed size n) ArrayView_<T> from a stream as a sequence of space- 
or comma-separated values of type T, optionally delimited by parentheses, 
square brackets, or curly braces. If there are no delimiters then we will read
size() values and then stop. Otherwise, there must be exactly size() values 
within the brackets. Each element is read in with its own operator ">>" so 
this won't work if no such operator is defined for type T.
@relates ArrayView_ **/
template <class T, class X> inline
std::istream& operator>>(std::istream& in, ArrayView_<T,X>& out) 
{   return fillArrayViewFromStream<T,X>(in, out); }

/*@}                       End of Array serialization. **/



/**@name                    Comparison operators
These operators permit lexicographical comparisons between two comparable
Array_ objects, possibly with differing element and index types, and between 
an Array_ object and a comparable std::vector object. @relates Array_ **/
/*@{*/

/** Two Array_ objects are equal if and only if they are the same size() and 
each element compares equal using an operator T1==T2. @relates Array_ **/
template <class T1, class X1, class T2, class X2> inline bool 
operator==(const ArrayViewConst_<T1,X1>& a1, const ArrayViewConst_<T2,X2>& a2) {
    // Avoid warnings in size comparison by using common type.
    const ptrdiff_t sz1 = a1.end()-a1.begin();
    const ptrdiff_t sz2 = a2.end()-a2.begin();
    if (sz1 != sz2) return false;
    const T1* p1 = a1.begin();
    const T2* p2 = a2.begin();
    while (p1 != a1.end())
        if (!(*p1++ == *p2++)) return false;
    return true;
}
/** The not equal operator is implemented using the equal operator.  
@relates Array_ **/
template <class T1, class X1, class T2, class X2> inline bool 
operator!=(const ArrayViewConst_<T1,X1>& a1, const ArrayViewConst_<T2,X2>& a2)
{   return !(a1 == a2); }

/** Array_ objects are ordered lexicographically; that is, by first differing 
element or by length if there are no differing elements up to the length of the
shorter array (in which case the shorter one is "less than" the longer). 
This depends on T1==T2 and T1<T2 operators working. @relates Array_ **/
template <class T1, class X1, class T2, class X2> inline bool 
operator<(const ArrayViewConst_<T1,X1>& a1, const ArrayViewConst_<T2,X2>& a2) {
    const T1* p1 = a1.begin();
    const T2* p2 = a2.begin();
    while (p1 != a1.end() && p2 != a2.end()) {
        if (!(*p1 == *p2))
            return *p1 < *p2; // otherwise p1 > p2
        ++p1; ++p2;
    }
    // All elements were equal until one or both arrays ran out of elements.
    // a1 is less than a2 only if a1 ran out and a2 didn't.
    return p1 == a1.end() && p2 != a2.end();
}
/** The greater than or equal operator is implemented using the less than 
operator. @relates Array_ **/
template <class T1, class X1, class T2, class X2> inline bool 
operator>=(const ArrayViewConst_<T1,X1>& a1, const ArrayViewConst_<T2,X2>& a2)
{   return !(a1 < a2); }
/** The greater than operator is implemented by using less than with the
arguments reversed, meaning the elements must have working comparison
operators of the form T2==T1 and T2<T1. @relates Array_ **/
template <class T1, class X1, class T2, class X2> inline bool 
operator>(const ArrayViewConst_<T1,X1>& a1, const ArrayViewConst_<T2,X2>& a2)
{   return a2 < a1; }
/** The less than or equal operator is implemented using the greater than 
operator. @relates Array_ **/
template <class T1, class X1, class T2, class X2> inline bool 
operator<=(const ArrayViewConst_<T1,X1>& a1, const ArrayViewConst_<T2,X2>& a2)
{   return !(a1 > a2); }

/** An Array_<T1> and an std::vector<T2> are equal if and only if they are the 
same size() and each element compares equal using an operator T1==T2.  
@relates Array_ **/
template <class T1, class X1, class T2, class A2> inline bool 
operator==(const ArrayViewConst_<T1,X1>& a1, const std::vector<T2,A2>& v2)
{   return a1 == ArrayViewConst_<T2,size_t>(v2); }

/** An std::vector<T1> and an Array_<T2> are equal if and only if they are the 
same size() and each element compares equal using an operator T2==T1.  
@relates Array_ **/
template <class T1, class A1, class T2, class X2> inline bool 
operator==(const std::vector<T1,A1>& v1, const ArrayViewConst_<T2,X2>& a2)
{   return a2 == v1; }

/** The not equal operator is implemented using the equal operator.  
@relates Array_ **/
template <class T1, class X1, class T2, class A2> inline bool 
operator!=(const ArrayViewConst_<T1,X1>& a1, const std::vector<T2,A2>& v2)
{   return !(a1 == v2); }
/** The not equal operator is implemented using the equal operator.  
@relates Array_ **/
template <class T1, class A1, class T2, class X2> inline bool 
operator!=(const std::vector<T1,A1>& v1, const ArrayViewConst_<T2,X2>& a2)
{   return !(a2 == v1); }

/** An Array_<T1> and std::vector<T2> are ordered lexicographically; that is, 
by first differing element or by length if there are no differing elements up 
to the length of the shorter container (in which case the shorter one is 
"less than" the longer). This depends on having working element operators 
T1==T2 and T1<T2. @relates Array_ **/
template <class T1, class X1, class T2, class A2> inline bool 
operator<(const ArrayViewConst_<T1,X1>& a1, const std::vector<T2,A2>& v2)
{   return a1 < ArrayViewConst_<T2,size_t>(v2); }

/** An std::vector<T1> and Array_<T2> are ordered lexicographically; that is, 
by first differing element or by length if there are no differing elements up 
to the length of the shorter container (in which case the shorter one is 
"less than" the longer). This depends on having working element operators 
T1==T2 and T1<T2. @relates Array_ **/
template <class T1, class A1, class T2, class X2> inline bool 
operator<(const std::vector<T1,A1>& v1, const ArrayViewConst_<T2,X2>& a2)
{   return ArrayViewConst_<T1,size_t>(v1) < a2; }

/** The greater than or equal operator is implemented using the less than 
operator. @relates Array_ **/
template <class T1, class X1, class T2, class A2> inline bool 
operator>=(const ArrayViewConst_<T1,X1>& a1, const std::vector<T2,A2>& v2)
{   return !(a1 < v2); }
/** The greater than or equal operator is implemented using the less than 
operator. @relates Array_ **/
template <class T1, class A1, class T2, class X2> inline bool 
operator>=(const std::vector<T1,A1>& v1, const ArrayViewConst_<T2,X2>& a2)
{   return !(v1 < a2); }

/** The greater than operator is implemented by using less than with the
arguments reversed, meaning the elements must have working comparison
operators of the form T2==T1 and T2<T1. @relates Array_  **/
template <class T1, class X1, class T2, class A2> inline bool 
operator>(const ArrayViewConst_<T1,X1>& a1, const std::vector<T2,A2>& v2)
{   return v2 < a1; }
/** The greater than operator is implemented by using less than with the
arguments reversed, meaning the elements must have working comparison
operators of the form T2==T1 and T2<T1. @relates Array_  **/
template <class T1, class A1, class T2, class X2> inline bool 
operator>(const std::vector<T1,A1>& v1, const ArrayViewConst_<T2,X2>& a2)
{   return a2 < v1; }

/** The less than or equal operator is implemented using the greater than 
operator. @relates Array_ **/
template <class T1, class X1, class T2, class A2> inline bool 
operator<=(const ArrayViewConst_<T1,X1>& a1, const std::vector<T2,A2>& v2)
{   return !(a1 > v2); }
/** The less than or equal operator is implemented using the greater than 
operator. @relates Array_ **/
template <class T1, class A1, class T2, class X2> inline bool 
operator<=(const std::vector<T1,A1>& v1, const ArrayViewConst_<T2,X2>& a2)
{   return !(v1 > a2); }

/*@}*/

} // namespace SimTK

namespace std {
/** This is a specialization of the STL std::swap() algorithm which uses the
constant time built-in swap() member of the Array_ class. 
@relates SimTK::Array_ **/
template <class T, class X> inline void
swap(SimTK::Array_<T,X>& a1, SimTK::Array_<T,X>& a2) {
    a1.swap(a2);
}

} // namespace std
  
#endif // SimTK_SimTKCOMMON_ARRAY_H_
