#ifndef SimTK_SimTKCOMMON_REFERENCE_PTR_H_
#define SimTK_SimTKCOMMON_REFERENCE_PTR_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2012-15 Stanford University and the Authors.        *
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

#include "SimTKcommon/internal/common.h"
#include <cassert>
#include <utility>

namespace SimTK {

/** This is a smart pointer that implements "cross reference" semantics where
a pointer data member of some object is intended to refer to some target
object in a larger data structure. Judicious use of this container will allow
you to use compiler-generated copy constructors and copy assignment operators
for classes which would otherwise have to implement their own in order to
properly initialize these pointer data members, which must not be copied.

The contained pointer is initialized to `nullptr` on construction, and it is 
reinitialized to null upon copy construction or copy assignment. That's 
because we are assuming this is part of copying the entire data structure and 
copying the old pointer would create a reference into the old data structure 
rather than the new copy. This pointer does not own the target to which it 
points, and there is no reference counting so it will become stale if the 
target is deleted.

The pointer *is* moved intact for move construction or move assignment. That
allows `std::vector<ReferencePtr<T>>` or `SimTK::Array_<ReferencePtr<T>>` to
behave properly when their contents have to be moved for expansion.

Whether you can write through the pointer is controlled by whether type T
is a const type. For example %ReferencePtr\<int> is equivalent to an int*, while
%ReferencePtr\<const int> is equivalent to a const int*.

This class is entirely inline and has no computational or space overhead; it
contains just a single pointer. 

@see ClonePtr, CloneOnWritePtr **/ 
template <class T> class ReferencePtr {
public:
    typedef T  element_type; ///< Type of the contained object.
    typedef T* pointer;      ///< Type of a pointer to the contained object.
    typedef T& reference;    ///< Type of a reference to the contained object.

    /** @name                    Constructors **/
    /**@{**/
    /** Default constructor creates an empty object. **/
    ReferencePtr() noexcept : p(nullptr) {}

    /** Constructor from `nullptr` is the same as the default constructor.
    This is an implicit conversion that allows `nullptr` to be used to
    initialize a %ReferencePtr. **/
    ReferencePtr(std::nullptr_t) noexcept : ReferencePtr() {}

    /** Construct from a given pointer stores the pointer. **/
    explicit ReferencePtr(T* tp) noexcept : p(tp) {}

    /** Construct from a reference stores the address of the supplied 
    object. **/
    explicit ReferencePtr(T& t) noexcept : ReferencePtr(&t) {}

    /** Copy constructor unconditionally sets the pointer to null; see class
    comments for why. **/
    ReferencePtr(const ReferencePtr&) noexcept : p(nullptr) {}

    /** Move constructor copies the pointer from the source and leaves the
    source empty. **/
    ReferencePtr(ReferencePtr&& src) noexcept : p(src.release()) {}

    /** <b>(Deprecated)</b> Use %ReferencePtr(nullptr) or just %ReferencePtr()
    instead. For backwards compatibility, this allows initialization 
    by "0" rather than `nullptr`. **/
    DEPRECATED_14("use ReferencePtr(nullptr) instead")
    ReferencePtr(int mustBeZero) noexcept : ReferencePtr()
    {   assert(mustBeZero==0); }
    /**@}**/

    /** @name                   Assignment **/
    /**@{**/

    /** Copy assignment sets the pointer to nullptr (except for a self-assign); 
    see class comments for why.  **/
    ReferencePtr& operator=(const ReferencePtr& src) noexcept
    {   if (&src != this) reset(); return *this; }

    /** Move assignment copies the pointer from the source and leaves the
    source empty. Nothing happens for self-assign. **/
    ReferencePtr& operator=(ReferencePtr&& src) noexcept {
        if (&src != this) 
            p = src.release();
        return *this; 
    }

    /** This form of assignment replaces the currently-referenced object by a 
    reference to the source object; no destruction occurs. **/    
    ReferencePtr& operator=(T& t) noexcept    
    {  reset(&t); return *this; }

    /** This form of assignment replaces the current pointer with the given
    one; no destruction occurs. **/
    ReferencePtr& operator=(T* tp) noexcept               
    {   reset(tp); return *this; }
    
    /** @name                    Destructor **/
    /**@{**/
    /** Destructor does nothing. **/
    ~ReferencePtr() noexcept {reset();} // just being tidy
    /**@}**/

    /** @name                     Accessors **/
    /**@{**/

    /** Return the contained pointer, or null if the container is empty. **/
    T* get() const noexcept {return p;}

    /** Return a reference to the target object. Fails if the pointer is
    null. **/
    T& getRef() const { 
        SimTK_ERRCHK(p!=nullptr, "ReferencePtr::getRef()", 
                    "An attempt was made to dereference a null pointer."); 
        return *p; 
    }

    /** Return the contained pointer. This will fail if the container is 
    empty. **/
    T* operator->() const {return &getRef();}

    /** The "dereference" operator returns a reference to the target object. 
    This will fail if the container is empty. **/
    T& operator*() const {return getRef();}

    /**@}**/

    /** @name                      Utility Methods **/
    /**@{**/

    /** Replace the stored pointer with a different one; no destruction
    occurs. **/
    void reset(T* tp=nullptr) noexcept {p=tp;}

    /** Swap the contents of this %ReferencePtr with another one. This is very
    fast. **/
    void swap(ReferencePtr& other) noexcept {
        std::swap(p, other.p);
    }

    /** Return true if this container is empty. **/
    bool empty() const noexcept {return !p;}

    /** This is a conversion to type bool that returns true if
    the container is non-null (that is, not empty). **/
    explicit operator bool() const noexcept {return !empty();}

    /** Extract the pointer from this container, leaving the container empty. 
    The pointer is returned. **/
    T* release() noexcept {
        T* x=p; 
        p=nullptr; 
        return x;
    }

    /** <b>(Deprecated)</b> Use reset() instead. **/
    DEPRECATED_14("use reset() instead")
    void clear() noexcept {reset();}

    /** <b>(Deprecated)</b> Use get() rather than implicit conversion from 
    %ReferencePtr\<T> to T*. **/
    DEPRECATED_14("use get() rather than implicit conversion to T*")
    operator T*() const noexcept {return p;}

    /**@}**/

private:
    T*    p;        // can be nullptr
};    



//==============================================================================
//                       SimTK namespace-scope functions
//==============================================================================
// These namespace-scope functions will be resolved by the compiler using
// "Koenig lookup" which examines the arguments' namespaces first.
// See Herb Sutter's discussion here: http://www.gotw.ca/publications/mill08.htm.

/** This is an overload of the STL std::swap() algorithm which uses the
cheap built-in swap() member of the ReferencePtr class. (This function
is defined in the `SimTK` namespace.)
@relates SimTK::ReferencePtr **/
template <class T> inline void
swap(ReferencePtr<T>& p1, ReferencePtr<T>& p2) noexcept {
    p1.swap(p2);
}

/** Output the system-dependent representation of the pointer contained
in a ReferencePtr object. This is equivalent to `os << p.get();`.
@relates ReferencePtr **/
template <class charT, class traits, class T>
inline std::basic_ostream<charT,traits>& 
operator<<(std::basic_ostream<charT,traits>& os, 
           const ReferencePtr<T>&  p) 
{   os << p.get(); return os; }

/** Compare for equality the managed pointers contained in two compatible
ReferencePtr containers. Returns `true` if the pointers refer to
the same object or if both are null. It must be possible for one of the
pointer types `T*` and `U*` to be implicitly converted to the other.
@relates ReferencePtr **/
template <class T, class U>
inline bool operator==(const ReferencePtr<T>& lhs,
                       const ReferencePtr<U>& rhs)
{   return lhs.get() == rhs.get(); }

/** Comparison against `nullptr`; same as `lhs.empty()`. 
@relates ReferencePtr **/
template <class T>
inline bool operator==(const ReferencePtr<T>& lhs, std::nullptr_t)
{   return lhs.empty(); }

/** Comparison against `nullptr`; same as `rhs.empty()`. 
@relates ReferencePtr **/
template <class T>
inline bool operator==(std::nullptr_t, const ReferencePtr<T>& rhs)
{   return rhs.empty(); }

/** Less-than operator for two compatible ReferencePtr containers, 
comparing the *pointers*, not the *objects* they point to. Returns `true` if the
lhs pointer tests less than the rhs pointer. A null pointer tests less than any
non-null pointer. It must be possible for one of the pointer types `T*` and 
`U*` to be implicitly converted to the other.
@relates ReferencePtr **/
template <class T, class U>
inline bool operator<(const ReferencePtr<T>& lhs,
                      const ReferencePtr<U>& rhs)
{   return lhs.get() < rhs.get(); }

/** Less-than comparison against a `nullptr`. A null pointer tests less than any
non-null pointer and equal to another null pointer, so this method always 
returns `false`.
@relates ReferencePtr **/
template <class T>
inline bool operator<(const ReferencePtr<T>& lhs, std::nullptr_t)
{   return false; }

/** Less-than comparison of a `nullptr` against this container. A null
pointer tests less than any non-null pointer and equal to another null pointer,
so this method returns `true` unless the container is empty.
@relates ReferencePtr **/
template <class T>
inline bool operator<(std::nullptr_t, const ReferencePtr<T>& rhs)
{   return !rhs.empty(); }


// These functions are derived from operator== and operator<.

/** Pointer inequality test defined as `!(lhs==rhs)`.
@relates ReferencePtr **/
template <class T, class U>
inline bool operator!=(const ReferencePtr<T>& lhs,
                       const ReferencePtr<U>& rhs)
{   return !(lhs==rhs); }

/** `nullptr` inequality test defined as `!(lhs==nullptr)`.
@relates ReferencePtr **/
template <class T>
inline bool operator!=(const ReferencePtr<T>& lhs, std::nullptr_t)
{   return !(lhs==nullptr); }

/** `nullptr` inequality test defined as `!(nullptr==rhs)`.
@relates ReferencePtr **/
template <class T>
inline bool operator!=(std::nullptr_t, const ReferencePtr<T>& rhs)
{   return !(nullptr==rhs); }

/** Pointer greater-than test defined as `rhs < lhs`.
@relates ReferencePtr **/
template <class T, class U>
inline bool operator>(const ReferencePtr<T>& lhs,
                      const ReferencePtr<U>& rhs)
{   return rhs < lhs; }
/** `nullptr` greater-than test defined as `nullptr < lhs`.
@relates ReferencePtr **/
template <class T>
inline bool operator>(const ReferencePtr<T>& lhs, std::nullptr_t)
{   return nullptr < lhs; }

/** `nullptr` greater-than test defined as `rhs < nullptr`.
@relates ReferencePtr **/
template <class T>
inline bool operator>(std::nullptr_t, const ReferencePtr<T>& rhs)
{   return rhs < nullptr; }


/** Pointer greater-or-equal test defined as `!(lhs < rhs)`.
@relates ReferencePtr **/
template <class T, class U>
inline bool operator>=(const ReferencePtr<T>& lhs,
                       const ReferencePtr<U>& rhs)
{   return !(lhs < rhs); }
/** `nullptr` greater-or-equal test defined as `!(lhs < nullptr)`.
@relates ReferencePtr **/
template <class T>
inline bool operator>=(const ReferencePtr<T>& lhs, std::nullptr_t)
{   return !(lhs < nullptr); }

/** `nullptr` greater-or-equal test defined as `!(nullptr < rhs)`.
@relates ReferencePtr **/
template <class T>
inline bool operator>=(std::nullptr_t, const ReferencePtr<T>& rhs)
{   return !(nullptr < rhs); }


/** Pointer less-or-equal test defined as `!(rhs < lhs)` (note reversed
arguments).
@relates ReferencePtr **/
template <class T, class U>
inline bool operator<=(const ReferencePtr<T>& lhs,
                       const ReferencePtr<U>& rhs)
{   return !(rhs < lhs); }
/** `nullptr` less-or-equal test defined as `!(nullptr < lhs)` (note reversed
arguments).
@relates ReferencePtr **/
template <class T>
inline bool operator<=(const ReferencePtr<T>& lhs, std::nullptr_t)
{   return !(nullptr < lhs); }
/** `nullptr` less-or-equal test defined as `!(rhs < nullptr)` (note reversed
arguments).
@relates ReferencePtr **/
template <class T>
inline bool operator<=(std::nullptr_t, const ReferencePtr<T>& rhs)
{   return !(rhs < nullptr); }


/** <b>(Deprecated)</b> Use `nullptr` instead of `0`.
@relates ReferencePtr **/
template <class T>
DEPRECATED_14("use refptr != nullptr instead")
inline bool operator!=(const ReferencePtr<T>& lhs, int mustBeZero)
{   assert(mustBeZero==0); return !(lhs==nullptr); }

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_REFERENCE_PTR_H_
