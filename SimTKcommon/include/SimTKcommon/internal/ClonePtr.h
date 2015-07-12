#ifndef SimTK_SimTKCOMMON_CLONE_PTR_H_
#define SimTK_SimTKCOMMON_CLONE_PTR_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-15 Stanford University and the Authors.        *
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

#include <memory>
#include <iosfwd>
#include <cassert>

namespace SimTK {

/** Smart pointer with deep copy semantics. This is similar to `std::unique_ptr` 
in that it does not permit shared ownership of the contained object. However, 
unlike `std::unique_ptr`, %ClonePtr supports copy and assignment operations, by 
insisting that the contained object have a clone() method that returns a pointer
to a heap-allocated deep copy of the *concrete* object. The API is modeled as 
closely as possible on the C++11 `std::unique_ptr` API. However, it always uses
a default deleter. Also, in keeping with Simbody's careful distinction between 
writable and const access, and for compatibility with CloneOnWritePtr, the get() 
method is modified to return only a const pointer, with upd() (update) added to 
return a writable pointer.

This class is entirely inline and has no computational or space overhead except
when cloning is required; it contains just a single pointer and does no 
reference counting. 

@tparam T   The type of the contained object, which *must* have a `clone()` 
            method. May be an abstract or concrete type.

@see CloneOnWritePtr, ReferencePtr **/ 
template <class T> class ClonePtr {
public:
    typedef T  element_type; ///< Type of the contained object.
    typedef T* pointer;      ///< Type of a pointer to the contained object.
    typedef T& reference;    ///< Type of a reference to the contained object.
    
    /** @name                    Constructors **/
    /**@{**/

    /** Default constructor stores a `nullptr`. No heap allocation is performed.
    The empty() method will return true when called on a default-constructed 
    %ClonePtr. **/
    ClonePtr() NOEXCEPT_11 {p=nullptr;}

    /** Constructor from `nullptr` is the same as the default constructor.
    This is an implicit conversion that allows `nullptr` to be used to
    initialize a %ClonePtr. **/
    ClonePtr(std::nullptr_t) NOEXCEPT_11 : ClonePtr() {}

    /** Given a pointer to a writable heap-allocated object, take over 
    ownership of that object. The `clone()` method is *not* invoked. **/
    explicit ClonePtr(T* x) NOEXCEPT_11 : p(x) {}

    /** Given a pointer to a read-only object, create a new heap-allocated copy 
    of that object via its `clone()` method and make this %ClonePtr the owner of
    the copy. Ownership of the original object is not affected. If the supplied 
    pointer is null, the resulting %ClonePtr will be empty. **/
    explicit ClonePtr(const T* x) : ClonePtr(cloneOrNull(x)) {}

    /** Given a read-only reference to an object, create a new heap-allocated 
    copy of that object via its `clone()` method and make this %ClonePtr
    object the owner of the copy. Ownership of the original object is not
    affected. **/
    explicit ClonePtr(const T& x) : ClonePtr(&x) {}

    /** Copy constructor is deep; the new %ClonePtr object contains a new
    copy of the object in the source, created via the source object's `clone()`
    method. If the source container is empty this one will be empty also. **/
    ClonePtr(const ClonePtr& src) : p(cloneOrNull(src.p)) {}

    /** Deep copy construction from a compatible %ClonePtr. Type `U*` must be 
    implicitly convertible to type `T*`. The new %ClonePtr object contains a new
    copy of the object in the source, created via the source object's `clone()`
    method. If the source container is empty this one will be empty also. **/
    template <class U>
    ClonePtr(const ClonePtr<U>& src) : p(cloneOrNull(src.p)) {} 

    /** Move constructor is very fast and leaves the source empty. Ownership
    is transferred from the source to the new %ClonePtr. If the source was empty
    this one will be empty also. No heap activity occurs. **/
    ClonePtr(ClonePtr&& src) NOEXCEPT_11 : p(src.release()) {} 

    /** Move construction from a compatible %ClonePtr. Type `U*` must be 
    implicitly convertible to type `T*`. Ownership is transferred from the 
    source to the new %ClonePtr. If the source was empty this one will be empty 
    also. No heap activity occurs. **/
    template <class U>
    ClonePtr(ClonePtr<U>&& src) NOEXCEPT_11 : p(src.release()) {
    }
    /**@}**/

    /** @name                   Assignment **/
    /**@{**/

    /** Copy assignment replaces the currently-held object by a copy of the 
    object held in the source container, created using the source object's 
    `clone()` method. The currently-held object (if any) is deleted. If the 
    source container is empty this one will be empty also after the assignment. 
    Nothing happens if the source and destination are the same container. **/
    ClonePtr& operator=(const ClonePtr& src) {
        if (&src != this) { 
            assert((p != src.p) || !p); // can't be same ptr unless null
            reset(cloneOrNull(src.p));
        }
        return *this;
    }

    /** Copy assignment from a compatible %ClonePtr. Type `U*` must be 
    implicitly convertible to type `T*`. The currently-held object is replaced
    by a copy of the object held in the source container, created using the 
    source object's `clone()` method. The currently-held object (if any) is 
    deleted. If the source container is empty this one will be empty also after 
    the assignment.**/
    template <class U>
    ClonePtr& operator=(const ClonePtr<U>& src) {
        // The source can't be the same container as this one since they are
        // different types. The managed pointers should never be the same either 
        // since ClonePtrs represent unique ownership. (OK if both nullptr.)
        assert((p != static_cast<const T*>(src.p)) || !p);

        reset(cloneOrNull(src.p));
        return *this;
    }

    /** Move assignment replaces the currently-held object by the source object, 
    leaving the source empty. The currently-held object (if any) is deleted. 
    The `clone()` method is *not* invoked. Nothing happens if the source and 
    destination are the same containers. **/
    ClonePtr& operator=(ClonePtr&& src) NOEXCEPT_11 { 
        if (&src != this) {   
            assert((p != src.p) || !p); // can't be same ptr unless null
            reset(src.p); src.p = nullptr; 
        }
        return *this;
    }

    /** Move assignment from a compatible %ClonePtr replaces the currently-held 
    object by the source object, leaving the source empty. Type U* must be 
    implicitly convertible to type T*. The currently-held object (if any) is 
    deleted. The `clone()` method is *not* invoked. **/
    template <class U>
    ClonePtr& operator=(ClonePtr<U>&& src) NOEXCEPT_11 {
        // The source can't be the same container as this one since they are
        // different types. The managed pointers should never be the same either 
        // since ClonePtrs represent unique ownership. (OK if both nullptr.)
        assert((p != static_cast<const T*>(src.p)) || !p);
        reset(src.p); src.p = nullptr;
        return *this;
    }

    /** This form of assignment replaces the currently-held object by a 
    heap-allocated copy of the source object, created using its `clone()`
    method. The currently-held object (if any) is deleted.  **/    
    ClonePtr& operator=(const T& x)          
    {   reset(cloneOrNull(&x)); return *this; }

    /** This form of assignment replaces the currently-held object by 
    the given source object and takes over ownership of the source object. The  
    currently-held object (if any) is deleted. **/ 
    ClonePtr& operator=(T* x) NOEXCEPT_11               
    {   reset(x); return *this; }
    /**@}**/
    
    /** @name                    Destructor **/
    /**@{**/    
    /** Destructor deletes the contained object. 
    @see reset() **/
    ~ClonePtr() NOEXCEPT_11 {reset();}
    /**@}**/

    /** @name                     Accessors **/
    /**@{**/

    /** Return a const pointer to the contained object if any, or `nullptr`.
    Note that this is different than `%get()` for the standard smart pointers 
    like `std::unique_ptr` which return a writable pointer. Use upd() here for 
    that purpose. 
    @see upd(), getRef() **/
    const T* get() const NOEXCEPT_11 {return p;}

    /** Return a writable pointer to the contained object if any, or `nullptr`. 
    Note that you need write access to this container in order to get write 
    access to the object it contains.
    @see get(), updRef() **/
    T* upd() NOEXCEPT_11 {return p;}

    /** Return a const reference to the contained object. Don't call this if 
    this container is empty. 
    @see get() **/
    const T& getRef() const { 
        SimTK_ERRCHK(!empty(), "ClonePtr::getRef()", 
                    "An attempt was made to dereference a null pointer."); 
        return *get(); 
    } 

    /** Return a writable reference to the contained object. Don't call this if 
    this container is empty. @see upd() **/
    T& updRef() { 
        SimTK_ERRCHK(!empty(), "ClonePtr::updRef()", 
                    "An attempt was made to dereference a null pointer.");
        return *upd(); 
    }

    /** Dereference a const pointer to the contained object. This will fail if 
    the container is empty. **/
    const T* operator->() const { return &getRef(); }

    /** Dereference a writable pointer to the contained object. This will fail 
    if the container is empty. **/
    T* operator->() { return &updRef(); }

    /** This "dereference" operator returns a const reference to the contained 
    object. This will fail if the container is empty. **/
    const T& operator*() const {return getRef();}

    /** Return a writable reference to the contained object. This will fail if 
    the container is empty. **/
    T& operator*() {return updRef();}
    /**@}**/

    /** @name                      Utility Methods **/
    /**@{**/

    /** Make this container empty if it isn't already, destructing the contained 
    object if there is one. The container is restored to its default-constructed
    state. 
    @see empty() **/
    void reset() NOEXCEPT_11 {
        delete p; 
        p = nullptr;
    }

    /** Replace the contents of this container with the supplied heap-allocated
    object, taking over ownership of that object and deleting the current one
    first if necessary. Nothing happens if the supplied pointer is the same
    as the one already being managed. **/
    void reset(T* x) NOEXCEPT_11 {
        if (x != p) {
            delete p;
            p = x;
        }
    }

    /** Swap the contents of this %ClonePtr with another one, with 
    ownership changing hands but no copying performed. This is very fast;
    no heap activity occurs. Both containers must have been instantiated with 
    the identical type. **/
    void swap(ClonePtr& other) NOEXCEPT_11 {
        std::swap(p, other.p);
    }
   
    /** Return true if this container is empty, which is the state the container
    is in immediately after default construction and various other 
    operations. **/
    bool empty() const NOEXCEPT_11 {return !p;} // count should be null also

    /** This is a conversion to type bool that returns true if the container is 
    non-null (that is, not empty). **/
    explicit operator bool() const NOEXCEPT_11 {return !empty();}

    /** Remove the contained object from management by this container and 
    transfer ownership to the caller. A writable pointer to the object is 
    returned. No object destruction occurs. This %ClonePtr is left empty. **/
    T* release() NOEXCEPT_11 {
        T* save = p;
        p = nullptr;
        return save;
    }

    /** <b>(Deprecated)</b> Same as `get()`. Use get() instead; it is more like 
    the API for `std::unique_ptr`. **/
    DEPRECATED_14("use get() instead")
    const T* getPtr() const NOEXCEPT_11 {return get();}
    /** <b>(Deprecated)</b> Same as `upd()`. Use upd() instead; it is a better 
    match for `get()` modeled after the API for `std::unique_ptr`. **/
    DEPRECATED_14("use upd() instead")
    T* updPtr() NOEXCEPT_11 {return upd();}  
    /** <b>(Deprecated)</b> Use reset() instead. **/
    DEPRECATED_14("use reset() instead")
    void clear() NOEXCEPT_11 {reset();}

    /**@}**/

private:
template <class U> friend class ClonePtr;

    // If src is non-null, clone it; otherwise return nullptr.
    static T* cloneOrNull(const T* src) {
        return src ? src->clone() : nullptr;
    }

    T*      p;          // this may be null
};



//==============================================================================
//                       SimTK namespace-scope functions
//==============================================================================
// These namespace-scope functions will be resolved by the compiler using
// "Koenig lookup" which examines the arguments' namespaces first.
// See Herb Sutter's discussion here: http://www.gotw.ca/publications/mill08.htm.

/** This is an overload of the STL std::swap() algorithm which uses the
cheap built-in swap() member of the ClonePtr class. (This function
is defined in the `SimTK` namespace.)
@relates ClonePtr **/
template <class T> inline void
swap(ClonePtr<T>& p1, ClonePtr<T>& p2) NOEXCEPT_11 {
    p1.swap(p2);
}

/** Output the system-dependent representation of the pointer contained
in a ClonePtr object. This is equivalent to `os << p.get();`.
@relates ClonePtr **/
template <class charT, class traits, class T>
inline std::basic_ostream<charT,traits>& 
operator<<(std::basic_ostream<charT,traits>& os, 
           const ClonePtr<T>&  p) 
{   os << p.get(); return os; }

/** Compare for equality the managed pointers contained in two compatible
ClonePtr containers. Returns `true` if the pointers refer to
the same object or if both are null. It must be possible for one of the
pointer types `T*` and `U*` to be implicitly converted to the other.
@relates ClonePtr **/
template <class T, class U>
inline bool operator==(const ClonePtr<T>& lhs,
                       const ClonePtr<U>& rhs)
{   return lhs.get() == rhs.get(); }

/** Comparison against `nullptr`; same as `lhs.empty()`. 
@relates ClonePtr **/
template <class T>
inline bool operator==(const ClonePtr<T>& lhs, std::nullptr_t)
{   return lhs.empty(); }

/** Comparison against `nullptr`; same as `rhs.empty()`. 
@relates ClonePtr **/
template <class T>
inline bool operator==(std::nullptr_t, const ClonePtr<T>& rhs)
{   return rhs.empty(); }

/** Less-than operator for two compatible ClonePtr containers, 
comparing the *pointers*, not the *objects* they point to. Returns `true` if the
lhs pointer tests less than the rhs pointer. A null pointer tests less than any
non-null pointer. It must be possible for one of the pointer types `T*` and 
`U*` to be implicitly converted to the other.
@relates ClonePtr **/
template <class T, class U>
inline bool operator<(const ClonePtr<T>& lhs,
                      const ClonePtr<U>& rhs)
{   return lhs.get() < rhs.get(); }

/** Less-than comparison against a `nullptr`. A null pointer tests less than any
non-null pointer and equal to another null pointer, so this method always 
returns `false`.
@relates ClonePtr **/
template <class T>
inline bool operator<(const ClonePtr<T>& lhs, std::nullptr_t)
{   return false; }

/** Less-than comparison of a `nullptr` against this container. A null
pointer tests less than any non-null pointer and equal to another null pointer,
so this method returns `true` unless the container is empty.
@relates ClonePtr **/
template <class T>
inline bool operator<(std::nullptr_t, const ClonePtr<T>& rhs)
{   return !rhs.empty(); }


// These functions are derived from operator== and operator<.

/** Pointer inequality test defined as `!(lhs==rhs)`.
@relates ClonePtr **/
template <class T, class U>
inline bool operator!=(const ClonePtr<T>& lhs,
                       const ClonePtr<U>& rhs)
{   return !(lhs==rhs); }
/** `nullptr` inequality test defined as `!(lhs==nullptr)`.
@relates ClonePtr **/
template <class T>
inline bool operator!=(const ClonePtr<T>& lhs, std::nullptr_t)
{   return !(lhs==nullptr); }
/** `nullptr` inequality test defined as `!(nullptr==rhs)`.
@relates ClonePtr **/
template <class T>
inline bool operator!=(std::nullptr_t, const ClonePtr<T>& rhs)
{   return !(nullptr==rhs); }

/** Pointer greater-than test defined as `rhs < lhs`.
@relates ClonePtr **/
template <class T, class U>
inline bool operator>(const ClonePtr<T>& lhs,
                      const ClonePtr<U>& rhs)
{   return rhs < lhs; }
/** `nullptr` greater-than test defined as `nullptr < lhs`.
@relates ClonePtr **/
template <class T>
inline bool operator>(const ClonePtr<T>& lhs, std::nullptr_t)
{   return nullptr < lhs; }

/** `nullptr` greater-than test defined as `rhs < nullptr`.
@relates ClonePtr **/
template <class T>
inline bool operator>(std::nullptr_t, const ClonePtr<T>& rhs)
{   return rhs < nullptr; }


/** Pointer greater-or-equal test defined as `!(lhs < rhs)`.
@relates ClonePtr **/
template <class T, class U>
inline bool operator>=(const ClonePtr<T>& lhs,
                       const ClonePtr<U>& rhs)
{   return !(lhs < rhs); }
/** `nullptr` greater-or-equal test defined as `!(lhs < nullptr)`.
@relates ClonePtr **/
template <class T>
inline bool operator>=(const ClonePtr<T>& lhs, std::nullptr_t)
{   return !(lhs < nullptr); }

/** `nullptr` greater-or-equal test defined as `!(nullptr < rhs)`.
@relates ClonePtr **/
template <class T>
inline bool operator>=(std::nullptr_t, const ClonePtr<T>& rhs)
{   return !(nullptr < rhs); }


/** Pointer less-or-equal test defined as `!(rhs < lhs)` (note reversed
arguments).
@relates ClonePtr **/
template <class T, class U>
inline bool operator<=(const ClonePtr<T>& lhs,
                       const ClonePtr<U>& rhs)
{   return !(rhs < lhs); }
/** `nullptr` less-or-equal test defined as `!(nullptr < lhs)` (note reversed
arguments).
@relates ClonePtr **/
template <class T>
inline bool operator<=(const ClonePtr<T>& lhs, std::nullptr_t)
{   return !(nullptr < lhs); }
/** `nullptr` less-or-equal test defined as `!(rhs < nullptr)` (note reversed
arguments).
@relates ClonePtr **/
template <class T>
inline bool operator<=(std::nullptr_t, const ClonePtr<T>& rhs)
{   return !(rhs < nullptr); }

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_CLONE_PTR_H_
