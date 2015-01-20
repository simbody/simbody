#ifndef SimTK_SimTKCOMMON_CLONE_ON_WRITE_PTR_H_
#define SimTK_SimTKCOMMON_CLONE_ON_WRITE_PTR_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2015 Stanford University and the Authors.           *
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

#include <memory>
#include <iosfwd>

namespace SimTK {

//==============================================================================
//                          CLONE ON WRITE PTR
//==============================================================================
/** Smart pointer with deep copy semantics but with the copying delayed until
an attempt is made to write on the contained object. 

This is like `std::shared_ptr` when an object is being read, but like 
SimTK::ClonePtr when the object is written. Like SimTK::ClonePtr, 
%CloneOnWritePtr supports copy and assigment operations, by insisting that the 
contained object have a `clone()` method that returns a pointer to a 
heap-allocated deep copy of the *concrete* object. The API is modeled as closely
as possible to the C++11 `std::shared_ptr` and `std::unique_ptr`. The get() 
method is modified to return a const pointer to avoid accidental copying, with 
upd() (update) added to return a writable pointer.

This class is entirely inline and has no computational or space overhead
beyond the cost of dealing with the reference count, except when a copy has
to be made due to a write attempt.

@tparam T   The type of the contained object, which *must* have a `clone()` 
            method. May be an abstract or concrete type.

@see ClonePtr, ReferencePtr **/ 
template <class T> class CloneOnWritePtr {
public:
    typedef T  element_type; ///< Type of the contained object.
    typedef T* pointer;      ///< Type of a pointer to the contained object.
    typedef T& reference;    ///< Type of a reference to the contained object.
    
    /** @name                    Constructors **/
    /**@{**/

    /** Default constructor stores a `nullptr` and sets use count to zero. No
    heap allocation is performed. The empty() method will return true when
    called on a default-constructed %CloneOnWritePtr. **/
    CloneOnWritePtr() {init();}

    /** Constructor from `nullptr` is the same as the default constructor.
    This is an implicit conversion that allows `nullptr` to be used to
    initialize a %CloneOnWritePtr. **/
    CloneOnWritePtr(std::nullptr_t) : CloneOnWritePtr() {}

    /** Given a pointer to a writable heap-allocated object, take over 
    ownership of that object. The use count will be one unless the pointer
    was null in which case it will be zero. **/
    explicit CloneOnWritePtr(T* x) : CloneOnWritePtr()
    {   if (x) {p=x; count=new long(1);} } 

    /** Given a pointer to a read-only object, create a new heap-allocated 
    copy of that object via its `clone()` method and make this %CloneOnWritePtr
    the owner of the copy. Ownership of the original object is not
    affected. If the supplied pointer is null, the resulting %CloneOnWritePtr 
    is as though default constructed, otherwise the use count will be one. **/
    explicit CloneOnWritePtr(const T* x) : CloneOnWritePtr(cloneOrNull(x)) {}

    /** Given a read-only reference to an object, create a new heap-allocated 
    copy of that object via its `clone()` method and make this %CloneOnWritePtr
    object the owner of the copy. Ownership of the original object is not
    affected. The use count will be one after construction. **/
    explicit CloneOnWritePtr(const T& x) : CloneOnWritePtr(&x) {}

    /** Copy constructor is deep but deferred so very fast here; the new 
    %CloneOnWritePtr object initially shares the source object, but if either
    source or destination are written to subsequently a deep copy is made and
    the objects become disconnected. If the source container is empty this one
    will be as though default constructed. **/
    CloneOnWritePtr(const CloneOnWritePtr& src) : CloneOnWritePtr() 
    {   shareWith(src); }

    /** Copy construction from a compatible %CloneOnWritePtr. Type `U*` must
    be implicitly convertible to type `T*`. **/
    template <class U>
    CloneOnWritePtr(const CloneOnWritePtr<U>& src) : CloneOnWritePtr() 
    {   shareWith<U>(src); }

    /** Move constructor is very fast and leaves the source empty. The use
    count is unchanged. If the source was empty this one will be as though
    default constructed. **/
    CloneOnWritePtr(CloneOnWritePtr&& src) : CloneOnWritePtr() 
    {   moveFrom(src); }

    /** Move construction from a compatible %CloneOnWritePtr. Type `U*` must
    be implicitly convertible to type `T*`. **/
    template <class U>
    CloneOnWritePtr(CloneOnWritePtr<U>&& src) : CloneOnWritePtr()
    {   moveFrom<U>(std::move(src)); } // std::move shouldn't be needed
    /**@}**/

    /** @name                   Assignment **/
    /**@{**/

    /** Copy assignment replaces the currently-held object by a deferred
    copy of the object held in the source container. The copy will be created 
    upon a subsequent write using the source object's `clone()` method.
    If the source container is empty this one will be empty after the 
    assignment. Nothing happens if the source and destination were already
    managing the same object. **/
    CloneOnWritePtr& operator=(const CloneOnWritePtr& src) { 
        if (src.p != p) 
        {   reset(); shareWith(src); }
        return *this;
    }

    /** Copy assignment from a compatible %CloneOnWritePtr. Type `U*` must
    be implicitly convertible to type `T*`. **/
    template <class U>
    CloneOnWritePtr& operator=(const CloneOnWritePtr<U>& src) { 
        if (static_cast<T*>(src.p) != p) 
        {   reset(); shareWith<U>(src); }
        return *this;
    }

    /** Move assignment replaces the currently-held object by the source
    object, leaving the source empty. The currently-held object's use count
    will be reduced by one and deleted if this was the last use. Nothing
    happens if the source and destination are the same containers. If they are
    different but are sharing the same object then the use count is reduced
    by one. **/
    CloneOnWritePtr& operator=(CloneOnWritePtr&& src) { 
        // The std::move here shouldn't be necessary but VS2013 needed it.
        if (&src != this) 
        {   reset(); moveFrom(std::move(src)); }
        return *this;
    }

    /** Move assignment from a compatible %CloneOnWritePtr. Type U* must
    be implicitly convertible to type T*. **/
    template <class U>
    CloneOnWritePtr& operator=(CloneOnWritePtr<U>&& src) {
        // Can't be the same container since the type is different.
        // The std::move here shouldn't be necessary but VS2013 needed it.
        reset(); moveFrom<U>(std::move(src));
        return *this;
    }

    /** This form of assignment replaces the currently-held object by a 
    heap-allocated copy of the source object, created using its `clone()`
    method. The use count of the currently-held object (if any) is decremented
    and the object is deleted if this was the last reference to it. On return
    the use count of this container will be one. **/    
    CloneOnWritePtr& operator=(const T& x)          
    {   reset(cloneOrNull(&x)); return *this; }

    /** This form of assignment replaces the currently-held object by the given
    source object and takes over ownership of the source object. The use count
    of the currently-held object is decremented and the object is deleted if 
    this was the last reference to it. **/ 
    CloneOnWritePtr& operator=(T* x)               
    {   reset(x); return *this; }
    /**@}**/
    
    /** @name                    Destructor **/
    /**@{**/    
    /** Destructor decrements the reference count and deletes the object
    if the count goes to zero. @see reset() **/
    ~CloneOnWritePtr() {reset();}
    /**@}**/

    /** @name                     Accessors **/
    /**@{**/

    /** Return a const pointer to the contained object if any, or `nullptr`. No
    cloning is needed so this is very fast. Note that this is different than 
    `%get()` for the standard smart pointers which return a writable pointer. 
    Use upd() here for that purpose. 
    @see upd(), getRef() **/
    const T* get() const {return p;}

    /** Clone if necessary to ensure the contained object is not shared, then 
    return a writable pointer to the contained (and now unshared) object if any,
    or `nullptr`. If you only need read access, use get() instead to avoid the
    potentially expensive cloning. Note that you need write access to this 
    container in order to get write access to the object it contains.
    @see get(), updRef(), detach() **/
    T* upd() {detach(); return p;}

    /** Return a const reference to the contained object. No cloning is needed
    so this is very fast. Don't call this if this container is empty. 
    @see get() **/
    const T& getRef() const { 
        SimTK_ERRCHK(!empty(), "CloneOnWritePtr::getRef()", 
                    "An attempt was made to dereference a null pointer."); 
        return *get(); 
    } 

    /** Clone if necessary to ensure the contained object is not shared, then 
    return a writable reference to the contained (and now unshared) object.
    Don't call this if this container is empty. @see upd() **/
    T& updRef() { 
        SimTK_ERRCHK(!empty(), "CloneOnWritePtr::updRef()", 
                    "An attempt was made to dereference a null pointer.");
        return *upd(); 
    }

    /** Dereference a const pointer to the contained object. This will fail if 
    the container is empty. 
    @warning This `const` operator will only be invoked if it is applied to a 
    `const` container; otherwise the non-const method will be called and an 
    unwanted copy may be performed. Use getRef() instead if you have a writable
    container but don't need write access to the contained object. **/
    const T* operator->() const { return &getRef(); }

    /** Clone if necessary, then dereference a writable pointer to the contained 
    object. This will fail if the container is empty. **/
    T* operator->() { return &updRef(); }

    /** This "dereference" operator returns a const reference to the contained 
    object. This will fail if the container is empty.
    @warning This `const` method will only be invoked if it is applied to a 
    `const` container; otherwise the non-const method will be called and an 
    unwanted copy may be performed. Use getRef() instead if you have a writable
    container but don't need write access to the contained object. **/
    const T& operator*() const {return getRef();}

    /** Clone if necessary, then return a writable reference to the 
    contained object. This will fail if the container is empty. **/
    T& operator*() {return updRef();}
    /**@}**/

    /** @name                      Utility Methods **/
    /**@{**/

    /** Make this container empty, decrementing the use count of the contained
    object (if any), and deleting it if this was the last use. The container
    is restored to its default-constructed state. 
    @see empty() **/
    void reset() {
        if (empty()) return;
        if (decr()==0) {delete p; delete count;} 
        init();
    }

    /** Replace the contents of this container with the supplied heap-allocated
    object, taking over ownership of that object and deleting the current one
    first if necessary. Nothing happens if the supplied pointer is the same
    as the one already being managed. Otherwise, on return the use count will be
    one if the supplied pointer was non-null or else zero. **/
    void reset(T* x) {
        if (x != p) {
            reset();
            if (x) {p=x; count=new long(1);}
        }
    }

    /** Swap the contents of this %CloneOnWritePtr with another one, with 
    ownership changing hands but no copying performed. This is very fast;
    no heap activity occurs. Both containers must have been instantiated with 
    the identical type. **/
    void swap(CloneOnWritePtr& other) {
        std::swap(p, other.p);
        std::swap(count, other.count);
    }

    /** Return count of how many %CloneOnWritePtr objects are currently
    sharing the referenced object. There is never more than
    one holding an object for writing. If the pointer is null the use 
    count is zero. **/
    long use_count() const {return count ? *count : 0;}

    /** Is this the only user of the referenced object? Note that this means
    there is exactly one; if the managed pointer is null `unique()` returns 
    `false`. **/
    bool unique() const {return use_count()==1;}
   
    /** Return true if this container is empty, which is the state the container
    is in immediately after default construction and various other 
    operations. **/
    bool empty() const {return !p;} // count should be null also

    /** This is a conversion to type bool that returns true if
    the container is non-null (that is, not empty). **/
    explicit operator bool() const {return !empty();}

    /** (Advanced) Remove the contained object from management by this 
    container and transfer ownership to the caller. Clone if necessary to 
    ensure that the contained object is not being shared with any other 
    container, then extract the object from this container, leaving the 
    container empty. A writable pointer to the object is returned. No object 
    destruction occurs. 
    @see detach() **/
    T* release() {
        detach(); // now use count is 1 or 0
        T* save = p; delete count; init();
        return save;
    }

    /** (Advanced) %Force the contained object to be unique, that is, not shared
    with any other container. (Normally this is done automatically when 
    necessary; this method is like useful mostly for debugging.) If the 
    referenced object is being shared (that is, use_count()>1), clone it and 
    replace the pointer with the newly cloned object. This %CloneOnWritePtr will
    have a use count of zero or one after this call; other sharers of the object 
    will see the shared use count reduced by one. If this is empty() or 
    unique() already then nothing happens. Note that you have to have write
    access to this container in order to detach it. **/
    void detach() {
        if (use_count() > 1) 
        {   decr(); p=p->clone(); count=new long(1); }
    }
    /**@}**/
     
private:
template <class U> friend class CloneOnWritePtr;

    // If src is non-null, clone it; otherwise return null.
    static T* cloneOrNull(const T* src) {
        return src ? src->clone() : nullptr;
    }

    // Set an empty pointer to share with the given object. Type U* must be
    // implicitly convertible to type T*.
    template <class U> void shareWith(const CloneOnWritePtr<U>& src) {
        assert(!(p||count)); 
        if (!src.empty()) {p=src.p; count=src.count; incr();}
    }

    // Steal the object and count from the source to initialize this *empty* 
    // pointer, leaving the source empty.
    template <class U> void moveFrom(CloneOnWritePtr<U>&& src) {
        assert(!(p||count)); 
        p=src.p; count=src.count; src.init();
    }

    // Increment/decrement use count and return the result.
    long incr() const {assert(count && *count>=0); return ++(*count);}
    long decr() const {assert(count && *count>=1); return --(*count);}

    void init() {p=nullptr; count=nullptr;}

    // Can't use std::shared_ptr here due to lack of release() method.
    T*      p;          // this may be null
    long*   count;      // if p is null so is count
};    
    
} // namespace SimTK



//==============================================================================
//                          std:: namespace
//==============================================================================
namespace std {
/** This is a specialization of the STL std::swap() algorithm which uses the
cheap built-in swap() member of the CloneOnWritePtr class. (This function
is defined in the `std` namespace.) 
@relates SimTK::CloneOnWritePtr **/
template <class T> inline void
swap(SimTK::CloneOnWritePtr<T>& p1, SimTK::CloneOnWritePtr<T>& p2) {
    p1.swap(p2);
}
} // namespace std



//==============================================================================
//                          global namespace
//==============================================================================
/** Output the system-dependent representation of the pointer contained
in a SimTK::CloneOnWritePtr object. This is equivalent to `os << p.get();`.
This operator is defined in the global namespace.
@relates SimTK::CloneOnWritePtr **/
template <class charT, class traits, class T>
inline std::basic_ostream<charT,traits>& 
operator<<(std::basic_ostream<charT,traits>& os, 
           const SimTK::CloneOnWritePtr<T>&  p) 
{   os << p.get(); }

/** Compare for equality the managed pointers contained in two compatible
SimTK::CloneOnWritePtr containers. Returns `true` if the pointers refer to
the same object or if both are null. It must be possible for one of the
pointer types `T*` and `U*` to be implicitly converted to the other.
This operator is defined in the global namespace.
@relates SimTK::CloneOnWritePtr **/
template <class T, class U>
inline bool operator==(const SimTK::CloneOnWritePtr<T>& lhs,
                       const SimTK::CloneOnWritePtr<U>& rhs)
{   return lhs.get() == rhs.get(); }

/** Comparison against `nullptr`; same as `lhs.empty()`. 
This operator is defined in the global namespace.
@relates SimTK::CloneOnWritePtr **/
template <class T>
inline bool operator==(const SimTK::CloneOnWritePtr<T>& lhs, std::nullptr_t)
{   return lhs.empty(); }

/** Comparison against `nullptr`; same as `rhs.empty()`. 
This operator is defined in the global namespace.
@relates SimTK::CloneOnWritePtr **/
template <class T>
inline bool operator==(std::nullptr_t, const SimTK::CloneOnWritePtr<T>& rhs)
{   return rhs.empty(); }

/** Less-than operator for two compatible SimTK::CloneOnWritePtr containers, 
comparing the *pointers*, not the *objects* they point to. Returns `true` if the
lhs pointer tests less than the rhs pointer. A null pointer tests less than any
non-null pointer. It must be possible for one of the pointer types `T*` and 
`U*` to be implicitly converted to the other. This operator is defined in the 
global namespace. @relates SimTK::CloneOnWritePtr **/
template <class T, class U>
inline bool operator<(const SimTK::CloneOnWritePtr<T>& lhs,
                      const SimTK::CloneOnWritePtr<U>& rhs)
{   return lhs.get() < rhs.get(); }

/** Less-than comparison against a `nullptr`. A null pointer tests less than any
non-null pointer and equal to another null pointer, so this method always 
returns `false`. This operator is defined in the global namespace.
@relates SimTK::CloneOnWritePtr **/
template <class T>
inline bool operator<(const SimTK::CloneOnWritePtr<T>& lhs, std::nullptr_t)
{   return false; }

/** Less-than comparison of a `nullptr` against this container. A null
pointer tests less than any non-null pointer and equal to another null pointer,
so this method returns `true` unless the container is empty. This operator is 
defined in the global namespace. @relates SimTK::CloneOnWritePtr **/
template <class T>
inline bool operator<(std::nullptr_t, const SimTK::CloneOnWritePtr<T>& rhs)
{   return !rhs.empty(); }


// These functions are derived from operator== and operator<.

/** Pointer inequality test defined as `!(lhs==rhs)`. This operator is defined 
in the global namespace. @relates SimTK::CloneOnWritePtr **/
template <class T, class U>
inline bool operator!=(const SimTK::CloneOnWritePtr<T>& lhs,
                       const SimTK::CloneOnWritePtr<U>& rhs)
{   return !(lhs==rhs); }
/** `nullptr` inequality test defined as `!(lhs==nullptr)`. This operator is 
defined in the global namespace. @relates SimTK::CloneOnWritePtr **/
template <class T>
inline bool operator!=(const SimTK::CloneOnWritePtr<T>& lhs, std::nullptr_t)
{   return !(lhs==nullptr); }
/** `nullptr` inequality test defined as `!(nullptr==rhs)`. This operator is 
defined in the global namespace. @relates SimTK::CloneOnWritePtr **/
template <class T>
inline bool operator!=(std::nullptr_t, const SimTK::CloneOnWritePtr<T>& rhs)
{   return !(nullptr==rhs); }

/** Pointer greater-than test defined as `rhs < lhs`. This operator is defined 
in the global namespace. @relates SimTK::CloneOnWritePtr **/
template <class T, class U>
inline bool operator>(const SimTK::CloneOnWritePtr<T>& lhs,
                      const SimTK::CloneOnWritePtr<U>& rhs)
{   return rhs < lhs; }
/** `nullptr` greater-than test defined as `nullptr < lhs`. This operator is 
defined in the global namespace. @relates SimTK::CloneOnWritePtr **/
template <class T>
inline bool operator>(const SimTK::CloneOnWritePtr<T>& lhs, std::nullptr_t)
{   return nullptr < lhs; }

/** `nullptr` greater-than test defined as `rhs < nullptr`. This operator is 
defined in the global namespace. @relates SimTK::CloneOnWritePtr **/
template <class T>
inline bool operator>(std::nullptr_t, const SimTK::CloneOnWritePtr<T>& rhs)
{   return rhs < nullptr; }


/** Pointer greater-or-equal test defined as `!(lhs < rhs)`. This operator is 
defined in the global namespace. @relates SimTK::CloneOnWritePtr **/
template <class T, class U>
inline bool operator>=(const SimTK::CloneOnWritePtr<T>& lhs,
                       const SimTK::CloneOnWritePtr<U>& rhs)
{   return !(lhs < rhs); }
/** `nullptr` greater-or-equal test defined as `!(lhs < nullptr)`. This operator
is defined in the global namespace. @relates SimTK::CloneOnWritePtr **/
template <class T>
inline bool operator>=(const SimTK::CloneOnWritePtr<T>& lhs, std::nullptr_t)
{   return !(lhs < nullptr); }

/** `nullptr` greater-or-equal test defined as `!(nullptr < rhs)`. This operator
is defined in the global namespace. @relates SimTK::CloneOnWritePtr **/
template <class T>
inline bool operator>=(std::nullptr_t, const SimTK::CloneOnWritePtr<T>& rhs)
{   return !(nullptr < rhs); }


/** Pointer less-or-equal test defined as `!(rhs < lhs)` (note reversed
arguments). This operator is defined in the global namespace.
@relates SimTK::CloneOnWritePtr **/
template <class T, class U>
inline bool operator<=(const SimTK::CloneOnWritePtr<T>& lhs,
                       const SimTK::CloneOnWritePtr<U>& rhs)
{   return !(rhs < lhs); }
/** `nullptr` less-or-equal test defined as `!(nullptr < lhs)` (note reversed
arguments). This operator is defined in the global namespace.
@relates SimTK::CloneOnWritePtr **/
template <class T>
inline bool operator<=(const SimTK::CloneOnWritePtr<T>& lhs, std::nullptr_t)
{   return !(nullptr < lhs); }
/** `nullptr` less-or-equal test defined as `!(rhs < nullptr)` (note reversed
arguments). This operator is defined in the global namespace.
@relates SimTK::CloneOnWritePtr **/
template <class T>
inline bool operator<=(std::nullptr_t, const SimTK::CloneOnWritePtr<T>& rhs)
{   return !(rhs < nullptr); }


#endif // SimTK_SimTKCOMMON_CLONE_ON_WRITE_PTR_H_
