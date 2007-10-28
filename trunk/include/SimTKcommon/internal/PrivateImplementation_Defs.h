#ifndef SimTK_PRIVATE_IMPLEMENTATION_DEFS_H_
#define SimTK_PRIVATE_IMPLEMENTATION_DEFS_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
 * Authors: Michael Sherman                                                   *
 * Contributors: Chris Bruns and Peter Eastman                                *
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

/**@file
 * This header provides the definitions of the PIMPLHandle template methods
 * as declared in PrivateImplementation.h, and also the
 * declaration and definition of the template classes used for creating well-
 * behaved implementation classes.
 *
 * Although this header is part of the SimTK Core installation, it is not automatically
 * included with SimTKcommon.h. It is intended to be included explicitly in
 * compilation units where the private implementations are being defined.
 * When used in SimTK Core software, this header will be included only
 * in library-side compilation units. 
 */

#include "SimTKcommon/internal/PrivateImplementation.h"

#include <cassert>
#include <iostream>
#include <typeinfo>

namespace SimTK {

    /////////////////////////
    // PIMPLImplementation //
    /////////////////////////

/**
 * This class provides some infrastructure useful in creating PIMPL 
 * Implementation classes (the ones referred to by Handles). Note that
 * this class is used by SimTK Core code ONLY on the library side; it never
 * appears in headers intended for use by clients. However it is 
 * generally useful enough that we include it here to assist people
 * who would like to make their own PIMPL classes. Consequently, there
 * are no binary compatibility issues raised by the exposure of data
 * members here, and no guarantees that they won't change from release to
 * relase of the SimTK Core; if the definition should change at some point
 * the library code will change but it will be using the updated definition
 * and does not have to coordinate in any way with client code.
 *
 * Other users of this class should be aware that if you include it in 
 * code you expose to your own users you may create binary compatibility 
 * problems for yourself. Better to restrict use of this class (and indeed
 * inclusion of this header file) to your private ".cpp" source code and
 * not in your API header files.
 *
 * The PIMPLImplementation base class keeps track of how many Handles
 * are referencing it, so that that the last handle to be deleted can
 * delete the implementation. One handle is designated as the "owner"
 * handle of this implementation. We keep a pointer to that handle here,
 * so special handling is required if the owner handle is deleted while
 * other references still exist.
 */
template <class HANDLE, class IMPL> 
class PIMPLImplementation {
    HANDLE*     ownerHandle;
    mutable int handleCount; // ref count determining when this is destructed 
public:
    /// This serves as a default constructor and as a way to construct
    /// an implementation class which already knows its owner handle.
    /// If the handle is supplied then the handle count is set to one.
    /// If not (default constructor) owner handle is null and the 
    /// handle count at 0.
    explicit PIMPLImplementation(HANDLE* h=0) 
      : ownerHandle(h), handleCount(h ? 1 : 0)
    {
    }

    /// Get the number of handles known to be referencing this implementation.
    int getHandleCount() const {return handleCount;}

    /// Register that a new handle is referencing this implementation so we
    /// won't delete the implementation prematurely.
    void incrementHandleCount() const {handleCount++;}

    /// Register the fact that one of the previously-referencing handles no
    /// longer references this implementation. The remaining number of references
    /// is returned; if it is zero the caller should delete the implementation.
    int decrementHandleCount() const {assert(handleCount>=1); return handleCount--;}

    /// Note that the base class destructor is non-virtual, although it is
    /// expected that derived classes will be abstract. Be sure to provide
    /// a virtual destructor in any abstract class which is derived from 
    /// this base, and be sure to delete a pointer to the abstract class
    /// <em>not</em> a pointer to this base class!
    ~PIMPLImplementation() {assert(handleCount==0); ownerHandle=0;}

    /// The copy constructor for the base class makes sure that the
    /// new object has a null owner handle. A derived class must set
    /// the appropriate owner handle after this is called, that is, in 
    /// the <em>body</em> (not the initializer list) of the derived
    /// class's copy constructor. Also the caller must make sure to
    /// increment the handle count.
    PIMPLImplementation(const PIMPLImplementation&) : ownerHandle(0), handleCount(0) { }

    /// Copy assignment for the base class just makes sure that the 
    /// owner handle is not copied, and that the handle count is zero
    /// for the copy. Caller is required to register a handle and increment
    /// the handle counter.
    PIMPLImplementation& operator=(const PIMPLImplementation& src) {
        if (&src != this)
            ownerHandle=0, handleCount=0;
        return *this;
    }

    /// Provide an owner handle for an implementation which currently does
    /// not have one. This can't be used to <em>replace</em> the owner handle.
    /// This will increment the handle count also.
    void setOwnerHandle(HANDLE& p) {
        assert(!hasOwnerHandle()); 
        ownerHandle=&p;
        incrementHandleCount();
    }

    /// Remove the owner reference from an implementation that currently has
    /// an owner. This decrements the handle count also. The number of remaining
    /// handles is returned.
    int removeOwnerHandle() {
        assert(hasOwnerHandle());
        ownerHandle=0;
        return decrementHandleCount();
    }

    /// Replace the current owner handle with another one. This can't be used to
    /// set the initial owner handle; just to replace an existing one with a
    /// new one. The handle count is not changed here.
    void replaceOwnerHandle(HANDLE& p) {
        assert(hasOwnerHandle()); 
        ownerHandle=&p;
    }

    /// Check whether this implementation currently has a reference to its
    /// owner handle.
    bool hasOwnerHandle() const {return ownerHandle != 0;}

    /// Check whether a given Handle of the appropriate type is the owner of
    /// this implementation.
    bool isOwnerHandle(const HANDLE& p) const {
        return hasOwnerHandle() && ownerHandle==&p;
    }

    /// Return a reference to the owner handle of this implementation. This will
    /// throw an exception if there is no owner handle currently known to this
    /// implementation.
    const HANDLE& getOwnerHandle() const {
        assert(hasOwnerHandle()); 
        return *ownerHandle;
    }
};

    /////////////////////////////
    // PIMPLHandle definitions //
    /////////////////////////////

// serves as default constructor
template <class HANDLE, class IMPL, bool PTR>
inline /*explicit*/ PIMPLHandle<HANDLE,IMPL,PTR>::
PIMPLHandle(IMPL* p) : impl(p) {
    // this bumps the reference count in the implementation
    if (impl) impl->setOwnerHandle(updDowncastToHandle());
}

// destructor
template <class HANDLE, class IMPL, bool PTR>
inline PIMPLHandle<HANDLE,IMPL,PTR>::~PIMPLHandle() {
    // reduces the implementation reference count and deletes it if it hits 0
    clearHandle();
} 

// copy constructor
template <class HANDLE, class IMPL, bool PTR>
PIMPLHandle<HANDLE,IMPL,PTR>::PIMPLHandle(const PIMPLHandle& src) : impl(0) {
    if (PTR) referenceAssign(src.downcastToHandle());
    else     copyAssign(src.downcastToHandle());
}                   

// copy assignment
template <class HANDLE, class IMPL, bool PTR>
PIMPLHandle<HANDLE,IMPL,PTR>& PIMPLHandle<HANDLE,IMPL,PTR>::
operator=(const PIMPLHandle& src) {
    if (PTR) referenceAssign(src.downcastToHandle());
    else     copyAssign(src.downcastToHandle());
    return *this;
} 

template <class HANDLE, class IMPL, bool PTR>
inline bool PIMPLHandle<HANDLE,IMPL,PTR>::isOwnerHandle() const {
    return impl && impl->hasOwnerHandle() &&
        static_cast<const PIMPLHandle*>(&impl->getOwnerHandle()) == this;
}

template <class HANDLE, class IMPL, bool PTR>
inline bool PIMPLHandle<HANDLE,IMPL,PTR>::isSameHandle(const HANDLE& other) const {
    return static_cast<const PIMPLHandle*>(&other) == this;
}


template <class HANDLE, class IMPL, bool PTR>
inline bool PIMPLHandle<HANDLE,IMPL,PTR>::hasSameImplementation(const HANDLE& other) const {
    return impl && (impl==other.impl);
}

// The current (this) handle is an owner. Here it transfers ownership to the supplied
// new empty handle, while retaining a reference to the implementation.
template <class HANDLE, class IMPL, bool PTR>
void PIMPLHandle<HANDLE,IMPL,PTR>::
disown(HANDLE& newOwner) {
    assert(!isSameHandle(newOwner));
    assert(!this->isEmptyHandle() && newOwner.isEmptyHandle());
    newOwner.impl = impl;
    impl->replaceOwnerHandle(newOwner);
    // since the old handle retains a reference, there is now one more handle
    impl->incrementHandleCount();
}

// Reference assignment:
//   - if target (this) is an owner handle, throw an exception; we don't allow that
//   - if source and target have same implementation, there is nothing to do
//   - otherwise, clear the handle, then set implementation and bump handle count
template <class HANDLE, class IMPL, bool PTR>
inline PIMPLHandle<HANDLE,IMPL,PTR>& PIMPLHandle<HANDLE,IMPL,PTR>::
referenceAssign(const HANDLE& src) {
    assert(!isOwnerHandle()); // owner can't be target of a reference assign
    if (!hasSameImplementation(src)) {
        clearHandle();
        impl = src.impl;
        if (impl)
            impl->incrementHandleCount();
    }
    return *this;
}

// Copy assignment:
//  - if same handle, nothing to do
//  - clear this handle, decrementing ref count and deleting implementation if necessary
//  - clone the source implementation, then reference the copy in this target handle
template <class HANDLE, class IMPL, bool PTR>
PIMPLHandle<HANDLE,IMPL,PTR>& PIMPLHandle<HANDLE,IMPL,PTR>::
copyAssign(const HANDLE& src) {
    if (isSameHandle(src)) return *this; // that was easy!
    clearHandle();
    if (src.impl) {
        impl = src.impl->clone(); // NOTE: instantiation requires definition of IMPL class
        impl->setOwnerHandle(updDowncastToHandle()); // bumps ref count (to 1)
        assert(impl->getHandleCount() == 1);
    }
    return *this;
}

// Provide an implementation for this empty handle, bumping the handle count.
// We do not assume this handle is the owner of the implementation; the caller
// must handle that separately.
template <class HANDLE, class IMPL, bool PTR>
inline void PIMPLHandle<HANDLE,IMPL,PTR>::
setImpl(IMPL* p){
    assert(isEmptyHandle());
    impl=p;
    impl->incrementHandleCount();
}

// Remove this handle from its current implementation (if any). If this was the
// owner handle, we clear the owner reference in the implementation. We decrement
// the implementation's handle count and delete the implementation if this
// was the last handle referencing it. 
template <class HANDLE, class IMPL, bool PTR>
void PIMPLHandle<HANDLE,IMPL,PTR>::
clearHandle() {
    if (isEmptyHandle()) return; // handle is already clear
    const int nHandlesLeft = 
        isOwnerHandle() ? impl->removeOwnerHandle() 
                        : impl->decrementHandleCount();
    if (nHandlesLeft == 0)
        delete impl;
    impl=0;
}

template <class HANDLE, class IMPL, bool PTR>
int PIMPLHandle<HANDLE,IMPL,PTR>::
getImplHandleCount() const {
    assert(!isEmptyHandle());
    return impl->getHandleCount(); 
}

    ////////////////////////////////////
    // PIMPLDerivedHandle definitions //
    ////////////////////////////////////

template <class DERIVED, class DERIVED_IMPL, class PARENT>
inline /*explicit*/ PIMPLDerivedHandle<DERIVED,DERIVED_IMPL,PARENT>::
PIMPLDerivedHandle(DERIVED_IMPL* pimpl) 
  : PARENT(pimpl) 
{ 
}

template <class DERIVED, class DERIVED_IMPL, class PARENT>
inline const DERIVED_IMPL& PIMPLDerivedHandle<DERIVED,DERIVED_IMPL,PARENT>::
getImpl() const {
    return dynamic_cast<const DERIVED_IMPL&>(PARENT::getImpl());
}

template <class DERIVED, class DERIVED_IMPL, class PARENT>
inline DERIVED_IMPL& PIMPLDerivedHandle<DERIVED,DERIVED_IMPL,PARENT>::
updImpl() {
    return dynamic_cast<DERIVED_IMPL&>(PARENT::updImpl());
}

template <class DERIVED, class DERIVED_IMPL, class PARENT>
inline const PARENT& PIMPLDerivedHandle<DERIVED,DERIVED_IMPL,PARENT>::
upcast() const {
    return static_cast<const PARENT&>(*this);
}

template <class DERIVED, class DERIVED_IMPL, class PARENT>
inline PARENT& PIMPLDerivedHandle<DERIVED,DERIVED_IMPL,PARENT>::
updUpcast() {
    return static_cast<PARENT&>(*this);
}

template <class DERIVED, class DERIVED_IMPL, class PARENT>
/*static*/ bool PIMPLDerivedHandle<DERIVED,DERIVED_IMPL,PARENT>::
isInstanceOf(const PARENT& p) {
    return dynamic_cast<const DERIVED_IMPL*>(&p.getImpl()) != 0;
}

template <class DERIVED, class DERIVED_IMPL, class PARENT>
/*static*/ const DERIVED& PIMPLDerivedHandle<DERIVED,DERIVED_IMPL,PARENT>::
downcast(const PARENT& p) {
    assert(isInstanceOf(p));
    return static_cast<const DERIVED&>(p);
}

template <class DERIVED, class DERIVED_IMPL, class PARENT>
/*static*/ DERIVED& PIMPLDerivedHandle<DERIVED,DERIVED_IMPL,PARENT>::
updDowncast(PARENT& p) {
    assert(isInstanceOf(p));
    return static_cast<DERIVED&>(p);
}

    ////////////////////////////////
    // TEMPLATIZED GLOBAL METHODS //
    ////////////////////////////////

template <class HANDLE, class IMPL, bool PTR>
std::ostream& operator<<(std::ostream& o, const PIMPLHandle<HANDLE,IMPL,PTR>& h) {
    o << "PIMPLHandle<" << typeid(HANDLE).name() << "," << typeid(IMPL).name() << "> @" << &h;
    if (h.isEmptyHandle())
        return o << " is EMPTY." << endl;

    if (h.isOwnerHandle()) o << " is OWNER of";
    else o << " is REFERENCE to";

    return o << " Implementation @" << &h.getImpl() << " (handle count=" << h.getImpl().getHandleCount() << ")" << endl;
}


} // namespace SimTK

#endif // SimTK_PRIVATE_IMPLEMENTATION_DEFS_H_
