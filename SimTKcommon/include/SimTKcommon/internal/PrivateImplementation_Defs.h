#ifndef SimTK_PRIVATE_IMPLEMENTATION_DEFS_H_
#define SimTK_PRIVATE_IMPLEMENTATION_DEFS_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors: Christopher Bruns, Peter Eastman                             *
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

    /////////////////////////////////////
    // PIMPLImplementation definitions //
    /////////////////////////////////////

template <class HANDLE, class IMPL>
PIMPLImplementation<HANDLE, IMPL>::PIMPLImplementation(HANDLE* h) : ownerHandle(h), handleCount(h ? 1 : 0) {
}

template <class HANDLE, class IMPL>
int PIMPLImplementation<HANDLE, IMPL>::getHandleCount() const {
    return handleCount;
}

template <class HANDLE, class IMPL>
void PIMPLImplementation<HANDLE, IMPL>::incrementHandleCount() const {
    handleCount++;
}

template <class HANDLE, class IMPL>
int PIMPLImplementation<HANDLE, IMPL>::decrementHandleCount() const {
    assert(handleCount>=1); return --handleCount;
}

template <class HANDLE, class IMPL>
PIMPLImplementation<HANDLE, IMPL>::~PIMPLImplementation() {
    assert(handleCount==0); ownerHandle=0;
}

template <class HANDLE, class IMPL>
PIMPLImplementation<HANDLE, IMPL>::PIMPLImplementation(const PIMPLImplementation&) : ownerHandle(0), handleCount(0) {
}

template <class HANDLE, class IMPL>
PIMPLImplementation<HANDLE, IMPL>& PIMPLImplementation<HANDLE, IMPL>::operator=(const PIMPLImplementation& src) {
    if (&src != this)
        ownerHandle=0, handleCount=0;
    return *this;
}

template <class HANDLE, class IMPL>
void PIMPLImplementation<HANDLE, IMPL>::setOwnerHandle(HANDLE& p) {
    assert(!hasOwnerHandle()); 
    ownerHandle=&p;
    incrementHandleCount();
}

template <class HANDLE, class IMPL>
int PIMPLImplementation<HANDLE, IMPL>::removeOwnerHandle() {
    assert(hasOwnerHandle());
    ownerHandle=0;
    return decrementHandleCount();
}

template <class HANDLE, class IMPL>
void PIMPLImplementation<HANDLE, IMPL>::replaceOwnerHandle(HANDLE& p) {
    assert(hasOwnerHandle()); 
    ownerHandle=&p;
}

template <class HANDLE, class IMPL>
bool PIMPLImplementation<HANDLE, IMPL>::hasOwnerHandle() const {
    return ownerHandle != 0;
}

template <class HANDLE, class IMPL>
bool PIMPLImplementation<HANDLE, IMPL>::isOwnerHandle(const HANDLE& p) const {
    return hasOwnerHandle() && ownerHandle==&p;
}

template <class HANDLE, class IMPL>
const HANDLE& PIMPLImplementation<HANDLE, IMPL>::getOwnerHandle() const {
    assert(hasOwnerHandle()); 
    return *ownerHandle;
}

    /////////////////////////////
    // PIMPLHandle definitions //
    /////////////////////////////

template <class HANDLE, class IMPL, bool PTR>
/*explicit*/ PIMPLHandle<HANDLE,IMPL,PTR>::
PIMPLHandle(IMPL* p) : impl(p) {
    // this bumps the reference count in the implementation
    if (impl) impl->setOwnerHandle(updDowncastToHandle());
}

// destructor
template <class HANDLE, class IMPL, bool PTR>
PIMPLHandle<HANDLE,IMPL,PTR>::~PIMPLHandle() {
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
bool PIMPLHandle<HANDLE,IMPL,PTR>::isOwnerHandle() const {
    return impl && impl->hasOwnerHandle() &&
        static_cast<const PIMPLHandle*>(&impl->getOwnerHandle()) == this;
}

template <class HANDLE, class IMPL, bool PTR>
bool PIMPLHandle<HANDLE,IMPL,PTR>::isSameHandle(const HANDLE& other) const {
    return static_cast<const PIMPLHandle*>(&other) == this;
}


template <class HANDLE, class IMPL, bool PTR>
bool PIMPLHandle<HANDLE,IMPL,PTR>::hasSameImplementation(const HANDLE& other) const {
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
PIMPLHandle<HANDLE,IMPL,PTR>& PIMPLHandle<HANDLE,IMPL,PTR>::
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
void PIMPLHandle<HANDLE,IMPL,PTR>::
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

    ////////////////////////////////
    // TEMPLATIZED GLOBAL METHODS //
    ////////////////////////////////

template <class HANDLE, class IMPL, bool PTR>
std::ostream& operator<<(std::ostream& o, const PIMPLHandle<HANDLE,IMPL,PTR>& h) {
    o << "PIMPLHandle<" << typeid(HANDLE).name() << "," << typeid(IMPL).name() << "> @" << &h;
    if (h.isEmptyHandle())
        return o << " is EMPTY." << std::endl;

    if (h.isOwnerHandle()) o << " is OWNER of";
    else o << " is REFERENCE to";

    return o << " Implementation @" << &h.getImpl() << " (handle count=" << h.getImpl().getHandleCount() << ")" << std::endl;
}


} // namespace SimTK

#endif // SimTK_PRIVATE_IMPLEMENTATION_DEFS_H_
