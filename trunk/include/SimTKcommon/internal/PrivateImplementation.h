#ifndef SimTK_PRIVATE_IMPLEMENTATION_H_
#define SimTK_PRIVATE_IMPLEMENTATION_H_

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
 * This header provides declarations of the user-visible portion of the 
 * PIMPLHandle template classes that are used in the SimTK Core to implement
 * the PIMPL (private implementation) design pattern.
 *
 * The definitions associated with these template method declarations are
 * separated into the companion header file PrivateImplementation_Defs.h 
 * with the intent that those definitions will be visible only in library-side
 * code where they need to be instantiated. The definition header file is
 * available for end users as part of the SimTK Core installation, but is
 * not automatically included with SimTKcomon.h as this file is.
 *
 * SimTK Core client-side headers include this declarations file in order to
 * define the various client side Handle classes, but SimTK Core client-side
 * code <em>never</em> includes the definitions; that is done exclusively in
 * private, library-side code.
 */

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/ExceptionMacros.h"

#include <cassert>
#include <iosfwd>

namespace SimTK {


/**
 * This class provides some infrastructure useful in making SimTK Private
 * Implementation (PIMPL) classes. These consist of a "handle" class and
 * an "implementation" class. The handle contains only a single pointer,
 * which points to the implementation class whose definition is unknown
 * to the SimTK client. The implementation class has a pointer back to 
 * *one* of the handles that points to it -- that one is called the "owner
 * handle" and is the only one which will delete the implementation object
 * when the handle is deleted or goes out of scope. All other handles are
 * merely references to the implementation object, and must not be used
 * after the owner handle is deleted.
 *
 * The methods defined below require a definition for the implementation
 * class, so can't be instantiated on the client side. Instead they are
 * instantiated on the library side when needed in the implementation of
 * the PIMPL handle class derived from the PIMPLHandle base.
 *
 * By the time of instantiation, we must have a definition for the IMPL
 * class supplied to the templatized base class. We expect that the IMPL
 * class will be derived from PIMPLImplementation declared below. We 
 * also expect to find certain methods defined, with these names and
 * meanings:
 * 
 *     IMPL* IMPL::clone() const
 *        This creates an implementation object identical to the one
 *        we have, except that its owner handle is set to null. We 
 *        expect the owner handle to be filled in by the derived 
 *        Handle class, which should have initiated the PIMPLHandle
 *        operation which had the need to clone().
 *
 *     const HANDLE& IMPL::getOwnerHandle() const
 *        If the IMPL object does not have an owner, this is expected
 *        to assert(); that would indicate that the derived Handle class
 *        did not properly register itself as the owner upon construction.
 *        Otherwise, this routine returns a reference to the *derived*
 *        Handle class, NOT to the PIMPLHandle parent class! We expect,
 *        however that the Handle class was derived from this PIMPLHandle
 *        so that we can static_cast up here and then compare with 'this'.
 *
 * Usage:
 *    class MySecretImpl;
 *    class MyHandle : public PIMPLHandle<MyHandle,MySecretImpl>
 *
 * There is an optional third template argument, a bool, which can be
 * set true if you want the handle to have pointer semantics rather than
 * the usual object ("value") semantics. Pointer semantics objects
 * have shallow copy constuctor and copy assignment implementations so
 * that they are normally references to objects rather than owners,
 * and pointer semantics owner handles can't be the target of an
 * assignment.
 */
template <class HANDLE, class IMPL, bool PTR=false> 
class PIMPLHandle {
private:
    /// This is the only data member allowed in a handle class. It is 
    /// guaranteed to stay this way for eternity. That is, a well-formed
    /// SimTK Core handle class is just a pointer; and this fact may be
    /// depended upon when necessary.
    IMPL *impl;

public:
    typedef PIMPLHandle<HANDLE, IMPL, PTR> HandleBase;
    typedef HandleBase                     ParentHandle;

    /// Returns true if this handle is empty, that is, does not refer
    /// to any implementation object.
    bool isEmptyHandle() const {return impl==0;}

    /// Returns true if this handle is the owner of the implementation
    /// object to which it refers. An empty handle is <em>not</em> 
    /// considered by this method to be an owner. You can check for an
    /// empty handle using isEmptyHandle().
    /// @see isEmptyHandle()
    bool isOwnerHandle() const;

    /// Determine whether the supplied handle is the same object as
    /// "this" PIMPLHandle.
    bool isSameHandle(const HANDLE& other) const;

    /// Give up ownership of the implementation to an empty handle. The
    /// current handle retains a reference to the implementation but is
    /// no longer its owner. This method requires the current handle to 
    /// be an owner, and the supplied handle to be empty.
    void disown(HANDLE& newOwner);


    /// "Copy" assignment but with shallow (pointer) semantics. As long as
    /// this is not an owner handle already, make it reference the source
    /// implementation. It is not allowed for an <em>owner</em> handle
    /// to be the target of a reference assignment; clear the handle explicitly
    /// first with clearHandle() if you want to do that.
    /// This is the default copy and assignment behavior for pointer semantics handle
    /// classes (that is, those which set the PTR template argument to true).
    /// Caution: although the PIMPLHandle is taken const here, we obtain
    /// a writable pointer to the implementation, meaning that it can
    /// be modified through the reference handle if that handle is non-const.
    /// @see copyAssign()
    /// @see operator=()
    /// @see clearHandle()
    PIMPLHandle& referenceAssign(const HANDLE& source);

    /// This is real copy assignment, with ordinary C++ object ("value") semantics.
    /// Deletes the current implementation if owned; then replaces with a new
    /// copy of the source implementation, of which this handle will be the owner.
    /// This is the default copy and assignment behavior for normal handle objects,
    /// that is, those that let the PTR template argument default or set it to
    /// false explicitly. Use referenceAssign() to make a handle refer to an
    /// existing implementation rather than creating a new copy.
    /// @see referenceAssign()
    PIMPLHandle& copyAssign(const HANDLE& source);

    /// Make this an empty handle, deleting the implementation object if
    /// this handle is the owner of it. A call to isEmptyHandle() will return
    /// true after this.
    /// @see isEmptyHandle()
    void clearHandle();

    /// Get a const reference to the implementation associated with this Handle.
    /// This will throw an exception if there is no implementation.
    const IMPL& getImpl() const {assert(!isEmptyHandle()); return *impl;}

    /// Get a writable reference to the implementation associated with this Handle.
    /// Note that this requires writable access to the handle also. This will
    /// throw an exception if there is no implementation.
    IMPL& updImpl() {assert(!isEmptyHandle()); return *impl;}

    /// Return the number of handles the implementation believes are referencing it.
    /// Throws an exception if there is no implementation.
    /// This is for degugging and consistency checking and shouldn't normally be used.
    int getImplHandleCount() const;

protected:
    /// This serves as the default constructor, which will construct the handle
    /// with an empty implementation, and as a way to construct a handle referencing
    /// an existing implementation object.
    explicit PIMPLHandle(IMPL* p=0);

    /// Note that the destructor is non-virtual. This is a concrete class and so
    /// should be all the handle classes derived from it. If this handle is the 
    /// owner of its implementation, the destructor will destroy the implementation
    /// object as well. Any other handles referencing the same implementation will
    /// then be invalid, although there will be automated detection of that. Be very
    /// careful to ensure that owner handles always outlive their reference handles.
    ~PIMPLHandle();

    /// The copy constructor makes either a deep (value) or shallow (reference) copy
    /// of the supplied source PIMPL object, based on whether this is a "pointer
    /// sematics" (PTR=true) or "object (value) semantics" (PTR=false, default)
    /// class.
    /// @see referenceAssign
    /// @see copyAssign
    PIMPLHandle(const PIMPLHandle& source);

    /// Copy assignment makes the current handle either a deep (value) or shallow
    /// (reference) copy of the supplied source PIMPL object, based on whether this
    /// is a "pointer sematics" (PTR=true) or "object (value) semantics" (PTR=false,
    /// default) class. In the case of a pointer semantics class, an owner handle
    /// can <em>not</em> be the target of an assignment. You can call copyAssign()
    /// directly if you want to make a fresh copy of the source, or you can clear
    /// this handle first with clearHandle() and then use operator=() or referenceAssign()
    /// to turn this handle into a mere reference to the source.
    /// @see referenceAssign
    /// @see copyAssign
    PIMPLHandle& operator=(const PIMPLHandle& source);

    /// Set the implementation for this empty handle. This may result in either
    /// an owner or reference handle, depending on the owner handle reference
    /// stored in the implementation object. This will throw an exception if the
    /// handle is already occupied; it <em>cannot</em> be used to replace one
    /// implementation with another.
    void setImpl(IMPL* p);

    /// Determine whether the supplied handle is a reference to the same implementation
    /// object as is referenced by "this" PIMPLHandle.
    bool hasSameImplementation(const HANDLE& other) const;

private:
    const HANDLE& downcastToHandle() const {return static_cast<const HANDLE&>(*this);}
    HANDLE& updDowncastToHandle() {return static_cast<HANDLE&>(*this);}
};

/**
 * This class is the parent of derived handle classes. We are assuming that the implementation
 * is a C++ abstract class, so we can dynamic_cast the handle's implementation pointer to the
 * indicated concrete implementation type DERIVED_IMPL.
 * 
 * Shallow or deep copy behavior is inherited from the PARENT class. Note that the PARENT
 * may itself be a derived handle from a still higher-level handle class in this hierarchy.
 */
template <class DERIVED, class DERIVED_IMPL, class PARENT>
class PIMPLDerivedHandle : public PARENT {
public:
    /// This is the PIMPLHandle type at the root of this Handle hierarchy.
    typedef typename PARENT::HandleBase HandleBase;


    /// This is the type of this derived Handle object's Parent.
    typedef PARENT ParentHandle;

    /// This serves as the default constructor.
    explicit PIMPLDerivedHandle(DERIVED_IMPL* pimpl=0);

    // default copy constructor, copy assignment, destructor

    /// Cast this derived handle to a const reference to its parent handle.
    /// This conversion is normally done automatically by C++; use this routine
    /// if you want to be explicit about it.
    /// @see updUpcast()
    const PARENT& upcast() const;

    /// Cast this writable derived handle to a writable reference to its parent handle.
    /// This conversion is normally done automatically by C++; use this routine
    /// if you want to be explicit about it.
    /// @see upcast()
    PARENT& updUpcast();

    /// Determine whether an object of this handle's parent type can be safely
    /// downcast to an object of the derived type.
    /// @see downcast()
    static bool isInstanceOf(const PARENT& p);

    /// Downcast a const parent-class object to a const reference to the derived
    /// class object. This will throw an exception of the supplied object is not
    /// actually of derived object type. If you're not sure, use isInstanceOf()
    /// first to check.
    /// @see isInstanceOf()
    /// @see updDowncast()
    static const DERIVED& downcast(const PARENT& p);

    /// Downcast a writable parent-class object to a writable reference to the derived
    /// class object. This will throw an exception of the supplied object is not
    /// actually of derived object type. If you're not sure, use isInstanceOf()
    /// first to check.
    /// @see isInstanceOf()
    /// @see downcast()
    static DERIVED& updDowncast(PARENT& p);

    /// Obtain a const reference to the implementation stored in the HandleBase, but
    /// dynamically down cast to the DERIVED_IMPL class. This will throw a C++
    /// "bad_cast" exception if the implementation object can't be appropriately
    /// downcast. If you aren't sure, use isInstanceOf() first to check.
    /// @see isInstanceOf()
    /// @see updImpl()
    const DERIVED_IMPL& getImpl() const;

    /// Obtain a writable reference to the implementation stored in the HandleBase.
    /// @see getImpl()
    DERIVED_IMPL& updImpl();

    typedef PIMPLDerivedHandle<DERIVED,DERIVED_IMPL,PARENT> PIMPLDerivedHandleBase;
};

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
    explicit PIMPLImplementation(HANDLE* h=0);

    /// Get the number of handles known to be referencing this implementation.
    int getHandleCount() const;

    /// Register that a new handle is referencing this implementation so we
    /// won't delete the implementation prematurely.
    void incrementHandleCount() const;

    /// Register the fact that one of the previously-referencing handles no
    /// longer references this implementation. The remaining number of references
    /// is returned; if it is zero the caller should delete the implementation.
    int decrementHandleCount() const;

    /// Note that the base class destructor is non-virtual, although it is
    /// expected that derived classes will be abstract. Be sure to provide
    /// a virtual destructor in any abstract class which is derived from 
    /// this base, and be sure to delete a pointer to the abstract class
    /// <em>not</em> a pointer to this base class!
    ~PIMPLImplementation();

    /// The copy constructor for the base class makes sure that the
    /// new object has a null owner handle. A derived class must set
    /// the appropriate owner handle after this is called, that is, in 
    /// the <em>body</em> (not the initializer list) of the derived
    /// class's copy constructor. Also the caller must make sure to
    /// increment the handle count.
    PIMPLImplementation(const PIMPLImplementation&);

    /// Copy assignment for the base class just makes sure that the 
    /// owner handle is not copied, and that the handle count is zero
    /// for the copy. Caller is required to register a handle and increment
    /// the handle counter.
    PIMPLImplementation& operator=(const PIMPLImplementation& src);

    /// Provide an owner handle for an implementation which currently does
    /// not have one. This can't be used to <em>replace</em> the owner handle.
    /// This will increment the handle count also.
    void setOwnerHandle(HANDLE& p);

    /// Remove the owner reference from an implementation that currently has
    /// an owner. This decrements the handle count also. The number of remaining
    /// handles is returned.
    int removeOwnerHandle();

    /// Replace the current owner handle with another one. This can't be used to
    /// set the initial owner handle; just to replace an existing one with a
    /// new one. The handle count is not changed here.
    void replaceOwnerHandle(HANDLE& p);

    /// Check whether this implementation currently has a reference to its
    /// owner handle.
    bool hasOwnerHandle() const;

    /// Check whether a given Handle of the appropriate type is the owner of
    /// this implementation.
    bool isOwnerHandle(const HANDLE& p) const;

    /// Return a reference to the owner handle of this implementation. This will
    /// throw an exception if there is no owner handle currently known to this
    /// implementation.
    const HANDLE& getOwnerHandle() const;
};

template <class H, class IMPL, bool PTR>
std::ostream& operator<<(std::ostream& o, const PIMPLHandle<H,IMPL,PTR>& h);

} // namespace SimTK

#endif // SimTK_PRIVATE_IMPLEMENTATION_H_
