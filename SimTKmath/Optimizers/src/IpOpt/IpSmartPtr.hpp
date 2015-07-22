// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpSmartPtr.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPSMARTPTR_HPP__
#define __IPSMARTPTR_HPP__

#include "IpReferenced.hpp"

//#define IP_DEBUG_SMARTPTR
#ifdef IP_DEBUG_SMARTPTR
# include "IpDebug.hpp"
#endif

namespace Ipopt
{

  /** Template class for Smart Pointers.
   * A SmartPtr behaves much like a raw pointer, but manages the lifetime
   * of an object, deleting the object automatically. This class implements
   * a reference-counting, intrusive smart pointer design, where all
   * objects pointed to must inherit off of ReferencedObject, which
   * stores the reference count. Although this is intrusive (native types
   * and externally authored classes require wrappers to be referenced
   * by smart pointers), it is a safer design. A more detailed discussion of
   * these issues follows after the usage information.
   *
   * Usage Example:
   * Note: to use the SmartPtr, all objects to which you point MUST
   * inherit off of ReferencedObject.
   *
   * \verbatim
   *
   * In MyClass.hpp...
   *
   * #include "IpReferenced.hpp"

   * namespace Ipopt {
   *
   *  class MyClass : public ReferencedObject // must derive from ReferencedObject
   *    {
   *      ...
   *    }
   * } // namespace Ipopt
   *
   *
   * In my_usage.cpp...
   *
   * #include "IpSmartPtr.hpp"
   * #include "MyClass.hpp"
   *
   * void func(AnyObject& obj)
   *  {
   *    SmartPtr<MyClass> ptr_to_myclass = new MyClass(...);
   *    // ptr_to_myclass now points to a new MyClass,
   *    // and the reference count is 1
   *
   *    ...
   *
   *    obj.SetMyClass(ptr_to_myclass);
   *    // Here, let's assume that AnyObject uses a
   *    // SmartPtr<MyClass> internally here.
   *    // Now, both ptr_to_myclass and the internal
   *    // SmartPtr in obj point to the same MyClass object
   *    // and its reference count is 2.
   *
   *    ...
   *
   *    // No need to delete ptr_to_myclass, this
   *    // will be done automatically when the
   *    // reference count drops to zero.
   *
   *  }
   *
   * \endverbatim
   *
   * It is not necessary to use SmartPtr's in all cases where an
   * object is used that has been allocated "into" a SmartPtr.  It is
   * possible to just pass objects by reference or regular pointers,
   * even if lower down in the stack a SmartPtr is to be held on to.
   * Everything should work fine as long as a pointer created by "new"
   * is immediately passed into a SmartPtr, and if SmartPtr's are used
   * to hold on to objects.
   *
   * Other Notes:
   *  The SmartPtr implements both dereference operators -> & *.
   *  The SmartPtr does NOT implement a conversion operator to
   *    the raw pointer. Use the GetRawPtr() method when this
   *    is necessary. Make sure that the raw pointer is NOT
   *    deleted.
   *  The SmartPtr implements the comparison operators == & !=
   *    for a variety of types. Use these instead of
   *    \verbatim
   *    if (GetRawPtr(smrt_ptr) == ptr) // Don't use this
   *    \endverbatim
   * SmartPtr's, as currently implemented, do NOT handle circular references.
   *    For example: consider a higher level object using SmartPtrs to point to
   *    A and B, but A and B also point to each other (i.e. A has a SmartPtr
   *    to B and B has a SmartPtr to A). In this scenario, when the higher
   *    level object is finished with A and B, their reference counts will
   *    never drop to zero (since they reference each other) and they
   *    will not be deleted. This can be detected by memory leak tools like
   *    valgrind. If the circular reference is necessary, the problem can be
   *    overcome by a number of techniques:
   *
   *    1) A and B can have a method that "releases" each other, that is
   *        they set their internal SmartPtrs to NULL.
   *        \verbatim
   *        void AClass::ReleaseCircularReferences()
   *          {
   *          smart_ptr_to_B = NULL;
   *          }
   *        \endverbatim
   *        Then, the higher level class can call these methods before
   *        it is done using A & B.
   *
   *    2) Raw pointers can be used in A and B to reference each other.
   *        Here, an implicit assumption is made that the lifetime is
   *        controlled by the higher level object and that A and B will
   *        both exist in a controlled manner. Although this seems
   *        dangerous, in many situations, this type of referencing
   *        is very controlled and this is reasonably safe.
   *
   *    3) This SmartPtr class could be redesigned with the Weak/Strong
   *        design concept. Here, the SmartPtr is identified as being
   *        Strong (controls lifetime of the object) or Weak (merely
   *        referencing the object). The Strong SmartPtr increments
   *        (and decrements) the reference count in ReferencedObject
   *        but the Weak SmartPtr does not. In the example above,
   *        the higher level object would have Strong SmartPtrs to
   *        A and B, but A and B would have Weak SmartPtrs to each
   *        other. Then, when the higher level object was done with
   *        A and B, they would be deleted. The Weak SmartPtrs in A
   *        and B would not decrement the reference count and would,
   *        of course, not delete the object. This idea is very similar
   *        to item (2), where it is implied that the sequence of events
   *        is controlled such that A and B will not call anything using
   *        their pointers following the higher level delete (i.e. in
   *        their destructors!). This is somehow safer, however, because
   *        code can be written (however expensive) to perform run-time
   *        detection of this situation. For example, the ReferencedObject
   *        could store pointers to all Weak SmartPtrs that are referencing
   *        it and, in its destructor, tell these pointers that it is
   *        dying. They could then set themselves to NULL, or set an
   *        internal flag to detect usage past this point.
   *
   * Comments on Non-Intrusive Design:
   * In a non-intrusive design, the reference count is stored somewhere other
   * than the object being referenced. This means, unless the reference
   * counting pointer is the first referencer, it must get a pointer to the
   * referenced object from another smart pointer (so it has access to the
   * reference count location). In this non-intrusive design, if we are
   * pointing to an object with a smart pointer (or a number of smart
   * pointers), and we then give another smart pointer the address through
   * a RAW pointer, we will have two independent, AND INCORRECT, reference
   * counts. To avoid this pitfall, we use an intrusive reference counting
   * technique where the reference count is stored in the object being
   * referenced.
   */
  template<class T>
  class SmartPtr : public Referencer
  {
  public:
#define dbg_smartptr_verbosity 0

    /**@name Constructors/Destructors */
    //@{
    /** Default constructor, initialized to NULL */
    SmartPtr();

    /** Copy constructor, initialized from copy */
    SmartPtr(const SmartPtr<T>& copy);

    /** Copy Constructor, initialized from copy of
     * a different type, will cause run-time error
     * if type is not valid.
     */
    /* Don't support this, use dynamic_cast(GetRawPtr(copy))
    template <class U>
    SmartPtr(const SmartPtr<U> &copy);
    */

    /** Constructor, initialized from T* ptr */
    SmartPtr(T* ptr);

    /** Destructor, automatically decrements the
     * reference count, deletes the object if
     * necessary.*/
    ~SmartPtr();
    //@}

    /**@name Overloaded operators. */
    //@{
    /** Overloaded arrow operator, allows the user to call
     * methods using the contained pointer. */
    T* operator->() const;

    /** Overloaded dereference operator, allows the user
     * to dereference the contained pointer. */
    T& operator*() const;

    /** Overloaded equals operator, allows the user to
     * set the value of the SmartPtr from a raw pointer */
    SmartPtr<T>& operator=(T* rhs);

    /** Overloaded equals operator, allows the user to
     * set the value of the SmartPtr from another
     * SmartPtr */
    SmartPtr<T>& operator=(const SmartPtr<T>& rhs);

    /** Overloaded equality comparison operator, allows the
     * user to compare the value of two SmartPtrs */
    template <class U1, class U2>
    friend
    bool operator==(const SmartPtr<U1>& lhs, const SmartPtr<U2>& rhs);

    /** Overloaded equality comparison operator, allows the
     * user to compare the value of a SmartPtr with a raw pointer. */
    template <class U1, class U2>
    friend
    bool operator==(const SmartPtr<U1>& lhs, U2* raw_rhs);

    /** Overloaded equality comparison operator, allows the
     * user to compare the value of a raw pointer with a SmartPtr. */
    template <class U1, class U2>
    friend
    bool operator==(U1* lhs, const SmartPtr<U2>& raw_rhs);

    /** Overloaded in-equality comparison operator, allows the
     * user to compare the value of two SmartPtrs */
    template <class U1, class U2>
    friend
    bool operator!=(const SmartPtr<U1>& lhs, const SmartPtr<U2>& rhs);

    /** Overloaded in-equality comparison operator, allows the
     * user to compare the value of a SmartPtr with a raw pointer. */
    template <class U1, class U2>
    friend
    bool operator!=(const SmartPtr<U1>& lhs, U2* raw_rhs);

    /** Overloaded in-equality comparison operator, allows the
     * user to compare the value of a SmartPtr with a raw pointer. */
    template <class U1, class U2>
    friend
    bool operator!=(U1* lhs, const SmartPtr<U2>& raw_rhs);
    //@}

    /**@name friend method declarations. */
    //@{
    /** Returns the raw pointer contained.
     * Use to get the value of
     * the raw ptr (i.e. to pass to other
     * methods/functions, etc.)
     * Note: This method does NOT copy,
     * therefore, modifications using this
     * value modify the underlying object
     * contained by the SmartPtr,
     * NEVER delete this returned value.
     */
    template <class U>
    friend
    U* GetRawPtr(const SmartPtr<U>& smart_ptr);

    /** Returns a const pointer */
    template <class U>
    friend
    SmartPtr<const U> ConstPtr(const SmartPtr<U>& smart_ptr);

    /** Returns true if the SmartPtr is NOT NULL.
     * Use this to check if the SmartPtr is not null
     * This is preferred to if(GetRawPtr(sp) != NULL)
     */
    template <class U>
    friend
    bool IsValid(const SmartPtr<U>& smart_ptr);

    /** Returns true if the SmartPtr is NULL.
     * Use this to check if the SmartPtr IsNull.
     * This is preferred to if(GetRawPtr(sp) == NULL)
     */
    template <class U>
    friend
    bool IsNull(const SmartPtr<U>& smart_ptr);
    //@}

  private:
    /**@name Private Data/Methods */
    //@{
    /** Actual raw pointer to the object. */
    T* ptr_;

    /** Set the value of the internal raw pointer
     * from another raw pointer, releasing the
     * previously referenced object if necessary. */
    SmartPtr<T>& SetFromRawPtr_(T* rhs);

    /** Set the value of the internal raw pointer
     * from a SmartPtr, releasing the previously referenced
     * object if necessary. */
    SmartPtr<T>& SetFromSmartPtr_(const SmartPtr<T>& rhs);

    /** Release the currently referenced object. */
    void ReleasePointer_();
    //@}
  };

  /**@name SmartPtr friend function declarations.*/
  //@{
  template <class U>
  U* GetRawPtr(const SmartPtr<U>& smart_ptr);

  template <class U>
  SmartPtr<const U> ConstPtr(const SmartPtr<U>& smart_ptr);

  template <class U>
  bool IsNull(const SmartPtr<U>& smart_ptr);

  template <class U>
  bool IsValid(const SmartPtr<U>& smart_ptr);

  template <class U1, class U2>
  bool operator==(const SmartPtr<U1>& lhs, const SmartPtr<U2>& rhs);

  template <class U1, class U2>
  bool operator==(const SmartPtr<U1>& lhs, U2* raw_rhs);

  template <class U1, class U2>
  bool operator==(U1* lhs, const SmartPtr<U2>& raw_rhs);

  template <class U1, class U2>
  bool operator!=(const SmartPtr<U1>& lhs, const SmartPtr<U2>& rhs);

  template <class U1, class U2>
  bool operator!=(const SmartPtr<U1>& lhs, U2* raw_rhs);

  template <class U1, class U2>
  bool operator!=(U1* lhs, const SmartPtr<U2>& raw_rhs);

  //@}


  template <class T>
  SmartPtr<T>::SmartPtr()
      :
      ptr_(NULL)
  {
#ifdef IP_DEBUG_SMARTPTR
    DBG_START_METH("SmartPtr<T>::SmartPtr()", dbg_smartptr_verbosity);
#endif

#ifdef CHECK_SMARTPTR

    const ReferencedObject* trying_to_use_SmartPtr_with_an_object_that_does_not_inherit_from_ReferencedObject_
    = ptr_;
    trying_to_use_SmartPtr_with_an_object_that_does_not_inherit_from_ReferencedObject_ = NULL;
#endif

  }


  template <class T>
  SmartPtr<T>::SmartPtr(const SmartPtr<T>& copy)
      :
      ptr_(NULL)
  {
#ifdef IP_DEBUG_SMARTPTR
    DBG_START_METH("SmartPtr<T>::SmartPtr(const SmartPtr<T>& copy)", dbg_smartptr_verbosity);
#endif

#ifdef CHECK_SMARTPTR

    const ReferencedObject* trying_to_use_SmartPtr_with_an_object_that_does_not_inherit_from_ReferencedObject_
    = ptr_;
    trying_to_use_SmartPtr_with_an_object_that_does_not_inherit_from_ReferencedObject_ = NULL;
#endif

    (void) SetFromSmartPtr_(copy);
  }


  /* Don't support this, use dynamic_cast(GetRawPtr(copy))
  template <class T>
  template <class U> SmartPtr<T>::SmartPtr(const SmartPtr<U>& copy)
   :
   ptr_(NULL)
   {
     DBG_START_METH("SmartPtr<T>::SmartPtr(const SmartPtr<U>& copy)", dbg_smartptr_verbosity);

     U* raw_U_ptr = GetRawPtr(copy);

     // try to cast to type T
     T* raw_T_ptr = dynamic_cast<T*>(raw_U_ptr);

     // Don't try to point to type U with a
     // type T SmartPtr if the types U & T
     // are not related
     DBG_ASSERT(raw_T_ptr);

     (void) SetFromRawPtr_(raw_T_ptr);
   }
   */


  template <class T>
  SmartPtr<T>::SmartPtr(T* ptr)
      :
      ptr_(NULL)
  {
#ifdef IP_DEBUG_SMARTPTR
    DBG_START_METH("SmartPtr<T>::SmartPtr(T* ptr)", dbg_smartptr_verbosity);
#endif

#ifdef CHECK_SMARTPTR

    const ReferencedObject* trying_to_use_SmartPtr_with_an_object_that_does_not_inherit_from_ReferencedObject_
    = ptr_;
    trying_to_use_SmartPtr_with_an_object_that_does_not_inherit_from_ReferencedObject_ = NULL;
#endif

    (void) SetFromRawPtr_(ptr);
  }

  template <class T>
  SmartPtr<T>::~SmartPtr()
  {
#ifdef IP_DEBUG_SMARTPTR
    DBG_START_METH("SmartPtr<T>::~SmartPtr(T* ptr)", dbg_smartptr_verbosity);
#endif

    ReleasePointer_();
  }


  template <class T>
  T* SmartPtr<T>::operator->() const
  {
#ifdef IP_DEBUG_SMARTPTR
    DBG_START_METH("T* SmartPtr<T>::operator->()", dbg_smartptr_verbosity);
#endif

    // cannot deref a null pointer
#ifdef IP_DEBUG

    assert(ptr_);
#endif

    return ptr_;
  }


  template <class T>
  T& SmartPtr<T>::operator*() const
  {
#ifdef IP_DEBUG_SMARTPTR
    DBG_START_METH("T& SmartPtr<T>::operator*()", dbg_smartptr_verbosity);
#endif

    // cannot dereference a null pointer
#ifdef IP_DEBUG

    assert(ptr_);
#endif

    return *ptr_;
  }


  template <class T>
  SmartPtr<T>& SmartPtr<T>::operator=(T* rhs)
  {
#ifdef IP_DEBUG_SMARTPTR
    DBG_START_METH("SmartPtr<T>& SmartPtr<T>::operator=(T* rhs)", dbg_smartptr_verbosity);
#endif

    return SetFromRawPtr_(rhs);
  }


  template <class T>
  SmartPtr<T>& SmartPtr<T>::operator=(const SmartPtr<T>& rhs)
  {
#ifdef IP_DEBUG_SMARTPTR
    DBG_START_METH(
      "SmartPtr<T>& SmartPtr<T>::operator=(const SmartPtr<T>& rhs)",
      dbg_smartptr_verbosity);
#endif

    return SetFromSmartPtr_(rhs);
  }


  template <class T>
  SmartPtr<T>& SmartPtr<T>::SetFromRawPtr_(T* rhs)
  {
#ifdef IP_DEBUG_SMARTPTR
    DBG_START_METH(
      "SmartPtr<T>& SmartPtr<T>::SetFromRawPtr_(T* rhs)", dbg_smartptr_verbosity);
#endif

    // Release any old pointer
    ReleasePointer_();

    if (rhs != NULL) {
      rhs->AddRef(this);
      ptr_ = rhs;
      //       const ReferencedObject* r_ptr =
      //         dynamic_cast<const ReferencedObject*>(rhs);

      //       // All objects pointed to by a SmartPtr MUST
      //       // inherit off of ReferencedObject
      //       DBG_ASSERT(r_ptr && "This is not inherited from ReferencedObject");

      //       if (r_ptr) {
      //         r_ptr->AddRef(this);
      //         ptr_ = rhs;
      //       }
    }

    return *this;
  }

  template <class T>
  SmartPtr<T>& SmartPtr<T>::SetFromSmartPtr_(const SmartPtr<T>& rhs)
  {
#ifdef IP_DEBUG_SMARTPTR
    DBG_START_METH(
      "SmartPtr<T>& SmartPtr<T>::SetFromSmartPtr_(const SmartPtr<T>& rhs)",
      dbg_smartptr_verbosity);
#endif

    T* ptr = GetRawPtr(rhs);
    /* AW: I changed this so that NULL is correctly copied from the
       right hand side */
    //     if (ptr != NULL) {
    //       SetFromRawPtr_(ptr);
    //     }
    SetFromRawPtr_(ptr);

    return (*this);
  }


  template <class T>
  void SmartPtr<T>::ReleasePointer_()
  {
#ifdef IP_DEBUG_SMARTPTR
    DBG_START_METH(
      "void SmartPtr<T>::ReleasePointer()",
      dbg_smartptr_verbosity);
#endif

    if (ptr_) {
      ptr_->ReleaseRef(this);
      if (ptr_->ReferenceCount() == 0) {
        delete ptr_;
      }
      ptr_ = NULL;
    }
    //     if (ptr_ != NULL) {
    //       const ReferencedObject* r_ptr =
    //         dynamic_cast<const ReferencedObject*>(ptr_);

    //       r_ptr->ReleaseRef(this);
    //       if (r_ptr->ReferenceCount() == 0) {
    //         delete ptr_;
    //       }
    //       ptr_ = NULL;
    //     }
  }


  template <class U>
  U* GetRawPtr(const SmartPtr<U>& smart_ptr)
  {
#ifdef IP_DEBUG_SMARTPTR
    DBG_START_FUN(
      "T* GetRawPtr(const SmartPtr<T>& smart_ptr)",
      0);
#endif

    return smart_ptr.ptr_;
  }

  template <class U>
  SmartPtr<const U> ConstPtr(const SmartPtr<U>& smart_ptr)
  {
    // compiler should implicitly cast
    return GetRawPtr(smart_ptr);
  }

  template <class U>
  bool IsValid(const SmartPtr<U>& smart_ptr)
  {
    return !IsNull(smart_ptr);
  }

  template <class U>
  bool IsNull(const SmartPtr<U>& smart_ptr)
  {
#ifdef IP_DEBUG_SMARTPTR
    DBG_START_FUN(
      "bool IsNull(const SmartPtr<T>& smart_ptr)",
      0);
#endif

    return (smart_ptr.ptr_ == NULL);
  }


  template <class U1, class U2>
  bool ComparePointers(const U1* lhs, const U2* rhs)
  {
#ifdef IP_DEBUG_SMARTPTR
    DBG_START_FUN(
      "bool ComparePtrs(const U1* lhs, const U2* rhs)",
      dbg_smartptr_verbosity);
#endif

    if (lhs == rhs) {
      return true;
    }

    // Even if lhs and rhs point to the same object
    // with different interfaces U1 and U2, we cannot guarantee that
    // the value of the pointers will be equivalent. We can
    // guarantee this if we convert to void*
    const void* v_lhs = dynamic_cast<const void*>(lhs);
    const void* v_rhs = dynamic_cast<const void*>(rhs);
    if (v_lhs == v_rhs) {
      return true;
    }

    // They must not be the same
    return false;
  }

  template <class U1, class U2>
  bool operator==(const SmartPtr<U1>& lhs, const SmartPtr<U2>& rhs)
  {
#ifdef IP_DEBUG_SMARTPTR
    DBG_START_FUN(
      "bool operator==(const SmartPtr<U1>& lhs, const SmartPtr<U2>& rhs)",
      dbg_smartptr_verbosity);
#endif

    U1* raw_lhs = GetRawPtr(lhs);
    U2* raw_rhs = GetRawPtr(rhs);
    return ComparePointers(raw_lhs, raw_rhs);
  }

  template <class U1, class U2>
  bool operator==(const SmartPtr<U1>& lhs, U2* raw_rhs)
  {
#ifdef IP_DEBUG_SMARTPTR
    DBG_START_FUN(
      "bool operator==(SmartPtr<U1>& lhs, U2* rhs)",
      dbg_smartptr_verbosity);
#endif

    U1* raw_lhs = GetRawPtr(lhs);
    return ComparePointers(raw_lhs, raw_rhs);
  }

  template <class U1, class U2>
  bool operator==(U1* raw_lhs, const SmartPtr<U2>& rhs)
  {
#ifdef IP_DEBUG_SMARTPTR
    DBG_START_FUN(
      "bool operator==(U1* raw_lhs, SmartPtr<U2>& rhs)",
      dbg_smartptr_verbosity);
#endif

    const U2* raw_rhs = GetRawPtr(rhs);
    return ComparePointers(raw_lhs, raw_rhs);
  }

  template <class U1, class U2>
  bool operator!=(const SmartPtr<U1>& lhs, const SmartPtr<U2>& rhs)
  {
#ifdef IP_DEBUG_SMARTPTR
    DBG_START_FUN(
      "bool operator!=(const SmartPtr<U1>& lhs, const SmartPtr<U2>& rhs)",
      dbg_smartptr_verbosity);
#endif

    bool retValue = operator==(lhs, rhs);
    return !retValue;
  }

  template <class U1, class U2>
  bool operator!=(const SmartPtr<U1>& lhs, U2* raw_rhs)
  {
#ifdef IP_DEBUG_SMARTPTR
    DBG_START_FUN(
      "bool operator!=(SmartPtr<U1>& lhs, U2* rhs)",
      dbg_smartptr_verbosity);
#endif

    bool retValue = operator==(lhs, raw_rhs);
    return !retValue;
  }

  template <class U1, class U2>
  bool operator!=(U1* raw_lhs, const SmartPtr<U2>& rhs)
  {
#ifdef IP_DEBUG_SMARTPTR
    DBG_START_FUN(
      "bool operator!=(U1* raw_lhs, SmartPtr<U2>& rhs)",
      dbg_smartptr_verbosity);
#endif

    bool retValue = operator==(raw_lhs, rhs);
    return !retValue;
  }

} // namespace Ipopt

#endif

