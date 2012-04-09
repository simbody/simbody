#ifndef SimTK_SimTKCOMMON_CLONE_PTR_H_
#define SimTK_SimTKCOMMON_CLONE_PTR_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Simbody: SimTKcommon                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-12 Stanford University and the Authors.        *
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

namespace SimTK {

/** Wrap a pointer to an abstract base class in a way that makes it behave like
a concrete class. This is similar to std::unique_ptr in that it does not permit
shared ownership of the object. However, unlike std::unique_ptr, %ClonePtr
supports copy and assigment operations, by insisting that the contained object
have a clone() method that returns a pointer to a heap-allocated deep copy of 
the <em>concrete</em> object.

We define operator==() and operator<() here that delegate to the contained 
object; if you want to use those the contained type must support that
operator.

This class is entirely inline and has no computational or space overhead; it
contains just a single pointer and does no reference counting. **/ 
template <class T> class ClonePtr {
public:
    typedef T  element_type;
    typedef T* pointer;
    typedef T& reference;

    /** Default constructor creates an empty object. **/
	ClonePtr() : p(0) { }
    /** Given a pointer to a writable heap-allocated object, take over 
    ownership of that object. **/
    explicit ClonePtr(T*  obj) : p(obj) { }
    /** Given a pointer to a writable heap-allocated object, take over 
    ownership of that object and set the original pointer to null. **/
    explicit ClonePtr(T** obj) : p(*obj) { *obj=0; }
    /** Given a pointer to a read-only object, create a new heap-allocated 
    copy of that object via its clone() method and make this %ClonePtr
    object the owner of the copy. Ownership of the original object is not
    affected. If the supplied pointer is null, the resulting %ClonePtr
    object is empty. **/
	explicit ClonePtr(const T* obj) : p(obj?obj->clone():0) { }
    /** Given a read-only reference to an object, create a new heap-allocated 
    copy of that object via its clone() method and make this %ClonePtr
    object the owner of the copy. Ownership of the original object is not
    affected. **/
    explicit ClonePtr(const T& obj) : p(&obj?obj.clone():0) { }
    /** Copy constructor is deep; the new %ClonePtr object contains a new
    copy of the object in the source, created via the source object's clone()
    method. If the source container is empty this one will be empty also. **/
	ClonePtr(const ClonePtr& c) : p(c.p?c.p->clone():0) { }
    /** Copy assignment replaces the currently-held object by a heap-allocated
    copy of the object held in the source container. The copy is created using 
    the source object's clone() method. The currently-held object is deleted. 
    If the source container is empty this one will be empty after the 
    assignment. **/
	ClonePtr& operator=(const ClonePtr& c) 
    {   reset(c.p?c.p->clone():0); return *this; }
    /** This form of assignment replaces the currently-held object by a 
    heap-allocated copy of the source object. The copy is created using the 
    source object's clone() method. The currently-held object is deleted. **/	
    ClonePtr& operator=(const T& t)          
    {   reset(&t ? t.clone()  :0); return *this; }
    /** This form of assignment replaces the currently-held object by the given
    source object and takes over ownership of the source object. The 
    currently-held object is deleted. **/ 
    ClonePtr& operator=(T* tp)               
    {   reset(tp); return *this; }
    
    /** Destructor deletes the referenced object. **/
    ~ClonePtr() { delete p; }

    /** Compare the contained objects for equality using the contained 
    objects' operator==() operator. If both containers are empty or both
    refer to the same object they will test equal; if only one is empty they 
    will test not equal. Otherwise the objects are compared. **/
    bool operator==(const ClonePtr& other) const {
        if (p == other.p) return true; // same object or both empty
        if (empty() || other.empty()) return false;
        return getRef()==other.getRef();
    }
    /** Compare the contained objects for inequality using operator==(). **/
    bool operator!=(const ClonePtr& other) const {return !((*this)==other);}

    /** Provide an ordering for use in sorted containers using the contained 
    objects' operator<(). If both containers are empty or both
    refer to the same object they will test equal; if only one is empty that 
    one is considered less than the other one. Otherwise the objects are 
    compared using T::operator<(). **/
    bool operator<(const ClonePtr& other) const {
        if (p == other.p)  return false;    // same object or both empty
        if (empty())       return true;     //  empty < !empty
        if (other.empty()) return false;    // !empty > empty
        return getRef() < other.getRef();
    }
    
    /** Dereference a const pointer to the contained object. This will fail if 
    the container is empty. **/
    const T* operator->() const { return &getRef(); }
    /** Dereference a writable pointer to the contained object. This will fail 
    if the container is empty. **/
    T*       operator->()       { return &updRef(); }

    /** This "dereference" operator returns a const reference to the contained 
    object. This will fail if the container is empty. **/
    const T& operator*() const { return getRef(); }
    /** This "dereference" operator returns a writable reference to the 
    contained object. This will fail if the container is empty. **/
    T&       operator*()       { return updRef(); }

    /** This "address of" operator returns a const pointer to the contained 
    object (or null if none). **/
    const T* operator&() const { return p; }
    /** This "address of" operator returns a writable pointer to the contained 
    object (or null if none). **/
    T*       operator&()       { return p; }
    
    /** This is an implicit conversion from %ClonePtr\<T> to a const 
    reference to the contained object. This will fail if the container is
    empty. **/
    operator const T&() const { return getRef(); }
    /** This is an implicit conversion from %ClonePtr\<T> to a writable 
    reference to the contained object. **/
    operator T&()             { return updRef(); } 

    /** Return a writable pointer to the contained object if any, or null. 
    You can use the "address of" operator\&() instead if you prefer. **/
	T* updPtr() { return *p; }
    /** Return a const pointer to the contained object if any, or null.  
    You can use the "address of" operator\&() instead if you prefer. **/
	const T* getPtr()  const  { return *p; }

    /** Return a writable reference to the contained object. Don't call this
    this container is empty. There is also an implicit conversion to reference
    that allows %ClonePtr\<T> to be used as though it were a T\&. **/
	T& updRef() { 
        SimTK_ERRCHK(p!=0, "ClonePtr::updRef()", 
                    "An attempt was made to dereference a null pointer."); 
        return *p; 
    }
    /** Return a const reference to the contained object. Don't call this if
    this container is empty. There is also an implicit conversion to reference
    that allows %ClonePtr\<T> to be used as though it were a T\&. **/
	const T& getRef() const { 
        SimTK_ERRCHK(p!=0, "ClonePtr::getRef()", 
                    "An attempt was made to dereference a null pointer."); 
        return *p; 
    }	

    /** Return true if this container is empty. **/
	bool     empty() const    { return p==0; }
    /** Make this container empty, deleting the currently contained object if
    there is one. **/
    void     clear()          { reset(0); }
    /** Extract the object from this container, leaving the container empty
    and transferring ownership to the caller. A pointer to the object is
    returned. **/
    T*       release()        { T* x=p; p=0; return x; }
    /** Replace the contents of this container with the supplied heap-allocated
    object, taking over ownership of that object and deleting the current one
    first if necessary. Nothing happens if the supplied pointer is the same
    as the one already stored. **/
    void     reset(T* tp)     { if (tp!=p) {delete p; p=tp;} }
    /** Replace the contents of this container with the supplied heap-allocated
    object, taking over ownership of that object, deleting the current one
    first if necessary, and setting the caller's supplied pointer to null. 
    If the supplied pointer is the same as the one we're already storing, then
    we just set the caller's pointer to null and return. **/
    void reset(T** tpp) { 
        if (*tpp!=p) {delete p; p=*tpp;} 
        *tpp=0; 
    }
    /** Swap the contents of this %ClonePtr with another one, with ownership
    changing hands but no copying performed. **/
    void swap(ClonePtr& other) {
        T* otherp = other.release();
        other.reset(p);
        reset(otherp);
    }
	 
private:
    // Warning: ClonePtr must be exactly the same size as type T*. That way
    // one can reinterpret_cast a T* to a ClonePtr<T> when needed.
    T*	p;  
};	
	
} // namespace SimTK

namespace std {
/** This is a specialization of the STL std::swap() algorithm which uses the
cheap built-in swap() member of the ClonePtr class. 
@relates SimTK::ClonePtr **/
template <class T> inline void
swap(SimTK::ClonePtr<T>& p1, SimTK::ClonePtr<T>& p2) {
    p1.swap(p2);
}

} // namespace std

#endif // SimTK_SimTKCOMMON_CONCRETIZE_H_
