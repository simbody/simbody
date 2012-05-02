#ifndef SimTK_SimTKCOMMON_REFERENCE_PTR_H_
#define SimTK_SimTKCOMMON_REFERENCE_PTR_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Simbody: SimTKcommon                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2012 Stanford University and the Authors.           *
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

/** This is a smart pointer that implements "cross reference" semantics where
a pointer data member of some object is intended to refer to some target
object in a larger data structure. Judicious use of this container will allow
you to use compiler-generated copy constructors and copy assignment operators
for classes which would otherwise have to implement their own in order to
property initialize these pointer data members, which must not be copied.

The contained pointer is initialized to null (0) on construction, and it is 
reinitialized to null upon copy construction or copy assignment. That's 
because we are assuming this is part of copying the entire data structure and 
copying the old pointer would create a reference into the old data structure 
rather than the new copy. This pointer does not own the target to which it 
points, and there is no reference counting so it will become stale if the 
target is deleted.

Whether you can write through the pointer is controlled by whether type T
is a const type. For example %ReferencePtr\<int> is equivalent to an int*, while
%ReferencePtr\<const int> is equivalent to a const int*.

This class is entirely inline and has no computational or space overhead; it
contains just a single pointer. 

@see ClonePtr **/ 
template <class T> class ReferencePtr {
public:
    typedef T  element_type;
    typedef T* pointer;
    typedef T& reference;

    /** Default constructor creates an empty object. **/
	ReferencePtr() : p(0) { }
    /** Construct from a given pointer stores the pointer. **/
    explicit ReferencePtr(T* tp) : p(tp) { }
    /** Construct from a reference stores the address of the supplied 
    object. **/
    explicit ReferencePtr(T& t) : p(&t) { }
    /** Copy constructor unconditionally sets the pointer to null; see class
    comments for why. **/
	ReferencePtr(const ReferencePtr&) : p(0) { }
    /** Copy assignment sets the pointer to null (except for a self-assign); 
    see class comments for why.  **/
	ReferencePtr& operator=(const ReferencePtr& r) 
    {   if (&r != this) clear(); return *this; }
    /** This form of assignment replaces the currently-referenced object by a 
    reference to the source object; no destruction occurs. **/	
    ReferencePtr& operator=(T& t)          
    {  reset(&t); return *this; }
    /** This form of assignment replaces the current pointer with the given
    one; no destruction occurs. **/
    ReferencePtr& operator=(T* tp)               
    {   reset(tp); return *this; }
    
    /** Destructor does nothing. **/
    ~ReferencePtr() {clear();} // just being tidy

    /** Return the contained pointer. This will fail if the container is 
    empty. **/
    T* operator->() const { return &getRef(); }

    /** The "dereference" operator returns a reference to the target object. 
    This will fail if the container is empty. **/
    T& operator*() const { return getRef(); }

    /** This is an implicit conversion from %ReferencePtr\<T> to T*. **/
    operator T*() const { return p; }

    /** This is an implicit conversion to type bool that returns true if
    the container is non-null (that is, not empty). **/
    operator bool() const { return !empty(); }

    /** Return the contained pointer, or null if the container is empty. **/
	T* get()  const  { return p; }

    /** Return a reference to the target object. Fails if the pointer is
    null. **/
	T& getRef() const { 
        SimTK_ERRCHK(p!=0, "ReferencePtr::getRef()", 
                    "An attempt was made to dereference a null pointer."); 
        return *p; 
    }

    /** Return true if this container is empty. **/
	bool     empty() const    { return p==0; }
    /** Make this container empty; no destruction occurs. **/
    void     clear()          { p=0; }
    /** Extract the pointer from this container, leaving the container empty. 
    The pointer is returned. **/
    T*       release()        { T* x=p; p=0; return x; }
    /** Replace the stored pointer with a different one; no destruction
    occurs. **/
    void     reset(T* tp)     { p=tp;}

    /** Swap the contents of this %ReferencePtr with another one. **/
    void swap(ReferencePtr& other) {
        T* otherp = other.release();
        other.reset(p);
        reset(otherp);
    }
	 
private:
    // Warning: ReferencePtr must be exactly the same size as type T*. That way
    // one can reinterpret_cast a T* to a ReferencePtr<T> when needed.
    T*	p;  
};	
	
} // namespace SimTK

namespace std {
/** This is a specialization of the STL std::swap() algorithm which uses the
cheap built-in swap() member of the ReferencePtr class. 
@relates SimTK::ReferencePtr **/
template <class T> inline void
swap(SimTK::ReferencePtr<T>& p1, SimTK::ReferencePtr<T>& p2) {
    p1.swap(p2);
}

} // namespace std

#endif // SimTK_SimTKCOMMON_REFERENCE_PTR_H_
