#ifndef SimTK_SimTKCOMMON_CONCRETIZE_H_
#define SimTK_SimTKCOMMON_CONCRETIZE_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-6 Stanford University and the Authors.         *
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

#include "SimTKcommon/internal/common.h"

namespace SimTK {

/** 
 * Wrap a pointer to an abstract base class in a way that makes
 * it behave like a concrete class (sometimes called a "ClonePtr").
 * 
 * The abstract base class must cooperate by containing a clone() method
 * which returns a pointer to a copy of the <em>concrete</em> object.
 *
 * The Concretize object normally makes a copy of the object passed to its 
 * constructor. However, if you pass it a non-const pointer  Concretize will 
 * steal it and (if possible) return your pointer null just to be tidy.
 */ 
template <class T> class Concretize {
public:
	Concretize() : p(0) { }
    explicit Concretize(T*  obj) : p(obj) { }           // steals object            
    explicit Concretize(T** obj) : p(*obj) { *obj=0; }  // steal & tidy up
	explicit Concretize(const T* obj) : p(obj?obj->clone():0) { }
	explicit Concretize(const T& obj) : p(&obj?obj.clone():0) { }

	Concretize(const Concretize& c) : p(c.p?c.p->clone():0) { }
	Concretize& operator=(const Concretize& c) { replace(c.p?c.p->clone():0); return *this; }
	Concretize& operator=(const T& t)          { replace(&t ? t.clone()  :0); return *this; }
    Concretize& operator=(T* tp)               { replace(tp); return *this; }
       
    ~Concretize() { delete p; }

    bool operator==(const Concretize& c) const {return getRef()==c.getRef();}
    bool operator!=(const Concretize& c) const {return !((*this)==c);}
    
    const T* operator->() const { return p; }
    T*       operator->()       { return p; }
    
    /// implicit conversions       
    operator const T&() const { return *p; }
    operator T&()             { return *p; }    
       
	T&       updRef()         { return *p; }
	const T& getRef()  const  { return *p; }	
	bool     isEmpty() const  { return p==0; }
    void     clear()          { replace(0); }
    T*       extract()        { T* x=p; p=0; return x; }
    void     replace(T* tp)   { delete p; p=tp; }
    void     replace(T** tpp) { delete p; p=*tpp; *tpp=0; }
	 
private:
    // Warning: Concretize must be exactly the same size as type T*. That way
    // one can reinterpret_cast a T* to a Concretize<T> when needed.
    T*	p;  
};	
	
} // namespace SimTK
#endif // SimTK_SimTKCOMMON_CONCRETIZE_H_
