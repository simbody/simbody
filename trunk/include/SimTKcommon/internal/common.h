#ifndef SimTK_SimTKCOMMON_COMMON_H_
#define SimTK_SimTKCOMMON_COMMON_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-7 Stanford University and the Authors.         *
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

/**@file
 * Mandatory first inclusion for any SimTK source or header file.
 * 
 * Every source and most header files using SimTK must include this 
 * header as its <em>first</em> inclusion. Declarations and definitions that 
 * must be available and compiler-and machine-specific issues are dealt
 * with here.
 *
 * This file must be includable from either C++ or ANSI C. It uses
 * the ANSI-C++ macro "__cplusplus" for any code that will compile
 * only under C++.
 */

/*****************************/
/* ANSI-C COMPATIBLE SECTION */
/*****************************/

/* Set up a few compile-time options that affect all SimTK Core headers. */

/* 
 * This compile-time constant determines the default precision used everywhere
 * in SimTK Core code. Wherever a SimTK::Real, SimTK::Vector, SimTK::Matrix,
 * etc. appears with no precision specified, it will have this underlying precision.
 * We use 1==float, 2==double, 4==long double. Any other value will cause
 * a compile time error. The default is 2.
 */
#ifndef SimTK_DEFAULT_PRECISION
#   define SimTK_DEFAULT_PRECISION 2
#endif

/* This type is for use in C. In C++ use SimTK::Real instead. */
#if   (SimTK_DEFAULT_PRECISION == 1)
    typedef float SimTK_Real;
#elif (SimTK_DEFAULT_PRECISION == 2)
    typedef double SimTK_Real;
#elif (SimTK_DEFAULT_PRECISION == 4)
    typedef long double SimTK_Real;
#else
    #error ILLEGAL VALUE FOR DEFAULT PRECISION
#endif

#ifndef NDEBUG
    #if defined(__cplusplus)
        #include <cstdio>
        #define SimTK_DEBUG(s) std::printf("DBG: " s)
        #define SimTK_DEBUG1(s,a1) std::printf("DBG: " s,a1)	
        #define SimTK_DEBUG2(s,a1,a2) std::printf("DBG: " s,a1,a2)	
        #define SimTK_DEBUG3(s,a1,a2,a3) std::printf("DBG: " s,a1,a2,a3)	
        #define SimTK_DEBUG4(s,a1,a2,a3,a4) std::printf("DBG: " s,a1,a2,a3,a4)
    #else
        #include <stdio.h>
        #define SimTK_DEBUG(s) printf("DBG: " s)
        #define SimTK_DEBUG1(s,a1) printf("DBG: " s,a1)	
        #define SimTK_DEBUG2(s,a1,a2) printf("DBG: " s,a1,a2)	
        #define SimTK_DEBUG3(s,a1,a2,a3) printf("DBG: " s,a1,a2,a3)	
        #define SimTK_DEBUG4(s,a1,a2,a3,a4) printf("DBG: " s,a1,a2,a3,a4)
    #endif
#else
    #define SimTK_DEBUG(s)
    #define SimTK_DEBUG1(s,a1)
    #define SimTK_DEBUG2(s,a1,a2)
    #define SimTK_DEBUG3(s,a1,a2,a3)	
    #define SimTK_DEBUG4(s,a1,a2,a3,a4)
#endif

/*
 * Shared libraries are messy in Visual Studio. We have to distinguish three
 * cases:
 *   (1) this header is being used to build the SimTKcommon shared library (dllexport)
 *   (2) this header is being used by a *client* of the SimTKcommon shared
 *       library (dllimport)
 *   (3) we are building the SimTKcommon static library, or the client is
 *       being compiled with the expectation of linking with the
 *       SimTKcommon static library (nothing special needed)
 * In the CMake script for building this library, we define one of the symbols
 *     SimTK_SimTKCOMMON_BUILDING_{SHARED|STATIC}_LIBRARY
 * Client code normally has no special symbol defined, in which case we'll
 * assume it wants to use the shared library. However, if the client defines
 * the symbol SimTK_USE_STATIC_LIBRARIES we'll suppress the dllimport so
 * that the client code can be linked with static libraries. Note that
 * the client symbol is not library dependent, while the library symbols
 * affect only the SimTKcommon library, meaning that other libraries can
 * be clients of this one. However, we are assuming all-static or all-shared.
 */

#ifdef WIN32
    #if defined(SimTK_SimTKCOMMON_BUILDING_SHARED_LIBRARY)
        #define SimTK_SimTKCOMMON_EXPORT __declspec(dllexport)
    #elif defined(SimTK_SimTKCOMMON_BUILDING_STATIC_LIBRARY) || defined(SimTK_USE_STATIC_LIBRARIES)
        #define SimTK_SimTKCOMMON_EXPORT
    #else
        #define SimTK_SimTKCOMMON_EXPORT __declspec(dllimport)   // i.e., a client of a shared library
    #endif
#else
    #define SimTK_SimTKCOMMON_EXPORT // Linux, Mac
#endif

/* Every SimTK Core library must provide these two routines, with the library
 * name appearing after the "version_" and "about_".
 */
#if defined(__cplusplus)
extern "C" {
#endif
    SimTK_SimTKCOMMON_EXPORT void SimTK_version_SimTKcommon(int* major, int* minor, int* build);
    SimTK_SimTKCOMMON_EXPORT void SimTK_about_SimTKcommon(const char* key, int maxlen, char* value);
#if defined(__cplusplus)
}
#endif

/************************************/
/* END OF ANSI-C COMPATIBLE SECTION */
/************************************/

#if defined(__cplusplus)

#include <cstddef>
#include <cassert>
#include <complex>
#include <limits>

/**
 * Use this macro to define a unique "Id" type which is just a type-safe
 * non-negative int, augmented with a "NaN" value given by the predefined
 * int constant SimTK::InvalidId. No namespace is assumed for the newly-
 * defined type; if you want the symbol in a namespace be sure to invoke
 * the macro within that namespace. Make sure that the system include
 * file <cassert> has been included by the point of invocation, because
 * the define Id type uses the assert() macro when in Debug mode.
 *
 * For most uses it will behave like an int, and it has an implicit
 * conversion *to* int. Importantly though, it has no implicit conversion
 * *from* int so you can't pass some other kind of number where a particular
 * kind of Id was expected. This is used to create Id types
 * which can be used as array indices but which prevent accidental mixing
 * of types. Examples: SubsystemId, ConstraintId.
 *
 * If you create a type "ThingId" you will also get a constant of
 * type ThingId named "ThingInvalidId" which will be the initial
 * value of any objects of type ThingId, and will have the same numerical
 * value as SimTK::InvalidId.
 */
namespace SimTK {
    static const int InvalidId = -1111111111;
}

#define SimTK_DEFINE_UNIQUE_ID_TYPE(NAME)   \
class NAME {                                \
    int id;                                 \
public:                                     \
    inline NAME();                          \
    inline explicit NAME(int i);            \
    inline explicit NAME(long i);           \
    inline explicit NAME(unsigned int u);   \
    inline explicit NAME(unsigned long u);  \
    operator int() const {return id;}       \
    bool isValid() const {return id>=0;}    \
    void invalidate();                      \
    const NAME& operator++() {assert(id>=0); ++id;return *this;}      /*prefix */   \
    NAME operator++(int)     {assert(id>=0); ++id; return NAME(id-1);}/*postfix*/   \
    const NAME& operator--() {assert(id>=1); --id;return *this;}      /*prefix */   \
    NAME operator--(int)     {assert(id>=1); --id; return NAME(id+1);}/*postfix*/   \
};                                                      \
static const NAME Invalid ## NAME(SimTK::InvalidId);    \
inline NAME::NAME() : id(Invalid ## NAME) { }           \
inline void NAME::invalidate() {id=Invalid ## NAME;}    \
inline NAME::NAME(unsigned int u) : id((int)u) {        \
    assert((int)u >= 0);                                \
}                                                       \
inline NAME::NAME(unsigned long u) : id((int)u) {       \
    assert((int)u >= 0);                                \
}                                                       \
inline NAME::NAME(int i) : id(i) {                      \
    assert(i>=0 || i==SimTK::InvalidId);                \
}                                                       \
inline NAME::NAME(long i) : id((int)i) {                \
    assert(i>=0 || i==SimTK::InvalidId);                \
}

/**
 * Add public static method declaration in class derived from an abstract
 * parent to assist in downcasting objects of the parent type to the 
 * derived type.
 */
#define SimTK_DOWNCAST(Derived,Parent) \
    static bool isA(const Parent& p)                        \
        { return dynamic_cast<const Derived*>(&p) != 0; }   \
    static const Derived& downcast(const Parent& p)         \
        { return dynamic_cast<const Derived&>(p); }         \
    static Derived& downcast(Parent& p)                     \
        { return dynamic_cast<Derived&>(p); }

/**
 * This is like SimTK_DOWNCAST except it allows for an intermediate
 * "helper" class between Derived and Parent.
 */
#define SimTK_DOWNCAST2(Derived,Helper,Parent) \
    static bool isA(const Parent& p)                                        \
        { return Helper::isA(p); }                                          \
    static const Derived& downcast(const Parent& p)                         \
        { return reinterpret_cast<const Derived&>(Helper::downcast(p)); }   \
    static Derived& downcast(Parent& p)                                     \
        { return reinterpret_cast<Derived&>(Helper::downcast(p)); }


// Similar to the above but for private implementation abstract classes, that
// is, abstract class hierarchies where the virtual function table is 
// hidden on the library side.
#define SimTK_PIMPL_DOWNCAST(Derived, Parent)           \
    static bool           isInstanceOf(const Parent&);  \
    static const Derived& downcast(const Parent&);      \
    static Derived&       updDowncast(Parent&)

namespace SimTK {
    
namespace Options { }
namespace Exception { }

typedef SimTK_Real              Real;
typedef std::complex<Real>      Complex;

struct Segment {
    Segment() : length(0), offset(0) { }
    explicit Segment(int l, int ofs=0) : length(l), offset(ofs) { 
        assert(l>=0 && ofs>=0);
    }
    // default copy, assignment, destructor
    int length;
    int offset;
};  

template <class T> class TypeInfo {
public:
    static const char* name() {return typeid(T).name();}
};

#define SimTK_TYPEINFO_SPECIALIZE(T)            \
template <> class TypeInfo< T > {               \
public:                                         \
    static const char* name() { return #T; }    \
};
// Built-in types
SimTK_TYPEINFO_SPECIALIZE(bool);            SimTK_TYPEINFO_SPECIALIZE(signed char); 
SimTK_TYPEINFO_SPECIALIZE(char);            SimTK_TYPEINFO_SPECIALIZE(unsigned char);
SimTK_TYPEINFO_SPECIALIZE(short);           SimTK_TYPEINFO_SPECIALIZE(int); 
SimTK_TYPEINFO_SPECIALIZE(long);
SimTK_TYPEINFO_SPECIALIZE(unsigned short);  SimTK_TYPEINFO_SPECIALIZE(unsigned int); 
SimTK_TYPEINFO_SPECIALIZE(unsigned long);
SimTK_TYPEINFO_SPECIALIZE(float);           SimTK_TYPEINFO_SPECIALIZE(double); 
SimTK_TYPEINFO_SPECIALIZE(long double);
SimTK_TYPEINFO_SPECIALIZE(std::complex<float>);
SimTK_TYPEINFO_SPECIALIZE(std::complex<double>); 
SimTK_TYPEINFO_SPECIALIZE(std::complex<long double>); 


} // namespace SimTK

namespace SimTKimpl {

///\{
/// Template-free signatures of TypeDescriptor methods	
typedef void*		(*IndexT)(void* tp, int n);
typedef const void*	(*IndexConstT)(const void* tp, int n);
typedef void*		(*CreateOneT)(const void* iptr);
typedef void		(*DestructOneT)(void*& tvptr);	
typedef void		(*AssignArrayOfT)(void* dest, const void* src, int n);
typedef void		(*SetT)(void* dest, const void* valuep, int n);
typedef void*		(*CreateArrayOfT)(int n, const void* iptr);
typedef void		(*DestructArrayOfT)(void*& tvptr);
///\}

struct TypeManipulatorT {
	TypeManipulatorT(int z, IndexT it, IndexConstT ict,
					 CreateOneT c1t, DestructOneT d1t, AssignArrayOfT aat, SetT st,
					 CreateArrayOfT cat, DestructArrayOfT dat)
		: sizeOfT(z), indexT(it), indexConstT(ict), 
		  createOneT(c1t), destructOneT(d1t), assignArrayOfT(aat), setT(st),
		  createArrayOfT(cat), destructArrayOfT(dat)
	{ }
    
    bool operator==(const TypeManipulatorT& t) const
      { return sizeOfT==t.sizeOfT && indexT==t.indexT; }

    // THESE MUST NEVER CHANGE ORDER! This is effectively a compiler-independent
    // virtual function table. The first entry is a size in bytes and the rest have
    // type "pointer to function" so all the data should have a very predictable
    // and stable layout in memory.
    //
    // This is the one place in the SimTK API where the client
    // side and library side must agree on the physical layout of the class. 
    // It is safe to add new entries to the end of this list, but inserting earlier,
    // deleting or reordering anything already here will break binary compatibility.
	const int               sizeOfT;
	const IndexT            indexT;
	const IndexConstT       indexConstT;
	const CreateOneT        createOneT;
	const DestructOneT      destructOneT;
	const AssignArrayOfT    assignArrayOfT;
	const SetT              setT;
	const CreateArrayOfT    createArrayOfT;
	const DestructArrayOfT  destructArrayOfT;
};

/** 
 * Templatized helper class builds a non-templatized descriptor for T.
 *
 * This can then be used to provide a hidden implementation of a visible
 * templatized class like List<T>. The hidden implementation uses services
 * from this class to manipulate the elements.
 * 
 * Note that all members are static, so there need be only one object
 * of class MakeTypeManipulator<T> for any type T. Also note that the whole
 * implementation of MakeTypeManipulator must be contained in this header file
 * so that its definition is always consistent with its declaration. This is
 * what allows us to use templatized classes in the SimTK API without risk
 * of source/binary incompatibilities.
 */
template <class T> class MakeTypeManipulator {
public:
	MakeTypeManipulator() { }
	// default copy, assignment, destructor
	
	static const TypeManipulatorT& getTypeManipulatorT() { return manipT; }
	static T& updAs(void* v) 
      { assert(v); return *reinterpret_cast<T*>(v); }
	static const T& getAs(const void* v) 
      { assert(v); return *reinterpret_cast<const T*>(v); }

private:	
	/// Return the number of bytes used to contain an object of type T.
	/// More specifically, this is the offset in bytes from one T to
	/// the next in an ordinary C++ array T[]. It does not matter if
	/// the object contains pointers to more memory in the heap; we
	/// expect the object to handle that itself.
	static int sizeT() { return (int)sizeof(T); }
	
	///\{
	/// These two routines perform an offset calculation on a void*
	/// pointer to an array of T to find the n'th following T and
	/// then return a pointer to it as a void*.
	static void* indexT(void* tp, int n)
		{ return (void*)(reinterpret_cast<T*>(tp) + n); }
	static const void* indexConstT(const void* tp, int n)
		{ return (const void*)(reinterpret_cast<const T*>(tp) + n); }
	///\}

	/// Heap allocate, default construct and optionally initialize a single
	/// object of type T. Return an obfuscated pointer to it as a void*.
	/// This object must be explicitly destructed with destructOneT.
	static void* createOneT(const void* iptr=0)
	{
		T* tptr = new T;
		if (iptr) *tptr = *reinterpret_cast<const T*>(iptr);
		return tptr;
	}
	
	/// Destruct an object that was created with createOneT. We zero out
	/// the passed-in pointer as a good hygiene measure.
	static void destructOneT(void*& tvptr)
	{
		T* tptr = reinterpret_cast<T*>(tvptr);
		delete tptr;
		tvptr = 0;
	}

	/// Copy an array of T's to another array of T's. The destination
	/// T's must have already been constructed for this to work properly
	/// although there is no way to check. 
	/// The n elements of source and destination must not overlap and
	/// both array pointers must be non-null, unless n is 0.		
	static void assignArrayOfT(void* dest, const void* src, int n)
	{   assert(n>=0);
		if (n==0) return;
		assert(dest && src);
		assert(indexT(dest,n) < src || indexConstT(src,n) < dest);		
		T*		 d = reinterpret_cast<T*>(dest);
		const T* s = reinterpret_cast<const T*>(src);
		for (int i=0; i < n; ++i) *d++ = *s++;
	}
		
	/// Assign a single value repeatedly to each element of an array of T's.
	/// we assume that the destination T's have already been constructed.
	/// It is an error to call this with null destination or value unless
	/// the number of elements is 0. 
	static void setT(void* dest, const void* valuep, int n=1)
	{   assert(n>=0);
		if (n==0) return;
		assert(dest && valuep);
		T*		 d = reinterpret_cast<T*>(dest);
		const T& v = *reinterpret_cast<const T*>(valuep);
		for (int i=0; i < n; ++i) *d++ = v; 
	}
	
	/// Allocate and default-construct an array of n T's, with 
	/// optional initialization.
	static void* createArrayOfT(int n, const void* iptr=0)
	{   assert(n>=0);
		if (n == 0) return 0;
		T* tptr = new T[n];
		if (iptr) {
			const T& init = *reinterpret_cast<const T*>(iptr);
			for (int i=0; i < n; ++i) tptr[i] = init;
		}		
		return tptr;
	}
	
	/// Destruct an array which was created with createArrayOfT. You
	/// <em>must</em> pass the original pointer (to the 0th element)
	/// here! You cannot use this to destruct part of an array, just
	/// the whole thing at once.
	static void destructArrayOfT(void*& tvptr)
	{
		T* tptr = reinterpret_cast<T*>(tvptr);
		delete[] tptr;
		tvptr = 0;
	}
	
private:	
	static const TypeManipulatorT manipT;
};

/*static*/ template <class T> const TypeManipulatorT
MakeTypeManipulator<T>::manipT = TypeManipulatorT(
									(int)sizeof(T),indexT,indexConstT,
				 					createOneT, destructOneT, assignArrayOfT, setT, 
				 					createArrayOfT, destructArrayOfT);
} // namespace SimTKimpl


#endif /* C++ stuff */

#endif /* SimTK_SimTKCOMMON_COMMON_H_ */
