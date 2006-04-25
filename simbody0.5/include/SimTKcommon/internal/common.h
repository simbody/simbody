#ifndef SimTK_SimTKCOMMON_COMMON_H_
#define SimTK_SimTKCOMMON_COMMON_H_

/* Copyright (c) 2005-6 Stanford University and Michael Sherman.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/**@file
 * Mandatory first inclusion for any SimTK source or header file.
 * 
 * Every source and most header files using SimTK must include this 
 * header as its <em>first</em> inclusion. Declarations and definitions that 
 * must be available and compiler-and machine-specific issues are dealt
 * with here.
 */

#include <cstddef>
#include <cassert>
#include <complex>

// Set up a few compile-time options that affect all SimTK Core headers.
#ifndef SimTK_DEFAULT_PRECISION
#   define SimTK_DEFAULT_PRECISION double
#endif

// Unless this symbol is set at compile time, we will not use dllimport on
// Windows. There *may* be a performance advantage to using it, although
// you can't then link with a static library.
// #define SimTK_OPTIMIZE_FOR_DYNAMIC_LIBRARY

#ifndef NDEBUG
#include <cstdio>
#define SimTK_DEBUG(s) printf("DBG: " s)
#define SimTK_DEBUG1(s,a1) printf("DBG: " s,a1)	
#define SimTK_DEBUG2(s,a1,a2) printf("DBG: " s,a1,a2)	
#define SimTK_DEBUG3(s,a1,a2,a3) printf("DBG: " s,a1,a2,a3)	
#define SimTK_DEBUG4(s,a1,a2,a3,a4) printf("DBG: " s,a1,a2,a3,a4)
#else
#define SimTK_DEBUG(s)
#define SimTK_DEBUG1(s,a1)
#define SimTK_DEBUG2(s,a1,a2)	
#define SimTK_DEBUG3(s,a1,a2,a3)	
#define SimTK_DEBUG4(s,a1,a2,a3,a4)
#endif

// When building a shared library 'xyz', CMake defines a symbol 'xyz_EXPORTS'
// for use in distinguishing builds from client use of a header. The following
// is specific for the current 'SimTK' library and doesn't affect other
// libraries even if they use this one.
#ifndef SimTK_SIMBODY_API
    #ifdef WIN32
        #ifdef simbody_EXPORTS
            #define SimTK_SIMBODY_API __declspec(dllexport)
        #elif defined(SimTK_OPTIMIZE_FOR_DYNAMIC_LIBRARY)
            #define SimTK_SIMBODY_API __declspec(dllimport)   // can't link with static lib now
        #else
            #define SimTK_SIMBODY_API // This works both for static & dynamic clients
        #endif
    #else
        #define SimTK_SIMBODY_API // Linux, Mac
    #endif
#endif

// Every SimTK Core library must provide these two routines, with the library
// name appearing after the "version_" and "about_".
extern "C" {
    SimTK_SIMBODY_API void SimTK_version_SimTKcommon(int* major, int* minor, int* build);
    SimTK_SIMBODY_API void SimTK_about_SimTKcommon(const char* key, int maxlen, char* value);
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

namespace SimTK {
    
namespace Options { }
namespace Exception { }
      
typedef float       FReal;
typedef double      DReal;
typedef long double LReal;

typedef std::complex<float>         FComplex;
typedef std::complex<double>        DComplex;
typedef std::complex<long double>   LComplex;

typedef SimTK_DEFAULT_PRECISION Real;
typedef std::complex<Real>      Complex;

struct Segment {
    Segment() : length(0), offset(0) { }
    explicit Segment(size_t l, size_t ofs=0) : length(l), offset(ofs) { }
    // default copy, assignment, destructor
    size_t length;
    size_t offset;
};  

template <class T> class TypeInfo {
public:
    static const char* name() { return typeid(T).name(); }
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
SimTK_TYPEINFO_SPECIALIZE(long);            SimTK_TYPEINFO_SPECIALIZE(long long); 
SimTK_TYPEINFO_SPECIALIZE(unsigned short);  SimTK_TYPEINFO_SPECIALIZE(unsigned int); 
SimTK_TYPEINFO_SPECIALIZE(unsigned long);   SimTK_TYPEINFO_SPECIALIZE(unsigned long long); 
SimTK_TYPEINFO_SPECIALIZE(float);           SimTK_TYPEINFO_SPECIALIZE(double); 
SimTK_TYPEINFO_SPECIALIZE(long double);
SimTK_TYPEINFO_SPECIALIZE(std::complex<float>);
SimTK_TYPEINFO_SPECIALIZE(std::complex<double>); 
SimTK_TYPEINFO_SPECIALIZE(std::complex<long double>); 

} // namespace SimTK

namespace SimTKimpl {

///\{
/// Template-free signatures of TypeDescriptor methods	
typedef void*		(*IndexT)(void* tp, size_t n);
typedef const void*	(*IndexConstT)(const void* tp, size_t n);
typedef void*		(*CreateOneT)(const void* iptr);
typedef void		(*DestructOneT)(void*& tvptr);	
typedef void		(*AssignArrayOfT)(void* dest, const void* src, size_t n);
typedef void		(*SetT)(void* dest, const void* valuep, size_t n);
typedef void*		(*CreateArrayOfT)(size_t n, const void* iptr);
typedef void		(*DestructArrayOfT)(void*& tvptr);
typedef bool        (*IsEqualT)(const void* t1, const void* t2);
///\}

struct TypeManipulatorT {
	TypeManipulatorT(size_t z, IndexT it, IndexConstT ict,
					 CreateOneT c1t, DestructOneT d1t, AssignArrayOfT aat, SetT st,
					 CreateArrayOfT cat, DestructArrayOfT dat, IsEqualT eq)
		: sizeOfT(z), indexT(it), indexConstT(ict), 
		  createOneT(c1t), destructOneT(d1t), assignArrayOfT(aat), setT(st),
		  createArrayOfT(cat), destructArrayOfT(dat), isEqualT(eq)
	{ }
    
    bool operator==(const TypeManipulatorT& t) const
      { return sizeOfT==t.sizeOfT && indexT==t.indexT; }

    // THESE MUST NEVER CHANGE ORDER! This is effectively a compiler-independent
    // virtual function table. The first entry is a size_t and the rest have
    // type "pointer to function" so all the data should have a very predictable
    // and stable layout in memory.
    //
    // This is the one place in the SimTK API where the client
    // side and library side must agree on the physical layout of the class. 
    // It is safe to add new entries to the end of this list, but inserting earlier,
    // deleting or reordering anything already here will break binary compatibility.
	const size_t			sizeOfT;
	const IndexT			indexT;
	const IndexConstT		indexConstT;
	const CreateOneT		createOneT;
	const DestructOneT		destructOneT;
	const AssignArrayOfT	assignArrayOfT;
	const SetT				setT;
	const CreateArrayOfT	createArrayOfT;
	const DestructArrayOfT	destructArrayOfT;
    const IsEqualT          isEqualT;
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
	static size_t sizeT() { return sizeof(T); }
	
	///\{
	/// These two routines perform an offset calculation on a void*
	/// pointer to an array of T to find the n'th following T and
	/// then return a pointer to it as a void*.
	static void* indexT(void* tp, size_t n)
		{ return (void*)(reinterpret_cast<T*>(tp) + n); }
	static const void* indexConstT(const void* tp, size_t n)
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
	static void assignArrayOfT(void* dest, const void* src, size_t n)
	{
		if (n==0) return;
		assert(dest && src);
		assert(indexT(dest,n) < src || indexConstT(src,n) < dest);		
		T*		 d = reinterpret_cast<T*>(dest);
		const T* s = reinterpret_cast<const T*>(src);
		for (size_t i=0; i < n; ++i) *d++ = *s++;
	}
		
	/// Assign a single value repeatedly to each element of an array of T's.
	/// we assume that the destination T's have already been constructed.
	/// It is an error to call this with null destination or value unless
	/// the number of elements is 0. 
	static void setT(void* dest, const void* valuep, size_t n=1)
	{
		if (n==0) return;
		assert(dest && valuep);
		T*		 d = reinterpret_cast<T*>(dest);
		const T& v = *reinterpret_cast<const T*>(valuep);
		for (size_t i=0; i < n; ++i) *d++ = v; 
	}
	
	/// Allocate and default-construct an array of n T's, with 
	/// optional initialization.
	static void* createArrayOfT(size_t n, const void* iptr=0)
	{
		if (n == 0) return 0;
		T* tptr = new T[n];
		if (iptr) {
			const T& init = *reinterpret_cast<const T*>(iptr);
			for (size_t i=0; i < n; ++i) tptr[i] = init;
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

    /// Compare two T's and return true if their values should be considered
    /// equal.
    static bool isEqualT(const void* v1, const void* v2) {
        const T& t1 = *reinterpret_cast<const T*>(v1);
        const T& t2 = *reinterpret_cast<const T*>(v2);
        return t1 == t2;
    }
	
private:	
	static const TypeManipulatorT manipT;
};

/*static*/ template <class T> const TypeManipulatorT
MakeTypeManipulator<T>::manipT = TypeManipulatorT(
									sizeof(T),indexT,indexConstT,
				 					createOneT, destructOneT, assignArrayOfT, setT, 
				 					createArrayOfT, destructArrayOfT, isEqualT);
} // namespace SimTKimpl

#endif // SimTK_SimTKCOMMON_COMMON_H_
