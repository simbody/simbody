#ifndef SimTK_SimTKCOMMON_COMMON_H_
#define SimTK_SimTKCOMMON_COMMON_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-15 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors: Chris Dembia                                                 *
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
Mandatory first inclusion for any Simbody source or header file.

Every source and most header files using %SimTK must include this 
header as its \e first inclusion. Declarations and definitions that 
must be available and compiler-and machine-specific issues are dealt
with here.

This file must be includable from either C++ or ANSI C. It uses
the ANSI-C++ macro "__cplusplus" for any code that will compile
only under C++. **/

// Provide doxygen documentation for the SimTK namespace.

/**@namespace SimTK
This is the top-level %SimTK namespace into which all %SimTK names are 
placed to avoid collision with other symbols. If you get tired of prefacing 
every symbol with "SimTK::", include the statement "using namespace SimTK;" 
at the beginning of your %SimTK-using compilation units. Any names which 
cannot be put in the namespace (macro names, for example) begin with the 
prefix "SimTK_" instead. **/

// Define shared doxygen "modules" and sub-modules here. We'll put things 
// in them at various places when appropriate.

/**@defgroup GlobalFunctions Global Functions in the SimTK namespace
These are functions at the top level of the SimTK namespace, meaning
that a function named funcName() is invoked as SimTK::funcName(), or
just funcName() if there is a "using namespace SimTK;" statement in effect. **/

/**@defgroup ScalarFunctions Scalar Functions
   @ingroup GlobalFunctions
These functions are overloaded to act on %SimTK scalar types and C++
built-in types, including integral types when appropriate. **/

/**@defgroup BitFunctions Bit-twiddling Functions
   @ingroup GlobalFunctions
These functions perform highly optimized bit-twiddling operations on
the built-in integral types, and sometimes on the representations of
floating point types as well. **/

/**@defgroup Serialization  Utilities for De/serializing
   @ingroup GlobalFunctions
These namespace-scope templatized utilities provide uniform serialization
and deserialization behavior for built-in and SimTK-defined types. See
SimTK::Xml for support of serialization to/from Xml files. **/
    
/**@defgroup UniqueIndexTypes    Type-Safe Integer Indices

It is common to store objects or information about them in randomly-indexable 
arrays, and then to support maximum-performance selection by allowing the
index to be used. We want these arrays indexable by simple ints for speed, but
this quickly leads to APIs in which there are multiple int arguments in a
function call, each intended to select a different kind of object. A common
error when there is a series of identical argument types is to put them in
the wrong order. To avoid that, we define unique index types here for 
accessing each category to help stay out of trouble.

A unique index type is just a type-safe non-negative int, augmented with a 
"NaN" value called InvalidBLAH where BLAH is the type name. For most uses it 
will behave like an int, and it has an implicit conversion *to* int. Importantly
though, it has no implicit conversion *from* int so you can't pass a plain int 
or any other Index type to an argument expecting a certain Index type. **/

/*****************************/
/* ANSI-C COMPATIBLE SECTION */
/*****************************/

/* Set up a few compile-time options that affect all SimTK Core headers. */

/**
 * This compile-time constant determines the default precision used everywhere
 * in %SimTK Core code. Wherever a SimTK::Real, SimTK::Vector, SimTK::Matrix,
 * etc. appears with no precision specified, it will have this underlying precision.
 * We use 1==float, 2==double. Any other value will cause
 * a compile time error. The default is 2, i.e., double precision.
 */
#ifndef SimTK_DEFAULT_PRECISION
#   define SimTK_DEFAULT_PRECISION 2
#endif

#if   (SimTK_DEFAULT_PRECISION == 1)
/** This type is for use in C; in C++ use SimTK::Real instead. */
    typedef float SimTK_Real;
#elif (SimTK_DEFAULT_PRECISION == 2)
/** This type is for use in C; in C++ use SimTK::Real instead. */
    typedef double SimTK_Real;
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

#ifdef _WIN32
    #ifdef _MSC_VER
    #pragma warning(disable:4231) /*need to use 'extern' template explicit instantiation*/
    #pragma warning(disable:4251) /*no DLL interface for type of member of exported class*/
    #pragma warning(disable:4275) /*no DLL interface for base class of exported class*/
    #pragma warning(disable:4345) /*warning about PODs being default-initialized*/


    /* Until VS2015 struct timespec was missing from <ctime> so is faked here 
    if needed. When Simbody used pthreads and provided its own pthread.h for
    Windows, we had to avoid a duplicate declaration with timespec in pthread.h
    via the HAVE_STRUCT_TIMESPEC guard. In 2018, we removed pthread.h, but we
    left in the HAVE_STRUCT_TIMESPEC guard in case a third party defines
    timespec.
    TODO: there is a potential problem here since VS2015's struct timespec 
    doesn't appear to match pthread's definition. */
    #ifndef HAVE_STRUCT_TIMESPEC
    #define HAVE_STRUCT_TIMESPEC 1
        #if _MSC_VER < 1900
        struct timespec {
            long tv_sec; /* TODO(sherm1,chrisdembia) should be time_t? */
            long tv_nsec;
        };
        #endif
    #endif /* HAVE_STRUCT_TIMESPEC */
    #endif
    #if defined(SimTK_SimTKCOMMON_BUILDING_SHARED_LIBRARY)
        #ifdef _MSC_VER
        #define SimTK_SimTKCOMMON_EXPORT __declspec(dllexport)
        /* Keep MS VC++ quiet when it tries to instantiate incomplete template classes in a DLL. */
        #pragma warning(disable:4661)
        #else
        #define SimTK_SimTKCOMMON_EXPORT
        #endif
    #elif defined(SimTK_SimTKCOMMON_BUILDING_STATIC_LIBRARY) || defined(SimTK_USE_STATIC_LIBRARIES)
        #define SimTK_SimTKCOMMON_EXPORT
    #else
        #ifdef _MSC_VER
        #define SimTK_SimTKCOMMON_EXPORT __declspec(dllimport) /*i.e., a client of a shared library*/
        #else
        #define SimTK_SimTKCOMMON_EXPORT
        #endif
    #endif
    /* VC++ tries to be secure by leaving bounds checking on for STL containers
     * even in Release mode. This macro exists to disable that feature and can
     * result in a considerable speedup.
     * CAUTION: every linked-together compilation unit must have this set the same
     * way. Everyone who properly includes this file first is fine; but as of this
     * writing Simmath's IpOpt doesn't do so.
     * NOTE: Microsoft corrected this problem with VC10 -- the feature is 
     * disabled by default in that compiler and later.
     */
    /* (sherm 081204 disabling for now: doesn't work on VC++ 8 and is 
     * tricky on VC++ 9 because all libraries, including 3rd party, must
     * be built the same way). Better to use the SimTK::Array_<T> class in
     * place of the std::vector<T> class to get better performance.
     #ifdef NDEBUG
         #undef _SECURE_SCL
         #define _SECURE_SCL 0
     #endif
     */
#else
    #define SimTK_SimTKCOMMON_EXPORT // Linux, Mac
#endif

/* Every SimTK Core library must provide these two routines, with the library
 * name appearing after the "version_" and "about_".
 */
#if defined(__cplusplus)
extern "C" {
#endif
    /** Obtain version information for the currently-loaded SimTKcommon library. */
    SimTK_SimTKCOMMON_EXPORT void SimTK_version_SimTKcommon(int* major, int* minor, int* build);
    /** 
     * Obtain "about" information for the currently-loaded SimTKcommon library.
     * Available keywords are "version" (major.minor.build), "library", 
     * "type" (shared or static), "copyright", "svn_revision", "authors", 
     * "debug" (debug or release).
     */
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
#include <cstring>
#include <cmath>
#include <cfloat>
#include <complex>
#include <limits>
#include <typeinfo>
#include <algorithm>

/* Be very careful with this macro -- don't use it unless you have measured
a performance improvement. You can end up with serious code bloat if you 
override the compiler's judgement about when to inline, and that can cause
cache misses which ultimately reduce performance. */
#ifdef _MSC_VER
    #define SimTK_FORCE_INLINE __forceinline
#else
    #define SimTK_FORCE_INLINE __attribute__((always_inline)) inline
#endif

/* Microsoft added noexcept in VS2015 */
#if defined(_MSC_VER) && _MSC_VER < 1900
    #define NOEXCEPT_11 throw()
#else
    #define NOEXCEPT_11 noexcept
#endif

/* C++14 introduces a standard way to mark deprecated declarations. Before
that we can use non-standard compiler hacks. */
#ifndef SWIG
    #if __cplusplus >= 201402L
        /* C++14 */
        #define DEPRECATED_14(MSG) [[deprecated(MSG)]]
    #elif _MSC_VER
        /* VC++ just says warning C4996 so add "DEPRECATED" to the message. */
        #define DEPRECATED_14(MSG) __declspec(deprecated("DEPRECATED: " MSG))
    #else /* gcc or clang */
        #define DEPRECATED_14(MSG) __attribute__((deprecated(MSG)))
    #endif
#else /* Swigging */
    #define DEPRECATED_14(MSG)
#endif

/* These macros are deprecated, leftover from before C++11 was available. 
Don't use them. Sorry, can't use the DEPRECATED_14 macro here! */
#define OVERRIDE_11 override
#define FINAL_11 final

namespace SimTK {


// This utility answers the question "if I put this integral value in an int and then
// get it back, will its value be the same?".
inline bool canStoreInInt(bool)            {return true;}
inline bool canStoreInInt(char)            {return true;}
inline bool canStoreInInt(unsigned char)   {return true;}
inline bool canStoreInInt(signed char)     {return true;}
inline bool canStoreInInt(short)           {return true;}
inline bool canStoreInInt(unsigned short)  {return true;}
inline bool canStoreInInt(int)             {return true;}
inline bool canStoreInInt(unsigned int  u) {return (unsigned int)(int(u)) == u;}
inline bool canStoreInInt(long i)          {return long(int(i)) == i;}
inline bool canStoreInInt(unsigned long u) {return (unsigned long)(int(u)) == u;}
inline bool canStoreInInt(long long i)          {return (long long)(int(i)) == i;}
inline bool canStoreInInt(unsigned long long u) {return (unsigned long long)(int(u)) == u;}

// This utility answers the question "is this integral value a nonnegative number
// that can be stored in an int?".
inline bool canStoreInNonnegativeInt(bool)             {return true;}
inline bool canStoreInNonnegativeInt(char c)           {return c >= 0;}
inline bool canStoreInNonnegativeInt(unsigned char)    {return true;}
inline bool canStoreInNonnegativeInt(signed char c)    {return c >= 0;}
inline bool canStoreInNonnegativeInt(short s)          {return s >= 0;}
inline bool canStoreInNonnegativeInt(unsigned short)   {return true;}
inline bool canStoreInNonnegativeInt(int  i)           {return i >= 0;}
inline bool canStoreInNonnegativeInt(long l)           {return canStoreInInt(l) && l >= 0;}
inline bool canStoreInNonnegativeInt(long long l)      {return canStoreInInt(l) && l >= 0;}
inline bool canStoreInNonnegativeInt(unsigned int  u)  {return canStoreInInt(u);}
inline bool canStoreInNonnegativeInt(unsigned long u)  {return canStoreInInt(u);}
inline bool canStoreInNonnegativeInt(unsigned long long u) {return canStoreInInt(u);}

// This utility answers the question of whether an integer is suitable as a size
// limited by the given maximum size. Signed types must be checked for being
// nonegative; doing that with unsigned types leads to compiler warnings.

// char can be signed or unsigned depending on the compiler; assume signed.
inline bool isSizeInRange(char           sz, char           mx){return 0<=sz&&sz<=mx;}
inline bool isSizeInRange(signed char    sz, signed char    mx){return 0<=sz&&sz<=mx;}
inline bool isSizeInRange(short          sz, short          mx){return 0<=sz&&sz<=mx;}
inline bool isSizeInRange(int            sz, int            mx){return 0<=sz&&sz<=mx;}
inline bool isSizeInRange(long           sz, long           mx){return 0<=sz&&sz<=mx;}
inline bool isSizeInRange(long long      sz, long long      mx){return 0<=sz&&sz<=mx;}
inline bool isSizeInRange(unsigned char  sz, unsigned char  mx){return sz<=mx;}
inline bool isSizeInRange(unsigned short sz, unsigned short mx){return sz<=mx;}
inline bool isSizeInRange(unsigned int   sz, unsigned int   mx){return sz<=mx;}
inline bool isSizeInRange(unsigned long  sz, unsigned long  mx){return sz<=mx;}
inline bool isSizeInRange(unsigned long long sz, unsigned long long mx){return sz<=mx;}

// This utility answers the question of whether an integer is suitable as an index
// for an array limited by the given maximum size. Signed types must be checked for being
// nonegative; doing that with unsigned types leads to compiler warnings. This is just
// like the "size in range" check above except the maximum value allowed for an index
// is one less that the size.

// char can be signed or unsigned depending on the compiler; assume signed.
inline bool isIndexInRange(char           ix, char           sz){return 0<=ix&&ix<sz;}
inline bool isIndexInRange(signed char    ix, signed char    sz){return 0<=ix&&ix<sz;}
inline bool isIndexInRange(short          ix, short          sz){return 0<=ix&&ix<sz;}
inline bool isIndexInRange(int            ix, int            sz){return 0<=ix&&ix<sz;}
inline bool isIndexInRange(long           ix, long           sz){return 0<=ix&&ix<sz;}
inline bool isIndexInRange(long long      ix, long long      sz){return 0<=ix&&ix<sz;}
inline bool isIndexInRange(unsigned char  ix, unsigned char  sz){return ix<sz;}
inline bool isIndexInRange(unsigned short ix, unsigned short sz){return ix<sz;}
inline bool isIndexInRange(unsigned int   ix, unsigned int   sz){return ix<sz;}
inline bool isIndexInRange(unsigned long  ix, unsigned long  sz){return ix<sz;}
inline bool isIndexInRange(unsigned long long ix, unsigned long long sz){return ix<sz;}

// This utility answers the question: is this integral value nonnegative? The answer
// is always true for unsigned types and you'll get a warning from some compilers if
// you check.

inline bool isNonnegative(bool)              {return true;}
// char can be signed or unsigned depending on the compiler; assume signed.
inline bool isNonnegative(char        n)     {return n>=0;}
inline bool isNonnegative(signed char n)     {return n>=0;}
inline bool isNonnegative(short       n)     {return n>=0;}
inline bool isNonnegative(int         n)     {return n>=0;}
inline bool isNonnegative(long        n)     {return n>=0;}
inline bool isNonnegative(long long   n)     {return n>=0;}
inline bool isNonnegative(unsigned char)     {return true;}
inline bool isNonnegative(unsigned short)    {return true;}
inline bool isNonnegative(unsigned int)      {return true;}
inline bool isNonnegative(unsigned long)     {return true;}
inline bool isNonnegative(unsigned long long){return true;}

// A NaN-like value for unique index types created using the macro
// SimTK_DEFINE_UNIQUE_INDEX_TYPE(). A unique, typed constant with
// this numerical value is created for each index type.
static const int InvalidIndex = -1111111111;
}



/**
 * Use this macro to define a unique "Index" type which is just a type-safe
 * non-negative int, augmented with a "NaN" value given by the predefined
 * int constant SimTK::InvalidIndex. We also allow the Index to take on
 * the value -1 if that is produced by a subtraction operation acting on a 
 * previously-valid Index, since that can occur during loops which are 
 * processed from the end towards the beginning. -1 is then allowed in 
 * comparison operators but not in any other operations, including further 
 * decrementing.
 *
 * No namespace is assumed for the newly-defined type; if you want the 
 * symbol in a namespace be sure to invoke the macro within that namespace. 
 * Make sure that the statement "#include <cassert>" appears somewhere before 
 * the point of invocation of this macro, because the defined Index type uses 
 * the assert() macro when in Debug mode.
 *
 * For most uses it will behave like an int, and it has an implicit
 * conversion \e to int. Importantly though, it has no implicit conversion
 * \e from int so you can't pass some other kind of number where a particular
 * kind of Index was expected. This is used to create Index types
 * which can be used as array indices but which prevent accidental mixing
 * of types. Examples: SubsystemIndex, ConstraintIndex.
 *
 * If you create a type "ThingIndex" you will also get a constant of
 * type ThingIndex named "InvalidThingIndex" which will be the initial
 * value of any objects of type ThingIndex, and will have the same numerical
 * value as SimTK::InvalidIndex.
 */

/** Define a global (that is, SimTK namespace level) Index class that is not 
exported in MS VC++ DLLs. **/
#define SimTK_DEFINE_UNIQUE_INDEX_TYPE(NAME)                   \
    SimTK_DEFINE_AND_EXPORT_UNIQUE_LOCAL_INDEX_TYPE(,,,NAME)   \
    static const NAME Invalid ## NAME;

/** Define a global (that is, SimTK namespace level) Index class with a MS VC++
"export" specification for DLLs. **/
#define SimTK_DEFINE_AND_EXPORT_UNIQUE_INDEX_TYPE(EXPORT,NAME)     \
    SimTK_DEFINE_AND_EXPORT_UNIQUE_LOCAL_INDEX_TYPE(EXPORT,,,NAME) \
    static const NAME Invalid ## NAME;

/** Define a local Index class within a Parent class. **/
#define SimTK_DEFINE_UNIQUE_LOCAL_INDEX_TYPE(PARENT,NAME) \
    SimTK_DEFINE_AND_EXPORT_UNIQUE_LOCAL_INDEX_TYPE(,PARENT,::,NAME)

/** The most general form allows a MS VC++ "export" specification for DLLs,
and a Parent class (with SEP=::) for local Index names. **/
#define SimTK_DEFINE_AND_EXPORT_UNIQUE_LOCAL_INDEX_TYPE(EXPORT,PARENT,SEP,NAME)   \
class EXPORT NAME {                         \
    int ix;                                 \
public:                                     \
    NAME() : ix(SimTK::InvalidIndex) { }       \
    explicit NAME(int i) : ix(i)      {assert(i>=0 || i==SimTK::InvalidIndex);} \
    explicit NAME(long l): ix((int)l) {assert(SimTK::canStoreInNonnegativeInt(l));}    \
    explicit NAME(long long l): ix((int)l) {assert(SimTK::canStoreInNonnegativeInt(l));}    \
    explicit NAME(unsigned int  u)  : ix((int)u)  {assert(SimTK::canStoreInInt(u));}   \
    explicit NAME(unsigned long ul) : ix((int)ul) {assert(SimTK::canStoreInInt(ul));}  \
    explicit NAME(unsigned long long ul) : ix((int)ul) {assert(SimTK::canStoreInInt(ul));}  \
    operator int() const {return ix;}               \
    bool isValid() const {return ix>=0;}            \
    bool isValidExtended() const {return ix>=-1;}   \
    void invalidate(){clear();}                     \
    void clear(){ix=SimTK::InvalidIndex;}           \
    \
    bool operator==(int  i) const {assert(isValidExtended() && isValidExtended(i)); return ix==i;}    \
    bool operator==(short s) const{assert(isValidExtended() && isValidExtended(s)); return ix==(int)s;}  \
    bool operator==(long l) const {assert(isValidExtended() && isValidExtended(l)); return ix==(int)l;}  \
    bool operator==(long long l) const {assert(isValidExtended() && isValidExtended(l)); return ix==(int)l;}  \
    bool operator==(unsigned int  u)  const {assert(isValidExtended() && isValid(u)); return ix==(int)u;}   \
    bool operator==(unsigned short us)const {assert(isValidExtended() && isValid(us)); return ix==(int)us;} \
    bool operator==(unsigned long ul) const {assert(isValidExtended() && isValid(ul)); return ix==(int)ul;} \
    bool operator==(unsigned long long ul) const {assert(isValidExtended() && isValid(ul)); return ix==(int)ul;} \
    bool operator!=(int  i)           const {return !operator==(i);}    \
    bool operator!=(short s)          const {return !operator==(s);}    \
    bool operator!=(long l)           const {return !operator==(l);}    \
    bool operator!=(long long l)      const {return !operator==(l);}    \
    bool operator!=(unsigned int  u)  const {return !operator==(u);}    \
    bool operator!=(unsigned long ul) const {return !operator==(ul);}   \
    bool operator!=(unsigned long long ul) const {return !operator==(ul);}   \
    \
    bool operator< (int  i) const {assert(isValidExtended() && isValidExtended(i)); return ix<i;}        \
    bool operator< (short s) const{assert(isValidExtended() && isValidExtended(s)); return ix<(int)s;}   \
    bool operator< (long l) const {assert(isValidExtended() && isValidExtended(l)); return ix<(int)l;}   \
    bool operator< (long long l) const {assert(isValidExtended() && isValidExtended(l)); return ix<(int)l;}   \
    bool operator< (unsigned int  u)  const {assert(isValidExtended() && isValid(u));  return ix<(int)u;}    \
    bool operator< (unsigned short us)const {assert(isValidExtended() && isValid(us)); return ix<(int)us;}   \
    bool operator< (unsigned long ul) const {assert(isValidExtended() && isValid(ul)); return ix<(int)ul;}   \
    bool operator< (unsigned long long ul) const {assert(isValidExtended() && isValid(ul)); return ix<(int)ul;}   \
    bool operator>=(int  i)           const {return !operator<(i);}    \
    bool operator>=(short s)          const {return !operator<(s);}    \
    bool operator>=(long l)           const {return !operator<(l);}    \
    bool operator>=(long long l)           const {return !operator<(l);}    \
    bool operator>=(unsigned int  u)  const {return !operator<(u);}    \
    bool operator>=(unsigned short us)const {return !operator<(us);}   \
    bool operator>=(unsigned long ul) const {return !operator<(ul);}   \
    bool operator>=(unsigned long long ul) const {return !operator<(ul);}   \
    \
    bool operator> (int  i) const {assert(isValidExtended() && isValidExtended(i)); return ix>i;}        \
    bool operator> (short s) const{assert(isValidExtended() && isValidExtended(s)); return ix>(int)s;}   \
    bool operator> (long l) const {assert(isValidExtended() && isValidExtended(l)); return ix>(int)l;}   \
    bool operator> (long long l) const {assert(isValidExtended() && isValidExtended(l)); return ix>(int)l;}   \
    bool operator> (unsigned int  u)  const {assert(isValidExtended() && isValid(u));  return ix>(int)u;}    \
    bool operator> (unsigned short us)const {assert(isValidExtended() && isValid(us)); return ix>(int)us;}   \
    bool operator> (unsigned long ul) const {assert(isValidExtended() && isValid(ul)); return ix>(int)ul;}   \
    bool operator> (unsigned long long ul) const {assert(isValidExtended() && isValid(ul)); return ix>(int)ul;}   \
    bool operator<=(int  i)           const {return !operator>(i);}    \
    bool operator<=(short s)          const {return !operator>(s);}    \
    bool operator<=(long l)           const {return !operator>(l);}    \
    bool operator<=(long long l)      const {return !operator>(l);}    \
    bool operator<=(unsigned int  u)  const {return !operator>(u);}    \
    bool operator<=(unsigned short us)const {return !operator>(us);}   \
    bool operator<=(unsigned long ul) const {return !operator>(ul);}   \
    bool operator<=(unsigned long long ul) const {return !operator>(ul);}   \
    \
    const NAME& operator++() {assert(isValid()); ++ix; return *this;}       /*prefix */   \
    NAME operator++(int)     {assert(isValid()); ++ix; return NAME(ix-1);}  /*postfix*/   \
    const NAME& operator--() {assert(isValid()); --ix; return *this;}       /*prefix */   \
    NAME operator--(int)     {assert(isValid()); --ix; return NAME(ix+1);}  /*postfix*/   \
    NAME next() const {assert(isValid()); return NAME(ix+1);}                             \
    NAME prev() const {assert(isValid()); return NAME(ix-1);} /*might return -1*/         \
    \
    NAME& operator+=(int i)  {assert(isValid() && isValidExtended(ix+i)); ix+=i; return *this;}     \
    NAME& operator-=(int i)  {assert(isValid() && isValidExtended(ix-i)); ix-=i; return *this;}     \
    NAME& operator+=(short s){assert(isValid() && SimTK::canStoreInInt(s) && isValidExtended(ix+(int)s)); ix+=(int)s; return *this;}     \
    NAME& operator-=(short s){assert(isValid() && SimTK::canStoreInInt(s) && isValidExtended(ix-(int)s)); ix-=(int)s; return *this;}     \
    NAME& operator+=(long l) {assert(isValid() && SimTK::canStoreInInt(l) && isValidExtended(ix+(int)l)); ix+=(int)l; return *this;}     \
    NAME& operator-=(long l) {assert(isValid() && SimTK::canStoreInInt(l) && isValidExtended(ix-(int)l)); ix-=(int)l; return *this;}     \
    NAME& operator+=(long long l) {assert(isValid() && SimTK::canStoreInInt(l) && isValidExtended(ix+(int)l)); ix+=(int)l; return *this;}     \
    NAME& operator-=(long long l) {assert(isValid() && SimTK::canStoreInInt(l) && isValidExtended(ix-(int)l)); ix-=(int)l; return *this;}     \
    NAME& operator+=(unsigned int  u)  {assert(isValid()&& SimTK::canStoreInInt(u)  && isValid(ix+(int)u));  ix+=(int)u;  return *this;}  \
    NAME& operator-=(unsigned int  u)  {assert(isValid()&& SimTK::canStoreInInt(u)  && isValidExtended(ix-(int)u));  ix-=(int)u;  return *this;}  \
    NAME& operator+=(unsigned short us){assert(isValid()&& SimTK::canStoreInInt(us) && isValid(ix+(int)us)); ix+=(int)us; return *this;}  \
    NAME& operator-=(unsigned short us){assert(isValid()&& SimTK::canStoreInInt(us) && isValidExtended(ix-(int)us)); ix-=(int)us; return *this;}  \
    NAME& operator+=(unsigned long ul) {assert(isValid()&& SimTK::canStoreInInt(ul) && isValid(ix+(int)ul)); ix+=(int)ul; return *this;}  \
    NAME& operator-=(unsigned long ul) {assert(isValid()&& SimTK::canStoreInInt(ul) && isValidExtended(ix-(int)ul)); ix-=(int)ul; return *this;}  \
    NAME& operator+=(unsigned long long ul) {assert(isValid()&& SimTK::canStoreInInt(ul) && isValid(ix+(int)ul)); ix+=(int)ul; return *this;}  \
    NAME& operator-=(unsigned long long ul) {assert(isValid()&& SimTK::canStoreInInt(ul) && isValidExtended(ix-(int)ul)); ix-=(int)ul; return *this;}  \
    \
    static const NAME& Invalid() {static const NAME invalid; return invalid;}       \
    static bool isValid(int  i) {return i>=0;}                                      \
    static bool isValid(short s){return s>=0;}                                      \
    static bool isValid(long l) {return SimTK::canStoreInNonnegativeInt(l);}        \
    static bool isValid(long long l) {return SimTK::canStoreInNonnegativeInt(l);}        \
    static bool isValid(unsigned int  u)  {return SimTK::canStoreInInt(u);}         \
    static bool isValid(unsigned short)   {return true;}                            \
    static bool isValid(unsigned long ul) {return SimTK::canStoreInInt(ul);}        \
    static bool isValid(unsigned long long ul) {return SimTK::canStoreInInt(ul);}        \
    static bool isValidExtended(int  i) {return i>=-1;}                             \
    static bool isValidExtended(short s){return s>=-1;}                             \
    static bool isValidExtended(long l) {return SimTK::canStoreInInt(l) && l>=-1;}  \
    static bool isValidExtended(long long l) {return SimTK::canStoreInInt(l) && l>=-1;}  \
    /* IndexTraits for use in Array_<T,X> with this as X; same as int */            \
    typedef int size_type;                                                  \
    typedef int difference_type;                                            \
    static size_type max_size() {return std::numeric_limits<int>::max();}   \
};

/** Use this macro to generate a cast that is dynamic_cast in Debug builds
but static_cast in Release builds, for uses where you don't want to
pay for the extra safety. Caution: these are not necessarily equivalent for
dynamic types that use multiple inheritance; don't use this macro in that
case, and don't use it where you are using dynamic_cast on a pointer to
check what type of derived object you're looking at. **/
#ifndef NDEBUG
    #define SimTK_DYNAMIC_CAST_DEBUG dynamic_cast   // safe but slow
#else
    #define SimTK_DYNAMIC_CAST_DEBUG static_cast    // unsafe but fast
#endif

/** Add public static method declaration in class derived from an abstract
parent to assist in downcasting objects of the parent type to the derived 
type. **/
#define SimTK_DOWNCAST(Derived,Parent)                          \
    static bool isA(const Parent& p)                            \
        { return dynamic_cast<const Derived*>(&p) != 0; }       \
    static const Derived& downcast(const Parent& p)             \
        { return SimTK_DYNAMIC_CAST_DEBUG<const Derived&>(p); } \
    static Derived& updDowncast(Parent& p)                      \
        { return SimTK_DYNAMIC_CAST_DEBUG<Derived&>(p); }        \
    static Derived& downcast(Parent& p)                         \
        { return SimTK_DYNAMIC_CAST_DEBUG<Derived&>(p); }

/** This is like SimTK_DOWNCAST except it allows for an intermediate "helper" 
class between Derived and Parent. **/
#define SimTK_DOWNCAST2(Derived,Helper,Parent)                          \
    static bool isA(const Parent& p)                                    \
        { return Helper::isA(p); }                                      \
    static const Derived& downcast(const Parent& p)                     \
        { return static_cast<const Derived&>(Helper::downcast(p)); }    \
    static Derived& updDowncast(Parent& p)                                \
        { return static_cast<Derived&>(Helper::downcast(p)); }            \
    static Derived& downcast(Parent& p)                                 \
        { return static_cast<Derived&>(Helper::downcast(p)); }


/** Similar to the above but for private implementation abstract classes, that
is, abstract class hierarchies where the virtual function table is hidden on 
the library side. **/
#define SimTK_PIMPL_DOWNCAST(Derived, Parent)           \
    static bool           isInstanceOf(const Parent&);  \
    static const Derived& downcast(const Parent&);      \
    static Derived&       updDowncast(Parent&)

namespace SimTK {

/** This sub-namespace of SimTK is used for the exception types that are
thrown by our error handing code. **/
namespace Exception { }

/** This is the default compiled-in floating point type for SimTK, either
float or double. @see SimTK_DEFAULT_PRECISION **/
typedef SimTK_Real              Real;
/** This is the default complex type for SimTK, with precision for the real 
and imaginary parts set to the compiled-in Real type. @see Real **/
typedef std::complex<Real>      Complex;
/** An abbreviation for std::complex<float> for consistency with others. **/
typedef std::complex<float>     fComplex;
/** An abbreviation for std::complex<double> for consistency with others. **/
typedef std::complex<double>    dComplex;


// Forward declaration giving template defaults must come before any
// other declarations.
template <int M, class ELT=Real, int STRIDE=1>              class Vec;
template <int N, class ELT=Real, int STRIDE=1>              class Row; 
template <int M, int N, class ELT=Real, int CS=M, int RS=1> class Mat;
template <int M, class ELT=Real, int RS=1>                  class SymMat;

/** A convenient struct for anything requiring an offset and length to specify
a segment of some larger sequence. **/
struct Segment {
    Segment() : length(0), offset(0) { }
    explicit Segment(int l, int ofs=0) : length(l), offset(ofs) { 
        assert(l>=0 && ofs>=0);
    }
    // default copy, assignment, destructor
    int length;
    int offset;
};  

// These next four methods supply the missing relational operators for any
// types L and R where L==R and L<R have been defined. This is like the
// operators in the std::rel_ops namespace, except that those require both
// types to be the same.

template<class L, class R> inline
bool operator!=(const L& left, const R& right)
{   // test for inequality, in terms of equality
    return !(left == right);
}

template<class L, class R> inline
bool operator>(const L& left, const R& right)
{   // test if left > right, in terms of operator<
    return right < left;
}

template<class L, class R> inline
bool operator<=(const L& left, const R& right)
{   // test if left <= right, in terms of operator<
    return !(right < left);
}

template<class L, class R> inline
bool operator>=(const L& left, const R& right)
{   // test if left >= right, in terms of operator<
    return !(left < right);
}


/** This is a special type used for causing invocation of a particular
constructor or method overload that will avoid making a copy of the source
(that is, perform a "shallow" copy rather than a "deep" copy). Typically these
methods will have some dangerous side effects so make sure you know what you're
doing. **/
struct DontCopy {};
/** This is a special type used for forcing invocation of a particularly
dangerous constructor or method overload; don't use this unless you are an
advanced user and know exactly what you're getting into. **/
struct TrustMe {};

/** This is a compile-time equivalent of "false", used in compile-time
condition checking in templatized implementations. **/
struct FalseType {};
/** This is a compile-time equivalent of "true", used in compile-time
condition checking in templatized implementations. **/
struct TrueType {};

/** This is an operator for and-ing compile-time truth types. */
template <class L, class R> struct AndOpType {};
template<> struct AndOpType<FalseType,FalseType> {typedef FalseType Result;};
template<> struct AndOpType<FalseType,TrueType>  {typedef FalseType Result;};
template<> struct AndOpType<TrueType, FalseType> {typedef FalseType Result;};
template<> struct AndOpType<TrueType, TrueType>  {typedef TrueType  Result;};

/** This is an operator for or-ing compile-time truth types. */
template <class L, class R> struct OrOpType {};
template<> struct OrOpType<FalseType,FalseType> {typedef FalseType Result;};
template<> struct OrOpType<FalseType,TrueType>  {typedef TrueType  Result;};
template<> struct OrOpType<TrueType, FalseType> {typedef TrueType  Result;};
template<> struct OrOpType<TrueType, TrueType>  {typedef TrueType  Result;};

/** This is an operator for exclusive or-ing compile-time truth types. */
template <class L, class R> struct XorOpType {};
template<> struct XorOpType<FalseType,FalseType> {typedef FalseType Result;};
template<> struct XorOpType<FalseType,TrueType>  {typedef TrueType  Result;};
template<> struct XorOpType<TrueType, FalseType> {typedef TrueType  Result;};
template<> struct XorOpType<TrueType, TrueType>  {typedef FalseType Result;};

/** Compile-time type test: is this one of the built-in integral types?. **/
template <class T> struct IsIntegralType {
    /** This typedef is TrueType if the template type T is an integral type;
    otherwise it is FalseType. **/
    typedef FalseType Result;
    /** This compile-time constant bool is true if the template type T is an
    integral type otherwise it is false. **/
    static const bool result = false;
};
/** This macro must be invoked once for each of the built-in integral types to
specialize the IsIntegralType struct template for those types. **/
#define SimTK_SPECIALIZE_INTEGRAL_TYPE(T)       \
    template<> struct IsIntegralType<T>         \
    {typedef TrueType Result; static const bool result = true;}

SimTK_SPECIALIZE_INTEGRAL_TYPE(bool); 
SimTK_SPECIALIZE_INTEGRAL_TYPE(char);
// This causes problems when used with Qt which for some crazy
// reason likes to make its own wchar_t rather than using the built in.
// SimTK_SPECIALIZE_INTEGRAL_TYPE(wchar_t);
SimTK_SPECIALIZE_INTEGRAL_TYPE(signed char);
SimTK_SPECIALIZE_INTEGRAL_TYPE(unsigned char);
SimTK_SPECIALIZE_INTEGRAL_TYPE(short);
SimTK_SPECIALIZE_INTEGRAL_TYPE(unsigned short);
SimTK_SPECIALIZE_INTEGRAL_TYPE(int);
SimTK_SPECIALIZE_INTEGRAL_TYPE(unsigned int); // a.k.a. "unsigned"
SimTK_SPECIALIZE_INTEGRAL_TYPE(long);
SimTK_SPECIALIZE_INTEGRAL_TYPE(unsigned long);
SimTK_SPECIALIZE_INTEGRAL_TYPE(long long);
SimTK_SPECIALIZE_INTEGRAL_TYPE(unsigned long long);

/** Compile-time type test: is this one of the built-in floating point types?. **/
template <class T> struct IsFloatingType {
    /** This typedef is TrueType if the template type T is a floating point type;
    otherwise it is FalseType. **/
    typedef FalseType Result;
    /** This compile-time constant bool is true if the template type T is a
    floating point type otherwise it is false. **/
    static const bool result = false;
};
/** This macro must be invoked once for each of the built-in floating point 
types to specialize the IsFloatingType struct template for those types. **/
#define SimTK_SPECIALIZE_FLOATING_TYPE(T)       \
    template<> struct IsFloatingType<T>         \
    {typedef TrueType Result; static const bool result = true;}

SimTK_SPECIALIZE_FLOATING_TYPE(float); 
SimTK_SPECIALIZE_FLOATING_TYPE(double); 

/** Compile-time type test: is this the void type?. **/
template <class T> struct IsVoidType {
    /** This typedef is TrueType if the template type T is "void";
    otherwise it is FalseType. **/
    typedef FalseType Result;
    /** This compile-time constant bool is true if the template type T is
    "void" otherwise it is false. **/
    static const bool result = false;
};
template<> struct IsVoidType<void> 
{typedef TrueType Result; static const bool result = true;};

/** Compile-time test: is this one of the built-in "arithmetic" types, meaning
an integral or floating type? **/
template <class T> struct IsArithmeticType {
    /** This typedef is TrueType if the template type T is one of the integral;
    or floating point types, otherwise it is FalseType. **/
    typedef OrOpType<typename IsIntegralType<T>::Result,
                     typename IsFloatingType<T>::Result>    Result;
    /** This compile-time constant bool is true if the template type T is
    one of the integral or floating point types, otherwise it is false. **/
    static const bool result = IsIntegralType<T>::result 
                            || IsFloatingType<T>::result;
};

// This struct's sole use is to allow us to define the typedef 
// Is64BitPlatformType as equivalent to either TrueType or FalseType.
template <bool is64Bit> struct Is64BitHelper {};
template<> struct Is64BitHelper<true>  
{typedef TrueType  Result; static const bool result = true;};
template<> struct Is64BitHelper<false> 
{typedef FalseType Result; static const bool result = false;};

/** Compile-time test: this typedef will be TrueType if this is a 64-bit 
platform, meaning that the size of a pointer is the same as the size of a 
long long; otherwise it will be FalseType and we have a 32-bit platform meaning
that the size of a pointer is the same as an int. **/
// We use a constexpr function to avoid a bug in SWIG.
constexpr bool detect64BitPlatform() { return (sizeof(size_t) > sizeof(int)); }
static const bool Is64BitPlatform = detect64BitPlatform();
typedef Is64BitHelper<Is64BitPlatform>::Result Is64BitPlatformType;


/** Attempt to demangle a type name as returned by typeid.name(), with the
result hopefully suitable for meaningful display to a human. Behavior is 
compiler-dependent. 
@relates SimTK::NiceTypeName **/
SimTK_SimTKCOMMON_EXPORT 
std::string demangle(const char* name);

/** Given a compiler-dependent demangled type name string as returned by 
SimTK::demangle(), attempt to form a canonicalized representation that will be
the same for any compiler. Unnecessary spaces and superfluous keywords like
"class" and "struct" are removed. The `namestr()` method of NiceTypeName\<T>
uses this function to produce a human-friendly type name that is the same on any
platform. The input argument is left empty. 
@relates SimTK::NiceTypeName **/
SimTK_SimTKCOMMON_EXPORT 
std::string canonicalizeTypeName(std::string&& demangledTypeName);

/** Same, but takes an lvalue reference so has to copy the input. 
@relates SimTK::NiceTypeName **/
inline std::string canonicalizeTypeName(const std::string& demangledTypeName)
{   return canonicalizeTypeName(std::string(demangledTypeName)); }

/** Given a canonicalized type name, produce a modified version that is 
better-suited to use as an XML attribute. This means replacing the angle
brackets with curly braces to avoid trouble. The input argument is left
empty. 
@relates SimTK::NiceTypeName **/
SimTK_SimTKCOMMON_EXPORT
std::string encodeTypeNameForXML(std::string&& canonicalizedTypeName);

/** Same, but takes an lvalue reference so has to copy the input. 
@relates SimTK::NiceTypeName **/
inline std::string encodeTypeNameForXML(const std::string& niceTypeName)
{   return encodeTypeNameForXML(std::string(niceTypeName)); }

/** Given a type name that was encoded for XML by SimTK::encodeTypeNameForXML(),
restore it to its canonicalized form. This means replacing curly braces by
angle brackets. The input argument is left empty. 
@relates SimTK::NiceTypeName **/
SimTK_SimTKCOMMON_EXPORT
std::string decodeXMLTypeName(std::string&& xmlTypeName);

/** Same, but takes an lvalue reference so has to copy the input. 
@relates SimTK::NiceTypeName **/
inline std::string decodeXMLTypeName(const std::string& xmlTypeName)
{   return decodeXMLTypeName(std::string(xmlTypeName)); }

/** Obtain human-readable and XML-usable names for arbitrarily-complicated
C++ types. Three methods `name()`, `namestr()`, and `xmlstr()` are provided
giving respectively the compiler-dependent output from `typeid(T).%name()`, 
a canonicalized human-readable string, and the canonicalized string with
XML-forbidden angle brackets replaced by curly braces. The default 
implementation is usable for most types, but if you don't like the result you 
can specialize to provide nicer names. For example, you may prefer SimTK::Vec3 
to SimTK::Vec\<3,double,1>.

@warning Don't expect usable names for types that are defined in an anonymous 
namespace or for function-local types. Names will still be produced but they 
won't be unique and won't necessarily be compiler-independent.

The output of `namestr()` is typically used for error messages and testing;
`xmlstr()` is used for type tags in XML for use in deserializing. **/
template <class T> struct NiceTypeName {
    /** The default implementation of name() here returns the raw result from
    `typeid(T).%name()` which will be fast but may be a mangled name in some 
    compilers (gcc and clang included). **/
    static const char* name() {return typeid(T).name();}
    /** The default implementation of namestr() attempts to return a nicely
    demangled and canonicalized type name on all platforms, using the 
    SimTK::demangle() and SimTK::canonicalizeTypeName() methods. This is an
    expensive operation but is only done once. **/
    static const std::string& namestr() {
        static const std::string canonical = 
            canonicalizeTypeName(demangle(name()));
        return canonical;
    }
    /** The default implementation of xmlstr() takes the output of namestr()
    and invokes SimTK::encodeTypeNameForXML() on it. **/
    static const std::string& xmlstr() {
        static const std::string xml = encodeTypeNameForXML(namestr());
        return xml;
    }
};

} // namespace SimTK

/** This specializes the name of a type to be exactly the text you use to
specify it, rather than whatever ugly thing might result on different platforms
from resolution of typedefs, default template arguments, etc. Note that this
macro generates a template specialization that must be done in the SimTK
namespace; consequently it opens and closes namespace SimTK and must not
be invoked if you already have that namespace open. **/
#define SimTK_NICETYPENAME_LITERAL(T)                                   \
namespace SimTK {                                                       \
template <> struct NiceTypeName< T > {                                  \
    static const char* name() { return #T; }                            \
    static const std::string& namestr() {                               \
        static const std::string str(#T);                               \
        return str;                                                     \
    }                                                                   \
    static const std::string& xmlstr() {                                \
        static const std::string xml = encodeTypeNameForXML(namestr()); \
        return xml;                                                     \
    }                                                                   \
};                                                                      \
}

// Some types for which we'd like to see nice type names.
SimTK_NICETYPENAME_LITERAL(bool);            
SimTK_NICETYPENAME_LITERAL(char); 
// This causes problems when used with Qt which for some crazy
// reason likes to make its own wchar_t rather than using the built in.
// SimTK_NICETYPENAME_LITERAL(wchar_t);            
SimTK_NICETYPENAME_LITERAL(signed char); 
SimTK_NICETYPENAME_LITERAL(unsigned char);
SimTK_NICETYPENAME_LITERAL(short);           
SimTK_NICETYPENAME_LITERAL(unsigned short);  
SimTK_NICETYPENAME_LITERAL(int); 
SimTK_NICETYPENAME_LITERAL(unsigned); // preferred to "unsigned int"
SimTK_NICETYPENAME_LITERAL(long);            
SimTK_NICETYPENAME_LITERAL(unsigned long);   
SimTK_NICETYPENAME_LITERAL(long long);
SimTK_NICETYPENAME_LITERAL(unsigned long long);
SimTK_NICETYPENAME_LITERAL(float);           
SimTK_NICETYPENAME_LITERAL(double); 
SimTK_NICETYPENAME_LITERAL(long double);
SimTK_NICETYPENAME_LITERAL(std::string);
SimTK_NICETYPENAME_LITERAL(std::complex<float>);
SimTK_NICETYPENAME_LITERAL(std::complex<double>); 
SimTK_NICETYPENAME_LITERAL(std::complex<long double>); 
SimTK_NICETYPENAME_LITERAL(SimTK::FalseType);
SimTK_NICETYPENAME_LITERAL(SimTK::TrueType); 


#endif /* C++ stuff */

#endif /* SimTK_SimTKCOMMON_COMMON_H_ */
