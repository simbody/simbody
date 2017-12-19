#ifndef SimTK_SIMMATRIX_SCALAR_H_
#define SimTK_SIMMATRIX_SCALAR_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
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

/** @file
 * This is a user-includable header which includes everything needed
 * to make use of SimMatrix Scalar code. More commonly, this will be
 * included from within Matrix code.
 */

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/Constants.h"
#include "SimTKcommon/internal/Exception.h"
#include "SimTKcommon/internal/ExceptionMacros.h"
#include "SimTKcommon/internal/String.h"
#include "SimTKcommon/internal/Array.h"

#include "SimTKcommon/internal/conjugate.h"
#include "SimTKcommon/internal/CompositeNumericalTypes.h"
#include "SimTKcommon/internal/NTraits.h"
#include "SimTKcommon/internal/negator.h"

#include <complex>
#include <cmath>
#include <climits>
#include <cassert>
#include <utility>

namespace SimTK {

    /////////////////////////////////////////
    // Handy default-precision definitions //
    /////////////////////////////////////////

typedef conjugate<Real> Conjugate;  // like Complex

/** @defgroup TypedNumConstants  Numerical Constants with Types
    @ingroup  PredefinedConstants

This is a set of predefined constants in the form of Real (%SimTK default 
precision) symbols that are important for writing precision-independent
numerical algorithms. These constants have memory addresses (that is, they
are not macros), so you can return references to them. These are global
external symbols rather than static members to avoid problems with static 
initialization order.

Constants defined here include common mathematical values like pi, e, and
sqrt(2) and also numerical constants related to the floating point
implementation of the Real type, such as NaN, Infinity, and the machine
precision Epsilon. For convenience we also provide several common numerical
values for which it is useful to have a precision-independent representation
(mostly to avoid warnings or casts to avoid them), and also where it is 
useful to have a referenceable memory location that contains those values.
These include small integers and common small fractions like 1/2 and 1/3.

Note that the Simbody convention for typed constants is to name them like 
ordinary variables except with an initial capital letter (like a class name). 
Typed constants are processed by the compiler rather than the preprocessor and 
do not require any special treatment when used; they behave just like variables 
of the same type and value would behave. So we don't feel the need to draw 
attention to them with ALL_CAPS like we do with preprocessor symbols. **/
/**@{**/

/** This is the IEEE "not a number" constant for this implementation of
the default-precision Real type; be very careful using this because it has
many strange properties such as not comparing equal to itself. You must use
the SimTK::isNaN() function instead to determine whether something contains
a NaN value. **/
extern SimTK_SimTKCOMMON_EXPORT const Real NaN; 
/** This is the IEEE positive infinity constant for this implementation of
the default-precision Real type; -Infinity will produce the negative infinity
constant. Infinity tests larger than any other Real value. **/
extern SimTK_SimTKCOMMON_EXPORT const Real Infinity;

/** Epsilon is the size of roundoff noise; it is the smallest positive number
of default-precision type Real such that 1+Eps != 1. If Real is double (the
normal case) then Eps ~= 1e-16; if Real is float then Eps ~= 1e-7. **/
extern SimTK_SimTKCOMMON_EXPORT const Real Eps;
/** This is the square root of Eps, ~1e-8 if Real is double, ~3e-4 if Real
is float. Many numerical algorithms are limited to accuracy of sqrt(Eps)
so this constant is useful in checking for termination of them. **/
extern SimTK_SimTKCOMMON_EXPORT const Real SqrtEps;
/** TinyReal is a floating point value smaller than the floating point
precision; it is defined as Eps^(5/4) which is ~1e-20 for Real==double and
~1e-9 for float. This is commonly used as a number to add to a computation in 
a denominator (such as a vector length) that might come out zero, just for
the purpose of avoiding a divide by zero. **/
extern SimTK_SimTKCOMMON_EXPORT const Real TinyReal; 
/** SignificantReal is the smallest value that we consider to be clearly 
distinct from roundoff error when it is the result of a computation;
it is defined as Eps^(7/8) which is ~1e-14 when Real==double, ~1e-6 when
Real==float. **/
extern SimTK_SimTKCOMMON_EXPORT const Real SignificantReal; 

/** This is the smallest positive real number that can be expressed in the
type Real; it is ~1e-308 when Real==double, ~1e-38 when Real==float. **/
extern SimTK_SimTKCOMMON_EXPORT const Real LeastPositiveReal; 
/** This is the largest finite positive real number that can be expressed in 
the Real type; ~1e+308 when Real==double, ~1e+38 when Real==float.\ Note
that there is also a value Infinity that will test larger than this one. **/
extern SimTK_SimTKCOMMON_EXPORT const Real MostPositiveReal;  
/** This is the largest negative real number (that is, closest to zero) that
can be expressed in values of type Real. **/
extern SimTK_SimTKCOMMON_EXPORT const Real LeastNegativeReal;
/** This is the smallest finite negative real number that
can be expressed in values of type Real.\ Note that -Infinity is a value
that will still test smaller than this one. **/
extern SimTK_SimTKCOMMON_EXPORT const Real MostNegativeReal;

/** This is the number of decimal digits that can be reliably stored and
retrieved in the default Real precision (typically log10(1/eps)-1), that is,
about 15 digits when Real==double and 6 digits when Real==float. **/
extern SimTK_SimTKCOMMON_EXPORT const int NumDigitsReal; 

/** This is the smallest number of decimal digits you should store in a text
file if you want to be able to get exactly the same bit pattern back when you 
read it back in and convert the text to a Real value. Typically, this is about
log10(1/tiny), which is about 20 digits when Real==double and 9 digits when
Real==float. **/
extern SimTK_SimTKCOMMON_EXPORT const int LosslessNumDigitsReal; // double ~20, 
                                                                 // float ~9

    // Carefully calculated constants, with convenient memory addresses.

extern SimTK_SimTKCOMMON_EXPORT const Real Zero;        ///< Real(0)
extern SimTK_SimTKCOMMON_EXPORT const Real One;         ///< Real(1)
extern SimTK_SimTKCOMMON_EXPORT const Real MinusOne;    ///< Real(-1)
extern SimTK_SimTKCOMMON_EXPORT const Real Two;         ///< Real(2)
extern SimTK_SimTKCOMMON_EXPORT const Real Three;       ///< Real(3)

extern SimTK_SimTKCOMMON_EXPORT const Real OneHalf;     ///< Real(1)/2
extern SimTK_SimTKCOMMON_EXPORT const Real OneThird;    ///< Real(1)/3
extern SimTK_SimTKCOMMON_EXPORT const Real OneFourth;   ///< Real(1)/4
extern SimTK_SimTKCOMMON_EXPORT const Real OneFifth;    ///< Real(1)/5
extern SimTK_SimTKCOMMON_EXPORT const Real OneSixth;    ///< Real(1)/6
extern SimTK_SimTKCOMMON_EXPORT const Real OneSeventh;  ///< Real(1)/7
extern SimTK_SimTKCOMMON_EXPORT const Real OneEighth;   ///< Real(1)/8
extern SimTK_SimTKCOMMON_EXPORT const Real OneNinth;    ///< Real(1)/9
extern SimTK_SimTKCOMMON_EXPORT const Real Pi;          ///< Real(pi)
extern SimTK_SimTKCOMMON_EXPORT const Real OneOverPi;   ///< 1/Real(pi)
extern SimTK_SimTKCOMMON_EXPORT const Real E;           ///< \e e = Real(exp(1))
extern SimTK_SimTKCOMMON_EXPORT const Real Log2E;       ///< Real(log2(\e e)) (log base 2)
extern SimTK_SimTKCOMMON_EXPORT const Real Log10E;      ///< Real(log10(\e e)) (log base 10)
extern SimTK_SimTKCOMMON_EXPORT const Real Sqrt2;       ///< Real(sqrt(2))
extern SimTK_SimTKCOMMON_EXPORT const Real OneOverSqrt2;///< 1/sqrt(2)==sqrt(2)/2 as Real
extern SimTK_SimTKCOMMON_EXPORT const Real Sqrt3;       ///< Real(sqrt(3))
extern SimTK_SimTKCOMMON_EXPORT const Real OneOverSqrt3;///< Real(1/sqrt(3))
extern SimTK_SimTKCOMMON_EXPORT const Real CubeRoot2;   ///< Real(2^(1/3)) (cube root of 2)
extern SimTK_SimTKCOMMON_EXPORT const Real CubeRoot3;   ///< Real(3^(1/3)) (cube root of 3)
extern SimTK_SimTKCOMMON_EXPORT const Real Ln2;         ///< Real(ln(2)) (natural log of 2)
extern SimTK_SimTKCOMMON_EXPORT const Real Ln10;        ///< Real(ln(10)) (natural log of 10)
       
/** We only need one complex constant, \e i = sqrt(-1).\ For the rest 
just multiply the real constant by \e i, or convert with 
Complex(the Real constant), or if you need an address you can use 
NTraits<Complex>::getPi(), etc. 
@see NTraits **/
extern SimTK_SimTKCOMMON_EXPORT const Complex I;
/**@}**/

    ///////////////////////////
    // SOME SCALAR UTILITIES //
    ///////////////////////////


/**
 * @defgroup atMostOneBitIsSet atMostOneBitIsSet()
 * @ingroup BitFunctions
 *
 * atMostOneBitIsSet(i) provides an extremely fast way to determine whether
 * an integral type is either zero or consists of a single set bit. This
 * question arises when using bits to represent set membership where one
 * may wish to verify that an integer represents a single element rather
 * than a set of elements.
 *
 * @see exactlyOneBitIsSet()
 */
//@{
// We use the trick that v & v-1 returns the value that is v with its
// rightmost bit cleared (if it has a rightmost bit set).
inline bool atMostOneBitIsSet(unsigned char v)      {return (v&(v-1))==0;}
inline bool atMostOneBitIsSet(unsigned short v)     {return (v&(v-1))==0;}
inline bool atMostOneBitIsSet(unsigned int v)       {return (v&(v-1))==0;}
inline bool atMostOneBitIsSet(unsigned long v)      {return (v&(v-1))==0;}
inline bool atMostOneBitIsSet(unsigned long long v) {return (v&(v-1))==0;}
inline bool atMostOneBitIsSet(signed char v)        {return (v&(v-1))==0;}
inline bool atMostOneBitIsSet(char v)               {return (v&(v-1))==0;}
inline bool atMostOneBitIsSet(short v)              {return (v&(v-1))==0;}
inline bool atMostOneBitIsSet(int v)                {return (v&(v-1))==0;}
inline bool atMostOneBitIsSet(long v)               {return (v&(v-1))==0;}
inline bool atMostOneBitIsSet(long long v)          {return (v&(v-1))==0;}
//@}

/**
 * @defgroup exactlyOneBitIsSet exactlyOneBitIsSet()
 * @ingroup BitFunctions
 *
 * exactlyOneBitIsSet(i) provides a very fast way to determine whether
 * an integral type has exactly one bit set. For unsigned and positive
 * signed values, this is equivalent to the value being a power of two.
 * Note that negative powers of two are <em>not</em> represented with
 * a single bit set -- negate it first if you want to use this routine
 * to determine if a signed value is a negative power of two.
 *
 * @see atMostOneBitIsSet()
 */
//@{
inline bool exactlyOneBitIsSet(unsigned char v)      {return v && atMostOneBitIsSet(v);}
inline bool exactlyOneBitIsSet(unsigned short v)     {return v && atMostOneBitIsSet(v);}
inline bool exactlyOneBitIsSet(unsigned int v)       {return v && atMostOneBitIsSet(v);}
inline bool exactlyOneBitIsSet(unsigned long v)      {return v && atMostOneBitIsSet(v);}
inline bool exactlyOneBitIsSet(unsigned long long v) {return v && atMostOneBitIsSet(v);}
inline bool exactlyOneBitIsSet(signed char v)        {return v && atMostOneBitIsSet(v);}
inline bool exactlyOneBitIsSet(char v)               {return v && atMostOneBitIsSet(v);}
inline bool exactlyOneBitIsSet(short v)              {return v && atMostOneBitIsSet(v);}
inline bool exactlyOneBitIsSet(int v)                {return v && atMostOneBitIsSet(v);}
inline bool exactlyOneBitIsSet(long v)               {return v && atMostOneBitIsSet(v);}
inline bool exactlyOneBitIsSet(long long v)          {return v && atMostOneBitIsSet(v);}
//@}

/**
 * @defgroup SignBitGroup signBit()
 * @ingroup BitFunctions
 *
 * signBit(i) provides a fast way to determine the value of the sign
 * bit (as a bool) for integral and floating types. Note that this is
 * significantly different than sign(x); be sure you know what you're doing
 * if you use this method. signBit() refers to the underlying representation
 * rather than the numerical value. For example, for floating types there
 * are two zeroes, +0 and -0 which have opposite sign bits but the same
 * sign (0). Also, unsigned types have sign() of 0 or 1, but they are 
 * considered here to have a sign bit of 0 always, since the stored high
 * bit does not indicate a negative value.
 *
 * \b Notes
 *  - signBit() is overloaded for 'signed char' and 'unsigned char', but not
 *    for plain 'char' because we don't know whether to interpret the high
 *    bit as a sign in that case (because the C++ standard leaves it unspecified).
 *  - complex and conjugate numbers do not have sign bits.
 *  - negator<float> and negator<double> have the \e same sign bit as the
 *    underlying representation -- it's up to you to realize that it is
 *    interpreted differently!
 *
 * @see sign()
 */
/*@{*/
// Don't remove these unused formal parameter names 'u'; doxygen barfs.
inline bool signBit(unsigned char u)      {return false;}
inline bool signBit(unsigned short u)     {return false;}
inline bool signBit(unsigned int u)       {return false;}
inline bool signBit(unsigned long u)      {return false;}
inline bool signBit(unsigned long long u) {return false;}

// Note that plain 'char' type is not overloaded -- see above.

// We're assuming sizeof(char)==1, short==2, int==4, long long==8 which is safe
// for all our anticipated platforms. But some 64 bit implementations have
// sizeof(long)==4 while others have sizeof(long)==8 so we'll use the ANSI
// standard value LONG_MIN which is a long value with just the high bit set.
// (We're assuming two's complement integers here; I haven't seen anything
// else in decades.)
inline bool signBit(signed char i) {return ((unsigned char)i      & 0x80U) != 0;}
inline bool signBit(short       i) {return ((unsigned short)i     & 0x8000U) != 0;}
inline bool signBit(int         i) {return ((unsigned int)i       & 0x80000000U) != 0;}
inline bool signBit(long long   i) {return ((unsigned long long)i & 0x8000000000000000ULL) != 0;}

inline bool signBit(long        i) {return ((unsigned long)i
                                            & (unsigned long)LONG_MIN) != 0;}

inline bool signBit(const float&  f) {return std::signbit(f);}
inline bool signBit(const double& d) {return std::signbit(d);}
inline bool signBit(const negator<float>&  nf) {return std::signbit(-nf);} // !!
inline bool signBit(const negator<double>& nd) {return std::signbit(-nd);} // !!

/*@}*/

/**
 * @defgroup SignGroup      sign()
 * @ingroup ScalarFunctions
 *
 * s=sign(n) returns int -1,0,1 according to n<0, n==0, n>0 for any integer
 * or real numeric type, unsigned 0 or 1 for any unsigned argument. This
 * routine is specialized for each of the int, unsigned, and real types
 * of all sizes. Sign is defined for "signed char" and "unsigned char" but
 * not plain "char" since the language leaves unspecified whether that is
 * a signed or unsigned type. Sign is not defined for complex or conjugate.
 */
//@{
inline unsigned int sign(unsigned char      u) {return u==0 ? 0 : 1;}
inline unsigned int sign(unsigned short     u) {return u==0 ? 0 : 1;}
inline unsigned int sign(unsigned int       u) {return u==0 ? 0 : 1;}
inline unsigned int sign(unsigned long      u) {return u==0 ? 0 : 1;}
inline unsigned int sign(unsigned long long u) {return u==0 ? 0 : 1;}

// Don't overload for plain "char" because it may be signed or unsigned
// depending on the compiler.

inline int sign(signed char i) {return i>0 ? 1 : (i<0 ? -1 : 0);}
inline int sign(short       i) {return i>0 ? 1 : (i<0 ? -1 : 0);}
inline int sign(int         i) {return i>0 ? 1 : (i<0 ? -1 : 0);}
inline int sign(long        i) {return i>0 ? 1 : (i<0 ? -1 : 0);}
inline int sign(long long   i) {return i>0 ? 1 : (i<0 ? -1 : 0);}

inline int sign(const float&       x) {return x>0 ? 1 : (x<0 ? -1 : 0);}
inline int sign(const double&      x) {return x>0 ? 1 : (x<0 ? -1 : 0);}

inline int sign(const negator<float>&       x) {return -sign(-x);} // -x is free
inline int sign(const negator<double>&      x) {return -sign(-x);}
//@}

/**
 * @defgroup square square()
 * @ingroup ScalarFunctions
 *
 * y=square(x) returns the square of the argument for any numeric type.
 * We promise to evaluate x only once. We assume is is acceptable for
 * the result type to be the same as the argument type; if it won't
 * fit caller must cast argument to a wider type first. This is an
 * inline routine which will run as fast as an explicit multiply
 * (x*x) in optimized code, and somewhat faster for complex and 
 * conjugate types (5 flops instead of the usual 6).
 *
 * Squaring a negated number loses the negator at no cost; squaring
 * a conjugate number returns a complex result at no additional cost.
 */
//@{
inline unsigned char  square(unsigned char  u) {return u*u;}
inline unsigned short square(unsigned short u) {return u*u;}
inline unsigned int   square(unsigned int   u) {return u*u;}
inline unsigned long  square(unsigned long  u) {return u*u;}
inline unsigned long long square(unsigned long long u) {return u*u;}

inline char        square(char c) {return c*c;}

inline signed char square(signed char i) {return i*i;}
inline short       square(short       i) {return i*i;}
inline int         square(int         i) {return i*i;}
inline long        square(long        i) {return i*i;}
inline long long   square(long long   i) {return i*i;}

inline float       square(const float&       x) {return x*x;}
inline double      square(const double&      x) {return x*x;}

// Negation is free for negators, so we can square them and clean
// them up at the same time at no extra cost.
inline float       square(const negator<float>&       x) {return square(-x);}
inline double      square(const negator<double>&      x) {return square(-x);}

// It is safer to templatize using complex classes, and doesn't make
// debugging any worse since complex is already templatized. 
// 5 flops vs. 6 for general complex multiply.
template <class P> inline 
std::complex<P> square(const std::complex<P>& x) {
    const P re=x.real(), im=x.imag();
    return std::complex<P>(re*re-im*im, 2*re*im);
}

// We can square a conjugate and clean it up back to complex at
// the same time at zero cost (or maybe 1 negation depending
// on how the compiler handles the "-2" below).
template <class P> inline 
std::complex<P> square(const conjugate<P>& x) {
    const P re=x.real(), nim=x.negImag();
    return std::complex<P>(re*re-nim*nim, -2*re*nim);
}

template <class P> inline
std::complex<P> square(const negator< std::complex<P> >& x) {
    return square(-x); // negation is free for negators
}

// Note return type here after squaring negator<conjugate>
// is complex, not conjugate.
template <class P> inline
std::complex<P> square(const negator< conjugate<P> >& x) {
    return square(-x); // negation is free for negators
}
//@}

/**
 * @defgroup cube cube()
 * @ingroup ScalarFunctions
 *
 * y=cube(x) returns the cube of the argument for any numeric type,
 * integral or floating point. We promise to evaluate x only once.
 * We assume is is acceptable for the result type to be the same
 * as the argument type; if it won't fit caller must cast argument
 * to a wider type first. This is an inline routine which will run
 * as fast as explicit multiplies (x*x*x) in optimized code, and
 * significantly faster for complex or conjugate types (8 flops
 * vs. 11).
 *
 * Cubing a negated real type returns a negated result. Cubing
 * a negated complex or conjugate returns a non-negated complex
 * result since that can be done with no additional cost.
 */
//@{
inline unsigned char  cube(unsigned char  u) {return u*u*u;}
inline unsigned short cube(unsigned short u) {return u*u*u;}
inline unsigned int   cube(unsigned int   u) {return u*u*u;}
inline unsigned long  cube(unsigned long  u) {return u*u*u;}
inline unsigned long long cube(unsigned long long u) {return u*u*u;}

inline char        cube(char c) {return c*c*c;}

inline signed char cube(signed char i) {return i*i*i;}
inline short       cube(short       i) {return i*i*i;}
inline int         cube(int         i) {return i*i*i;}
inline long        cube(long        i) {return i*i*i;}
inline long long   cube(long long   i) {return i*i*i;}

inline float       cube(const float&       x) {return x*x*x;}
inline double      cube(const double&      x) {return x*x*x;}

// To keep this cheap we'll defer getting rid of the negator<> until
// some other operation. We cube -x and then recast that to negator<>
// on return, for a total cost of 2 flops.
inline negator<float> cube(const negator<float>& x) {
    return negator<float>::recast(cube(-x));
}
inline negator<double> cube(const negator<double>& x) {
    return negator<double>::recast(cube(-x));
}

// Cubing a complex this way is cheaper than doing it by
// multiplication. Cost here is 8 flops vs. 11 for a square
// followed by a multiply.
template <class P> inline
std::complex<P> cube(const std::complex<P>& x) {
    const P re=x.real(), im=x.imag(), rr=re*re, ii=im*im;
    return std::complex<P>(re*(rr-3*ii), im*(3*rr-ii));
}

// Cubing a negated complex allows us to cube and eliminate the
// negator at the same time for zero extra cost. Compare the 
// expressions here to the normal cube above to see the free
// sign changes in both parts. 8 flops here.
template <class P> inline
std::complex<P> cube(const negator< std::complex<P> >& x) {
    // -x is free for a negator
    const P nre=(-x).real(), nim=(-x).imag(), rr=nre*nre, ii=nim*nim;
    return std::complex<P>(nre*(3*ii-rr), nim*(ii-3*rr));
}

// Cubing a conjugate this way saves a lot over multiplying, and
// also lets us convert the result to complex for free.
template <class P> inline
std::complex<P> cube(const conjugate<P>& x) {
    const P re=x.real(), nim=x.negImag(), rr=re*re, ii=nim*nim;
    return std::complex<P>(re*(rr-3*ii), nim*(ii-3*rr));
}


// Cubing a negated conjugate this way saves a lot over multiplying, and
// also lets us convert the result to non-negated complex for free.
template <class P> inline
std::complex<P> cube(const negator< conjugate<P> >& x) {
    // -x is free for a negator
    const P nre=(-x).real(), im=(-x).negImag(), rr=nre*nre, ii=im*im;
    return std::complex<P>(nre*(3*ii-rr), im*(3*rr-ii));
}
//@}

/** @defgroup ClampingGroup     clamp(), clampInPlace()
    @ingroup ScalarFunctions

Limit a numerical value so that it does not go outside a given range.
There are two functions (plus overloads) defined here: clamp() and 
clampInPlace(). The first one is the most commonly used and simply 
calculates an in-range value from an input value and a given range
[low,high]. clampInPlace() is given a reference to a variable and if
necessary modifies that variable so that its value is in the given range
[low,high]. Both functions are overloaded for all the integral and real 
types but are not defined for complex or conjugate types. 

The following examples shows how clamp() and clampInPlace() can be defined 
in terms of std::min() and std::max(). 

clamp():
@code
    const double low, high; // from somewhere
    const double v;
    double clampedValue;
    clampedValue = clamp(low,v,high); // clampedValue is in [low,high]
    clampedValue = std::min(std::max(low,v), high); // equivalent
@endcode
clampInPlace():
@code
    const double low, high; // from somewhere
    double v;
    clampInPlace(low,v,high); // now v is in [low,high]
    v = std::min(std::max(low,v), high); // equivalent
@endcode **/
/*@{*/

/** Check that low <= v <= high and modify v in place if necessary to 
bring it into that range.
@param[in]      low     The lower bound; must be <= high.
@param[in,out]  v       The variable whose value is changed if necessary.
@param[in]      high    The upper bound; must be >= low.
@return A writable reference to the now-possibly-modified input variable v.

This method is overloaded for all standard integral and floating point
types. All the arguments and the return type will be the same type; there
are no explicit overloads for mixed types so you may find you need to
cast the bounds in some cases to ensure you are calling the correct
overload.

Cost: These are very fast inline methods; the floating point ones use
just two flops. **/
inline double& clampInPlace(double low, double& v, double high) 
{   assert(low<=high); if (v<low) v=low; else if (v>high) v=high; return v; }
/** @copydoc SimTK::clampInPlace(double,double&,double) **/
inline float& clampInPlace(float low, float& v, float high) 
{   assert(low<=high); if (v<low) v=low; else if (v>high) v=high; return v; }

    // Floating point clamps with integer bounds; without these
    // explicit casts are required.

/** @copydoc SimTK::clampInPlace(double,double&,double) 
Takes integer bounds to avoid need for explicit casts. **/
inline double& clampInPlace(int low, double& v, int high) 
{   return clampInPlace((double)low,v,(double)high); }
/** @copydoc SimTK::clampInPlace(double,double&,double) 
Takes integer bounds to avoid need for explicit casts. **/
inline float& clampInPlace(int low, float& v, int high) 
{   return clampInPlace((float)low,v,(float)high); }

/** @copydoc SimTK::clampInPlace(double,double&,double) 
Takes an integer bound to avoid need for explicit casts. **/
inline double& clampInPlace(int low, double& v, double high) 
{   return clampInPlace((double)low,v,high); }
/** @copydoc SimTK::clampInPlace(double,double&,double) 
Takes an integer bound to avoid need for explicit casts. **/
inline float& clampInPlace(int low, float& v, float high) 
{   return clampInPlace((float)low,v,high); }

/** @copydoc SimTK::clampInPlace(double,double&,double) 
Takes an integer bound to avoid need for explicit casts. **/
inline double& clampInPlace(double low, double& v, int high) 
{   return clampInPlace(low,v,(double)high); }
/** @copydoc SimTK::clampInPlace(double,double&,double) 
Takes an integer bound to avoid need for explicit casts. **/
inline float& clampInPlace(float low, float& v, int high) 
{   return clampInPlace(low,v,(float)high); }

/** @copydoc SimTK::clampInPlace(double,double&,double) **/
inline unsigned char& clampInPlace(unsigned char low, unsigned char& v, unsigned char high) 
{   assert(low<=high); if (v<low) v=low; else if (v>high) v=high; return v; }
/** @copydoc SimTK::clampInPlace(double,double&,double) **/
inline unsigned short& clampInPlace(unsigned short low, unsigned short& v, unsigned short high) 
{   assert(low<=high); if (v<low) v=low; else if (v>high) v=high; return v; }
/** @copydoc SimTK::clampInPlace(double,double&,double) **/
inline unsigned int& clampInPlace(unsigned int low, unsigned int& v, unsigned int high) 
{   assert(low<=high); if (v<low) v=low; else if (v>high) v=high; return v; }
/** @copydoc SimTK::clampInPlace(double,double&,double) **/
inline unsigned long& clampInPlace(unsigned long low, unsigned long& v, unsigned long high) 
{   assert(low<=high); if (v<low) v=low; else if (v>high) v=high; return v; }
/** @copydoc SimTK::clampInPlace(double,double&,double) **/
inline unsigned long long& clampInPlace(unsigned long long low, unsigned long long& v, unsigned long long high) 
{   assert(low<=high); if (v<low) v=low; else if (v>high) v=high; return v; }

/** @copydoc SimTK::clampInPlace(double,double&,double) **/
inline char& clampInPlace(char low, char& v, char high) 
{   assert(low<=high); if (v<low) v=low; else if (v>high) v=high; return v; }
/** @copydoc SimTK::clampInPlace(double,double&,double) **/
inline signed char& clampInPlace(signed char low, signed char& v, signed char high) 
{   assert(low<=high); if (v<low) v=low; else if (v>high) v=high; return v; }

/** @copydoc SimTK::clampInPlace(double,double&,double) **/
inline short& clampInPlace(short low, short& v, short high) 
{   assert(low<=high); if (v<low) v=low; else if (v>high) v=high; return v; }
/** @copydoc SimTK::clampInPlace(double,double&,double) **/
inline int& clampInPlace(int low, int& v, int high) 
{   assert(low<=high); if (v<low) v=low; else if (v>high) v=high; return v; }
/** @copydoc SimTK::clampInPlace(double,double&,double) **/
inline long& clampInPlace(long low, long& v, long high) 
{   assert(low<=high); if (v<low) v=low; else if (v>high) v=high; return v; }
/** @copydoc SimTK::clampInPlace(double,double&,double) **/
inline long long& clampInPlace(long long low, long long& v, long long high) 
{   assert(low<=high); if (v<low) v=low; else if (v>high) v=high; return v; }

/** @copydoc SimTK::clampInPlace(double,double&,double) **/
inline negator<float>& clampInPlace(float low, negator<float>& v, float high) 
{   assert(low<=high); if (v<low) v=low; else if (v>high) v=high; return v; }
/** @copydoc SimTK::clampInPlace(double,double&,double) **/
inline negator<double>& clampInPlace(double low, negator<double>& v, double high) 
{   assert(low<=high); if (v<low) v=low; else if (v>high) v=high; return v; }



/** If v is in range low <= v <= high then return v, otherwise return
the nearest bound; this function does not modify the input variable v.
@param[in]      low     The lower bound; must be <= high.
@param[in]      v       The value to be put into range [low,high].
@param[in]      high    The upper bound; must be >= low.
@return Either the value v or one of the bounds.

This method is overloaded for all standard integral and floating point
types. All the arguments and the return type will be the same type; there
are no explicit overloads for mixed types so you may find you need to
cast the bounds in some cases to ensure you are calling the correct
overload.

@warning Viewed as a function of the input value v, clamp(v) is not
very smooth -- it is continuous (C0) but even its first derivative is
infinite at the bounds.

Cost: These are very fast, inline methods; the floating point ones use
just two flops.
@see clampInPlace() **/
inline double clamp(double low, double v, double high) 
{   return clampInPlace(low,v,high); }
/** @copydoc SimTK::clamp(double,double,double) **/
inline float clamp(float low, float v, float high) 
{   return clampInPlace(low,v,high); }

/** @copydoc SimTK::clamp(double,double,double) 
Takes integer bounds to avoid need for explicit casts. **/
inline double clamp(int low, double v, int high) 
{   return clampInPlace(low,v,high); }
/** @copydoc SimTK::clamp(double,double,double) 
Takes integer bounds to avoid need for explicit casts. **/
inline float clamp(int low, float v, int high) 
{   return clampInPlace(low,v,high); }

/** @copydoc SimTK::clamp(double,double,double) 
Takes an integer bound to avoid need for explicit casts. **/
inline double clamp(int low, double v, double high) 
{   return clampInPlace(low,v,high); }
/** @copydoc SimTK::clamp(double,double,double) 
Takes an integer bound to avoid need for explicit casts. **/
inline float clamp(int low, float v, float high) 
{   return clampInPlace(low,v,high); }

/** @copydoc SimTK::clamp(double,double,double) 
Takes an integer bound to avoid need for explicit casts. **/
inline double clamp(double low, double v, int high) 
{   return clampInPlace(low,v,high); }
/** @copydoc SimTK::clamp(double,double,double) 
Takes an integer bound to avoid need for explicit casts. **/
inline float clamp(float low, float v, int high) 
{   return clampInPlace(low,v,high); }

/** @copydoc SimTK::clamp(double,double,double) **/
inline unsigned char clamp(unsigned char low, unsigned char v, unsigned char high) 
{   return clampInPlace(low,v,high); }
/** @copydoc SimTK::clamp(double,double,double) **/
inline unsigned short clamp(unsigned short low, unsigned short v, unsigned short high) 
{   return clampInPlace(low,v,high); }
/** @copydoc SimTK::clamp(double,double,double) **/
inline unsigned int clamp(unsigned int low, unsigned int v, unsigned int high) 
{   return clampInPlace(low,v,high); }
/** @copydoc SimTK::clamp(double,double,double) **/
inline unsigned long clamp(unsigned long low, unsigned long v, unsigned long high) 
{   return clampInPlace(low,v,high); }
/** @copydoc SimTK::clamp(double,double,double) **/
inline unsigned long long clamp(unsigned long long low, unsigned long long v, unsigned long long high) 
{   return clampInPlace(low,v,high); }

/** @copydoc SimTK::clamp(double,double,double) **/
inline char clamp(char low, char v, char high) 
{   return clampInPlace(low,v,high); }
/** @copydoc SimTK::clamp(double,double,double) **/
inline signed char clamp(signed char low, signed char v, signed char high) 
{   return clampInPlace(low,v,high); }

/** @copydoc SimTK::clamp(double,double,double) **/
inline short clamp(short low, short v, short high) 
{   return clampInPlace(low,v,high); }
/** @copydoc SimTK::clamp(double,double,double) **/
inline int clamp(int low, int v, int high) 
{   return clampInPlace(low,v,high); }
/** @copydoc SimTK::clamp(double,double,double) **/
inline long clamp(long low, long v, long high) 
{   return clampInPlace(low,v,high); }
/** @copydoc SimTK::clamp(double,double,double) **/
inline long long clamp(long long low, long long v, long long high) 
{   return clampInPlace(low,v,high); }


// These aren't strictly necessary but are here to help the
// compiler find the right overload to call. These cost an
// extra flop because the negator<T> has to be cast to T which
// requires that the pending negation be performed. Note that
// the result types are not negated.

/** @copydoc SimTK::clamp(double,double,double)
Explicitly takes a negator<float> argument to help the compiler find the
right overload, but the negation is performed (1 extra flop) and the result 
type is an ordinary float. **/
inline float clamp(float low, negator<float> v, float high)
{   return clamp(low,(float)v,high); }
/** @copydoc SimTK::clamp(double,double,double)
Explicitly takes a negator<double> argument to help the compiler find the
right overload, but the negation is performed (1 extra flop) and the result 
type is an ordinary double. **/
inline double clamp(double low, negator<double> v, double high) 
{   return clamp(low,(double)v,high); }
/*@}*/



/** @defgroup SmoothedStepFunctions   Smoothed step functions
    @ingroup ScalarFunctions

Functions stepUp(), stepDown() and stepAny() provide smooth, S-shaped 
step functions that are useful for "softening" abrupt transitions between 
two values. Functions that return the first three derivatives of these
step functions are also provided.

y=stepUp(x) for x=0:1 returns a smooth, S-shaped step function y(x),
symmetric about the midpoint x=0.5, such that y(0)=0, y(1)=1, 
y'(0)=y''(0)=y'(1)=y''(1)=0, where the primes indicate differentiation 
with respect to x. No guarantees are made about the behavior of higher 
derivatives except to say that y''' does exist.

y=stepDown(x) for x=0:1 is a mirror image S curve that takes y
down from 1 to 0 smoothly as x goes from 0 to 1. As with stepUp()
y' and y'' are 0 at both ends, but here y(0)=1 and y(1)=0.

The stepAny() function is also available to make a step function
that takes y smoothly from y0 to y1 as x goes from x0 to x1.

We also provide functions that calculate the first three derivatives
of stepUp(), stepDown(), and stepAny(), namely: dstepUp(), d2stepUp(), 
d3stepUp(), and similarly dstepDown(), dstepAny() and so on. Note again 
that the third derivative is not guaranteed to be well behaved.

Costs of the current implementations of these inline functions are:
    - stepUp() 7 flops
    - dstepUp() 4 flops
    - d2stepUp() 6 flops
    - d3stepUp() 4 flops

The corresponding stepDown() methods cost one extra flop; the
stepAny() methods cost 6 extra flops (provided you do some precalculation;
see the stepAny() documentation). **/
/*@{*/

/** Interpolate smoothly from 0 up to 1 as the input argument goes from 0 to 1,
with first and second derivatives zero at either end of the interval.

@param[in]  x       The control parameter, in range [0,1].
@return The smoothed output value, in range [0,1] but with first and second 
derivatives smoothly approaching zero as x approaches either end of its 
interval.

See the documentation for stepAny() for a discussion about how to shift and
scale this function to produce arbitrary steps.

This function is overloaded for all the floating point precisions.
Cost is 7 flops.
@see stepDown(), stepAny() **/
inline double stepUp(double x)
{   assert(0 <= x && x <= 1);
    return x*x*x*(10+x*(6*x-15)); }  //10x^3-15x^4+6x^5


/** Interpolate smoothly from 1 down to 0 as the input argument goes from 
0 to 1, with first and second derivatives zero at either end of the interval.

@param[in]  x       The control parameter, in range [0,1].
@return The smoothed output value, in range [1,0] but with first and second 
derivatives smoothly approaching zero as x approaches either end of its 
interval.

See the documentation for stepAny() for a discussion about how to shift and
scale this function to produce arbitrary steps.

This function is overloaded for all the floating point precisions.
Cost is 8 flops.
@see stepUp(), stepAny() **/
inline double stepDown(double x) {return 1.0 -stepUp(x);}


/** Interpolate smoothly from y0 to y1 as the input argument goes from x0 to x1,
with first and second derivatives zero at either end of the interval.

@param[in]          y0
    The output value when x=x0.
@param[in]          yRange  
    The amount by which the output can change over the full interval, that is, 
    yRange=(y1-y0) where y1 is the value of the output when x=x1.
@param[in]          x0      
    The minimum allowable value for x.
@param[in]          oneOverXRange
    1/xRange, where xRange is the amount by which the input variable x can 
    change over its full interval, that is, xRange=(x1-x0) where x1 is the
    maximum allowable value for x.
@param[in]          x       
    The control parameter, in range [x0,x1]. This is often a time over
    which the output transition from y0 to y1 is to occur.
@return 
    The smoothed output value, in range [y0,y1] (where y1=y0+yRange) but with 
    first and second derivatives smoothly approaching zero as x approaches 
    either end of its interval.

Note that the desired curve is defined in terms of y0 and (y1-y0), and
x0 and 1/(x1-x0), rather than y0,y1,x0,x1 which would make for a nicer
calling signature. This is a concession to efficiency since the two ranges
are likely to be unchanged during many calls to this function and can thus
be precalculated. It wouldn't matter except that division is so expensive
(equivalent to 15-20 floating point operations). If you aren't concerned
about that (and in most cases it won't matter), you can call this
function using y0,y1,x0,x1 like this: @code
    y = stepAny(y0,y1-y0,x0,1/(x1-x0), x);
@endcode

Not counting the cost of calculating the ranges, each call to this function
requires 13 flops. Calculating the ranges in the argument list as shown
above raises the per-call cost to about 30 flops. However, there are many
common special cases that are much simpler. For example, if y is to go from
-1 to 1 while x goes from 0 to 1, then you can write: @code
    y = stepAny(-1,2,0,1, x);
@endcode
which is still only 13 flops despite the lack of pre-calculation.

@par Theory

stepUp() and stepDown() and their derivatives are defined only for arguments 
in the range 0 to 1. To create a general step function that smoothly 
interpolates from y=y0 to y1 while x goes from x0 to x1, and three derivatives,
we use the stepUp() functions like this:
@code
     ooxr = 1/(x1-x0);      // one over x range (signed)
     yr   = y1-y0;          // y range (signed)
     xadj = (x-x0)*ooxr     // x adjusted into 0:1 range (watch for roundoff)
     y    = y0 + yr * stepUp(xadj);
     dy   = yr*ooxr   * dstepUp(xadj);
     d2y  = yr*ooxr^2 * d2stepUp(xadj);
     d3y  = yr*ooxr^3 * d3stepUp(xadj);
@endcode

As a common special case, note that when y has a general range but x is
still in [0,1], the above simplifies considerably and you can save a few 
flops if you want by working with stepUp() or stepDown() directly. For 
example, in the common case where you want y to go from -1 to 1 as x goes 
from 0 to 1:
@code
    y   = stepAny(-1,2,0,1, x);  // y in [-1,1]; 13 flops
    y   = -1 + 2*stepUp(x);      // equivalent, saves 4 flops
    dy  = 2*dstepUp(x);          // these save 5 flops over dstepAny(), etc.
    d2y = 2*d2stepUp(x);
    d3y = 2*d3stepUp(x);
@endcode

It would be extremely rare for these few flops to matter at all; you should
almost always choose based on what looks better and/or is less error prone
instead. **/
inline double stepAny(double y0, double yRange,
                      double x0, double oneOverXRange,
                      double x) 
{   double xadj = (x-x0)*oneOverXRange;    
    assert(-NTraits<double>::getSignificant() <= xadj
           && xadj <= 1 + NTraits<double>::getSignificant());
    clampInPlace(0.0,xadj,1.0);
    return y0 + yRange*stepUp(xadj); }

/** First derivative of stepUp(): d/dx stepUp(x). 
@param[in] x    Control parameter in range [0,1]. 
@return First derivative of stepUp() at x. **/
inline double dstepUp(double x) {
    assert(0 <= x && x <= 1);
    const double xxm1=x*(x-1);
    return 30*xxm1*xxm1;        //30x^2-60x^3+30x^4
}

/** First derivative of stepDown(): d/dx stepDown(x). 
@param[in] x    Control parameter in range [0,1]. 
@return First derivative of stepDown() at x. **/
inline double dstepDown(double x) {return -dstepUp(x);}

/** First derivative of stepAny(): d/dx stepAny(x). 
See stepAny() for parameter documentation. 
@return First derivative of stepAny() at x. **/
inline double dstepAny(double yRange,
                       double x0, double oneOverXRange,
                       double x) 
{   double xadj = (x-x0)*oneOverXRange;    
    assert(-NTraits<double>::getSignificant() <= xadj
           && xadj <= 1 + NTraits<double>::getSignificant());
    clampInPlace(0.0,xadj,1.0);
    return yRange*oneOverXRange*dstepUp(xadj); }

/** Second derivative of stepUp(): d^2/dx^2 stepUp(x). 
@param[in] x    Control parameter in range [0,1]. 
@return Second derivative of stepUp() at x. **/
inline double d2stepUp(double x) {
    assert(0 <= x && x <= 1);
    return 60*x*(1+x*(2*x-3));  //60x-180x^2+120x^3
}

/** Second derivative of stepDown(): d^2/dx^2 stepDown(x). 
@param[in] x    Control parameter in range [0,1]. 
@return Second derivative of stepDown() at x. **/
inline double d2stepDown(double x) {return -d2stepUp(x);}

/** Second derivative of stepAny(): d^2/dx^2 stepAny(x). 
See stepAny() for parameter documentation. 
@return Second derivative of stepAny() at x. **/
inline double d2stepAny(double yRange,
                        double x0, double oneOverXRange,
                        double x) 
{   double xadj = (x-x0)*oneOverXRange;    
    assert(-NTraits<double>::getSignificant() <= xadj
           && xadj <= 1 + NTraits<double>::getSignificant());
    clampInPlace(0.0,xadj,1.0);
    return yRange*square(oneOverXRange)*d2stepUp(xadj); }

/** Third derivative of stepUp(): d^3/dx^3 stepUp(x). 
@param[in] x    Control parameter in range [0,1]. 
@return Third derivative of stepUp() at x. **/
inline double d3stepUp(double x) {
    assert(0 <= x && x <= 1);
    return 60+360*x*(x-1);      //60-360*x+360*x^2
}

/** Third derivative of stepDown(): d^3/dx^3 stepDown(x). 
@param[in] x    Control parameter in range [0,1]. 
@return Third derivative of stepDown() at x. **/
inline double d3stepDown(double x) {return -d3stepUp(x);}

/** Third derivative of stepAny(): d^3/dx^3 stepAny(x).
See stepAny() for parameter documentation. 
@return Third derivative of stepAny() at x. **/
inline double d3stepAny(double yRange,
                        double x0, double oneOverXRange,
                        double x) 
{   double xadj = (x-x0)*oneOverXRange;    
    assert(-NTraits<double>::getSignificant() <= xadj
           && xadj <= 1 + NTraits<double>::getSignificant());
    clampInPlace(0.0,xadj,1.0);
    return yRange*cube(oneOverXRange)*d3stepUp(xadj); }

            // float

/** @copydoc SimTK::stepUp(double) **/
inline float stepUp(float x) 
{   assert(0 <= x && x <= 1);
    return x*x*x*(10+x*(6*x-15)); }  //10x^3-15x^4+6x^5
/** @copydoc SimTK::stepDown(double) **/
inline float stepDown(float x) {return 1.0f-stepUp(x);}
/** @copydoc SimTK::stepAny(double,double,double,double,double) **/
inline float stepAny(float y0, float yRange,
                     float x0, float oneOverXRange,
                     float x) 
{   float xadj = (x-x0)*oneOverXRange;    
    assert(-NTraits<float>::getSignificant() <= xadj
           && xadj <= 1 + NTraits<float>::getSignificant());
    clampInPlace(0.0f,xadj,1.0f);
    return y0 + yRange*stepUp(xadj); }

/** @copydoc SimTK::dstepUp(double) **/
inline float dstepUp(float x) 
{   assert(0 <= x && x <= 1);
    const float xxm1=x*(x-1);
    return 30*xxm1*xxm1; }  //30x^2-60x^3+30x^4
/** @copydoc SimTK::dstepDown(double) **/
inline float dstepDown(float x) {return -dstepUp(x);}
/** @copydoc SimTK::dstepAny(double,double,double,double) **/
inline float dstepAny(float yRange,
                      float x0, float oneOverXRange,
                      float x) 
{   float xadj = (x-x0)*oneOverXRange;    
    assert(-NTraits<float>::getSignificant() <= xadj
           && xadj <= 1 + NTraits<float>::getSignificant());
    clampInPlace(0.0f,xadj,1.0f);
    return yRange*oneOverXRange*dstepUp(xadj); }

/** @copydoc SimTK::d2stepUp(double) **/
inline float d2stepUp(float x) 
{   assert(0 <= x && x <= 1);
    return 60*x*(1+x*(2*x-3)); }    //60x-180x^2+120x^3
/** @copydoc SimTK::d2stepDown(double) **/
inline float d2stepDown(float x) {return -d2stepUp(x);}
/** @copydoc SimTK::d2stepAny(double,double,double,double) **/
inline float d2stepAny(float yRange,
                       float x0, float oneOverXRange,
                       float x) 
{   float xadj = (x-x0)*oneOverXRange;    
    assert(-NTraits<float>::getSignificant() <= xadj
           && xadj <= 1 + NTraits<float>::getSignificant());
    clampInPlace(0.0f,xadj,1.0f);
    return yRange*square(oneOverXRange)*d2stepUp(xadj); }

/** @copydoc SimTK::d3stepUp(double) **/
inline float d3stepUp(float x) 
{   assert(0 <= x && x <= 1);
    return 60+360*x*(x-1); }    //60-360*x+360*x^2
/** @copydoc SimTK::d3stepDown(double) **/
inline float d3stepDown(float x) {return -d3stepUp(x);}
/** @copydoc SimTK::d3stepAny(double,double,double,double) **/
inline float d3stepAny(float yRange,
                       float x0, float oneOverXRange,
                       float x) 
{   float xadj = (x-x0)*oneOverXRange;    
    assert(-NTraits<float>::getSignificant() <= xadj
           && xadj <= 1 + NTraits<float>::getSignificant());
    clampInPlace(0.0f,xadj,1.0f);
    return yRange*cube(oneOverXRange)*d3stepUp(xadj); }

        // int converts to double; only supplied for stepUp(), stepDown()
/** @copydoc SimTK::stepUp(double)
Treats int argument as a double (avoids ambiguity). **/
inline double stepUp(int x) {return stepUp((double)x);}
/** @copydoc SimTK::stepDown(double)
Treats int argument as a double (avoids ambiguity). **/
inline double stepDown(int x) {return stepDown((double)x);}


/*@}*/

/** @defgroup EllipticIntegralsGroup   Elliptic integrals
    @ingroup ScalarFunctions

Elliptic integrals arise occasionally in contexts relevant to us, particularly
in geometric calculations involving ellipses or ellipsoids. Here we provide
functions for calculating the complete elliptic integrals of the first and
second kinds, which arise in Hertz contact theory for point contacts where
the contact area is elliptical. We use the following definitions for these
two integrals:
<pre>
    K(m) = integ[0,Pi/2] {1 / sqrt(1 - m sin^2(t)) dt}     1st kind
    E(m) = integ[0,Pi/2] {    sqrt(1 - m sin^2(t)) dt}     2nd kind
    0 <= m <= 1
</pre>
Elliptic integrals are defined only for arguments in range [0,1] inclusive.

We provide a function that calculates K(m) and E(m) to machine precision
(float or double) with a fast-converging iterative method adapted from 
ref. 1, which was in turn adapted from ref. 2. A much faster approximate 
version is also available, using the higher-precision approximation of the
two provided in ref. 2. The approximate version provides a smooth function 
that gives at least 7 digits of accuracy (in either float or double 
precision) across the full range at about 1/4 the cost of the machine 
precision version. For many applications, including engineering- or 
scientific-quality contact, 7 digits is more than adequate and in float 
precision that's all you can expect anyway.

@warning In the literature there are two different definitions for
elliptic integrals. The other definition (call them K'(k) and E'(k)) uses 
the argument k=sqrt(m) so that K'(k)=K(m^2), E'(k)=E(m^2). Our definition
is used by Matlab's ellipke() function and much engineering literature while
the K',E' definitions seem to be preferred by mathematicians.

For an ellipse with semimajor axis \a a and semiminor axis \a b (so a >= b),
the eccentricity e=sqrt(1-(b/a)^2). Our argument to the elliptic integrals in
that case is m = e^2 = 1-(b/a)^2. In constrast, K.L. Johnson uses the
mathematicians' definition in Chapter 4, pg. 95 of his book (ref. 4)
where he writes K(e) and E(e), where e is eccentricity as defined above,
that is e=sqrt(m), so we would call his K'(e)=K(e^2) and E'=E(e^2).

@author Michael Sherman

<h3>References</h3>
(1) Dyson, Evans, Snidle. "A simple accurate method for calculation of stresses
and deformations in elliptical Hertzian contacts", J. Mech. Eng. Sci. Proc.
IMechE part C 206:139-141, 1992.

(2) Abramovitz, Stegun, eds. Handbook of Mathematical Functions with
Formulas, Graphs, and Mathematical Tables, Dover, NY, 1972.

(3) Antoine, Visa, Sauvey, Abba. "Approximate analytical model for 
Hertzian Elliptical Contact Problems", ASME J. Tribology 128:660, 2006.

(4) Johnson, K.L. %Contact Mechanics. Cambridge University Press 1987 
(corrected edition). **/
/*@{*/

/** @cond **/ // Doxygen should skip this helper template function
template <class T> // float or double 
static inline std::pair<T,T> approxCompleteEllipticIntegralsKE_T(T m) {
    static const T a[] =
    {   (T)1.38629436112, (T)0.09666344259, (T)0.03590092383,
        (T)0.03742563713, (T)0.01451196212 };
    static const T b[] = 
    {   (T)0.5,           (T)0.12498593597, (T)0.06880248576,
        (T)0.03328355346, (T)0.00441787012 };
    static const T c[] =
    {   (T)1,             (T)0.44325141463, (T)0.06260601220,
        (T)0.04757383546, (T)0.01736506451 };
    static const T d[] =
    {   (T)0,             (T)0.24998368310, (T)0.09200180037,
        (T)0.04069697526, (T)0.00526449639 };

    const T SignificantReal = NTraits<T>::getSignificant();
    const T PiOver2         = NTraits<T>::getPi()/2;
    const T Infinity        = NTraits<T>::getInfinity();

    SimTK_ERRCHK1_ALWAYS(-SignificantReal < m && m < 1+SignificantReal,
        "approxCompleteEllipticIntegralsKE()", 
        "Argument m (%g) is outside the legal range [0,1].", (double)m);
    if (m >= 1) return std::make_pair(Infinity, (T)1);
    if (m <= 0) return std::make_pair(PiOver2, PiOver2);

    const T m1=1-m, m2=m1*m1, m3=m1*m2, m4=m2*m2; // 4 flops
    const T lnm2 = std::log(m1);  // ~50 flops

    // The rest is 35 flops.
    const T K = (a[0] + a[1]*m1 + a[2]*m2 + a[3]*m3 + a[4]*m4) 
         - lnm2*(b[0] + b[1]*m1 + b[2]*m2 + b[3]*m3 + b[4]*m4);
    const T E = (c[0] + c[1]*m1 + c[2]*m2 + c[3]*m3 + c[4]*m4) 
         - lnm2*(       d[1]*m1 + d[2]*m2 + d[3]*m3 + d[4]*m4);

    return std::make_pair(K,E);
}
/** @endcond **/

/** Given 0<=m<=1, return complete elliptic integrals of the first and 
second kinds, K(m) and E(m), approximated but with a maximum error of
2e-8 so at least 7 digits are correct (same in float or double 
precision).\ See @ref EllipticIntegralsGroup "Elliptic integrals" 
for a discussion.

Note that a full-precision computation of these integrals is also 
available; see completeEllipticIntegralsKE(). However, if you are using
float precision there is no point in using the exact computation since
this approximation is just as good and much faster.

@param[in] m  The argument to the elliptic integrals. Must be in range [0,1]
              although we allow for a very small amount of numerical error
              (see @ref SimTK::SignificantReal "SignificantReal") that might
              put m outside that range.
@return A std::pair p from which K(m)=p.first and E(m)=p.second.

You can find the approximating formula we're using in ref. 2, sections
17.3.34 and 17.3.36, pp. 591-2, or you can look at the code for this function 
(in the header file). The formulas are accurate to 2e-8 over the full [0-1] 
argument range, according to the authors. I checked our implementation against
Matlab's ellipke() function and the results are very good, providing at 
least 7 correct digits for a range of sample values.

The cost is about 90 flops.

@see completeEllipticIntegralsKE()
@see @ref EllipticIntegralsGroup "Elliptic integrals"

<h3>References</h3>
(1) Antoine, Visa, Sauvey, Abba. "Approximate analytical model for 
Hertzian Elliptical Contact Problems", ASME J. Tribology 128:660, 2006.

(2) Abramovitz, Stegun, eds. Handbook of Mathematical Functions with
Formulas, Graphs, and Mathematical Tables, Dover, NY, 1972.
**/
inline std::pair<double,double> 
approxCompleteEllipticIntegralsKE(double m)
{   return approxCompleteEllipticIntegralsKE_T<double>(m); }
/** This is the single precision (float) version of the approximate calculation
of elliptic integrals, still yielding about 7 digits of accuracy even 
though all calculations are done in float precision.
@see approxCompleteEllipticIntegralsKE(double) 
@see @ref EllipticIntegralsGroup "Elliptic integrals" **/
inline std::pair<float,float> 
approxCompleteEllipticIntegralsKE(float m)
{   return approxCompleteEllipticIntegralsKE_T<float>(m); }
/** This integer overload is present to prevent ambiguity; it converts its
argument to double precision and then calls 
approxCompleteEllipticIntegralsKE(double). Note that the only legal values
here are 0 and 1. 
@see approxCompleteEllipticIntegralsKE(double) 
@see @ref EllipticIntegralsGroup "Elliptic integrals" **/
inline std::pair<double,double> 
approxCompleteEllipticIntegralsKE(int m)
{   return approxCompleteEllipticIntegralsKE_T<double>((double)m); }


/** @cond **/ // Doxygen should skip this template helper function
template <class T> // float or double
static inline std::pair<T,T> completeEllipticIntegralsKE_T(T m) {
    const T SignificantReal = NTraits<T>::getSignificant();
    const T TenEps          = 10*NTraits<T>::getEps();
    const T Infinity        = NTraits<T>::getInfinity();
    const T PiOver2         = NTraits<T>::getPi() / 2;

    // Allow a little slop in the argument since it may have resulted
    // from a numerical operation that gave 0 or 1 plus or minus
    // roundoff noise.
    SimTK_ERRCHK1_ALWAYS(-SignificantReal < m && m < 1+SignificantReal,
        "completeEllipticIntegralsKE()", 
        "Argument m (%g) is outside the legal range [0,1].", (double)m);
    if (m >= 1) return std::make_pair(Infinity, (T)1);
    if (m <= 0) return std::make_pair(PiOver2, PiOver2);

    const T k = std::sqrt(1-m); // ~25 flops
    T v1=1, w1=k, c1=1, d1=k*k; // initialize iteration
    do { // ~50 flops per iteration
        T v2 = (v1+w1)/2;
        T w2 = std::sqrt(v1*w1);
        T c2 = (c1+d1)/2;
        T d2 = (w1*c1+v1*d1)/(v1+w1);
        v1=v2; w1=w2; c1=c2; d1=d2;
    } while(std::abs(v1-w1) >= TenEps);

    const T K = PiOver2/v1; // ~20 flops
    const T E = K*c1;

    return std::make_pair(K,E);
}
/** @endcond **/

/** Given 0<=m<=1, return complete elliptic integrals of the first and
second kinds, K(m) and E(m), calculated to (roughly) machine precision 
(float or double).\ See @ref EllipticIntegralsGroup "Elliptic integrals" 
for a discussion.

Note that we also provide a faster approximate method for calculating these
functions (see approxCompleteEllipticIntegralsKE()). The approximate method
is good enough for many scientific and engineering applications and is always
preferred if you are using float precision.

@param[in] m  The argument to the elliptic integrals. Must be in range [0,1]
              although we allow for a very small amount of numerical error
              (see @ref SimTK::SignificantReal "SignificantReal") that 
              might put m outside that range.
@return A std::pair p from which K(m)=p.first and E(m)=p.second.

Cost here is about 50 + 50*n flops, where n is the number of iterations
required. The number of iterations n you can expect to get a 
double-precision result is 4 for a:b < 2, 5 for a:b < 20, 6 for a:b < 1000,
7 after that; for single-precision it will take one fewer iterations.
In flops that's 250, 300, 350, 400 -- 300 will be typical for double,
250 for float. 
@see approxCompleteEllipticIntegralsKE() 
@see @ref EllipticIntegralsGroup "Elliptic integrals" **/
inline std::pair<double,double> completeEllipticIntegralsKE(double m)
{   return completeEllipticIntegralsKE_T<double>(m); }
/** This is the single precision (float) version of the machine-precision
calculation of elliptic integrals, providing accuracy to float precision
(about 7 digits) which is no better than you'll get with the much faster
approximate version, so use that instead!
@see approxCompleteEllipticIntegralsKE()
@see completeEllipticIntegralsKE(double) 
@see @ref EllipticIntegralsGroup "Elliptic integrals" **/
inline std::pair<float,float> completeEllipticIntegralsKE(float m)
{   return completeEllipticIntegralsKE_T<float>(m); }
/** This integer overload is present to prevent ambiguity; it converts its
argument to double precision and then calls 
completeEllipticIntegralsKE(double). Note that the only legal values
here are 0 and 1. 
@see completeEllipticIntegralsKE(double) 
@see @ref EllipticIntegralsGroup "Elliptic integrals" **/
inline std::pair<double,double> completeEllipticIntegralsKE(int m)
{   return completeEllipticIntegralsKE_T<double>((double)m); }

/*@}*/

} // namespace SimTK

#endif //SimTK_SIMMATRIX_SCALAR_H_
