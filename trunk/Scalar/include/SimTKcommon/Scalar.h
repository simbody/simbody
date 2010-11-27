#ifndef SimTK_SIMMATRIX_SCALAR_H_
#define SimTK_SIMMATRIX_SCALAR_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmatrix(tm)                       *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-10 Stanford University and the Authors.        *
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

namespace SimTK {

    /////////////////////////////////////////
    // Handy default-precision definitions //
    /////////////////////////////////////////

typedef conjugate<Real> Conjugate;  // like Complex

    // Note that these constants have memory addresses, so you can
    // return references to them.
    // These are static variables rather than static members to avoid
    // problems with static initialization order.

    // Properties of the floaing point representation.

static const Real& NaN               = NTraits<Real>::getNaN();      // "not a number"
static const Real& Infinity          = NTraits<Real>::getInfinity();

// Epsilon is the size of roundoff noise; it is the smallest positive number
// such that 1+Eps != 1.
static const Real& Eps               = NTraits<Real>::getEps();         // double ~1e-16, float ~1e-7
static const Real& SqrtEps           = NTraits<Real>::getSqrtEps();     // eps^(1/2): double ~1e- 8, float ~3e-4
static const Real& TinyReal          = NTraits<Real>::getTiny();        // eps^(5/4): double ~1e-20, float ~1e-9
static const Real& SignificantReal   = NTraits<Real>::getSignificant(); // eps^(7/8): double ~1e-14, float ~1e-6

static const Real& LeastPositiveReal = NTraits<Real>::getLeastPositive(); // double ~1e-308, float ~1e-38
static const Real& MostPositiveReal  = NTraits<Real>::getMostPositive();  // double ~1e+308, float ~1e+38
static const Real& LeastNegativeReal = NTraits<Real>::getLeastNegative();
static const Real& MostNegativeReal  = NTraits<Real>::getMostNegative();

// This is the number of decimal digits that can be reliably stored and
// retrieved in the default Real precision (typically log10(1/eps)-1).
static const int NumDigitsReal = NTraits<Real>::getNumDigits(); // double ~15, float ~6

// This is the smallest number of decimal digits you should store in a file
// if you want to be able to get exactly the same bit pattern back when you
// read it in. Typically, this is about log10(1/tiny).
static const int LosslessNumDigitsReal = NTraits<Real>::getLosslessNumDigits(); // double ~20, float ~9

    // Carefully calculated constants, with convenient memory addresses.

static const Real& Zero         = NTraits<Real>::getZero();
static const Real& One          = NTraits<Real>::getOne();
static const Real& MinusOne     = NTraits<Real>::getMinusOne();
static const Real& Two          = NTraits<Real>::getTwo();
static const Real& Three        = NTraits<Real>::getThree();

static const Real& OneHalf      = NTraits<Real>::getOneHalf();
static const Real& OneThird     = NTraits<Real>::getOneThird();
static const Real& OneFourth    = NTraits<Real>::getOneFourth();
static const Real& OneFifth     = NTraits<Real>::getOneFifth();
static const Real& OneSixth     = NTraits<Real>::getOneSixth();
static const Real& OneSeventh   = NTraits<Real>::getOneSeventh();
static const Real& OneEighth    = NTraits<Real>::getOneEighth();
static const Real& OneNinth     = NTraits<Real>::getOneNinth();
static const Real& Pi           = NTraits<Real>::getPi();
static const Real& OneOverPi    = NTraits<Real>::getOneOverPi();
static const Real& E            = NTraits<Real>::getE();
static const Real& Log2E        = NTraits<Real>::getLog2E();
static const Real& Log10E       = NTraits<Real>::getLog10E();
static const Real& Sqrt2        = NTraits<Real>::getSqrt2();
static const Real& OneOverSqrt2 = NTraits<Real>::getOneOverSqrt2();  // also sqrt(2)/2
static const Real& Sqrt3        = NTraits<Real>::getSqrt3();
static const Real& OneOverSqrt3 = NTraits<Real>::getOneOverSqrt3();
static const Real& CubeRoot2    = NTraits<Real>::getCubeRoot2();
static const Real& CubeRoot3    = NTraits<Real>::getCubeRoot3();
static const Real& Ln2          = NTraits<Real>::getLn2();
static const Real& Ln10         = NTraits<Real>::getLn10();

// We only need one complex constant. For the rest just use
// Complex(the Real constant), or if you need an address use
// NTraits<Complex>::getPi(), etc.
static const Complex& I = NTraits<Complex>::getI();


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
inline bool atMostOneBitIsSet(unsigned char v)      {return (v&v-1)==0;}
inline bool atMostOneBitIsSet(unsigned short v)     {return (v&v-1)==0;}
inline bool atMostOneBitIsSet(unsigned int v)       {return (v&v-1)==0;}
inline bool atMostOneBitIsSet(unsigned long v)      {return (v&v-1)==0;}
inline bool atMostOneBitIsSet(unsigned long long v) {return (v&v-1)==0;}
inline bool atMostOneBitIsSet(signed char v)        {return (v&v-1)==0;}
inline bool atMostOneBitIsSet(char v)               {return (v&v-1)==0;}
inline bool atMostOneBitIsSet(short v)              {return (v&v-1)==0;}
inline bool atMostOneBitIsSet(int v)                {return (v&v-1)==0;}
inline bool atMostOneBitIsSet(long v)               {return (v&v-1)==0;}
inline bool atMostOneBitIsSet(long long v)          {return (v&v-1)==0;}
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

inline bool signBit(unsigned char i)      {return false;}
inline bool signBit(unsigned short i)     {return false;}
inline bool signBit(unsigned int i)       {return false;}
inline bool signBit(unsigned long i)      {return false;}
inline bool signBit(unsigned long long i) {return false;}

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
inline int sign(const long double& x) {return x>0 ? 1 : (x<0 ? -1 : 0);}

inline int sign(const negator<float>&       x) {return -sign(-x);} // -x is free
inline int sign(const negator<double>&      x) {return -sign(-x);}
inline int sign(const negator<long double>& x) {return -sign(-x);}
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
inline long double square(const long double& x) {return x*x;}

// Negation is free for negators, so we can square them and clean
// them up at the same time at no extra cost.
inline float       square(const negator<float>&       x) {return square(-x);}
inline double      square(const negator<double>&      x) {return square(-x);}
inline long double square(const negator<long double>& x) {return square(-x);}

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
inline long double cube(const long double& x) {return x*x*x;}

// To keep this cheap we'll defer getting rid of the negator<> until
// some other operation. We cube -x and then recast that to negator<>
// on return, for a total cost of 2 flops.
inline negator<float> cube(const negator<float>& x) {
    return negator<float>::recast(cube(-x));
}
inline negator<double> cube(const negator<double>& x) {
    return negator<double>::recast(cube(-x));
}
inline negator<long double> cube(const negator<long double>& x) {
    return negator<long double>::recast(cube(-x));
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
[los,high]. Both functions are overloaded for all the integral and real 
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
/** @copydoc SimTK::clampInPlace(double,double&,double) **/
inline long double& clampInPlace(long double low, long double& v, long double high) 
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
Takes integer bounds to avoid need for explicit casts. **/
inline long double& clampInPlace(int low, long double& v, int high) 
{   return clampInPlace((long double)low,v,(long double)high); }

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
inline long double& clampInPlace(int low, long double& v, long double high) 
{   return clampInPlace((long double)low,v,high); }

/** @copydoc SimTK::clampInPlace(double,double&,double) 
Takes an integer bound to avoid need for explicit casts. **/
inline double& clampInPlace(double low, double& v, int high) 
{   return clampInPlace(low,v,(double)high); }
/** @copydoc SimTK::clampInPlace(double,double&,double) 
Takes an integer bound to avoid need for explicit casts. **/
inline float& clampInPlace(float low, float& v, int high) 
{   return clampInPlace(low,v,(float)high); }
/** @copydoc SimTK::clampInPlace(double,double&,double) 
Takes an integer bound to avoid need for explicit casts. **/
inline long double& clampInPlace(long double low, long double& v, int high) 
{   return clampInPlace(low,v,(long double)high); }

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
/** @copydoc SimTK::clampInPlace(double,double&,double) **/
inline negator<long double>& clampInPlace(long double low, negator<long double>& v, long double high) 
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
/** @copydoc SimTK::clamp(double,double,double) **/
inline long double clamp(long double low, long double v, long double high) 
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
Takes integer bounds to avoid need for explicit casts. **/
inline long double clamp(int low, long double v, int high) 
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
inline long double clamp(int low, long double v, long double high) 
{   return clampInPlace(low,v,high); }

/** @copydoc SimTK::clamp(double,double,double) 
Takes an integer bound to avoid need for explicit casts. **/
inline double clamp(double low, double v, int high) 
{   return clampInPlace(low,v,high); }
/** @copydoc SimTK::clamp(double,double,double) 
Takes an integer bound to avoid need for explicit casts. **/
inline float clamp(float low, float v, int high) 
{   return clampInPlace(low,v,high); }
/** @copydoc SimTK::clamp(double,double,double)
Takes an integer bound to avoid need for explicit casts. **/
inline long double clamp(long double low, long double v, int high) 
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
/** @copydoc SimTK::clamp(double,double,double)
Explicitly takes a negator<long double> argument to help the compiler find the
right overload, but the negation is performed (1 extra flop) and the result 
type is an ordinary long double. **/
inline long double clamp(long double low, negator<long double> v, long double high) 
{   return clamp(low,(long double)v,high); }
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

            // long double

/** @copydoc SimTK::stepUp(double) **/
inline long double stepUp(long double x) 
{   assert(0 <= x && x <= 1);
    return x*x*x*(10+x*(6*x-15)); }  //10x^3-15x^4+6x^5
/** @copydoc SimTK::stepDown(double) **/
inline long double stepDown(long double x) {return 1.0L-stepUp(x);}
/** @copydoc SimTK::stepAny(double,double,double,double,double) **/
inline long double stepAny(long double y0, long double yRange,
                           long double x0, long double oneOverXRange,
                           long double x) 
{   long double xadj = (x-x0)*oneOverXRange;    
    assert(-NTraits<long double>::getSignificant() <= xadj
           && xadj <= 1 + NTraits<long double>::getSignificant());
    clampInPlace(0.0L,xadj,1.0L);
    return y0 + yRange*stepUp(xadj); }


/** @copydoc SimTK::dstepUp(double) **/
inline long double dstepUp(long double x) 
{   assert(0 <= x && x <= 1);
    const long double xxm1=x*(x-1);
    return 30*xxm1*xxm1; }          //30x^2-60x^3+30x^4
/** @copydoc SimTK::dstepDown(double) **/
inline long double dstepDown(long double x) {return -dstepUp(x);}
/** @copydoc SimTK::dstepAny(double,double,double,double) **/
inline long double dstepAny(long double yRange,
                            long double x0, long double oneOverXRange,
                            long double x) 
{   long double xadj = (x-x0)*oneOverXRange;    
    assert(-NTraits<long double>::getSignificant() <= xadj
           && xadj <= 1 + NTraits<long double>::getSignificant());
    clampInPlace(0.0L,xadj,1.0L);
    return yRange*oneOverXRange*dstepUp(xadj); }

/** @copydoc SimTK::d2stepUp(double) **/
inline long double d2stepUp(long double x) 
{   assert(0 <= x && x <= 1);
    return 60*x*(1+x*(2*x-3)); }    //60x-180x^2+120x^3
/** @copydoc SimTK::d2stepDown(double) **/
inline long double d2stepDown(long double x) {return -d2stepUp(x);}
/** @copydoc SimTK::d2stepAny(double,double,double,double) **/
inline long double d2stepAny(long double yRange,
                             long double x0, long double oneOverXRange,
                             long double x) 
{   long double xadj = (x-x0)*oneOverXRange;    
    assert(-NTraits<long double>::getSignificant() <= xadj
           && xadj <= 1 + NTraits<long double>::getSignificant());
    clampInPlace(0.0L,xadj,1.0L);
    return yRange*square(oneOverXRange)*d2stepUp(xadj); }


/** @copydoc SimTK::d3stepUp(double) **/
inline long double d3stepUp(long double x) 
{   assert(0 <= x && x <= 1);
    return 60+360*x*(x-1); }        //60-360*x+360*x^2
/** @copydoc SimTK::d3stepDown(double) **/
inline long double d3stepDown(long double x) {return -d3stepUp(x);}
/** @copydoc SimTK::d3stepAny(double,double,double,double) **/
inline long double d3stepAny(long double yRange,
                             long double x0, long double oneOverXRange,
                             long double x) 
{   long double xadj = (x-x0)*oneOverXRange;    
    assert(-NTraits<long double>::getSignificant() <= xadj
           && xadj <= 1 + NTraits<long double>::getSignificant());
    clampInPlace(0.0L,xadj,1.0L);
    return yRange*cube(oneOverXRange)*d3stepUp(xadj); }

        // int converts to double; only supplied for stepUp(), stepDown()
/** @copydoc SimTK::stepUp(double)
Treats int argument as a double (avoids ambiguity). **/
inline double stepUp(int x) {return stepUp((double)x);}
/** @copydoc SimTK::stepDown(double)
Treats int argument as a double (avoids ambiguity). **/
inline double stepDown(int x) {return stepDown((double)x);}


/*@}*/

} // namespace SimTK

#endif //SimTK_SIMMATRIX_SCALAR_H_
