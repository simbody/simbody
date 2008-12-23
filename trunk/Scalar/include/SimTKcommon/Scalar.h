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

#include "SimTKcommon/internal/conjugate.h"
#include "SimTKcommon/internal/CompositeNumericalTypes.h"
#include "SimTKcommon/internal/NTraits.h"
#include "SimTKcommon/internal/negator.h"

#include <complex>
#include <cmath>
#include <climits>

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

// This utility answers the question "if I put this integral value in an int and then
// get it back, will its value be the same?".
inline bool canStoreInInt(char)            {return true;}
inline bool canStoreInInt(unsigned char)   {return true;}
inline bool canStoreInInt(signed char)     {return true;}
inline bool canStoreInInt(short)           {return true;}
inline bool canStoreInInt(unsigned short)  {return true;}
inline bool canStoreInInt(int)             {return true;}
inline bool canStoreInInt(unsigned int  u) {return (unsigned int)(int(u)) == u;}
inline bool canStoreInInt(long i)          {return long(int(i)) == i;}
inline bool canStoreInInt(unsigned long u) {return (unsigned long)(int(u)) == u;}

// This utility answers the question "is this integral value a nonnegative number
// that can be stored in an int?".
inline bool canStoreInNonnegativeInt(char c)           {return c >= 0;}
inline bool canStoreInNonnegativeInt(unsigned char c)  {return true;}
inline bool canStoreInNonnegativeInt(signed char c)    {return c >= 0;}
inline bool canStoreInNonnegativeInt(short s)          {return s >= 0;}
inline bool canStoreInNonnegativeInt(unsigned short s) {return true;}
inline bool canStoreInNonnegativeInt(int  i)           {return i >= 0;}
inline bool canStoreInNonnegativeInt(long l)           {return canStoreInInt(l) && l >= 0;}
inline bool canStoreInNonnegativeInt(unsigned int  u)  {return canStoreInInt(u);}
inline bool canStoreInNonnegativeInt(unsigned long u)  {return canStoreInInt(u);}

/**
 * isNaN(x) provides a reliable way to determine if x is one of the
 * "not a number" floating point forms. Comparing x==NaN does not
 * work because any relational operation involving NaN always
 * return false, even (NaN==NaN)! 
 * This routine is specialized for all SimTK scalar types:
 * float, double std::complex<P>, SimTK::conjugate<P>,
 * and SimTK::negator<T>, where T is any of the above. For
 * complex and conjugate types, isNaN() returns true if either
 * the real or imaginary part or both are NaN.
 */
inline bool isNaN(const float& x)  {return std::isnan(x);}
inline bool isNaN(const double& x) {return std::isnan(x);}

inline bool isNaN(const negator<float>&  x) {return isNaN(-x);}
inline bool isNaN(const negator<double>& x) {return isNaN(-x);}

template <class P> inline bool
isNaN(const std::complex<P>& x)
{   return isNaN(x.real()) || isNaN(x.imag());}

template <class P> inline bool
isNaN(const conjugate<P>& x)
{   return isNaN(x.real()) || isNaN(x.negImag());}

template <class P> inline bool
isNaN(const negator< std::complex<P> >& x) {return isNaN(-x);}
template <class P> inline bool
isNaN(const negator< conjugate<P> >&    x) {return isNaN(-x);}


/**
 * isFinite(x) provides a reliable way to determine if x is a "normal"
 * floating point number, meaning not a NaN or +/- Infinity.
 * This routine is specialized for all SimTK scalar types:
 * float, double, std::complex<P>, SimTK::conjugate<P>,
 * and SimTK::negator<T>, where T is any of the above. For
 * complex and conjugate types, isFinite() returns true if
 * the real and imaginary parts are both finite.
 */
inline bool isFinite(const float&  x) {return std::isfinite(x);}
inline bool isFinite(const double& x) {return std::isfinite(x);}

template <class P> inline bool
isFinite(const std::complex<P>& x)
{   return isFinite(x.real()) && isFinite(x.imag());}

template <class P> inline bool
isFinite(const conjugate<P>& x)
{   return isFinite(x.real()) && isFinite(x.negImag());}

inline bool isFinite(const negator<float>&  x) {return isFinite(-x);}
inline bool isFinite(const negator<double>& x) {return isFinite(-x);}

template <class P> inline bool
isFinite(const negator< std::complex<P> >& x) {return isFinite(-x);}
template <class P> inline bool
isFinite(const negator< conjugate<P> >&    x) {return isFinite(-x);}

/**
 * isInf(x) provides a reliable way to determine if x is one of
 * the two infinities (either negative or positive).
 * This routine is specialized for all SimTK scalar types:
 * float, double std::complex<P>, SimTK::conjugate<P>,
 * and SimTK::negator<T>, where T is any of the above. For
 * complex and conjugate types, isInf() returns true if both
 * components are infinite, or one is infinite and the other
 * finite. That is, isInf() will never return true if one
 * component is NaN.
 */
inline bool isInf(const float&  x) {return std::isinf(x);}
inline bool isInf(const double& x) {return std::isinf(x);}

inline bool isInf(const negator<float>&  x) {return isInf(-x);}
inline bool isInf(const negator<double>& x) {return isInf(-x);}

template <class P> inline bool
isInf(const std::complex<P>& x) {
    return (isInf(x.real()) && !isNaN(x.imag()))
        || (isInf(x.imag()) && !isNaN(x.real()));
}

template <class P> inline bool
isInf(const conjugate<P>& x) {
    return (isInf(x.real())    && !isNaN(x.negImag()))
        || (isInf(x.negImag()) && !isNaN(x.real()));
}

template <class P> inline bool
isInf(const negator< std::complex<P> >& x) {return isInf(-x);}
template <class P> inline bool
isInf(const negator< conjugate<P> >&    x) {return isInf(-x);}

/**
 * signBit(x) provides a fast way to determine the value of the sign
 * bit (as a bool) for integral and floating types. Note that this is
 * significantly different than sign(x); be sure you know what you're doing
 * if you use this method. signBit() refers to the underlying representation
 * rather than the numerical value. For example, for floating types there
 * are two zeroes, +0 and -0 which have opposite sign bits but the same
 * sign (0). Also, unsigned types have sign() of 0 or 1, but they are 
 * considered here to have a sign bit of 0 always, since the stored high
 * bit does not indicate a negative value.
 *
 * NOTES:
 *  - signBit() is overloaded for 'signed char' and 'unsigned char', but not
 *    for plain 'char' because we don't know whether to interpret the high
 *    bit as a sign in that case (because the C++ standard leaves it unspecified).
 *  - complex and conjugate numbers do not have sign bits.
 *  - negator<float> and negator<double> have the *same* sign bit as the
 *    underlying representation -- it's up to you to realize that it is
 *    interpreted differently!
 */
inline bool signBit(unsigned char)      {return 0;}
inline bool signBit(unsigned short)     {return 0;}
inline bool signBit(unsigned int)       {return 0;}
inline bool signBit(unsigned long)      {return 0;}
inline bool signBit(unsigned long long) {return 0;}

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

/**
 * s=sign(n) returns int -1,0,1 according to n<0, n==0, n>0 for any integer
 * or real numeric type, unsigned 0 or 1 for any unsigned argument. This
 * routine is specialized for each of the int, unsigned, and real types
 * of all sizes. Sign is defined for "signed char" and "unsigned char" but
 * not plain "char" since the language leaves unspecified whether that is
 * a signed or unsigned type. Sign is not defined for complex or conjugate.
 */
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

/**
 * y=square(x) returns the square of the argument for any numeric type.
 * We promise to evaluate x only once. We assume is is acceptable for
 * the result type to be the same as the argument type; if it won't
 * fit caller must cast argument to a wider type first. This is an
 * inline routine which will run as fast as an explicit multiply
 * (x*x) in optimized code, and somewhat faster for complex types.
 */
inline unsigned char  square(unsigned char  u) {return u*u;}
inline unsigned short square(unsigned short u) {return u*u;}
inline unsigned int   square(unsigned int   u) {return u*u;}
inline unsigned long  square(unsigned long  u) {return u*u;}

inline char        square(char c) {return c*c;}

inline signed char square(signed char i) {return i*i;}
inline short       square(short       i) {return i*i;}
inline int         square(int         i) {return i*i;}
inline long        square(long        i) {return i*i;}

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

/**
 * y=cube(x) returns the cube of the argument for any numeric type,
 * integral or floating point. We promise to evaluate x only once.
 * We assume is is acceptable for the result type to be the same
 * as the argument type; if it won't fit caller must cast argument
 * to a wider type first. This is an inline routine which will run
 * as fast as explicit multiplies (x*x*x) in optimized code, and
 * even faster for complex types.
 */
inline unsigned char  cube(unsigned char  u) {return u*u*u;}
inline unsigned short cube(unsigned short u) {return u*u*u;}
inline unsigned int   cube(unsigned int   u) {return u*u*u;}
inline unsigned long  cube(unsigned long  u) {return u*u*u;}

inline char        cube(char c) {return c*c*c;}

inline signed char cube(signed char i) {return i*i*i;}
inline short       cube(short       i) {return i*i*i;}
inline int         cube(int         i) {return i*i*i;}
inline long        cube(long        i) {return i*i*i;}

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

} // namespace SimTK

#endif //SimTK_SIMMATRIX_SCALAR_H_
