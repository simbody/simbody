#ifndef SimTK_SIMMATRIX_NTRAITS_H_
#define SimTK_SIMMATRIX_NTRAITS_H_

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
 * This file contains classes and typedefs needed to provide uniform handling of 
 * floating point numeric values. There are three numeric types: real, complex, 
 * conjugate and each comes in float and double precision. Each of 
 * these may be modified by a negator, which does not change the in-memory 
 * \e representation but negates the \e interpretation. Thus there are 12 
 * distinct scalar types: 2 precisions each of real, complex, and conjugate and 
 * their negators.
 * @verbatim
 *      The Scalar Types
 *      ----------------
 *      Here is a complete taxonomy of the scalar types we support.
 *
 *      <scalar>    ::= <number> | negator< <number> >
 *      <number>    ::= <standard> | <conjugate>
 *      <standard>  ::= <real> | <complex>
 *
 *      <real>      ::= float | double
 *      <complex>   ::= complex< <real> >
 *      <conjugate> ::= conjugate< <real> >
 * @endverbatim
 */

#include "SimTKcommon/Constants.h"
#include "SimTKcommon/internal/CompositeNumericalTypes.h"

#include <cstddef>
#include <cassert>
#include <complex>
#include <iostream>
#include <limits>

using std::complex;
    
namespace SimTK {

// This is the 3rd type of number, conjugate. It is like complex except that
// the represented value is the conjugate of the value represented by a complex
// number containing the same bit pattern. That is, complex and conjugate both
// contain two real numbers, re and im, with complex(re,im) meaning
// (re + im*i) while conjugate(re,im) means (re - im*i). It is guaranteed that
// our conjugate type and complex type have identical sizes and representations.
// Together, these defininitions and guarantees permit conjugation
// to be done by reinterpretation rather than be computation.
template <class R> class conjugate; // Only defined for float, double

// Specializations of this class provide information about Composite Numerical 
// Types in the style of std::numeric_limits<T>. It is specialized for the 
// numeric types but can be invoked on any composite numerical type as well.
template <class T> class CNT;

// NTraits provides information like CNT, but for numeric types only.
template <class N> class NTraits; // allowed only for N=<number>
template <class R> class NTraits< complex<R> >;
template <class R> class NTraits< conjugate<R> >;
template <> class NTraits<float>;
template <> class NTraits<double>;

// This is an adaptor for numeric types which negates the apparent values. A
// negator<N> has exactly the same internal representation as a numeric value N, 
// but it is to be interpreted has having the negative of the value it would 
// have if interpreted as an N. This permits negation to be done by 
// reinterpretation rather than computation. A full set of arithmetic operators
// are provided involving negator<N>'s and N's. Sometimes we can save an op or
// two this way. For example negator<N>*negator<N> can be performed as an N*N
// since the negations cancel, and we saved two floating point negations.
template <class N> class negator;      // Only defined for numbers

/// This class is specialized for all 16 combinations of standard types
/// (that is, real and complex types in each of two precisions)
/// and has typedefs "Type" which is the appropriate "widened"
/// type for use when R1 & R2 appear in an operation together, and
/// "Precision" which is the wider precision (float,double). 
/// For example, if R1=complex< float > and R2=double, Widest<R1,R2>::Type is
/// complex<double> and Widest<R1,R2>::Precision is double.
template <class R1, class R2> struct Widest {/* Only defined for built-ins. */};
template <> struct Widest<float,float>              {typedef float       Type;  typedef float       Precision;};
template <> struct Widest<float,double>             {typedef double      Type;  typedef double      Precision;};
template <> struct Widest<double,float>             {typedef double      Type;  typedef double      Precision;};
template <> struct Widest<double,double>            {typedef double      Type;  typedef double      Precision;};
template <class R1, class R2> struct Widest< complex<R1>,complex<R2> > { 
    typedef complex< typename Widest<R1,R2>::Type > Type; 
    typedef typename Widest<R1,R2>::Precision       Precision; 
};
template <class R1, class R2> struct Widest< complex<R1>,R2 > { 
    typedef complex< typename Widest<R1,R2>::Type > Type; 
    typedef typename Widest<R1,R2>::Precision       Precision;  
};
template <class R1, class R2> struct Widest< R1,complex<R2> > { 
    typedef complex< typename Widest<R1,R2>::Type > Type; 
    typedef typename Widest<R1,R2>::Precision       Precision; 
};

/// This class is specialized for all 16 combinations of standard types
/// (that is, real and complex types in each of two precisions)
/// and has typedefs "Type" which is the appropriate "narrowed"
/// type for use when R1 & R2 appear in an operation together where the
/// result must be of the narrower precision, and "Precision" which is
/// the expected precision of the result (float,
/// double). For example, if R1=complex< double > and R2=float, 
/// Narrowest<R1,R2>::Type is complex< float > and Narrowest<R1,R2>::Precision
/// is float.
template <class R1, class R2> struct Narrowest {/* Only defined for built-ins. */};
template <> struct Narrowest<float,float>              {typedef float  Type; typedef float Precision;};
template <> struct Narrowest<float,double>             {typedef float  Type; typedef float Precision;};
template <> struct Narrowest<double,float>             {typedef float  Type; typedef float Precision;};
template <> struct Narrowest<double,double>            {typedef double Type; typedef double Precision;};
template <class R1, class R2> struct Narrowest< complex<R1>,complex<R2> > { 
    typedef complex< typename Narrowest<R1,R2>::Type >  Type; 
    typedef typename Narrowest<R1,R2>::Precision        Precision;
};
template <class R1, class R2> struct Narrowest< complex<R1>,R2 > { 
    typedef complex< typename Narrowest<R1,R2>::Type >  Type; 
    typedef typename Narrowest<R1,R2>::Precision        Precision; 
};
template <class R1, class R2> struct Narrowest< R1,complex<R2> > { 
    typedef complex< typename Narrowest<R1,R2>::Type >  Type; 
    typedef typename Narrowest<R1,R2>::Precision        Precision; 
};

/// RTraits is a helper class for NTraits.
template <class R> class RTraits {/* Only defined for real types */};
template <> class RTraits<float> {
public:
    /// Attainable accuracy at this precision.
    static const float& getEps()         {static const float c=std::numeric_limits<float>::epsilon(); return c;}
    /// What multiple of attainable accuracy do we consider significant? 
    static const float& getSignificant() {static const float c=std::pow(getEps(), 0.875f); return c;}
    /// The default numerical error tolerance is always given in double precision.
    static double getDefaultTolerance()  {return (double)getSignificant();}
};
template <> class RTraits<double> {
public:
    static const double& getEps()         {static const double c=std::numeric_limits<double>::epsilon(); return c;}
    static const double& getSignificant() {static const double c=std::pow(getEps(), 0.875); return c;}
    static double getDefaultTolerance()   {return getSignificant();}
};

/**
 * @defgroup isNaN isNaN()
 * @ingroup ScalarFunctions
 *
 * isNaN(x) provides a reliable way to determine if x is one of the "not a 
 * number" floating point forms. Comparing x==NaN does not work because any 
 * relational operation involving NaN always return false, even (NaN==NaN)! 
 * This routine is specialized for all SimTK scalar types:
 *  - float, double
 *  - std::complex<P>        (P is one of the above precisions)
 *  - SimTK::conjugate<P>
 *  - SimTK::negator<T>      (T is any of the above)
 *
 * For complex and conjugate types, isNaN() returns true if either the real or 
 * imaginary part or both are NaN.
 */
// See negator.h for isNaN() applied to negated scalars.
//@{
inline bool isNaN(const float& x)  {return std::isnan(x);}
inline bool isNaN(const double& x) {return std::isnan(x);}

template <class P> inline bool
isNaN(const std::complex<P>& x)
{   return isNaN(x.real()) || isNaN(x.imag());}

template <class P> inline bool
isNaN(const conjugate<P>& x)
{   return isNaN(x.real()) || isNaN(x.negImag());}
//@}

/**
 * @defgroup isFinite isFinite()
 * @ingroup ScalarFunctions
 *
 * isFinite(x) provides a reliable way to determine if x is a "normal"
 * floating point number, meaning not a NaN or +/- Infinity.
 * This routine is specialized for all SimTK scalar types:
 * float, double, std::complex<P>, SimTK::conjugate<P>,
 * and SimTK::negator<T>, where T is any of the above. For
 * complex and conjugate types, isFinite() returns true if
 * the real and imaginary parts are both finite.
 */
// See negator.h for isFinite() applied to negated scalars.
//@{
inline bool isFinite(const float&  x) {return std::isfinite(x);}
inline bool isFinite(const double& x) {return std::isfinite(x);}

template <class P> inline bool
isFinite(const std::complex<P>& x)
{   return isFinite(x.real()) && isFinite(x.imag());}

template <class P> inline bool
isFinite(const conjugate<P>& x)
{   return isFinite(x.real()) && isFinite(x.negImag());}
//@}

/**
 * @defgroup isInf isInf()
 * @ingroup ScalarFunctions
 *
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
// See negator.h for isInf() applied to negated scalars.
//@{
inline bool isInf(const float&  x) {return std::isinf(x);}
inline bool isInf(const double& x) {return std::isinf(x);}

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
//@}

/**
 * @defgroup isNumericallyEqual isNumericallyEqual()
 * @ingroup ScalarFunctions
 *
 * isNumericallyEqual(x,y) compares two scalar types using a tolerance (default
 * or explicitly specified) and returns true if they are close enough.
 *
 * The default tolerance used is the NTraits<P>::getSignificant() value (about 
 * 1e-14 in double precision, 1e-6 in float) for the narrower of the types being 
 * compared but you can override that. The tolerance is both a relative and 
 * absolute tolerance; for two numbers a and b and tolerance tol we compute the 
 * following condition:
 * <pre>
 *      scale = max(|a|,|b|,1)
 *      isNumericallyEqual = |a-b| <= scale*tol
 * </pre>
 * For complex or conjugate numbers we insist that both the real and
 * imaginary parts independently satisfy the above condition.
 *
 * \par Mixed precision
 * We support mixed argument types here in which case the default
 * tolerance used is the one appropriate to the \e lower-precision
 * argument. When one argument is an integer, the default tolerance
 * used is that of the floating point argument. Comparisons may be performed
 * at a higher precision than the tolerance to avoid incorrect truncation; for
 * example an int cannot generally be contained in a float so an int/float
 * comparison is done as double/double, but to float tolerance.
 *
 * \par Treatment of NaN
 * When both arguments are NaN they are considered equal here, which is different
 * than the behavior of the IEEE-sanctioned "==" comparison for which 
 * NaN==NaN returns false. If only one argument is NaN we return false. When
 * comparing complex or conjugate numbers the real and imaginary parts are
 * tested separately, so (NaN,0) and (0,NaN) don't test equal despite the
 * fact that isNaN() would return true for both of them. We don't distinguish
 * among *types* of NaNs, though, so NaN and -NaN (if there is such a thing)
 * will test numerically equal.
 */
//@{
/// Compare two floats for approximate equality.
inline bool isNumericallyEqual(const float& a, const float& b, 
                               double tol = RTraits<float>::getDefaultTolerance())
{   if (isNaN(a)) return isNaN(b); else if (isNaN(b)) return false;
    const float scale = std::max(std::max(std::abs(a),std::abs(b)), 1.f);
    return std::abs(a-b) <= scale*(float)tol; }
/// Compare two doubles for approximate equality.
inline bool isNumericallyEqual(const double& a, const double& b, 
                               double tol = RTraits<double>::getDefaultTolerance())
{   if (isNaN(a)) return isNaN(b); else if (isNaN(b)) return false;
    const double scale = std::max(std::max(std::abs(a),std::abs(b)), 1.);
    return std::abs(a-b) <= scale*tol; }

/// Compare a float and a double for approximate equality at float precision.
inline bool isNumericallyEqual(const float& a, const double& b, 
                               double tol = RTraits<float>::getDefaultTolerance())
{   return isNumericallyEqual((double)a, b, tol); }
/// Compare a float and a double for approximate equality at float precision.
inline bool isNumericallyEqual(const double& a, const float& b, 
                               double tol = RTraits<float>::getDefaultTolerance())
{   return isNumericallyEqual(a, (double)b, tol); }

/// %Test a float for approximate equality to an integer.
inline bool isNumericallyEqual(const float& a, int b,
                               double tol = RTraits<float>::getDefaultTolerance())
{   return isNumericallyEqual(a, (double)b, tol); }
/// %Test a float for approximate equality to an integer.
inline bool isNumericallyEqual(int a, const float& b,
                               double tol = RTraits<float>::getDefaultTolerance())
{   return isNumericallyEqual((double)a, b, tol); }
/// %Test a double for approximate equality to an integer.
inline bool isNumericallyEqual(const double& a, int b,
                               double tol = RTraits<double>::getDefaultTolerance())
{   return isNumericallyEqual(a, (double)b, tol); }
/// %Test a double for approximate equality to an integer.
inline bool isNumericallyEqual(int a, const double& b,
                               double tol = RTraits<double>::getDefaultTolerance())
{   return isNumericallyEqual((double)a, b, tol); }

/// Compare two complex numbers for approximate equality, using the numerical 
/// accuracy expectation of the narrower of the two precisions in the case of mixed 
/// precision.
template <class P, class Q>
inline bool isNumericallyEqual
  ( const std::complex<P>& a, const std::complex<Q>& b, 
    double tol = RTraits<typename Narrowest<P,Q>::Precision>::getDefaultTolerance())
{   return isNumericallyEqual(a.real(),b.real(),tol)
        && isNumericallyEqual(a.imag(),b.imag(),tol); }

/// Compare two conjugate numbers for approximate equality, using the numerical 
/// accuracy expectation of the narrower of the two precisions in the case of mixed 
/// precision.
template <class P, class Q>
inline bool isNumericallyEqual
  ( const conjugate<P>& a, const conjugate<Q>& b, 
    double tol = RTraits<typename Narrowest<P,Q>::Precision>::getDefaultTolerance())
{   return isNumericallyEqual(a.real(),b.real(),tol)
        && isNumericallyEqual(a.imag(),b.imag(),tol); }

/// Compare a complex and a conjugate number for approximate equality, using the 
/// numerical accuracy expectation of the narrower of the two precisions in the 
/// case of mixed precision.
template <class P, class Q>
inline bool isNumericallyEqual
  ( const std::complex<P>& a, const conjugate<Q>& b, 
    double tol = RTraits<typename Narrowest<P,Q>::Precision>::getDefaultTolerance())
{   return isNumericallyEqual(a.real(),b.real(),tol)
        && isNumericallyEqual(a.imag(),b.imag(),tol); }

/// Compare a complex and a conjugate number for approximate equality, using the 
/// numerical accuracy expectation of the narrower of the two precisions in the 
/// case of mixed precision.
template <class P, class Q>
inline bool isNumericallyEqual
  ( const conjugate<P>& a, const std::complex<Q>& b, 
    double tol = RTraits<typename Narrowest<P,Q>::Precision>::getDefaultTolerance())
{   return isNumericallyEqual(a.real(),b.real(),tol)
        && isNumericallyEqual(a.imag(),b.imag(),tol); }

/// %Test whether a complex number is approximately equal to a particular real float.
template <class P> inline bool 
isNumericallyEqual(const std::complex<P>& a, const float& b, 
                   double tol = RTraits<float>::getDefaultTolerance())
{   return isNumericallyEqual(a.real(),b,tol) && isNumericallyEqual(a.imag(),0.f,tol); }
/// %Test whether a complex number is approximately equal to a particular real float.
template <class P> inline bool 
isNumericallyEqual(const float& a, const std::complex<P>& b,
                   double tol = RTraits<float>::getDefaultTolerance())
{   return isNumericallyEqual(b,a,tol); }
/// %Test whether a complex number is approximately equal to a particular real double.
template <class P> inline bool 
isNumericallyEqual(const std::complex<P>& a, const double& b, 
                   double tol = RTraits<typename Narrowest<P,double>::Precision>::getDefaultTolerance())
{   return isNumericallyEqual(a.real(),b,tol) && isNumericallyEqual(a.imag(),0.,tol); }
/// %Test whether a complex number is approximately equal to a particular real double.
template <class P> inline bool 
isNumericallyEqual(const double& a, const std::complex<P>& b,
                   double tol = RTraits<typename Narrowest<P,double>::Precision>::getDefaultTolerance())
{   return isNumericallyEqual(b,a,tol); }
/// %Test whether a complex number is approximately equal to a particular integer.
template <class P> inline bool 
isNumericallyEqual(const std::complex<P>& a, int b, 
                   double tol = RTraits<P>::getDefaultTolerance())
{   typedef typename Widest<P,double>::Precision W; return isNumericallyEqual(a,(W)b,tol); }
/// %Test whether a complex number is approximately equal to a particular integer.
template <class P> inline bool 
isNumericallyEqual(int a, const std::complex<P>& b, 
                   double tol = RTraits<P>::getDefaultTolerance())
{   return isNumericallyEqual(b,a,tol); }

/// %Test whether a conjugate number is approximately equal to a particular real float.
template <class P> inline bool 
isNumericallyEqual(const conjugate<P>& a, const float& b, 
                   double tol = RTraits<float>::getDefaultTolerance())
{   return isNumericallyEqual(a.real(),b,tol) && isNumericallyEqual(a.imag(),0.f,tol); }
/// %Test whether a conjugate number is approximately equal to a particular real float.
template <class P> inline bool 
isNumericallyEqual(const float& a, const conjugate<P>& b,
                   double tol = RTraits<float>::getDefaultTolerance())
{   return isNumericallyEqual(b,a,tol); }
/// %Test whether a conjugate number is approximately equal to a particular real double.
template <class P> inline bool 
isNumericallyEqual(const conjugate<P>& a, const double& b, 
                   double tol = RTraits<typename Narrowest<P,double>::Precision>::getDefaultTolerance())
{   return isNumericallyEqual(a.real(),b,tol) && isNumericallyEqual(a.imag(),0.,tol); }
/// %Test whether a conjugate number is approximately equal to a particular real double.
template <class P> inline bool 
isNumericallyEqual(const double& a, const conjugate<P>& b,
                   double tol = RTraits<typename Narrowest<P,double>::Precision>::getDefaultTolerance())
{   return isNumericallyEqual(b,a,tol); }
/// %Test whether a conjugate number is approximately equal to a particular integer.
template <class P> inline bool 
isNumericallyEqual(const conjugate<P>& a, int b, 
                   double tol = RTraits<P>::getDefaultTolerance())
{   typedef typename Widest<P,double>::Precision W; return isNumericallyEqual(a,(W)b,tol); }
/// %Test whether a conjugate number is approximately equal to a particular integer.
template <class P> inline bool 
isNumericallyEqual(int a, const conjugate<P>& b, 
                   double tol = RTraits<P>::getDefaultTolerance())
{   return isNumericallyEqual(b,a,tol); }

//@}


template <class N> class NTraits { 
    // only the specializations below are allowed 
};

/// Partial specialization for complex numbers -- underlying real R is
/// still a template parameter.
template <class R> class NTraits< complex<R> > {
    typedef complex<R>  C;
public:
    typedef C                T;
    typedef negator<C>       TNeg;            // type of this after *cast* to its negative
    typedef C                TWithoutNegator; // type of this ignoring negator (there isn't one!)

    typedef R                TReal;
    typedef R                TImag;
    typedef C                TComplex;    
    typedef conjugate<R>     THerm;     // this is a recast
    typedef C                TPosTrans;
    typedef R                TSqHermT;  // ~C*C
    typedef R                TSqTHerm;  // C*~C (same)
    typedef C                TElement;
    typedef C                TRow;
    typedef C                TCol;

    typedef C                TSqrt;
    typedef R                TAbs;
    typedef C                TStandard; // complex is a standard type
    typedef C                TInvert;   // this is a calculation, so use standard number
    typedef C                TNormalize;

    typedef C                Scalar;
    typedef C                ULessScalar;
    typedef C                Number;
    typedef C                StdNumber;
    typedef R                Precision;
    typedef R                ScalarNormSq;

    // For complex scalar C, op result types are:
    //   Typeof(C*P) = Typeof(P*C)
    //   Typeof(C/P) = Typeof(inv(P)*C)
    //   Typeof(C+P) = Typeof(P+C)
    //   typeof(C-P) = Typeof(P::TNeg + C)
    // These must be specialized for P=real, complex, conjugate.
    template <class P> struct Result { 
        typedef typename CNT<P>::template Result<C>::Mul Mul;
        typedef typename CNT< typename CNT<P>::THerm >::template Result<C>::Mul Dvd;
        typedef typename CNT<P>::template Result<C>::Add Add;
        typedef typename CNT< typename CNT<P>::TNeg >::template Result<C>::Add Sub;
    };

    // Shape-preserving element substitution (easy for scalars!)
    template <class P> struct Substitute {
        typedef P Type;
    };

    enum {
        NRows               = 1,
        NCols               = 1,
        RowSpacing          = 1,
        ColSpacing          = 1,
        NPackedElements     = 1,      // not two!
        NActualElements     = 1,
        NActualScalars      = 1,
        ImagOffset          = 1,
        RealStrideFactor    = 2,      // double stride when casting to real or imaginary
        ArgDepth            = SCALAR_DEPTH,
        IsScalar            = 1,
        IsULessScalar       = 1,
        IsNumber            = 1,
        IsStdNumber         = 1,
        IsPrecision         = 0,
        SignInterpretation  = 1       // after cast to Number, what is the sign?
    }; 
    static const T* getData(const T& t) { return &t; } 
    static T*       updData(T& t)       { return &t; }
    static const R& real(const T& t) { return (reinterpret_cast<const R*>(&t))[0]; }
    static R&       real(T& t)       { return (reinterpret_cast<R*>(&t))[0]; }
    static const R& imag(const T& t) { return (reinterpret_cast<const R*>(&t))[1]; }
    static R&       imag(T& t)       { return (reinterpret_cast<R*>(&t))[1]; }

    static const TNeg& negate(const T& t) {return reinterpret_cast<const TNeg&>(t);}
    static       TNeg& negate(T& t)       {return reinterpret_cast<TNeg&>(t);}

    static const THerm& transpose(const T& t) {return reinterpret_cast<const THerm&>(t);}
    static       THerm& transpose(T& t)       {return reinterpret_cast<THerm&>(t);}

    static const TPosTrans& positionalTranspose(const T& t)
        {return reinterpret_cast<const TPosTrans&>(t);}
    static       TPosTrans& positionalTranspose(T& t)
        {return reinterpret_cast<TPosTrans&>(t);} 

    static const TWithoutNegator& castAwayNegatorIfAny(const T& t)
        {return reinterpret_cast<const TWithoutNegator&>(t);}
    static       TWithoutNegator& updCastAwayNegatorIfAny(T& t)
        {return reinterpret_cast<TWithoutNegator&>(t);}

    static ScalarNormSq scalarNormSqr(const T& t)
        { return t.real()*t.real() + t.imag()*t.imag(); }
    static TSqrt    sqrt(const T& t)
        { return std::sqrt(t); }
    static TAbs     abs(const T& t)
        { return std::abs(t); } // no, not just sqrt of scalarNormSqr()!
    static const TStandard& standardize(const T& t) {return t;} // already standard
    static TNormalize normalize(const T& t) {return t/abs(t);}
    static TInvert    invert(const T& t)    {return TReal(1)/t;}

    static const T& getNaN() {
        static const T c=T(NTraits<R>::getNaN(), NTraits<R>::getNaN());
        return c;
    }
    static const T& getInfinity() {
        static const T c=T(NTraits<R>::getInfinity(),NTraits<R>::getInfinity());
        return c;
    }

    static const T& getI() {
        static const T c = T(0,1);
        return c;
    }

    static bool isFinite(const T& t) {return SimTK::isFinite(t);}
    static bool isNaN(const T& t) {return SimTK::isNaN(t);}
    static bool isInf(const T& t) {return SimTK::isInf(t);}

    static double getDefaultTolerance() {return RTraits<R>::getDefaultTolerance();}

    template <class R2> static bool isNumericallyEqual(const T& a, const complex<R2>& b)
    {   return SimTK::isNumericallyEqual(a,b); }
    template <class R2> static bool isNumericallyEqual(const T& a, const complex<R2>& b, double tol)
    {   return SimTK::isNumericallyEqual(a,b,tol); }
    template <class R2> static bool isNumericallyEqual(const T& a, const conjugate<R2>& b)
    {   return SimTK::isNumericallyEqual(a,b); }
    template <class R2> static bool isNumericallyEqual(const T& a, const conjugate<R2>& b, double tol)
    {   return SimTK::isNumericallyEqual(a,b,tol); }

    static bool isNumericallyEqual(const T& a, const float& b) {return SimTK::isNumericallyEqual(a,b);}
    static bool isNumericallyEqual(const T& a, const float& b, double tol) {return SimTK::isNumericallyEqual(a,b,tol);}
    static bool isNumericallyEqual(const T& a, const double& b) {return SimTK::isNumericallyEqual(a,b);}
    static bool isNumericallyEqual(const T& a, const double& b, double tol) {return SimTK::isNumericallyEqual(a,b,tol);}
    static bool isNumericallyEqual(const T& a, int b) {return SimTK::isNumericallyEqual(a,b);}
    static bool isNumericallyEqual(const T& a, int b, double tol) {return SimTK::isNumericallyEqual(a,b,tol);}

    // The rest are the same as the real equivalents, with zero imaginary part.              
    static const T& getZero()         {static const T c(NTraits<R>::getZero());         return c;}
    static const T& getOne()          {static const T c(NTraits<R>::getOne());          return c;}
    static const T& getMinusOne()     {static const T c(NTraits<R>::getMinusOne());     return c;}
    static const T& getTwo()          {static const T c(NTraits<R>::getTwo());          return c;}
    static const T& getThree()        {static const T c(NTraits<R>::getThree());        return c;}
    static const T& getOneHalf()      {static const T c(NTraits<R>::getOneHalf());      return c;}
    static const T& getOneThird()     {static const T c(NTraits<R>::getOneThird());     return c;}
    static const T& getOneFourth()    {static const T c(NTraits<R>::getOneFourth());    return c;}
    static const T& getOneFifth()     {static const T c(NTraits<R>::getOneFifth());     return c;}
    static const T& getOneSixth()     {static const T c(NTraits<R>::getOneSixth());     return c;}
    static const T& getOneSeventh()   {static const T c(NTraits<R>::getOneSeventh());   return c;}
    static const T& getOneEighth()    {static const T c(NTraits<R>::getOneEighth());    return c;}
    static const T& getOneNinth()     {static const T c(NTraits<R>::getOneNinth());     return c;}
    static const T& getPi()           {static const T c(NTraits<R>::getPi());           return c;}
    static const T& getOneOverPi()    {static const T c(NTraits<R>::getOneOverPi());    return c;}
    static const T& getE()            {static const T c(NTraits<R>::getE());            return c;}
    static const T& getLog2E()        {static const T c(NTraits<R>::getLog2E());        return c;}
    static const T& getLog10E()       {static const T c(NTraits<R>::getLog10E());       return c;}
    static const T& getSqrt2()        {static const T c(NTraits<R>::getSqrt2());        return c;}
    static const T& getOneOverSqrt2() {static const T c(NTraits<R>::getOneOverSqrt2()); return c;}
    static const T& getSqrt3()        {static const T c(NTraits<R>::getSqrt3());        return c;}
    static const T& getOneOverSqrt3() {static const T c(NTraits<R>::getOneOverSqrt3()); return c;}
    static const T& getCubeRoot2()    {static const T c(NTraits<R>::getCubeRoot2());    return c;}
    static const T& getCubeRoot3()    {static const T c(NTraits<R>::getCubeRoot3());    return c;}
    static const T& getLn2()          {static const T c(NTraits<R>::getLn2());          return c;}
    static const T& getLn10()         {static const T c(NTraits<R>::getLn10());         return c;}
};


// Specialize NTraits<complex>::Results for <complex> OP <scalar>. Result type is
// always just the complex type of sufficient precision.
#define SimTK_BNTCMPLX_SPEC(T1,T2)  \
template<> template<> struct NTraits< complex<T1> >::Result<T2> {      \
    typedef Widest< complex<T1>,T2 >::Type W;                      \
    typedef W Mul; typedef W Dvd; typedef W Add; typedef W Sub;         \
};                                                                      \
template<> template<> struct NTraits< complex<T1> >::Result< complex<T2> > {  \
    typedef Widest< complex<T1>,complex<T2> >::Type W;        \
    typedef W Mul; typedef W Dvd; typedef W Add; typedef W Sub;         \
};                                                                      \
template<> template<> struct NTraits< complex<T1> >::Result< conjugate<T2> > {  \
    typedef Widest< complex<T1>,complex<T2> >::Type W;        \
    typedef W Mul; typedef W Dvd; typedef W Add; typedef W Sub;         \
}
SimTK_BNTCMPLX_SPEC(float,float);SimTK_BNTCMPLX_SPEC(float,double);
SimTK_BNTCMPLX_SPEC(double,float);SimTK_BNTCMPLX_SPEC(double,double);
#undef SimTK_BNTCMPLX_SPEC


// conjugate -- should be instantiated only for float, double.
template <class R> class NTraits< conjugate<R> > {
    typedef complex<R>          C;
public:
    typedef conjugate<R>        T;
    typedef negator<T>          TNeg;            // type of this after *cast* to its negative
    typedef conjugate<R>        TWithoutNegator; // type of this ignoring negator (there isn't one!)
    typedef R                   TReal;
    typedef negator<R>          TImag;
    typedef conjugate<R>        TComplex;     
    typedef complex<R>          THerm;      // conjugate evaporates
    typedef conjugate<R>        TPosTrans;  // Positional transpose of scalar does nothing
    typedef R                   TSqHermT;   // C*C~
    typedef R                   TSqTHerm;   // ~C*C (same)
    typedef conjugate<R>        TElement;
    typedef conjugate<R>        TRow;
    typedef conjugate<R>        TCol;

    typedef complex<R>          TSqrt;
    typedef R                   TAbs;
    typedef complex<R>          TStandard;
    typedef conjugate<R>        TInvert;
    typedef conjugate<R>        TNormalize;

    typedef conjugate<R>        Scalar;
    typedef conjugate<R>        ULessScalar;
    typedef conjugate<R>        Number;
    typedef complex<R>          StdNumber;
    typedef R                   Precision;
    typedef R                   ScalarNormSq;

    // Typeof( Conj<S>*P ) is Typeof(P*Conj<S>)
    // Typeof( Conj<S>/P ) is Typeof(inv(P)*Conj<S>)
    // Typeof( Conj<S>+P ) is Typeof(P+Conj<S>)
    // Typeof( Conj<S>-P ) is Typeof(P::TNeg+Conj<S>)
    // Must specialize for P=real or P=complex or P=conjugate
    template <class P> struct Result {
        typedef typename CNT<P>::template Result<T>::Mul Mul;
        typedef typename CNT<typename CNT<P>::THerm>::template Result<T>::Mul Dvd;
        typedef typename CNT<P>::template Result<T>::Add Add;
        typedef typename CNT<typename CNT<P>::TNeg>::template Result<T>::Add Sub;
    };

    // Shape-preserving element substitution (easy for scalars!)
    template <class P> struct Substitute {
        typedef P Type;
    };

    enum {
        NRows               = 1,
        NCols               = 1,
        RowSpacing          = 1,
        ColSpacing          = 1,
        NPackedElements     = 1,      // not two!
        NActualElements     = 1,
        NActualScalars      = 1,
        ImagOffset          = 1,
        RealStrideFactor    = 2,      // double stride when casting to real or imaginary
        ArgDepth            = SCALAR_DEPTH,
        IsScalar            = 1,
        IsULessScalar       = 1,
        IsNumber            = 1,
        IsStdNumber         = 0,
        IsPrecision         = 0,
        SignInterpretation  = 1       // after cast to Number, what is the sign?
    }; 

    static const T*     getData(const T& t) { return &t; } 
    static T*           updData(T& t)       { return &t; }
    static const TReal& real(const T& t) { return t.real(); }
    static TReal&       real(T& t)       { return t.real(); }
    static const TImag& imag(const T& t) { return t.imag(); }
    static TImag&       imag(T& t)       { return t.imag(); }

    static const TNeg& negate(const T& t) {return reinterpret_cast<const TNeg&>(t);}
    static       TNeg& negate(T& t)       {return reinterpret_cast<TNeg&>(t);}

    static const THerm& transpose(const T& t) {return t.conj();}
    static       THerm& transpose(T& t)       {return t.conj();}

    static const TPosTrans& positionalTranspose(const T& t)
        {return reinterpret_cast<const TPosTrans&>(t);}
    static       TPosTrans& positionalTranspose(T& t)
        {return reinterpret_cast<TPosTrans&>(t);} 

    static const TWithoutNegator& castAwayNegatorIfAny(const T& t)
        {return reinterpret_cast<const TWithoutNegator&>(t);}
    static       TWithoutNegator& updCastAwayNegatorIfAny(T& t)
        {return reinterpret_cast<TWithoutNegator&>(t);}

    static ScalarNormSq scalarNormSqr(const T& t)
        { return t.real()*t.real() + t.negImag()*t.negImag(); }
    static TSqrt    sqrt(const T& t)
        { return std::sqrt(C(t)); } // cast to complex (one negation)
    static TAbs     abs(const T& t)
        { return std::abs(t.conj()); }  // no, not just sqrt of scalarNormSqr()!
    static TStandard standardize(const T& t)
        { return TStandard(t); }        // i.e., convert to complex
    static TNormalize normalize(const T& t) {return TNormalize(t/abs(t));}

    // 1/conj(z) = conj(1/z), for complex z.
    static TInvert invert(const T& t)    
    {   const typename NTraits<THerm>::TInvert cmplx(NTraits<THerm>::invert(t.conj()));
        return reinterpret_cast<const TInvert&>(cmplx); } // recast complex to conjugate it

    // We want a "conjugate NaN", NaN - NaN*i, meaning both reals should
    // be positive NaN.
    static const T& getNaN() { 
        static const T c=T(NTraits<R>::getNaN(),NTraits<R>::getNaN());
        return c;
    }
    // We want a "conjugate infinity", Inf - Inf*i, meaning both stored reals
    // are positive Inf.
    static const T& getInfinity() {
        static const T c=T(NTraits<R>::getInfinity(),NTraits<R>::getInfinity());
        return c;
    }
    // But we want the constant i (=sqrt(-1)) to be the same however we represent it,
    // so for conjugate i = 0 - (-1)i.
    static const T& getI() {
        static const T c = T(0,-1);
        return c;
    }

    static bool isFinite(const T& t) {return SimTK::isFinite(t);}
    static bool isNaN(const T& t) {return SimTK::isNaN(t);}
    static bool isInf(const T& t) {return SimTK::isInf(t);}

    static double getDefaultTolerance() {return RTraits<R>::getDefaultTolerance();}

    template <class R2> static bool isNumericallyEqual(const T& a, const conjugate<R2>& b)
    {   return SimTK::isNumericallyEqual(a,b); }
    template <class R2> static bool isNumericallyEqual(const T& a, const conjugate<R2>& b, double tol)
    {   return SimTK::isNumericallyEqual(a,b,tol); }
    template <class R2> static bool isNumericallyEqual(const T& a, const complex<R2>& b)
    {   return SimTK::isNumericallyEqual(a,b); }
    template <class R2> static bool isNumericallyEqual(const T& a, const complex<R2>& b, double tol)
    {   return SimTK::isNumericallyEqual(a,b,tol); }

    static bool isNumericallyEqual(const T& a, const float& b) {return SimTK::isNumericallyEqual(a,b);}
    static bool isNumericallyEqual(const T& a, const float& b, double tol) {return SimTK::isNumericallyEqual(a,b,tol);}
    static bool isNumericallyEqual(const T& a, const double& b) {return SimTK::isNumericallyEqual(a,b);}
    static bool isNumericallyEqual(const T& a, const double& b, double tol) {return SimTK::isNumericallyEqual(a,b,tol);}
    static bool isNumericallyEqual(const T& a, int b) {return SimTK::isNumericallyEqual(a,b);}
    static bool isNumericallyEqual(const T& a, int b, double tol) {return SimTK::isNumericallyEqual(a,b,tol);}

    // The rest are the same as the real equivalents, with zero imaginary part.              
    static const T& getZero()         {static const T c(NTraits<R>::getZero());         return c;}
    static const T& getOne()          {static const T c(NTraits<R>::getOne());          return c;}
    static const T& getMinusOne()     {static const T c(NTraits<R>::getMinusOne());     return c;}
    static const T& getTwo()          {static const T c(NTraits<R>::getTwo());          return c;}
    static const T& getThree()        {static const T c(NTraits<R>::getThree());        return c;}
    static const T& getOneHalf()      {static const T c(NTraits<R>::getOneHalf());      return c;}
    static const T& getOneThird()     {static const T c(NTraits<R>::getOneThird());     return c;}
    static const T& getOneFourth()    {static const T c(NTraits<R>::getOneFourth());    return c;}
    static const T& getOneFifth()     {static const T c(NTraits<R>::getOneFifth());     return c;}
    static const T& getOneSixth()     {static const T c(NTraits<R>::getOneSixth());     return c;}
    static const T& getOneSeventh()   {static const T c(NTraits<R>::getOneSeventh());   return c;}
    static const T& getOneEighth()    {static const T c(NTraits<R>::getOneEighth());    return c;}
    static const T& getOneNinth()     {static const T c(NTraits<R>::getOneNinth());     return c;}
    static const T& getPi()           {static const T c(NTraits<R>::getPi());           return c;}
    static const T& getOneOverPi()    {static const T c(NTraits<R>::getOneOverPi());    return c;}
    static const T& getE()            {static const T c(NTraits<R>::getE());            return c;}
    static const T& getLog2E()        {static const T c(NTraits<R>::getLog2E());        return c;}
    static const T& getLog10E()       {static const T c(NTraits<R>::getLog10E());       return c;}
    static const T& getSqrt2()        {static const T c(NTraits<R>::getSqrt2());        return c;}
    static const T& getOneOverSqrt2() {static const T c(NTraits<R>::getOneOverSqrt2()); return c;}
    static const T& getSqrt3()        {static const T c(NTraits<R>::getSqrt3());        return c;}
    static const T& getOneOverSqrt3() {static const T c(NTraits<R>::getOneOverSqrt3()); return c;}
    static const T& getCubeRoot2()    {static const T c(NTraits<R>::getCubeRoot2());    return c;}
    static const T& getCubeRoot3()    {static const T c(NTraits<R>::getCubeRoot3());    return c;}
    static const T& getLn2()          {static const T c(NTraits<R>::getLn2());          return c;}
    static const T& getLn10()         {static const T c(NTraits<R>::getLn10());         return c;}
};

// Any op involving conjugate & a real is best left as a conjugate. However,
// an op involving conjugate & a complex or conjugate can lose the conjugate at zero cost
// and return just a complex in some cases. Also, we prefer negator<complex> to conjugate.
//
// Conj op complex
//   a-bi * r+si = (ar+bs) + (as-br)i               (complex)
//   a-bi / r+si = hairy and slow anyway; we'll convert to complex
//   a-bi + r+si = (a+r) + (s-b)i                   (complex)
//   a-bi - r+si = (a-r) - (b+s)i = -[(r-a)+(b+s)i] (negator<complex>)
//
// Conj op Conj
//   a-bi * r-si = (ar-bs) - (as+br)i = -[(bs-ar)+(as+br)i] (negator<complex>)
//   a-bi / r-si = hairy and slow anyway; we'll convert to complex
//   a-bi + r-si = (a+r) - (b+s)i  (conjugate)
//   a-bi - r-si = (a-r) + (s-b)i  (complex)

#define SimTK_NTRAITS_CONJ_SPEC(T1,T2)                                      \
template<> template<> struct NTraits< conjugate<T1> >::Result<T2> {         \
  typedef conjugate<Widest<T1,T2>::Type> W;                                 \
  typedef W Mul; typedef W Dvd; typedef W Add; typedef W Sub;               \
};                                                                          \
template<> template<> struct NTraits< conjugate<T1> >::Result<complex<T2> >{\
  typedef Widest<complex<T1>,complex<T2> >::Type W;               \
  typedef W Mul; typedef W Dvd; typedef W Add; typedef negator<W> Sub;      \
};                                                                          \
template<> template<> struct NTraits< conjugate<T1> >::Result<conjugate<T2> >{\
    typedef Widest<T1,T2>::Type W; typedef complex<W> WC;              \
    typedef negator<WC> Mul; typedef WC Dvd; typedef conjugate<W> Add; typedef WC Sub;\
}
SimTK_NTRAITS_CONJ_SPEC(float,float);SimTK_NTRAITS_CONJ_SPEC(float,double);
SimTK_NTRAITS_CONJ_SPEC(double,float);SimTK_NTRAITS_CONJ_SPEC(double,double);
#undef SimTK_NTRAITS_CONJ_SPEC 


// Specializations for real numbers.
// For real scalar R, op result types are:
//   Typeof(R*P) = Typeof(P*R)
//   Typeof(R/P) = Typeof(inv(P)*R)
//   Typeof(R+P) = Typeof(P+R)
//   typeof(R-P) = Typeof(P::TNeg + R)
// These must be specialized for P=Real and P=Complex.
#define SimTK_DEFINE_REAL_NTRAITS(R)            \
template <> class NTraits<R> {                  \
public:                                         \
    typedef R                T;                 \
    typedef negator<T>       TNeg;              \
    typedef T                TWithoutNegator;   \
    typedef T                TReal;             \
    typedef T                TImag;             \
    typedef complex<T>       TComplex;          \
    typedef T                THerm;             \
    typedef T                TPosTrans;         \
    typedef T                TSqHermT;          \
    typedef T                TSqTHerm;          \
    typedef T                TElement;          \
    typedef T                TRow;              \
    typedef T                TCol;              \
    typedef T                TSqrt;             \
    typedef T                TAbs;              \
    typedef T                TStandard;         \
    typedef T                TInvert;           \
    typedef T                TNormalize;        \
    typedef T                Scalar;            \
    typedef T                ULessScalar;       \
    typedef T                Number;            \
    typedef T                StdNumber;         \
    typedef T                Precision;         \
    typedef T                ScalarNormSq;      \
    template <class P> struct Result {          \
        typedef typename CNT<P>::template Result<R>::Mul Mul;   \
        typedef typename CNT< typename CNT<P>::THerm >::template Result<R>::Mul Dvd;    \
        typedef typename CNT<P>::template Result<R>::Add Add;   \
        typedef typename CNT< typename CNT<P>::TNeg >::template Result<R>::Add Sub;     \
    };                                          \
    template <class P> struct Substitute {      \
        typedef P Type;                         \
    };                                          \
    enum {                                      \
        NRows               = 1,                \
        NCols               = 1,                \
        RowSpacing          = 1,                \
        ColSpacing          = 1,                \
        NPackedElements     = 1,                \
        NActualElements     = 1,                \
        NActualScalars      = 1,                \
        ImagOffset          = 0,                \
        RealStrideFactor    = 1,                \
        ArgDepth            = SCALAR_DEPTH,     \
        IsScalar            = 1,                \
        IsULessScalar       = 1,                \
        IsNumber            = 1,                \
        IsStdNumber         = 1,                \
        IsPrecision         = 1,                \
        SignInterpretation  = 1                 \
    };                                          \
    static const T* getData(const T& t) { return &t; }  \
    static T*       updData(T& t)       { return &t; }  \
    static const T& real(const T& t) { return t; }      \
    static T&       real(T& t)       { return t; }      \
    static const T& imag(const T&)   { return getZero(); }   \
    static T&       imag(T&)         { assert(false); return *reinterpret_cast<T*>(0); } \
    static const TNeg& negate(const T& t) {return reinterpret_cast<const TNeg&>(t);}        \
    static       TNeg& negate(T& t) {return reinterpret_cast<TNeg&>(t);}                    \
    static const THerm& transpose(const T& t) {return reinterpret_cast<const THerm&>(t);}   \
    static       THerm& transpose(T& t) {return reinterpret_cast<THerm&>(t);}               \
    static const TPosTrans& positionalTranspose(const T& t)                 \
        {return reinterpret_cast<const TPosTrans&>(t);}                     \
    static       TPosTrans& positionalTranspose(T& t)                       \
        {return reinterpret_cast<TPosTrans&>(t);}                           \
    static const TWithoutNegator& castAwayNegatorIfAny(const T& t)          \
        {return reinterpret_cast<const TWithoutNegator&>(t);}               \
    static       TWithoutNegator& updCastAwayNegatorIfAny(T& t)             \
        {return reinterpret_cast<TWithoutNegator&>(t);}                     \
    static ScalarNormSq scalarNormSqr(const T& t) {return t*t;}             \
    static TSqrt        sqrt(const T& t) {return std::sqrt(t);}             \
    static TAbs         abs(const T& t) {return std::abs(t);}               \
    static const TStandard& standardize(const T& t) {return t;}             \
    static TNormalize normalize(const T& t) {return (t>0?T(1):(t<0?T(-1):getNaN()));} \
    static TInvert invert(const T& t) {return T(1)/t;}                      \
    /* properties of this floating point representation, with memory addresses */     \
    static const T& getEps()          {return RTraits<T>::getEps();}                                    \
    static const T& getSignificant()  {return RTraits<T>::getSignificant();}                            \
    static const T& getNaN()          {static const T c=std::numeric_limits<T>::quiet_NaN(); return c;} \
    static const T& getInfinity()     {static const T c=std::numeric_limits<T>::infinity();  return c;} \
    static const T& getLeastPositive(){static const T c=std::numeric_limits<T>::min();       return c;} \
    static const T& getMostPositive() {static const T c=std::numeric_limits<T>::max();       return c;} \
    static const T& getLeastNegative(){static const T c=-std::numeric_limits<T>::min();      return c;} \
    static const T& getMostNegative() {static const T c=-std::numeric_limits<T>::max();      return c;} \
    static const T& getSqrtEps()      {static const T c=std::sqrt(getEps());                 return c;} \
    static const T& getTiny()         {static const T c=std::pow(getEps(), (T)1.25L);        return c;} \
    static bool isFinite(const T& t) {return SimTK::isFinite(t);}   \
    static bool isNaN   (const T& t) {return SimTK::isNaN(t);}      \
    static bool isInf   (const T& t) {return SimTK::isInf(t);}      \
    /* Methods to use for approximate comparisons. Perform comparison in the wider of the two */                \
    /* precisions, using the default tolerance from the narrower of the two precisions.       */                \
    static double getDefaultTolerance() {return RTraits<T>::getDefaultTolerance();}                             \
    static bool isNumericallyEqual(const T& t, const float& f) {return SimTK::isNumericallyEqual(t,f);}         \
    static bool isNumericallyEqual(const T& t, const double& d) {return SimTK::isNumericallyEqual(t,d);}        \
    static bool isNumericallyEqual(const T& t, int i) {return SimTK::isNumericallyEqual(t,i);}                  \
    /* Here the tolerance is given so we don't have to figure it out. */                                                        \
    static bool isNumericallyEqual(const T& t, const float& f, double tol){return SimTK::isNumericallyEqual(t,f,tol);}          \
    static bool isNumericallyEqual(const T& t, const double& d, double tol){return SimTK::isNumericallyEqual(t,d,tol);}         \
    static bool isNumericallyEqual(const T& t, int i, double tol){return SimTK::isNumericallyEqual(t,i,tol);}                   \
    /* Carefully calculated constants with convenient memory addresses. */               \
    static const T& getZero()         {static const T c=(T)(0);               return c;} \
    static const T& getOne()          {static const T c=(T)(1);               return c;} \
    static const T& getMinusOne()     {static const T c=(T)(-1);              return c;} \
    static const T& getTwo()          {static const T c=(T)(2);               return c;} \
    static const T& getThree()        {static const T c=(T)(3);               return c;} \
    static const T& getOneHalf()      {static const T c=(T)(0.5L);            return c;} \
    static const T& getOneThird()     {static const T c=(T)(1.L/3.L);         return c;} \
    static const T& getOneFourth()    {static const T c=(T)(0.25L);           return c;} \
    static const T& getOneFifth()     {static const T c=(T)(0.2L);            return c;} \
    static const T& getOneSixth()     {static const T c=(T)(1.L/6.L);         return c;} \
    static const T& getOneSeventh()   {static const T c=(T)(1.L/7.L);         return c;} \
    static const T& getOneEighth()    {static const T c=(T)(0.125L);          return c;} \
    static const T& getOneNinth()     {static const T c=(T)(1.L/9.L);         return c;} \
    static const T& getPi()           {static const T c=(T)(SimTK_PI);        return c;} \
    static const T& getOneOverPi()    {static const T c=(T)(1.L/SimTK_PI);    return c;} \
    static const T& getE()            {static const T c=(T)(SimTK_E);         return c;} \
    static const T& getLog2E()        {static const T c=(T)(SimTK_LOG2E);     return c;} \
    static const T& getLog10E()       {static const T c=(T)(SimTK_LOG10E);    return c;} \
    static const T& getSqrt2()        {static const T c=(T)(SimTK_SQRT2);     return c;} \
    static const T& getOneOverSqrt2() {static const T c=(T)(1.L/SimTK_SQRT2); return c;} \
    static const T& getSqrt3()        {static const T c=(T)(SimTK_SQRT3);     return c;} \
    static const T& getOneOverSqrt3() {static const T c=(T)(1.L/SimTK_SQRT3); return c;} \
    static const T& getCubeRoot2()    {static const T c=(T)(SimTK_CBRT2);     return c;} \
    static const T& getCubeRoot3()    {static const T c=(T)(SimTK_CBRT3);     return c;} \
    static const T& getLn2()          {static const T c=(T)(SimTK_LN2);       return c;} \
    static const T& getLn10()         {static const T c=(T)(SimTK_LN10);      return c;} \
    /* integer digit counts useful for formatted input and output */                     \
    static int getNumDigits()         {static const int c=(int)(std::log10(1/getEps()) -0.5); return c;} \
    static int getLosslessNumDigits() {static const int c=(int)(std::log10(1/getTiny())+0.5); return c;} \
}; \
template<> struct NTraits<R>::Result<float> \
  {typedef Widest<R,float>::Type Mul;typedef Mul Dvd;typedef Mul Add;typedef Mul Sub;};    \
template<> struct NTraits<R>::Result<double> \
  {typedef Widest<R,double>::Type Mul;typedef Mul Dvd;typedef Mul Add;typedef Mul Sub;};    \
template<> struct NTraits<R>::Result<complex<float> > \
  {typedef Widest<R,complex<float> >::Type Mul;typedef Mul Dvd;typedef Mul Add;typedef Mul Sub;}; \
template<> struct NTraits<R>::Result<complex<double> > \
  {typedef Widest<R,complex<double> >::Type Mul;typedef Mul Dvd;typedef Mul Add;typedef Mul Sub;}; \
template<> struct NTraits<R>::Result<conjugate<float> > \
  {typedef conjugate<Widest<R,float>::Type> Mul;typedef Mul Dvd;typedef Mul Add;typedef Mul Sub;}; \
template<> struct NTraits<R>::Result<conjugate<double> > \
  {typedef conjugate<Widest<R,double>::Type> Mul;typedef Mul Dvd;typedef Mul Add;typedef Mul Sub;}

#if defined(__clang__)
#pragma clang diagnostic push
// The function `T& imag(T&)` generates a null-dereference warning.
#pragma clang diagnostic ignored "-Wnull-dereference"
#endif
SimTK_DEFINE_REAL_NTRAITS(float);
SimTK_DEFINE_REAL_NTRAITS(double);
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

#undef SimTK_DEFINE_REAL_NTRAITS

/// Specializations of CNT for numeric types.
template <class R> class CNT< complex<R> > : public NTraits< complex<R> > { };
template <class R> class CNT< conjugate<R> > : public NTraits< conjugate<R> > { };
template <> class CNT<float> : public NTraits<float> { };
template <> class CNT<double> : public NTraits<double> { };


} // namespace SimTK

#endif //SimTK_SIMMATRIX_NTRAITS_H_
