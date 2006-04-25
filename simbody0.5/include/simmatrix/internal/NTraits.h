#ifndef SimTK_SIMMATRIX_NTRAITS_H_
#define SimTK_SIMMATRIX_NTRAITS_H_

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

/** @file
 * This file contains classes and typedefs needed to provide uniform
 * handling of floating point numeric values. There are three numeric
 * classes: Real, Complex, conjugate and each comes in float, double,
 * and long double precision. Each of these may be modified by a
 * negator, which does not change the in-memory *representation* but 
 * negates the *interpretation*. Thus there are 18 distinct scalar
 * types: 3 precisions each of real, complex, and conjugate and their
 * Negators.
 *
 *      The Scalar Types
 *      ----------------
 *      Here is a complete taxonomy of the scalar types we support.
 *
 *      <scalar>    ::= <number> | negator< <number> >
 *      <number>    ::= <standard> | <conjugate>
 *      <standard>  ::= <real> | <complex>
 *
 *      <real>      ::= float | double | long double
 *      <complex>   ::= complex< <real> >
 *      <conjugate> ::= SimTK::conjugate< <real> >
 *
 *   
 */

#include "SimTKcommon.h"
#include "simmatrix/internal/common.h"
#include "simmatrix/internal/CompositeNumericalTypes.h"

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
template <class R> class conjugate;    // Only defined for float, double, long double

// Specializations of this class provide information about Composite Numerical Types
// in the style of std::numeric_limits<T>. It is specialized for the numeric types
// but can be invoked on any composite Basic Numerical Type as well.
template <class T> class CNT;

// NTraits provides information like CNT, but for numeric types only.
template <class N> class NTraits; // allowed only for N=<number>
template <class R> class NTraits< complex<R> >;
template <class R> class NTraits< conjugate<R> >;
template <> class NTraits<float>;
template <> class NTraits<double>;
template <> class NTraits<long double>;

// This is an adaptor for numeric types which negates the apparent values. A
// negator<N> has exactly the same internal representation as a numeric
// value N, but it is to be interpreted has having the negative of the value
// it would have if interpreted as an N. This permits negation to be done
// by reinterpretation rather than compuation. A full set of arithmetic operators
// are provided involving negator<N>'s and N's. Sometimes we can save an op or
// two this way. For example negator<N>*negator<N> can be performed as an N*N
// since the negations cancel, and we saved two floating point negations.
template <class N> class negator;      // Only defined for numbers

// This is here so we can provide references to 0's when needed, e.g.
// when returning the imaginary part of a real number. These are local to
// the compilation unit, so the returned addresses will differ in different
// files. There are enough zeroes here for the widest number, complex<long double>
// (or conjugate<long double>).
static const complex<long double> zeroes(0);

// This class is specialized for all 36 combinations of standard types
// (that is, real and complex types in each of three precisions)
// and has a single typedef "Type" which is the appropriate "widened"
// type for use when R1 & R2 appear in an operation together. For example,
// if R1=complex<float> and R2=long double, Widest<R1,R2>::Type is
// complex<long double>.
template <class R1, class R2> struct Widest {/* Only defined for built-ins. */};
template <> struct Widest<float,float>              {typedef float Type;};
template <> struct Widest<float,double>             {typedef double Type;};
template <> struct Widest<float,long double>        {typedef long double Type;};
template <> struct Widest<double,float>             {typedef double Type;};
template <> struct Widest<double,double>            {typedef double Type;};
template <> struct Widest<double,long double>       {typedef long double Type;};
template <> struct Widest<long double,float>        {typedef long double Type;};
template <> struct Widest<long double,double>       {typedef long double Type;};
template <> struct Widest<long double,long double>  {typedef long double Type;};
template <class R1, class R2> struct Widest< complex<R1>,complex<R2> >
  { typedef complex< typename Widest<R1,R2>::Type > Type; };
template <class R1, class R2> struct Widest< complex<R1>,R2 >
  { typedef complex< typename Widest<R1,R2>::Type > Type; };
template <class R1, class R2> struct Widest< R1,complex<R2> >
  { typedef complex< typename Widest<R1,R2>::Type > Type; };

template <class N> class SimTK_SIMBODY_API NTraits { /*only the specializations below are allowed*/ };

/// Partial specialization for complex numbers -- underlying real R is
/// still a template parameter.
template <class R> class NTraits< complex<R> > {
    typedef complex<R>  C;
public:
    typedef C                T;
    typedef negator<C>       TNeg;
    typedef R                TAbs;
    typedef R                TReal;
    typedef R                TImag;
    typedef C                TComplex;    
    typedef conjugate<R>     THerm;     // this is a recast
    typedef C                TInvert;   // this is a calculation, so use standard number
    typedef C                TPosTrans;
    typedef R                TSqHermT;  // ~C*C
    typedef R                TSqTHerm;  // C*~C (same)
    typedef C                TElement;
    typedef C                TRow;
    typedef C                TCol;

    typedef C                Scalar;
    typedef C                Number;
    typedef C                StdNumber;
    typedef R                Precision;
    typedef R                ScalarSq;

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
        IsNumber            = 1,
        IsStdNumber         = 1,
        IsPrecision         = 0
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

    static ScalarSq scalarNormSqr(const T& t)
        { return t.real()*t.real() + t.imag()*t.imag(); }
    static TAbs     abs(const T& t)
        { return std::abs(t); } // no, not just sqrt of scalarNormSqr()!

    static T getNaN()      { return T(CNT<R>::getNaN(),     CNT<R>::getNaN()); }
    static T getInfinity() { return T(CNT<R>::getInfinity(),CNT<R>::getInfinity()); }
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
SimTK_BNTCMPLX_SPEC(float,float);SimTK_BNTCMPLX_SPEC(float,double);SimTK_BNTCMPLX_SPEC(float,long double);
SimTK_BNTCMPLX_SPEC(double,float);SimTK_BNTCMPLX_SPEC(double,double);SimTK_BNTCMPLX_SPEC(double,long double);
SimTK_BNTCMPLX_SPEC(long double,float);SimTK_BNTCMPLX_SPEC(long double,double);SimTK_BNTCMPLX_SPEC(long double,long double);
#undef SimTK_BNTCMPLX_SPEC


// conjugate -- should be instantiated only for float, double, long double.
template <class R> class NTraits< conjugate<R> > {
    typedef complex<R>          C;
public:
    typedef conjugate<R>        T;
    typedef negator<T>          TNeg;
    typedef R                   TAbs;
    typedef R                   TReal;
    typedef negator<R>          TImag;
    typedef conjugate<R>        TComplex;     
    typedef complex<R>          THerm;      // conjugate evaporates
    typedef complex<R>          TInvert;
    typedef conjugate<R>        TPosTrans;  // Positional transpose of scalar does nothing
    typedef R                   TSqHermT;   // C*C~
    typedef R                   TSqTHerm;   // ~C*C (same)
    typedef conjugate<R>        TElement;
    typedef conjugate<R>        TRow;
    typedef conjugate<R>        TCol;

    typedef conjugate<R>        Scalar;
    typedef conjugate<R>        Number;
    typedef C                   StdNumber;
    typedef R                   Precision;
    typedef R                   ScalarSq;

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
        IsNumber            = 1,
        IsStdNumber         = 0,
        IsPrecision         = 0
    }; 

    static const T*     getData(const T& t) { return &t; } 
    static T*           updData(T& t)       { return &t; }
    static const TReal& real(const T& t) { return t.real(); }
    static TReal&       real(T& t)       { return t.real(); }
    static const TImag& imag(const T& t) { return t.imag(); }
    static TImag&       imag(T& t)       { return t.imag(); }

    static const TNeg& negate(const T& t) {return reinterpret_cast<const TNeg&>(t);}
    static       TNeg& negate(T& t)       {return reinterpret_cast<TNeg&>(t);}

    static const THerm& transpose(const T& t) {return reinterpret_cast<const THerm&>(t);}
    static       THerm& transpose(T& t)       {return reinterpret_cast<THerm&>(t);}

    static const TPosTrans& positionalTranspose(const T& t)
        {return reinterpret_cast<const TPosTrans&>(t);}
    static       TPosTrans& positionalTranspose(T& t)
        {return reinterpret_cast<TPosTrans&>(t);} 

    static ScalarSq scalarNormSqr(const T& t)
        { return t.real()*t.real() + t.negImag()*t.negImag(); }
    static TAbs     abs(const T& t)
        { return std::abs(t.conj()); } // no, not just sqrt of scalarNormSqr()!

    static conjugate<R> getNaN()      { return conjugate<R>(CNT<C>::getNaN()); }
    static conjugate<R> getInfinity() { 
        const C inf = CNT<C>::getInfinity();        // this is    Inf + Inf*i
        return conjugate<R>(inf.real(),inf.imag()); // so this is Inf - Inf*i
    }
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
SimTK_NTRAITS_CONJ_SPEC(float,long double);
SimTK_NTRAITS_CONJ_SPEC(double,float);SimTK_NTRAITS_CONJ_SPEC(double,double);
SimTK_NTRAITS_CONJ_SPEC(double,long double);
SimTK_NTRAITS_CONJ_SPEC(long double,float);SimTK_NTRAITS_CONJ_SPEC(long double,double);
SimTK_NTRAITS_CONJ_SPEC(long double,long double);
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
    typedef T                TAbs;              \
    typedef T                TReal;             \
    typedef T                TImag;             \
    typedef complex<T>       TComplex;          \
    typedef T                THerm;             \
    typedef T                TInvert;           \
    typedef T                TPosTrans;         \
    typedef T                TSqHermT;          \
    typedef T                TSqTHerm;          \
    typedef T                TElement;          \
    typedef T                TRow;              \
    typedef T                TCol;              \
    typedef T                Scalar;            \
    typedef T                Number;            \
    typedef T                StdNumber;         \
    typedef T                Precision;         \
    typedef T                ScalarSq;          \
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
        IsNumber            = 1,                \
        IsStdNumber         = 1,                \
        IsPrecision         = 1                 \
    };                                          \
    static const T* getData(const T& t) { return &t; }  \
    static T*       updData(T& t)       { return &t; }  \
    static const T& real(const T& t) { return t; }      \
    static T&       real(T& t)       { return t; }      \
    static const T& imag(const T&)   { return reinterpret_cast<const T&>(zeroes); }   \
    static T&       imag(T&)         { assert(false); return *reinterpret_cast<T*>(0); } \
    static const TNeg& negate(const T& t) {return reinterpret_cast<const TNeg&>(t);}        \
    static       TNeg& negate(T& t) {return reinterpret_cast<TNeg&>(t);}                    \
    static const THerm& transpose(const T& t) {return reinterpret_cast<const THerm&>(t);}   \
    static       THerm& transpose(T& t) {return reinterpret_cast<THerm&>(t);}               \
    static const TPosTrans& positionalTranspose(const T& t)                 \
        {return reinterpret_cast<const TPosTrans&>(t);}                     \
    static       TPosTrans& positionalTranspose(T& t)                       \
        {return reinterpret_cast<TPosTrans&>(t);}                           \
    static ScalarSq scalarNormSqr(const T& t) {return t*t;}                 \
    static TAbs     abs(const T& t) {return std::abs(t);}                   \
    static T getNaN()      { return std::numeric_limits<T>::quiet_NaN(); }  \
    static T getInfinity() { return std::numeric_limits<T>::infinity();  }  \
}; \
template<> struct NTraits<R>::Result<float> \
  {typedef Widest<R,float>::Type Mul;typedef Mul Dvd;typedef Mul Add;typedef Mul Sub;};    \
template<> struct NTraits<R>::Result<double> \
  {typedef Widest<R,double>::Type Mul;typedef Mul Dvd;typedef Mul Add;typedef Mul Sub;};    \
template<> struct NTraits<R>::Result<long double> \
  {typedef Widest<R,long double>::Type Mul;typedef Mul Dvd;typedef Mul Add;typedef Mul Sub;};    \
template<> struct NTraits<R>::Result<complex<float> > \
  {typedef Widest<R,complex<float> >::Type Mul;typedef Mul Dvd;typedef Mul Add;typedef Mul Sub;}; \
template<> struct NTraits<R>::Result<complex<double> > \
  {typedef Widest<R,complex<double> >::Type Mul;typedef Mul Dvd;typedef Mul Add;typedef Mul Sub;}; \
template<> struct NTraits<R>::Result<complex<long double> > \
  {typedef Widest<R,complex<long double> >::Type Mul;typedef Mul Dvd;typedef Mul Add;typedef Mul Sub;}; \
template<> struct NTraits<R>::Result<conjugate<float> > \
  {typedef conjugate<Widest<R,float>::Type> Mul;typedef Mul Dvd;typedef Mul Add;typedef Mul Sub;}; \
template<> struct NTraits<R>::Result<conjugate<double> > \
  {typedef conjugate<Widest<R,double>::Type> Mul;typedef Mul Dvd;typedef Mul Add;typedef Mul Sub;}; \
template<> struct NTraits<R>::Result<conjugate<long double> > \
  {typedef conjugate<Widest<R,long double>::Type> Mul;typedef Mul Dvd;typedef Mul Add;typedef Mul Sub;}
SimTK_DEFINE_REAL_NTRAITS(float);
SimTK_DEFINE_REAL_NTRAITS(double);
SimTK_DEFINE_REAL_NTRAITS(long double);
#undef SimTK_DEFINE_REAL_NTRAITS

/// Specializations of CNT for numeric types.
template <class R> class CNT< complex<R> > : public NTraits< complex<R> > { };
template <class R> class CNT< conjugate<R> > : public NTraits< conjugate<R> > { };
template <> class CNT<float> : public NTraits<float> { };
template <> class CNT<double> : public NTraits<double> { };
template <> class CNT<long double> : public NTraits<long double> { };


} // namespace SimTK

#endif //SimTK_SIMMATRIX_NTRAITS_H_
