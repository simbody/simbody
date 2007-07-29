#ifndef SimTK_SIMMATRIX_NEGATOR_H_
#define SimTK_SIMMATRIX_NEGATOR_H_

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
 * This file defines the negator<N> template which is an adaptor for
 * the numeric types N (Real, Complex, conjugate). negator must NOT
 * be instantiated for anything other than these three types (each
 * of which comes in three precisions). negator<N> is guaranteed to
 * have the same memory layout as N, except that the stored values
 * represent the *negative* value of that number.
 *
 * This is part of the SimTK Scalar package, which forms the basis
 * for composite numerical types like vectors and matrices. The negator
 * class allows negation to be performed by reinterpretation rather
 * than by computation.
 *
 * @verbatim
 *
 * The Scalar Types
 * ----------------
 * Here is a complete taxonomy of the scalar types we support.
 *
 * <scalar>    ::= <number> | negator< <number> >
 * <number>    ::= <standard> | <conjugate>
 * <standard>  ::= <real> | <complex>
 *
 * <real>      ::= float | double | long double
 * <complex>   ::= std::complex< <real> >
 * <conjugate> ::= SimTK::conjugate< <real> >
 *
 * @endverbatim
 */

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/NTraits.h"
#include "SimTKcommon/internal/CompositeNumericalTypes.h"

#include <iostream>
    
namespace SimTK {

// Specializations of this class provide information about Composite Numerical Types
// (i.e. composite types) in the style of std::numeric_limits<T>.
template <class T> class CNT;
template <class N> class NTraits;   // Same, but only defined for numeric types.
template <class N> class negator;   // negator is only defined for numbers.


/**
 * negator<N>, where N is a number type (real, complex, conjugate), is represented in 
 * memory identically to N, but behaves as though multiplied by -1, though at zero
 * cost. Only negators instantiated with the nine number types (real, complex, 
 * conjugate) are allowed.
 */ 
template <class NUMBER> 
class SimTK_SimTKCOMMON_EXPORT negator {
    typedef NUMBER N;
    typedef typename NTraits<N>::TReal      NReal;
    typedef typename NTraits<N>::TImag      NImag;
    typedef typename NTraits<N>::TComplex   NComplex;
    typedef typename NTraits<N>::THerm      NHerm;
public:
    typedef negator<N>                                          T;
    typedef NUMBER                                              TNeg;   // negator evaporates
    typedef NUMBER                                              TWithoutNegator; // "
    typedef typename CNT<NReal>::TNeg                           TReal;
    typedef typename CNT<NImag>::TNeg                           TImag;
    typedef typename CNT<NComplex>::TNeg                        TComplex;
    typedef typename CNT<NHerm>::TNeg                           THerm;
    typedef negator<N>                                          TPosTrans;
    typedef typename NTraits<N>::TSqHermT                       TSqHermT;
    typedef typename NTraits<N>::TSqTHerm                       TSqTHerm;
    typedef negator<N>                                          TElement;
    typedef negator<N>                                          TRow;
    typedef negator<N>                                          TCol;

    typedef typename NTraits<N>::TAbs                           TAbs;
    typedef typename NTraits<N>::TStandard                      TStandard;
    typedef typename NTraits<N>::TInvert                        TInvert;
    typedef typename NTraits<N>::TStandard                      TNormalize; // neg<conj> -> complex


    typedef negator<N>                                          Scalar;
    typedef NUMBER                                              Number;
    typedef typename NTraits<N>::StdNumber                      StdNumber;
    typedef typename NTraits<N>::Precision                      Precision;
    typedef typename NTraits<N>::ScalarSq                       ScalarSq;

    // negator may be used in combination with any composite numerical type, not just
    // numbers. Hence we must use CNT<P> here rather than NTraits<P> (they are 
    // the same when P turns out to be a numeric type).
    //      Typeof( Neg<N>*P ) is Typeof(P*N)::TNeg
    //      Typeof( Neg<N>/P ) is Typeof(N/P)::TNeg
    //      Typeof( Neg<N>+P ) is Typeof(P-N)
    //      Typeof( Neg<N>-P ) is Typeof(P+N)::TNeg
    // No need to specialize because we strip off the "negator" here.
    template <class P> struct Result {
    private:
        // These are the type of the calculations we actually perform.
        typedef typename CNT<P>::template Result<N>::Mul     PMul;
        typedef typename NTraits<N>::template Result<P>::Dvd PDvd;
        typedef typename CNT<P>::template Result<N>::Sub     PAdd;
        typedef typename CNT<P>::template Result<N>::Add     PSub;
    public:
        // These are the types to which we must recast the results.
        typedef typename CNT<PMul>::TNeg    Mul;
        typedef typename CNT<PDvd>::TNeg    Dvd;
        typedef PAdd                        Add;
        typedef typename CNT<PSub>::TNeg    Sub;
    };

    // Shape-preserving element substitution (easy for scalars!)
    template <class P> struct Substitute {
        typedef P Type;
    };

    enum {
        NRows               = 1, // Negators are always scalars
        NCols               = 1,
        RowSpacing          = 1,
        ColSpacing          = 1,
        NPackedElements     = 1,
        NActualElements     = 1,
        NActualScalars      = 1,
        ImagOffset          = NTraits<N>::ImagOffset,
        RealStrideFactor    = NTraits<N>::RealStrideFactor,
        ArgDepth            = SCALAR_DEPTH,
        IsScalar            = 1,
        IsNumber            = 0,
        IsStdNumber         = 0,
        IsPrecision         = 0,
        SignInterpretation  = -1 // if you cast away the negator, don't forget this!
    };
    static negator<N> getNaN()      {return recast(NTraits<N>::getNaN());}
	static negator<N> getInfinity() {return recast(NTraits<N>::getInfinity());}

    const negator<N>* getData() const {return this;}
    negator<N>*       updData()       {return this;}

    const TReal& real() const {return reinterpret_cast<const TReal&>(NTraits<N>::real(v));}
    TReal&       real()       {return reinterpret_cast<      TReal&>(NTraits<N>::real(v));}
    const TImag& imag() const {return reinterpret_cast<const TImag&>(NTraits<N>::imag(v));}
    TImag&       imag()       {return reinterpret_cast<      TImag&>(NTraits<N>::imag(v));}

    ScalarSq   scalarNormSqr() const {return NTraits<N>::scalarNormSqr(v);}
    TAbs       abs()           const {return NTraits<N>::abs(v);}
    TStandard  standardize()   const {return -NTraits<N>::standardize(v);}
    TNormalize normalize()     const {return -NTraits<N>::normalize(v);}

    negator() {
    #ifndef NDEBUG
        v = NTraits<N>::getNaN();
    #endif
    }
    negator(const negator& n) : v(n.v) { }
    negator& operator=(const negator& n) { v=n.v; return *this; }

    // These are explicit conversion from numeric type NN to negator<N>. The value must
    // be unchanged, so we must negate. Note that NN==N is a certainty for one of these cases.

    explicit negator(int t) {v = -N((typename NTraits<N>::Precision)t);}
    explicit negator(const float& t) {v = -N(t);}
    explicit negator(const double& t) {v = -N(t);}
    explicit negator(const long double& t) {v = -N(t);}
    explicit negator(const std::complex<float>& t) {v = -N(t);}
    explicit negator(const std::complex<double>& t) {v = -N(t);}
    explicit negator(const std::complex<long double>& t) {v = -N(t);}
    explicit negator(const conjugate<float>& t) {v = -N(t);}
    explicit negator(const conjugate<double>& t) {v = -N(t);}
    explicit negator(const conjugate<long double>& t) {v = -N(t);}

    // This can be used to negate a value of type N at zero cost. It is typically
    // used for recasting temporary expressions to apply a final negation. Note that
    // this is *not* the same as constructing a negator<N> from an N, which actually
    // peforms a floating point negation.
    static const negator<N>& recast(const N& val)
        { return reinterpret_cast<const negator<N>&>(val); }

    const N& operator-() const { return v;  }
    N&       operator-()       { return v;  } // an lvalue!
    N        operator+() const { return N(-v); } // performs the actual negation (expensive)

    operator N() const { return N(-v); } // implicit conversion to N (expensive)

    template <class P> negator& operator =(const P& t) { v = -t; return *this; }
    template <class P> negator& operator+=(const P& t) { v -= t; return *this; } //swap sign
    template <class P> negator& operator-=(const P& t) { v += t; return *this; }
    template <class P> negator& operator*=(const P& t) { v *= t; return *this; } //don't swap!
    template <class P> negator& operator/=(const P& t) { v /= t; return *this; }

    // If we know we've got a negator as an argument, get rid of its negation and change
    // signs as necessary. We're guaranteed to get rid of at least one negator<> this way.
    // Nothing to gain for multiplication or division, though.
    template <class NN> negator& operator =(const negator<NN>& t) { v =  -t; return *this; }
    template <class NN> negator& operator+=(const negator<NN>& t) { v += -t; return *this; } //swap sign
    template <class NN> negator& operator-=(const negator<NN>& t) { v -= -t; return *this; }

private:
    N v;
};


// Handle all binary numerical operators involving a negator<A> and a B, or negator<A>
// and negator<B>, obtaining results by stripping away the negator<>s and fiddling
// with signs appropriately.
// To appreciate the beauty of these operators, remember that -x is free when x
// is a negator.

template <class DEST, class SRC> static inline const DEST&
negRecast(const SRC& s) { return reinterpret_cast<const DEST&>(s); }

// Binary '+' with a negator as one or both arguments
template <class A, class B> inline typename negator<A>::template Result<B>::Add
operator+(const negator<A>& l, const B& r)
  {return negRecast<typename negator<A>::template Result<B>::Add>(r-(-l));}
template <class A, class B> inline typename CNT<A>::template Result<negator<B> >::Add
operator+(const A& l, const negator<B>& r)
  {return negRecast<typename CNT<A>::template Result<negator<B> >::Add>(l-(-r));}
template <class A, class B> inline typename negator<A>::template Result<negator<B> >::Add
operator+(const negator<A>& l, const negator<B>& r) 
  {return negRecast<typename negator<A>::template Result<negator<B> >::Add>((-l)+(-r)); }

// Binary '-' with a negator as one or both arguments
template <class A, class B> inline typename negator<A>::template Result<B>::Sub
operator-(const negator<A>& l, const B& r)
  {return negRecast<typename negator<A>::template Result<B>::Sub>((-l)+r);}
template <class A, class B> inline typename CNT<A>::template Result<negator<B> >::Sub
operator-(const A& l, const negator<B>& r)
  {return negRecast<typename CNT<A>::template Result<negator<B> >::Sub>(l+(-r));}
template <class A, class B> inline typename negator<A>::template Result<negator<B> >::Sub
operator-(const negator<A>& l, const negator<B>& r) 
  {return negRecast<typename negator<A>::template Result<negator<B> >::Sub>((-r)-(-l));}

// Binary '*' with a negator as one or both arguments
template <class A, class B> inline typename negator<A>::template Result<B>::Mul
operator*(const negator<A>& l, const B& r) 
  {return negRecast<typename negator<A>::template Result<B>::Mul>((-l)*r);}
template <class A, class B> inline typename CNT<A>::template Result<negator<B> >::Mul
operator*(const A& l, const negator<B>& r)
  {return negRecast<typename CNT<A>::template Result<negator<B> >::Mul>(l*(-r));}
template <class A, class B> inline typename negator<A>::template Result<negator<B> >::Mul
operator*(const negator<A>& l, const negator<B>& r)
  {return negRecast<typename negator<A>::template Result<negator<B> >::Mul>((-l)*(-r));}

// Binary '/' with a negator as one or both arguments
template <class A, class B> inline typename negator<A>::template Result<B>::Dvd
operator/(const negator<A>& l, const B& r) 
  {return negRecast<typename negator<A>::template Result<B>::Dvd>((-l)/r);}
template <class A, class B> inline typename CNT<A>::template Result<negator<B> >::Dvd
operator/(const A& l, const negator<B>& r)
  {return negRecast<typename CNT<A>::template Result<negator<B> >::Dvd>(l/(-r));}
template <class A, class B> inline typename negator<A>::template Result<negator<B> >::Dvd
operator/(const negator<A>& l, const negator<B>& r)
  {return negRecast<typename negator<A>::template Result<negator<B> >::Dvd>((-l)/(-r));}

// Binary '==' with a negator as one or both arguments
template <class A, class B> inline bool
operator==(const negator<A>& l, const B& r) { return (A)l == r; }
template <class A, class B> inline bool
operator==(const A& l, const negator<B>& r) { return l == (B)r; }
template <class A, class B> inline bool
operator==(const negator<A>& l, const negator<B>& r) { return (-l) == (-r); }   // cheap!

// The lazy man's '!=' operator.
template <class A, class B> inline bool
operator!=(const negator<A>& l, const B& r) { return !(l==r); }
template <class A, class B> inline bool
operator!=(const A& l, const negator<B>& r) { return !(l==r); }
template <class A, class B> inline bool
operator!=(const negator<A>& l, const negator<B>& r) { return !(l==r); }

// And a final touch of elegance allowing smooth interoperability with iostream.
template <class NUM, class CHAR, class TRAITS> inline std::basic_istream<CHAR,TRAITS>&
operator>>(std::basic_istream<CHAR,TRAITS>& is, negator<NUM>& nn) {
    NUM z; is >> z; nn=z;
    return is;
}
template <class NUM, class CHAR, class TRAITS> inline std::basic_ostream<CHAR,TRAITS>&
operator<<(std::basic_ostream<CHAR,TRAITS>& os, const negator<NUM>& nn) {
    return os << NUM(nn);
}

} // namespace SimTK

#endif //SimTK_SIMMATRIX_NEGATOR_H_
