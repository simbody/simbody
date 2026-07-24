#ifndef SimTK_SIMMATRIX_NEGATOR_H_
#define SimTK_SIMMATRIX_NEGATOR_H_

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
 * <real>      ::= float | double
 * <complex>   ::= std::complex< <real> >
 * <conjugate> ::= SimTK::conjugate< <real> >
 *
 * @endverbatim
 */

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
    typedef typename NTraits<N>::TInvert    NInvert;
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

    typedef typename NTraits<N>::TSqrt                          TSqrt;
    typedef typename NTraits<N>::TAbs                           TAbs;
    typedef typename NTraits<N>::TStandard                      TStandard;
    typedef typename CNT<NInvert>::TNeg                         TInvert;
    typedef typename NTraits<N>::TStandard                      TNormalize; // neg<conj> -> complex


    typedef negator<N>                                          Scalar;
    typedef negator<N>                                          ULessScalar;
    typedef NUMBER                                              Number;
    typedef typename NTraits<N>::StdNumber                      StdNumber;
    typedef typename NTraits<N>::Precision                      Precision;
    typedef typename NTraits<N>::ScalarNormSq                   ScalarNormSq;

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
        IsULessScalar       = 1,
        IsNumber            = 0,
        IsStdNumber         = 0,
        IsPrecision         = 0,
        SignInterpretation  = -1 // if you cast away the negator, don't forget this!
    };
    const negator<N>* getData() const {return this;}
    negator<N>*       updData()       {return this;}

    const TReal& real() const {return reinterpret_cast<const TReal&>(NTraits<N>::real(v));}
    TReal&       real()       {return reinterpret_cast<      TReal&>(NTraits<N>::real(v));}
    const TImag& imag() const {return reinterpret_cast<const TImag&>(NTraits<N>::imag(v));}
    TImag&       imag()       {return reinterpret_cast<      TImag&>(NTraits<N>::imag(v));}

    ScalarNormSq    scalarNormSqr() const {return NTraits<N>::scalarNormSqr(v);}
    TSqrt           sqrt()          const {return NTraits<N>::sqrt(N(v));}
    TAbs            abs()           const {return NTraits<N>::abs(v);}
    TStandard       standardize()   const {return -NTraits<N>::standardize(v);}
    TNormalize      normalize()     const {return -NTraits<N>::normalize(v);}

    // Inverse (1/x) of a non-negated type N will return a non-negated type, so we
    // can cast it to a negated type here to save a flop. The return type might
    // not be N (a negated conjugate comes back as a complex), so there may be
    // a flop done in the final conversion to TInvert.
    TInvert invert() const {return recast(NTraits<N>::invert(v));}

    static negator<N> getNaN()      {return recast(NTraits<N>::getNaN());}
    static negator<N> getInfinity() {return recast(NTraits<N>::getInfinity());}

    /// Returns true if the negated value is finite (i.e., not NaN or Inf).
    inline bool isFinite() const;
    /// Returns true if the negated value contains a NaN.
    inline bool isNaN() const;
    /// Returns true if the negated value contains an Inf or -Inf and does not
    /// contain a NaN.
    inline bool isInf() const;

    static double getDefaultTolerance() {return NTraits<N>::getDefaultTolerance();}

    /// In the generic case we'll perform the negation here to get a number, 
    /// and then delegate to the other type which can be any CNT.
    template <class T2> bool isNumericallyEqual(const T2& t2) const
    {   return CNT<T2>::isNumericallyEqual(t2, -v); } // perform negation

    /// In this partial specialization we know that both types have negators so we
    /// can just compare the underlying numbers, each of which has the reversed sign,
    /// using the global SimTK method available for comparing numbers.
    template <class N2> bool isNumericallyEqual(const negator<N2>& t2) const 
    {   return SimTK::isNumericallyEqual(v, t2.v); }

    /// This is the generic case (see above) but with an explicitly-provided tolerance.
    template <class T2> bool isNumericallyEqual(const T2& t2, double tol) const
    {   return CNT<T2>::isNumericallyEqual(t2, -v, tol); } // perform negation

    /// This is the partially specialized case again (see above) but with an explicitly-provided
    /// tolerance.
    template <class N2> bool isNumericallyEqual(const negator<N2>& t2, double tol) const
    {   return SimTK::isNumericallyEqual(v, t2.v, tol); }


    negator() {
    #ifndef NDEBUG
        v = NTraits<N>::getNaN();
    #endif
    }
    negator(const negator& n) : v(n.v) { }
    negator& operator=(const negator& n) { v=n.v; return *this; }

    // These are implicit conversions from numeric type NN to negator<N>. The
    // value must be unchanged, so we must negate. Note that NN==N is a 
    // certainty for one of these cases.
    negator(int                t) {v = -N((typename NTraits<N>::Precision)t);}
    negator(const float&       t) {v = -N((typename NTraits<N>::Precision)t);}
    negator(const double&      t) {v = -N((typename NTraits<N>::Precision)t);}

    // Some of these may not compile if instantiated -- you can't cast a complex
    // to a float, for example.
    template <class P> negator(const std::complex<P>& t) {v = -N(t);}
    template <class P> negator(const conjugate<P>&    t) {v = -N(t);}

    // This can be used to negate a value of type N at zero cost. It is typically
    // used for recasting temporary expressions to apply a final negation. Note that
    // this is *not* the same as constructing a negator<N> from an N, which actually
    // peforms a floating point negation.
    static const negator<N>& recast(const N& val)
    {   return reinterpret_cast<const negator<N>&>(val); }

    const N& operator-() const { return v;  }
    N&       operator-()       { return v;  } // an lvalue!
    N        operator+() const { return N(-v); } // performs the actual negation (expensive)

    operator N() const { return N(-v); } // implicit conversion to N (expensive)

    template <class P> negator& operator =(const P& t) { v = -t; return *this; }
    template <class P> negator& operator+=(const P& t) { v -= t; return *this; } //swap sign
    template <class P> negator& operator-=(const P& t) { v += t; return *this; }
    template <class P> negator& operator*=(const P& t) { v *= t; return *this; } //don't swap!
    template <class P> negator& operator/=(const P& t) { v /= t; return *this; }

    // If we know we've got a negator as an argument, get rid of its negation 
    // and change signs as necessary. We're guaranteed to get rid of at least 
    // one negator<> this way. Nothing to gain for multiplication or division,
    // though.
    template <class NN> negator& operator =(const negator<NN>& t) 
    {   v =  -t; return *this; }
    template <class NN> negator& operator+=(const negator<NN>& t) 
    {   v += -t; return *this; } //swap sign
    template <class NN> negator& operator-=(const negator<NN>& t) 
    {   v -= -t; return *this; }

private:
    N v;

template <class N2> friend class negator;
};

// isNaN() for real, complex, and conjugate numbers is provided in
// NTraits. Here we add isNaN() for negated scalar types.

/// @addtogroup isNaN
//@{
inline bool isNaN(const negator<float>&  x) {return isNaN(-x);}
inline bool isNaN(const negator<double>& x) {return isNaN(-x);}
template <class P> inline bool
isNaN(const negator< std::complex<P> >& x) {return isNaN(-x);}
template <class P> inline bool
isNaN(const negator< conjugate<P> >&    x) {return isNaN(-x);}
//@}

// isFinite() for real, complex, and conjugate numbers is provided in
// NTraits. Here we add isFinite() for negated scalar types.

/// @addtogroup isFinite
//@{
inline bool isFinite(const negator<float>&  x) {return isFinite(-x);}
inline bool isFinite(const negator<double>& x) {return isFinite(-x);}
template <class P> inline bool
isFinite(const negator< std::complex<P> >& x) {return isFinite(-x);}
template <class P> inline bool
isFinite(const negator< conjugate<P> >&    x) {return isFinite(-x);}
//@}

// isInf(x) for real, complex, and conjugate numbers is provided in
// NTraits. Here we add isInf() for negated scalar types.

/// @addtogroup isInf
//@{
inline bool isInf(const negator<float>&  x) {return isInf(-x);}
inline bool isInf(const negator<double>& x) {return isInf(-x);}
template <class P> inline bool
isInf(const negator< std::complex<P> >& x) {return isInf(-x);}
template <class P> inline bool
isInf(const negator< conjugate<P> >&    x) {return isInf(-x);}
//@}

// The member functions call the global ones just defined.
template <class N> inline bool
negator<N>::isFinite() const {return SimTK::isFinite(*this);}
template <class N> inline bool
negator<N>::isNaN()    const {return SimTK::isNaN(*this);}
template <class N> inline bool
negator<N>::isInf()    const {return SimTK::isInf(*this);}

// Handle all binary numerical operators involving a negator<A> and a B, or negator<A>
// and negator<B>, obtaining results by stripping away the negator<>s and fiddling
// with signs appropriately.
// Careful: don't remove both negators in one step because Result isn't specialized
// for negators so it might not predict the same result type as you actually get.
// Be patient and let it strip one negator at a time -- in Release mode that won't
// add any code since all this stuff drops away.
//
// To appreciate the beauty of these operators, remember that -x is free when x
// is a negator.

template <class DEST, class SRC> static inline const DEST&
negRecast(const SRC& s) { return reinterpret_cast<const DEST&>(s); }

    // Binary '+' with a negator as one or both arguments //
template <class A, class B> inline typename negator<A>::template Result<B>::Add
operator+(const negator<A>& l, const B& r)
  {return negRecast<typename negator<A>::template Result<B>::Add>(r-(-l));}
template <class A, class B> inline typename CNT<A>::template Result<negator<B> >::Add
operator+(const A& l, const negator<B>& r)
  {return negRecast<typename CNT<A>::template Result<negator<B> >::Add>(l-(-r));}
// One step at a time here (see above).
template <class A, class B> inline typename negator<A>::template Result<negator<B> >::Add
operator+(const negator<A>& l, const negator<B>& r) 
  {return negRecast<typename negator<A>::template Result<negator<B> >::Add>(r-(-l)); }

    // Binary '-' with a negator as one or both arguments //
template <class A, class B> inline typename negator<A>::template Result<B>::Sub
operator-(const negator<A>& l, const B& r)
  {return negRecast<typename negator<A>::template Result<B>::Sub>(r+(-l));}
template <class A, class B> inline typename CNT<A>::template Result<negator<B> >::Sub
operator-(const A& l, const negator<B>& r)
  {return negRecast<typename CNT<A>::template Result<negator<B> >::Sub>(l+(-r));}
// One step at a time here (see above).
template <class A, class B> inline typename negator<A>::template Result<negator<B> >::Sub
operator-(const negator<A>& l, const negator<B>& r) 
  {return negRecast<typename negator<A>::template Result<negator<B> >::Sub>(r+(-l));}

// Binary '*' with a negator as one or both arguments
template <class A, class B> inline typename negator<A>::template Result<B>::Mul
operator*(const negator<A>& l, const B& r) 
  {return negRecast<typename negator<A>::template Result<B>::Mul>((-l)*r);}
template <class A, class B> inline typename CNT<A>::template Result<negator<B> >::Mul
operator*(const A& l, const negator<B>& r)
  {return negRecast<typename CNT<A>::template Result<negator<B> >::Mul>(l*(-r));}
// One step at a time here (see above).
template <class A, class B> inline typename negator<A>::template Result<negator<B> >::Mul
operator*(const negator<A>& l, const negator<B>& r)
  {return negRecast<typename negator<A>::template Result<negator<B> >::Mul>((-l)*r);}

// Binary '/' with a negator as one or both arguments
template <class A, class B> inline typename negator<A>::template Result<B>::Dvd
operator/(const negator<A>& l, const B& r) 
  {return negRecast<typename negator<A>::template Result<B>::Dvd>((-l)/r);}
template <class A, class B> inline typename CNT<A>::template Result<negator<B> >::Dvd
operator/(const A& l, const negator<B>& r)
  {return negRecast<typename CNT<A>::template Result<negator<B> >::Dvd>(l/(-r));}
// One step at a time here (see above).
template <class A, class B> inline typename negator<A>::template Result<negator<B> >::Dvd
operator/(const negator<A>& l, const negator<B>& r)
  {return negRecast<typename negator<A>::template Result<negator<B> >::Dvd>((-l)/r);}

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
