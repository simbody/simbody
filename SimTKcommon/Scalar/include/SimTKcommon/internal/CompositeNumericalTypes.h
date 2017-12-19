#ifndef SimTK_SIMMATRIX_COMPOSITE_NUMERICAL_TYPES_H_
#define SimTK_SIMMATRIX_COMPOSITE_NUMERICAL_TYPES_H_

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
 * The purpose of the CNT<T> class is to hide the differences between
 * built-in numerical types and composite ones like Vec<3>. We can
 * decorate all the composite types with whatever information we
 * need, but we cannot add to the built-in types in the same way.
 * So we define a templatized class CNT<T>, where the template
 * parameter can be any composite numerical type whether built in or
 * composite. Then CNT members are used to access information about
 * class T. When T is a built-in, that information comes from 
 * specializations of CNT<T> for CNT<float>, CNT<std::complex<double>>,
 * etc. When T is composite, CNT<T> acts as a pass-through to allow
 * the composite type to provide its own information.
 *
 * Here is the information that must be provided by a SimTK Composite
 * Numerical Type T (or faked up for built-ins). Some of these are
 * given friendly names and documented since they are useful in
 * user code.
 * <pre>
 *
 *     Type             Meaning
 *     ---------------- --------------------------------------------------------
 *                      
 *     TNeg             same shape as T, but elements are negated
 *     TReal            same shape as T, but with real elements
 *     TImag            same shape as T, with real elements from imaginary part
 *     TComplex         same shape as T, but with Complex or conjugate elements
 *     THerm            transpose of T, with Hermitian transposed elements
 *     TPosTrans        positional transpose of T; i.e., elements not transposed
 *     TSqHermT         type of ~T*T (default vector & matrix square; symmetric)
 *     TSqTHerm         type of T*~T (row square; symmetric)
 *
 *     Scalar           the underlying scalar type (see below)
 *     ScalarNormSq     type of the "conjugate square" ~s*s of underlying scalar 
 *                        (always real)
 *
 *     Substitute<E>::Type  
 *                      A CNT of the same shape and container type as this one,
 *                        but with elements of type E instead of ElementType.
 *                        Special case: if this CNT is a scalar then 
 *                        Substitute<E>::Type just returns E.
 *     Result<RHS>::Mul (Dvd,Add,Sub)
 *                      The type of the result of T op RHS, where RHS is *any* CNT
 *
 *          ENUMS (all sizes are in units of T's elements)
 *
 *     NRows            logical num rows in type T (i.e., num elements in a column)
 *     NCols            logical number of columns in type T
 *     RowSpacing       num elements from one row to the next (default 1)
 *     ColSpacing       num elements from one col to the next (default NRows for Mat)
 *     NPackedElements  minimum num elements it would take to store this data
 *     NActualElements  num elements covered by T due to element spacing
 *     NActualScalars   NActualElements * CNT<ElementType>::NActualScalars. This 
 *                        should be the physical spacing between array elements 
 *                        in an array containing this kind of CNT. Our big 
 *                        Matrix/Vector types guarantee this packing.
 * </pre>
 *
 * @verbatim
 *
 * The Scalar Types
 * ----------------
 * Here is a complete taxonomy of the scalar types we support.
 *
 * <scalar>         ::= quantity< unitType, <unitlessScalar> >
 * <unitlessScalar> ::= <number> | negator< <number> >
 * <number>         ::= <standard> | <conjugate>
 * <standard>       ::= <real> | <complex>
 *
 * <real>           ::= float | double
 * <complex>        ::= std::complex< <real> >
 * <conjugate>      ::= SimTK::conjugate< <real> >
 *
 * @endverbatim
 *
 * With this in hand, we can build a clean facility in which scalars,
 * vectors, matrices, vectors of vectors, matrices of vectors of
 * matrices, etc. can all be treated uniformly.
 */

#include "SimTKcommon/internal/common.h"
    
namespace SimTK {

// These are CNT "depths". 0 means the corresponding CNT is a scalar,
// 1 means it is a composite with scalar elements, 2 means a composite
// with composite elements, and 3 means a composite with depth-2
// composite elements. Beyond that the user will have to diambiguate
// operations by using named routines rather than operators.
enum  {
    SCALAR_DEPTH              = 0,
    SCALAR_COMPOSITE_DEPTH    = 1,
    COMPOSITE_COMPOSITE_DEPTH = 2,
    COMPOSITE_3_DEPTH         = 3,
    MAX_RESOLVED_DEPTH        = COMPOSITE_3_DEPTH
};

/** 
 * Specialized information about Composite Numerical Types which
 * allows us to define appropriate templatized classes using them.
 * Transpose is particularly tricky -- we insist on Hermitian 
 * transpose meaning the elements must also be transposed and
 * complex subelements must be conjugated.
 * 
 * This class exists because the built-in scalar types don't
 * have the members we need. CNT<> is specialized for those types
 * only; it is just a pass-through for the rest. The idea
 * is to capture everything that has to be specialized here rather
 * than in the template classes which use these types.
 */
template <class K> class CNT : private K {
public:
    typedef K                        T;
    typedef typename K::TNeg         TNeg;
    typedef typename K::TWithoutNegator TWithoutNegator;
    typedef typename K::TReal        TReal;
    typedef typename K::TImag        TImag;
    typedef typename K::TComplex     TComplex;
    typedef typename K::THerm        THerm;
    typedef typename K::TPosTrans    TPosTrans;
    typedef typename K::TSqHermT     TSqHermT;
    typedef typename K::TSqTHerm     TSqTHerm;
    typedef typename K::TElement     TElement;
    typedef typename K::TRow         TRow;          // type of a row or column
    typedef typename K::TCol         TCol;

    // These are the results of calculations and should be packed regardless
    // of the spacing of this CNT.
    typedef typename K::TSqrt        TSqrt;         // also turns unit^2 to unit
    typedef typename K::TAbs         TAbs;
    typedef typename K::TStandard    TStandard;     // packed, StdNumbers
    typedef typename K::TInvert      TInvert;       // also turns units into 1/units
    typedef typename K::TNormalize   TNormalize;    // TODO: what effect on units?

    typedef typename K::Scalar       Scalar;        // quantity< units, <unitlessScalar> >
    typedef typename K::ULessScalar  ULessScalar;   // <number> or negator<number>
    typedef typename K::Number       Number;        // <real>, <complex> or <conjugate>
    typedef typename K::StdNumber    StdNumber;     // <real>, <complex>
    typedef typename K::Precision    Precision;     // float, double

    typedef typename K::ScalarNormSq ScalarNormSq;  // type of conjugate square of underlying scalar or
                                                    //   numeric value (squares the units too)

    template <class P> struct Result {
        typedef typename K::template Result<P>::Mul Mul;
        typedef typename K::template Result<P>::Dvd Dvd;
        typedef typename K::template Result<P>::Add Add;
        typedef typename K::template Result<P>::Sub Sub;
    };

    // Shape-preserving element substitution
    template <class P> struct Substitute {
        typedef typename K::template Substitute<P>::Type Type;
    };

    enum {
        NRows               = K::NRows,
        NCols               = K::NCols,
        RowSpacing          = K::RowSpacing,
        ColSpacing          = K::ColSpacing,
        NPackedElements     = K::NPackedElements,
        NActualElements     = K::NActualElements,
        NActualScalars      = K::NActualScalars,
        ImagOffset          = K::ImagOffset,
        RealStrideFactor    = K::RealStrideFactor,
        ArgDepth            = K::ArgDepth,
        IsScalar            = K::IsScalar,          // scalar with units, real, complex, conjugate, negator
        IsULessScalar       = K::IsULessScalar,     // real, complex, conjugate, negator
        IsNumber            = K::IsNumber,          // real, complex, conjugate
        IsStdNumber         = K::IsStdNumber,       // real, complex
        IsPrecision         = K::IsPrecision,       // real (float, double)
        SignInterpretation  = K::SignInterpretation // 1 normally, -1 if elements are negated
    };

    static const Scalar* getData(const T& t) { return t.getData(); }
    static       Scalar* updData(T& t)       { return t.updData(); }

    static const TReal& real(const T& t) { return t.real(); }
    static       TReal& real(T& t)       { return t.real(); }
    static const TImag& imag(const T& t) { return t.imag(); }
    static       TImag& imag(T& t)       { return t.imag(); }

    // We expect to be able to negate and transpose (hermitian or
    // positional) with just type casting; no need for help from class
    // K except to tell us the appropriate types.
    static const TNeg& negate(const T& t)
      { return reinterpret_cast<const TNeg&>(t); }
    static       TNeg& negate(T& t)
      { return reinterpret_cast<TNeg&>(t); }

    static const THerm& transpose(const K& t)
      { return reinterpret_cast<const THerm&>(t); }
    static       THerm& transpose(K& t)
      { return reinterpret_cast<THerm&>(t); }

    static const TPosTrans& positionalTranspose(const K& t)
      { return reinterpret_cast<const TPosTrans&>(t); }
    static       TPosTrans& positionalTranspose(K& t)
      { return reinterpret_cast<TPosTrans&>(t); }

    // If the underlying scalars of this CNT are negator<N> for some numeric type N,
    // this method removes the negator<>, effectively negating the entire CNT. You
    // can still deal with the sign correctly by using the above enum SignInterpretation
    // which will be -1 in that case, 1 if there was no negator<> to remove. Note:
    // I'm not talking about TWithoutNegator::SignInterpretation -- that one is guaranteed
    // to be 1! T::SignInterpretation is the one you want.
    static const TWithoutNegator& castAwayNegatorIfAny(const T& t)
        {return reinterpret_cast<const TWithoutNegator&>(t);}
    static       TWithoutNegator& updCastAwayNegatorIfAny(T& t)
        {return reinterpret_cast<TWithoutNegator&>(t);}

    static ScalarNormSq scalarNormSqr(const K& t) {return t.scalarNormSqr();}

    static TSqrt      sqrt(const K& t)          {return t.sqrt();}
    static TAbs       abs(const K& t)           {return t.abs();}
    static TStandard  standardize(const K& t)   {return t.standardize();}
    static TNormalize normalize(const K& t)     {return t.normalize();}
    static TInvert    invert(const K& t)        {return t.invert();}

    static K getInfinity() {return K::getInfinity();}
    static K getNaN()      {return K::getNaN();}

    /// This is true if any element contains a NaN anywhere.
    static bool isNaN(const K& t) {return t.isNaN();}
    /// This is true if at least one element contains a +Infinity or -Infinity
    /// and no element contains a NaN.
    static bool isInf(const K& t) {return t.isInf();}
    /// This is true only if no element has any entry that it NaN or Infinity.
    static bool isFinite(const K& t) {return t.isFinite();}

    /// CNTs are expected to support an "==" operator for exact, bitwise equality.
    /// This method implements approximate, numerical equality. For scalar types,
    /// this should boil down to the isNumericallyEqual() scalar method. For 2D composite
    /// types, the default tolerance should be loosened from the element's default
    /// tolerance by the shorter of the two dimensions. For example, if element E's
    /// default numerical tolerance is tol, then a Mat<3,5,E>'s default tolerance
    /// should be 3*tol.
    template <class K2> static bool 
    isNumericallyEqual(const K& t1, const K2& t2) 
    {   return t1.isNumericallyEqual(t2);}
    template <class K2> static bool 
    isNumericallyEqual(const K& t1, const K2& t2, double tol)
    {   return t1.isNumericallyEqual(t2,tol);}
    static double getDefaultTolerance() {return K::getDefaultTolerance();}

};

} // namespace SimTK

#endif // SimTK_SIMMATRIX_COMPOSITE_NUMERICAL_TYPES_H_
