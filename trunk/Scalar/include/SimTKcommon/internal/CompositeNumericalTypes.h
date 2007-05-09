#ifndef SimTK_SIMMATRIX_COMPOSITE_NUMERICAL_TYPES_H_
#define SimTK_SIMMATRIX_COMPOSITE_NUMERICAL_TYPES_H_

/* Portions copyright (c) 2005-6 Stanford University and Michael Sherman.
 * Contributors:
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
 *
 *     Type           Meaning
 *     -------------  ------------------------------------------------
 *
 *          REFERENCE TYPES
 *
 *     RefType           type of self
 *     RowRefType        type of a row of this CNT
 *     ColumnRefType     type of a column of this CNT
 *     ElementRefType    type of elements
 *     ScalarRefType     type of stored data
 *
 *     TransposeRefType  type returned by (Hermitian) transpose() or the
 *                         operator~; this is a *cast*, not a copy so will
 *                         have weird spacing if Type does.
 *     PositionalTransposeRefType
 *                       type returned by positionalTranspose(). This is
 *                         a *cast*, not a copy, so will have weird
 *                         spacing if Type does.
 *     NegateRefType     a type *cast*, which negates the interpretation
 *                         of the CNT's data
 *     RealPartRefType   the type *cast* which is used by real() to 
 *                         extract the real part of this CNT if it is
 *                         complex (or conjugate).
 *     ImagPartRefType   the type *cast* which is used by imag() to 
 *                         extract the imaginary part of this CNT if
 *                         it is complex (or conjugate). Note that this
 *                         is not necessarily the same type as RealPartType;
 *                         they can differ by a negator<>.
 *
 *          RESULTS (PACKED) TYPES -- these are always packed & use std numbers
 *
 *     PackedType        least weird type that can hold a copy of this
 *                         CNT object. The shape, CNT type, and numerical values
 *                         will be unchanged, but elements will be packed
 *                         together in columns and will use standard numbers (real &
 *                         complex). ElementRefType::PackedType (ElementType) is used
 *                         recursively for the elements.
 *     RowPackedType     same as above but packing by rows instead of columns
 *
 *     RowType           RowRefType::PackedType    
 *     ColumnType        ColumnRefType::PackedType
 *     ElementType       ElementRefType::PackedType
 *     ScalarType        ScalarRefType::PackedType
 *     TransposeType     TransposeRefType::PackedType
 *     PositionalTransposeType
 *                       PositionalTransposeRefType::PackedType
 *
 *     RealType          same shape as CNT, but elements are real
 *                         same as RealPartRefType::PackedType and ImagPartRefType::PackedType
 *     ComplexType       same shape as CNT, but elements are complex
 *
 *     InverseType       result of an invert() applied to this CNT. Looks like
 *                         a cleaned up version of TransposeType.
 *     SquareType        type returned by square()=~T*T; symmetric, scalar
 *                         if T=Vec
 *     RowSquareType     type returned by rowSquare()=T*~T; symmetric, scalar
 *                         if T=Row
 *     ScalarNormType    the type of ~s*s where s is the ScalarType of this CNT.
 *                         This is always a standard real number. This is the
 *                         result type of scalar norms, and abs(s).
 *     AbsType           type returned by elementwise abs() when applied to this
 *                         CNT. It has the same shape but each element is replaced
 *                         by abs() of that element. The result is always real.
 *                      
 *
 *     TNeg         same shape as T, but elements are negated
 *     TReal        same shape as T, but with real elements
 *     TImag        same shape as T, with real elements from the imaginary part
 *     TComplex     same shape as T, but with Complex or conjugate elements
 *     THerm        transpose of T, with Hermitian transposed elements
 *     TPosTrans    positional transpose of T, that is, elements not transposed
 *     TSqHermT     type of ~T*T (default vector and matrix square; symmetric)
 *     TSqTHerm     type of T*~T (row square; symmetric)
 *
 *     Scalar       the underlying <scalar> type (see below)
 *     ScalarSq     type of square of underlying scalar (always real)
 *
 *     Substitute<E>::Type
 *                    a CNT of the same shape and container type as this one,
 *                      but with elements of type E instead of ElementType.
 *                      Special case: if this CNT is a scalar then 
 *                      Substitute<E>::Type just returns E.
 *     Result<RHS>::Mul ::Dvd
 *                ::Add ::Sub
 *                  the type of the result of T op RHS, where RHS is *any* CNT
 *
 *          ENUMS (all sizes are in units of T's elements)
 *
 *     NRows           logical number of rows in type T (i.e., # elements in a column)
 *     NCols           logical number of columns in type T
 *     RowSpacing      # elements from one row to the next (default 1)
 *     ColSpacing      # elements from one col to the next (default NRows for Mat)
 *     NPackedElements minimum #elements it would take to store this data
 *     NActualElements #elements covered by T due to element spacing
 *     NActualScalars  NActualElements * CNT<ElementType>::NActualScalars. This should
 *                       be the physical spacing between array elements in an array
 *                       containing this kind of CNT. Our big Matrix/Vector types guarantee
 *                       this packing.
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
    typedef typename K::TAbs         TAbs;
    typedef typename K::TStandard    TStandard;     // packed, StdNumbers
    typedef typename K::TInvert      TInvert;
    typedef typename K::TNormalize   TNormalize;

    typedef typename K::Scalar       Scalar;        // <number> or negator<number>
    typedef typename K::Number       Number;        // <real>, <complex> or <conjugate>
    typedef typename K::StdNumber    StdNumber;     // <real>, <complex>
    typedef typename K::Precision    Precision;     // float, double, long double

    typedef typename K::ScalarSq     ScalarSq;      // type of square of scalar or
                                                    //   numeric value

    template <class P> struct Result {
        typedef typename K::template Result<P>::Mul Mul;
        typedef typename K::template Result<P>::Dvd Dvd;;
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
        IsScalar            = K::IsScalar,          // real, complex, conjugate, negator
        IsNumber            = K::IsNumber,          // real, complex, conjugate
        IsStdNumber         = K::IsStdNumber,       // real, complex
        IsPrecision         = K::IsPrecision,       // real (float, double, long double)
        SignInterpretation  = K::SignInterpretation // 1 normally, -1 if elements are negated
    };

    static const Scalar* getData(const T& t) { return t.getData(); }
    static       Scalar* updData(T& t)       { return t.updData(); }

    static const TReal& real(const T& t) { return t.real(); }
    static       TReal& real(T& t)       { return t.real(); }
    static const TImag& imag(const T& t) { return t.imag(); }
    static       TImag& imag(T& t)       { return t.imag(); }

    static K getInfinity() { return K::getInfinity(); }
    static K getNaN()      { return K::getNaN();      }

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

    static ScalarSq  scalarNormSqr(const K& t) {return t.scalarNormSqr();}
    static TAbs      abs(const K& t)           {return t.abs();}
    static TStandard standardize(const K& t)   {return t.standardize();}
    static TNormalize normalize(const K& t)    {return t.normalize();}
};

} // namespace SimTK

#endif // SimTK_SIMMATRIX_COMPOSITE_NUMERICAL_TYPES_H_
