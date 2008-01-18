#ifndef SimTK_SIMMATRIX_SMALLMATRIX_SYMMAT_H_
#define SimTK_SIMMATRIX_SMALLMATRIX_SYMMAT_H_

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

/**@file
 * This file declares class SymMat<M, ELEMENT_TYPE, ROW_SPACING>.
 * This is logically a square MxM Hermitian matrix, but only the diagonal
 * and *lower* triangle are stored. Elements above the diagonal are
 * the Hermitan transpose of their mirrored elements below the diagonal.
 *
 * The storage is packed by column, with an optional stride RS when going element-
 * to-element. We use an unconventional storage scheme here which
 * provides substantial conveniences and some performance advantage over the
 * conventional (LAPACK) scheme which stores the matrix in packed
 * columns. Instead, we store the diagonal first, as a Vec<M> and then store
 * the lower triangle separately in the conventional columnwise form. This
 * allows easy handling of the diagonal elements separately, which is very
 * common in a symmetric matrix. The diagonals are "special"; for example
 * they must be "self-Hermitian" meaning their imaginary parts must be zero.
 * We could just store real elements for the diagonals, but we don't for
 * a number of good reasons which I'll leave as an exercise for the reader.
 * However, we will ensure that their imaginary parts are zero, although
 * a determined user could screw this up.
 *
 * Here is an example. Imagine a full 4x4 Hermitian matrix with complex
 * elements. We'll use a' to mean conjugate(a). Then the elements are
 *     w  a' b' c'
 *     a  x  d' e'
 *     b  d  y  f'
 *     c  e  f  z
 * Note that it must be the case that w=w', x=x', y=y', and z=z' so the
 * diagonal elements must be real. Here's how we store that as a SymMat:
 *       w x y z a b c d e f
 * These are stored in consecutive memory locations (possibly separated
 * by meaningless elements as indicated by the row spacing). Mere 
 * recasting allows us to view this as a 4x4 matrix, or a 4-element
 * diagonal and 6-element "lower" vector, or as these vectors 
 * apparently containing the Hermitian types (i.e. one may cast
 * the "lower" portion into what looks like the "upper" portion).
 *
 * This is a "composite numerical type".
 */

#include "SimTKcommon/internal/common.h"

namespace SimTK {

/// RS is total spacing between rows in memory (default 1) 
template <int M, class ELT, int RS> class SymMat {
    typedef ELT                                 E;
    typedef typename CNT<E>::TNeg               ENeg;
    typedef typename CNT<E>::TWithoutNegator    EWithoutNegator;
    typedef typename CNT<E>::TReal              EReal;
    typedef typename CNT<E>::TImag              EImag;
    typedef typename CNT<E>::TComplex           EComplex;
    typedef typename CNT<E>::THerm              EHerm;
    typedef typename CNT<E>::TPosTrans          EPosTrans;
    typedef typename CNT<E>::TSqHermT           ESqHermT;
    typedef typename CNT<E>::TSqTHerm           ESqTHerm;

    typedef typename CNT<E>::TAbs               EAbs;
    typedef typename CNT<E>::TStandard          EStandard;
    typedef typename CNT<E>::TInvert            EInvert;
    typedef typename CNT<E>::TNormalize         ENormalize;

    typedef typename CNT<E>::Scalar             EScalar;
    typedef typename CNT<E>::Number             ENumber;
    typedef typename CNT<E>::StdNumber          EStdNumber;
    typedef typename CNT<E>::Precision          EPrecision;
    typedef typename CNT<E>::ScalarSq           EScalarSq;

public:

    enum {
        NRows               = M,
        NCols               = M,
        NPackedElements     = (M*(M+1))/2,
        NActualElements     = RS * NPackedElements,
        NActualScalars      = CNT<E>::NActualScalars * NActualElements,
        RowSpacing          = RS,
        ColSpacing          = NActualElements,
        ImagOffset          = NTraits<ENumber>::ImagOffset,
        RealStrideFactor    = 1, // composite types don't change size when
                                 // cast from complex to real or imaginary
        ArgDepth            = ((int)CNT<E>::ArgDepth < (int)MAX_RESOLVED_DEPTH 
                                ? CNT<E>::ArgDepth + 1 
                                : MAX_RESOLVED_DEPTH),
        IsScalar            = 0,
        IsNumber            = 0,
        IsStdNumber         = 0,
        IsPrecision         = 0,
        SignInterpretation  = CNT<E>::SignInterpretation
    };

    typedef SymMat<M,E,RS>                  T;
    typedef SymMat<M,ENeg,RS>               TNeg;
    typedef SymMat<M,EWithoutNegator,RS>    TWithoutNegator;

    typedef SymMat<M,EReal,RS*CNT<E>::RealStrideFactor>          
                                            TReal;
    typedef SymMat<M,EImag,RS*CNT<E>::RealStrideFactor>          
                                            TImag;
    typedef SymMat<M,EComplex,RS>           TComplex;
    typedef T                               THerm;   // These two are opposite of what you might expect
    typedef SymMat<M,EHerm,RS>              TPosTrans;
    typedef E                               TElement;
    typedef Vec<M,E,RS>                     TDiag;   // storage type for the diagonal elements
    typedef Vec<(M*(M-1))/2,E,RS>           TLower;  // storage type for the below-diag elements
    typedef Vec<(M*(M-1))/2,EHerm,RS>       TUpper;  // cast TLower to this for upper elements
    typedef Vec<(M*(M+1))/2,E,RS>           TAsVec;  // the whole SymMat as a single Vec

    // These are the results of calculations, so are returned in new, packed
    // memory. Be sure to refer to element types here which are also packed.
    typedef Row<M,E,1>                  TRow; // packed since we have to copy
    typedef Vec<M,E,1>                  TCol;
    typedef SymMat<M,EAbs,1>            TAbs;
    typedef SymMat<M,EStandard,1>       TStandard;
    typedef SymMat<M,EInvert,1>         TInvert;
    typedef SymMat<M,ENormalize,1>      TNormalize;
    typedef SymMat<M,ESqHermT,1>        TSqHermT; // ~Mat*Mat
    typedef SymMat<M,ESqTHerm,1>        TSqTHerm; // Mat*~Mat
    typedef SymMat<M,E,1>               TPacked;  // no extra row spacing for new data
    
    typedef EScalar                     Scalar;
    typedef ENumber                     Number;
    typedef EStdNumber                  StdNumber;
    typedef EPrecision                  Precision;
    typedef EScalarSq                   ScalarSq;

    int size() const { return (M*(M+1))/2; }
    int nrow() const { return M; }
    int ncol() const { return M; }

    // Scalar norm square is sum( squares of all scalars ). The off-diagonals
    // come up twice.
    ScalarSq scalarNormSqr() const { 
        return getDiag().scalarNormSqr() + 2.*getLower().scalarNormSqr();
    }

    // abs() is elementwise absolute value; that is, the return value has the same
    // dimension as this Mat but with each element replaced by whatever it thinks
    // its absolute value is.
    TAbs abs() const { 
        return TAbs(getAsVec().abs());
    }

    TStandard standardize() const {
        return TStandard(getAsVec().standardize());
    }

    StdNumber trace() const {return getDiag().sum();}

    // This gives the resulting SymMat type when (m[i] op P) is applied to each element of m.
    // It is a SymMat of dimension M, spacing 1, and element types which are the regular
    // composite result of E op P. Typically P is a scalar type but it doesn't have to be.
    template <class P> struct EltResult { 
        typedef SymMat<M, typename CNT<E>::template Result<P>::Mul, 1> Mul;
        typedef SymMat<M, typename CNT<E>::template Result<P>::Dvd, 1> Dvd;
        typedef SymMat<M, typename CNT<E>::template Result<P>::Add, 1> Add;
        typedef SymMat<M, typename CNT<E>::template Result<P>::Sub, 1> Sub;
    };

    // This is the composite result for m op P where P is some kind of appropriately shaped
    // non-scalar type.
    template <class P> struct Result { 
        typedef MulCNTs<M,M,ArgDepth,SymMat,ColSpacing,RowSpacing,
            CNT<P>::NRows, CNT<P>::NCols, CNT<P>::ArgDepth,
            P, CNT<P>::ColSpacing, CNT<P>::RowSpacing> MulOp;
        typedef typename MulOp::Type Mul;

        typedef MulCNTsNonConforming<M,M,ArgDepth,SymMat,ColSpacing,RowSpacing,
            CNT<P>::NRows, CNT<P>::NCols, CNT<P>::ArgDepth,
            P, CNT<P>::ColSpacing, CNT<P>::RowSpacing> MulOpNonConforming;
        typedef typename MulOpNonConforming::Type MulNon;


        typedef DvdCNTs<M,M,ArgDepth,SymMat,ColSpacing,RowSpacing,
            CNT<P>::NRows, CNT<P>::NCols, CNT<P>::ArgDepth,
            P, CNT<P>::ColSpacing, CNT<P>::RowSpacing> DvdOp;
        typedef typename DvdOp::Type Dvd;

        typedef AddCNTs<M,M,ArgDepth,SymMat,ColSpacing,RowSpacing,
            CNT<P>::NRows, CNT<P>::NCols, CNT<P>::ArgDepth,
            P, CNT<P>::ColSpacing, CNT<P>::RowSpacing> AddOp;
        typedef typename AddOp::Type Add;

        typedef SubCNTs<M,M,ArgDepth,SymMat,ColSpacing,RowSpacing,
            CNT<P>::NRows, CNT<P>::NCols, CNT<P>::ArgDepth,
            P, CNT<P>::ColSpacing, CNT<P>::RowSpacing> SubOp;
        typedef typename SubOp::Type Sub;
    };

    // Shape-preserving element substitution (always packed)
    template <class P> struct Substitute {
        typedef SymMat<M,P> Type;
    };

    // Default construction initializes to NaN when debugging but
    // is left uninitialized otherwise.
	SymMat(){ 
    #ifndef NDEBUG
        setToNaN();
    #endif
    }

    SymMat(const SymMat& src) {
        updAsVec() = src.getAsVec();
    }

    SymMat& operator=(const SymMat& src) {    // no harm if src and 'this' are the same
        updAsVec() = src.getAsVec();
        return *this;
    }

    // Allow an explicit conversion from square Mat of right size, looking only at lower
    // elements and real part of diagonal elements.
    template <class EE, int CSS, int RSS>
    explicit SymMat(const Mat<M,M,EE,CSS,RSS>& m) {
        updDiag() = m.diag().real();
        for (int j=0; j<M; ++j)
            for (int i=j+1; i<M; ++i)
                updEltLower(i,j) = m(i,j);
    }

    // We want an implicit conversion from a SymMat of the same length
    // and element type but with different spacings.
    template <int RSS> SymMat(const SymMat<M,E,RSS>& src) 
      { updAsVec() = src.getAsVec(); }

    // We want an implicit conversion from a SymMat of the same length
    // and *negated* element type, possibly with different spacings.
    template <int RSS> SymMat(const SymMat<M,ENeg,RSS>& src)
      { updAsVec() = src.getAsVec(); }

    // Construct a SymMat from a SymMat of the same dimensions, with any
    // spacings. Works as long as the element types are assignment compatible.
    template <class EE, int RSS> explicit SymMat(const SymMat<M,EE,RSS>& src)
      { updAsVec() = src.getAsVec(); }

    // Construction using an element repeats that element on the diagonal
    // but sets the rest of the matrix to zero.
    // TODO: diag should just use real part
    explicit SymMat(const E& e)
      { updDiag() = e; updLower() = E(0); }


    // Construction from a pointer to anything assumes we're pointing
    // at a packed element list of the right length, providing the
    // lower triangle in row order, so a b c d e f means
    //      a
    //      b c
    //      d e f
    //      g h i j
    // This has to be mapped to our diagonal/lower layout, which in
    // the above example will be:
    //      [a c f j][b d g e h i]
    //
    // In the input layout, the i'th row begins at element i(i+1)/2,
    // so diagonals are at i(i+1)/2 + i, while lower
    // elements (i,j; i>j) are at i(i+1)/2 + j.
    template <class EE> explicit SymMat(const EE* p) {
        assert(p);
        for (int i=0; i<M; ++i) {
            const int rowStart = (i*(i+1))/2;
            updDiag()[i] = p[rowStart + i];
            for (int j=0; j<i; ++j)
                updEltLower(i,j) = p[rowStart + j];
        }
    }

    // This is the same thing except as an assignment to pointer rather
    // than a constructor from pointer.
    template <class EE> SymMat& operator=(const EE* p) {
        assert(p);
        for (int i=0; i<M; ++i) {
            const int rowStart = (i*(i+1))/2;
            updDiag()[i] = p[rowStart + i];
            for (int j=0; j<i; ++j)
                updEltLower(i,j) = p[rowStart + j];
        }
        return *this;
    }

    // Assignment works similarly to copy -- if the lengths match,
    // go element-by-element. Otherwise, zero and then copy to each 
    // diagonal element.
    template <class EE, int RSS> SymMat& operator=(const SymMat<M,EE,RSS>& mm) {
        updAsVec() = mm.getAsVec();
        return *this;
    }


    // Conforming assignment ops
    template <class EE, int RSS> SymMat& 
    operator+=(const SymMat<M,EE,RSS>& mm) {
        updAsVec() += mm.getAsVec();
        return *this;
    }
    template <class EE, int RSS> SymMat&
    operator+=(const SymMat<M,negator<EE>,RSS>& mm) {
        updAsVec() -= -mm.getAsVec();   // negation of negator is free
        return *this;
    }

    template <class EE, int RSS> SymMat&
    operator-=(const SymMat<M,EE,RSS>& mm) {
        updAsVec() -= mm.getAsVec();
        return *this;
    }
    template <class EE, int RSS> SymMat&
    operator-=(const SymMat<M,negator<EE>,RSS>& mm) {
        updAsVec() += -mm.getAsVec();   // negation of negator is free
        return *this;
    }

    // In place matrix multiply can only be done when the RHS matrix is the same
    // size and the elements are also *= compatible.
    template <class EE, int RSS> SymMat&
    operator*=(const SymMat<M,EE,RSS>& mm) {
        assert(false); // TODO
        return *this;
    }

    // Conforming binary ops with 'this' on left, producing new packed result.
    // Cases: sy=sy+-sy, m=sy+-m, m=sy*sy, m=sy*m, v=sy*v

    // sy= this + sy
    template <class E2, int RS2> 
    typename Result<SymMat<M,E2,RS2> >::Add
    conformingAdd(const SymMat<M,E2,RS2>& r) const {
        return typename Result<SymMat<M,E2,RS2> >::Add
            (getAsVec().conformingAdd(r.getAsVec()));
    }
    // m= this - m
    template <class E2, int RS2> 
    typename Result<SymMat<M,E2,RS2> >::Sub
    conformingSubtract(const SymMat<M,E2,RS2>& r) const {
        return typename Result<SymMat<M,E2,RS2> >::Sub
            (getAsVec().conformingSubtract(r.getAsVec()));
    }

    // TODO: need the rest of the SymMat operators
    
    // Must be i >= j.
    const E& operator()(int i,int j) const 
      { return i==j ? getDiag()[i] : getEltLower(i,j); }
    E& operator()(int i,int j)
      { return i==j ? updDiag()[i] : updEltLower(i,j); }

    // This is the scalar Frobenius norm.
    ScalarSq normSqr() const { return scalarNormSqr(); }
    ScalarSq norm()    const { return std::sqrt(scalarNormSqr()); }

    // There is no conventional meaning for normalize() applied to a matrix. We
    // choose to define it as follows:
    // If the elements of this SymMat are scalars, the result is what you get by
    // dividing each element by the Frobenius norm() calculated above. If the elements are
    // *not* scalars, then the elements are *separately* normalized.
    //
    // Normalize returns a matrix of the same dimension but in new, packed storage
    // and with a return type that does not include negator<> even if the original
    // SymMat<> does, because we can eliminate the negation here almost for free.
    // But we can't standardize (change conjugate to complex) for free, so we'll retain
    // conjugates if there are any.
    TNormalize normalize() const {
        if (CNT<E>::IsScalar) {
            return castAwayNegatorIfAny() / (SignInterpretation*norm());
        } else {
            TNormalize elementwiseNormalized;
            // punt to the equivalent Vec to get the elements normalized
            elementwiseNormalized.updAsVec() = getAsVec().normalize();
            return elementwiseNormalized;
        }
    }

    TInvert invert() const {assert(false); return TInvert();} // TODO default inversion

    const SymMat& operator+() const { return *this; }
    const TNeg&   operator-() const { return negate(); }
    TNeg&         operator-()       { return updNegate(); }
    const THerm&  operator~() const { return transpose(); }
    THerm&        operator~()       { return updTranspose(); }

    const TNeg&  negate() const { return *reinterpret_cast<const TNeg*>(this); }
    TNeg&        updNegate()    { return *reinterpret_cast<TNeg*>(this); }

    const THerm& transpose()    const { return *reinterpret_cast<const THerm*>(this); }
    THerm&       updTranspose()       { return *reinterpret_cast<THerm*>(this); }

    const TPosTrans& positionalTranspose() const
        { return *reinterpret_cast<const TPosTrans*>(this); }
    TPosTrans&       updPositionalTranspose()
        { return *reinterpret_cast<TPosTrans*>(this); }

    const TReal& real() const { return *reinterpret_cast<const TReal*>(this); }
    TReal&       real()       { return *reinterpret_cast<      TReal*>(this); }

    // Had to contort these routines to get them through VC++ 7.net
    const TImag& imag()    const { 
        const int offs = ImagOffset;
        const EImag* p = reinterpret_cast<const EImag*>(this);
        return *reinterpret_cast<const TImag*>(p+offs);
    }
    TImag& imag() { 
        const int offs = ImagOffset;
        EImag* p = reinterpret_cast<EImag*>(this);
        return *reinterpret_cast<TImag*>(p+offs);
    }

    const TWithoutNegator& castAwayNegatorIfAny() const {return *reinterpret_cast<const TWithoutNegator*>(this);}
    TWithoutNegator&       updCastAwayNegatorIfAny()    {return *reinterpret_cast<TWithoutNegator*>(this);}

    // These are elementwise binary operators, (this op ee) by default but (ee op this) if
    // 'FromLeft' appears in the name. The result is a packed SymMat<M> but the element type
    // may change. These are mostly used to implement global operators.
    // We call these "scalar" operators but actually the "scalar" can be a composite type.

    //TODO: consider converting 'e' to Standard Numbers as precalculation and changing
    // return type appropriately.
    template <class EE> SymMat<M, typename CNT<E>::template Result<EE>::Mul>
    scalarMultiply(const EE& e) const {
        SymMat<M, typename CNT<E>::template Result<EE>::Mul> result;
        result.updAsVec() = getAsVec().scalarMultiply(e);
        return result;
    }
    template <class EE> SymMat<M, typename CNT<EE>::template Result<E>::Mul>
    scalarMultiplyFromLeft(const EE& e) const {
        SymMat<M, typename CNT<EE>::template Result<E>::Mul> result;
        result.updAsVec() = getAsVec().scalarMultiplyFromLeft(e);
        return result;
    }

    // TODO: should precalculate and store 1/e, while converting to Standard Numbers. Note
    // that return type should change appropriately.
    template <class EE> SymMat<M, typename CNT<E>::template Result<EE>::Dvd>
    scalarDivide(const EE& e) const {
        SymMat<M, typename CNT<E>::template Result<EE>::Dvd> result;
        result.updAsVec() = getAsVec().scalarDivide(e);
        return result;
    }
    template <class EE> SymMat<M, typename CNT<EE>::template Result<E>::Dvd>
    scalarDivideFromLeft(const EE& e) const {
        SymMat<M, typename CNT<EE>::template Result<E>::Dvd> result;
        result.updAsVec() = getAsVec().scalarDivideFromLeft(e);
        return result;
    }

    // Additive ops involving a scalar update only the diagonal
    template <class EE> SymMat<M, typename CNT<E>::template Result<EE>::Add>
    scalarAdd(const EE& e) const {
        SymMat<M, typename CNT<E>::template Result<EE>::Add> result(*this);
        result.updDiag() += e;
        return result;
    }
    // Add is commutative, so no 'FromLeft'.

    template <class EE> SymMat<M, typename CNT<E>::template Result<EE>::Sub>
    scalarSubtract(const EE& e) const {
        SymMat<M, typename CNT<E>::template Result<EE>::Sub> result(*this);
        result.diag() -= e;
        return result;
    }
    // This is s-m; negate m and add s to diagonal
    // TODO: Should do something clever with negation here.
    template <class EE> SymMat<M, typename CNT<EE>::template Result<E>::Sub>
    scalarSubtractFromLeft(const EE& e) const {
        SymMat<M, typename CNT<EE>::template Result<E>::Sub> result(-(*this));
        result.diag() += e;
        return result;
    }

    // Generic assignments for any element type not listed explicitly, including scalars.
    // These are done repeatedly for each element and only work if the operation can
    // be performed leaving the original element type.
    template <class EE> SymMat& operator =(const EE& e) {return scalarEq(e);}
    template <class EE> SymMat& operator+=(const EE& e) {return scalarPlusEq(e);}
    template <class EE> SymMat& operator-=(const EE& e) {return scalarMinusEq(e);}
    template <class EE> SymMat& operator*=(const EE& e) {return scalarTimesEq(e);}
    template <class EE> SymMat& operator/=(const EE& e) {return scalarDivideEq(e);}

    // Generalized scalar assignment & computed assignment methods. These will work
    // for any assignment-compatible element, not just scalars.
    template <class EE> SymMat& scalarEq(const EE& ee)
      { updDiag() = ee; updLower() = E(0); return *this; }
    template <class EE> SymMat& scalarPlusEq(const EE& ee)
      { updDiag() += ee; return *this; }
    template <class EE> SymMat& scalarMinusEq(const EE& ee)
      { updDiag() -= ee; return *this; }

    // this is m = s-m; negate off diagonal, do d=s-d for each diagonal element d
    template <class EE> SymMat& scalarMinusEqFromLeft(const EE& ee)
      { updLower() *= E(-1); updDiag().scalarMinusEqFromLeft(ee); return *this; }

    template <class EE> SymMat& scalarTimesEq(const EE& ee)
      { updAsVec() *= ee; return *this; }
    template <class EE> SymMat& scalarTimesEqFromLeft(const EE& ee)
      { updAsVec().scalarTimesEqFromLeft(ee); return *this; }
    template <class EE> SymMat& scalarDivideEq(const EE& ee)
      { updAsVec() /= ee; return *this; }
    template <class EE> SymMat& scalarDivideEqFromLeft(const EE& ee)
      { updAsVec().scalarDivideEqFromLeft(ee); return *this; } 

    void setToNaN() { updAsVec().setToNaN(); }

    // These assume we are given a pointer to d[0] of a SymMat<M,E,RS> like this one.
    static const SymMat& getAs(const ELT* p)  {return *reinterpret_cast<const SymMat*>(p);}
    static SymMat&       updAs(ELT* p)        {return *reinterpret_cast<SymMat*>(p);}

    // Note packed spacing
    static TPacked getNaN() {
        return TPacked(CNT<typename TPacked::TDiag>::getNaN(),
                       CNT<typename TPacked::TLower>::getNaN());
    }

    const TDiag&  getDiag()  const {return TDiag::getAs(d);}
    TDiag&        updDiag()        {return TDiag::updAs(d);}

    // Conventional synonym
    const TDiag& diag() const {return getDiag();}
    TDiag&       diag()       {return updDiag();}

    const TLower& getLower() const {return TLower::getAs(&d[M*RS]);}
    TLower&       updLower()       {return TLower::updAs(&d[M*RS]);}

    const TUpper& getUpper() const {return reinterpret_cast<const TUpper&>(getLower());}
    TUpper&       updUpper()       {return reinterpret_cast<      TUpper&>(updLower());}

    const TAsVec& getAsVec() const {return TAsVec::getAs(d);}
    TAsVec&       updAsVec()       {return TAsVec::updAs(d);}

    // must be i > j
    const E& getEltLower(int i, int j) const {return getLower()[lowerIx(i,j)];}
    E&       updEltLower(int i, int j)       {return updLower()[lowerIx(i,j)];}

    // must be i < j
    const EHerm& getEltUpper(int i, int j) const {return getUpper()[lowerIx(j,i)];}
    EHerm&       updEltUpper(int i, int j)       {return updUpper()[lowerIx(j,i)];}
    
    TRow sum() const {
        TRow temp(~getDiag());
        for (int i = 1; i < M; ++i)
            for (int j = 0; j < i; ++j) {
                E value = getEltLower(i, j);;
                temp[i] += value;
                temp[j] += value;
            }
        return temp;
    }

private:
    E d[NActualElements];

    SymMat(const TDiag& di, const TLower& low) {
        updDiag() = di; updLower() = low;
    }

    explicit SymMat(const TAsVec& v) {
        TAsVec::updAs(d) = v;
    }
    
    // Convert i,j with i>j to an index in "lower". This also gives the
    // right index for "upper" with i and j reversed.
    static int lowerIx(int i, int j) {
        assert(0 <= j && j < i && i < M);
        return (i-j-1) + j*(M-1) - (j*(j-1))/2;
    }

    template <int MM, class EE, int RSS> friend class SymMat;
};

////////////////////////////////////////////////////////
// Global operators involving two symmetric matrices. //
//   sy=sy+sy, sy=sy-sy, m=sy*sy, sy==sy, sy!=sy      //
////////////////////////////////////////////////////////

// m3 = m1 + m2 where all m's have the same dimension M. 
template <int M, class E1, int S1, class E2, int S2> inline
typename SymMat<M,E1,S1>::template Result< SymMat<M,E2,S2> >::Add
operator+(const SymMat<M,E1,S1>& l, const SymMat<M,E2,S2>& r) {
    return SymMat<M,E1,S1>::template Result< SymMat<M,E2,S2> >
        ::AddOp::perform(l,r);
}

// m3 = m1 - m2 where all m's have the same dimension M. 
template <int M, class E1, int S1, class E2, int S2> inline
typename SymMat<M,E1,S1>::template Result< SymMat<M,E2,S2> >::Sub
operator-(const SymMat<M,E1,S1>& l, const SymMat<M,E2,S2>& r) {
    return SymMat<M,E1,S1>::template Result< SymMat<M,E2,S2> >
        ::SubOp::perform(l,r);
}

// m3 = m1 * m2 where all m's have the same dimension M. 
// The result will not be symmetric.
template <int M, class E1, int S1, class E2, int S2> inline
typename SymMat<M,E1,S1>::template Result< SymMat<M,E2,S2> >::Mul
operator-(const SymMat<M,E1,S1>& l, const SymMat<M,E2,S2>& r) {
    return SymMat<M,E1,S1>::template Result< SymMat<M,E2,S2> >
        ::MulOp::perform(l,r);
}

// bool = m1 == m2, m1 and m2 have the same dimension M
template <int M, class E1, int S1, class E2, int S2> inline bool
operator==(const SymMat<M,E1,S1>& l, const SymMat<M,E2,S2>& r) {
    return l.getAsVec() == r.getAsVec();
}

// bool = m1 == m2, m1 and m2 have the same dimension M
template <int M, class E1, int S1, class E2, int S2> inline bool
operator!=(const SymMat<M,E1,S1>& l, const SymMat<M,E2,S2>& r) {return !(l==r);} 

///////////////////////////////////////////////////////
// Global operators involving a SymMat and a scalar. //
///////////////////////////////////////////////////////

// SCALAR MULTIPLY

// m = m*real, real*m 
template <int M, class E, int S> inline
typename SymMat<M,E,S>::template Result<float>::Mul
operator*(const SymMat<M,E,S>& l, const float& r)
  { return SymMat<M,E,S>::template Result<float>::MulOp::perform(l,r); }
template <int M, class E, int S> inline
typename SymMat<M,E,S>::template Result<float>::Mul
operator*(const float& l, const SymMat<M,E,S>& r) {return r*l;}

template <int M, class E, int S> inline
typename SymMat<M,E,S>::template Result<double>::Mul
operator*(const SymMat<M,E,S>& l, const double& r)
  { return SymMat<M,E,S>::template Result<double>::MulOp::perform(l,r); }
template <int M, class E, int S> inline
typename SymMat<M,E,S>::template Result<double>::Mul
operator*(const double& l, const SymMat<M,E,S>& r) {return r*l;}

template <int M, class E, int S> inline
typename SymMat<M,E,S>::template Result<long double>::Mul
operator*(const SymMat<M,E,S>& l, const long double& r)
  { return SymMat<M,E,S>::template Result<long double>::MulOp::perform(l,r); }
template <int M, class E, int S> inline
typename SymMat<M,E,S>::template Result<long double>::Mul
operator*(const long double& l, const SymMat<M,E,S>& r) {return r*l;}

// m = m*int, int*m -- just convert int to m's precision float
template <int M, class E, int S> inline
typename SymMat<M,E,S>::template Result<typename CNT<E>::Precision>::Mul
operator*(const SymMat<M,E,S>& l, int r) {return l * (typename CNT<E>::Precision)r;}
template <int M, class E, int S> inline
typename SymMat<M,E,S>::template Result<typename CNT<E>::Precision>::Mul
operator*(int l, const SymMat<M,E,S>& r) {return r * (typename CNT<E>::Precision)l;}

// Complex, conjugate, and negator are all easy to templatize.

// m = m*complex, complex*m
template <int M, class E, int S, class R> inline
typename SymMat<M,E,S>::template Result<std::complex<R> >::Mul
operator*(const SymMat<M,E,S>& l, const std::complex<R>& r)
  { return SymMat<M,E,S>::template Result<std::complex<R> >::MulOp::perform(l,r); }
template <int M, class E, int S, class R> inline
typename SymMat<M,E,S>::template Result<std::complex<R> >::Mul
operator*(const std::complex<R>& l, const SymMat<M,E,S>& r) {return r*l;}

// m = m*conjugate, conjugate*m (convert conjugate->complex)
template <int M, class E, int S, class R> inline
typename SymMat<M,E,S>::template Result<std::complex<R> >::Mul
operator*(const SymMat<M,E,S>& l, const conjugate<R>& r) {return l*(std::complex<R>)r;}
template <int M, class E, int S, class R> inline
typename SymMat<M,E,S>::template Result<std::complex<R> >::Mul
operator*(const conjugate<R>& l, const SymMat<M,E,S>& r) {return r*(std::complex<R>)l;}

// m = m*negator, negator*m: convert negator to standard number
template <int M, class E, int S, class R> inline
typename SymMat<M,E,S>::template Result<typename negator<R>::StdNumber>::Mul
operator*(const SymMat<M,E,S>& l, const negator<R>& r) {return l * (typename negator<R>::StdNumber)(R)r;}
template <int M, class E, int S, class R> inline
typename SymMat<M,E,S>::template Result<typename negator<R>::StdNumber>::Mul
operator*(const negator<R>& l, const SymMat<M,E,S>& r) {return r * (typename negator<R>::StdNumber)(R)l;}


// SCALAR DIVIDE. This is a scalar operation when the scalar is on the right,
// but when it is on the left it means scalar * pseudoInverse(mat), which is 
// a similar symmetric matrix.

// m = m/real, real/m 
template <int M, class E, int S> inline
typename SymMat<M,E,S>::template Result<float>::Dvd
operator/(const SymMat<M,E,S>& l, const float& r)
  { return SymMat<M,E,S>::template Result<float>::DvdOp::perform(l,r); }
template <int M, class E, int S> inline
typename CNT<float>::template Result<SymMat<M,E,S> >::Dvd
operator/(const float& l, const SymMat<M,E,S>& r)
  { return CNT<float>::template Result<SymMat<M,E,S> >::DvdOp::perform(l,r); }

template <int M, class E, int S> inline
typename SymMat<M,E,S>::template Result<double>::Dvd
operator/(const SymMat<M,E,S>& l, const double& r)
  { return SymMat<M,E,S>::template Result<double>::DvdOp::perform(l,r); }
template <int M, class E, int S> inline
typename CNT<double>::template Result<SymMat<M,E,S> >::Dvd
operator/(const double& l, const SymMat<M,E,S>& r)
  { return CNT<double>::template Result<SymMat<M,E,S> >::DvdOp::perform(l,r); }

template <int M, class E, int S> inline
typename SymMat<M,E,S>::template Result<long double>::Dvd
operator/(const SymMat<M,E,S>& l, const long double& r)
  { return SymMat<M,E,S>::template Result<long double>::DvdOp::perform(l,r); }
template <int M, class E, int S> inline
typename CNT<long double>::template Result<SymMat<M,E,S> >::Dvd
operator/(const long double& l, const SymMat<M,E,S>& r)
  { return CNT<long double>::template Result<SymMat<M,E,S> >::DvdOp::perform(l,r); }

// m = m/int, int/m -- just convert int to m's precision float
template <int M, class E, int S> inline
typename SymMat<M,E,S>::template Result<typename CNT<E>::Precision>::Dvd
operator/(const SymMat<M,E,S>& l, int r) {return l / (typename CNT<E>::Precision)r;}
template <int M, class E, int S> inline
typename CNT<typename CNT<E>::Precision>::template Result<SymMat<M,E,S> >::Dvd
operator/(int l, const SymMat<M,E,S>& r) {return (typename CNT<E>::Precision)l / r;}


// Complex, conjugate, and negator are all easy to templatize.

// m = m/complex, complex/m
template <int M, class E, int S, class R> inline
typename SymMat<M,E,S>::template Result<std::complex<R> >::Dvd
operator/(const SymMat<M,E,S>& l, const std::complex<R>& r)
  { return SymMat<M,E,S>::template Result<std::complex<R> >::DvdOp::perform(l,r); }
template <int M, class E, int S, class R> inline
typename CNT<std::complex<R> >::template Result<SymMat<M,E,S> >::Dvd
operator/(const std::complex<R>& l, const SymMat<M,E,S>& r)
  { return CNT<std::complex<R> >::template Result<SymMat<M,E,S> >::DvdOp::perform(l,r); }

// m = m/conjugate, conjugate/m (convert conjugate->complex)
template <int M, class E, int S, class R> inline
typename SymMat<M,E,S>::template Result<std::complex<R> >::Dvd
operator/(const SymMat<M,E,S>& l, const conjugate<R>& r) {return l/(std::complex<R>)r;}
template <int M, class E, int S, class R> inline
typename CNT<std::complex<R> >::template Result<SymMat<M,E,S> >::Dvd
operator/(const conjugate<R>& l, const SymMat<M,E,S>& r) {return (std::complex<R>)l/r;}

// m = m/negator, negator/m: convert negator to number
template <int M, class E, int S, class R> inline
typename SymMat<M,E,S>::template Result<typename negator<R>::StdNumber>::Dvd
operator/(const SymMat<M,E,S>& l, const negator<R>& r) {return l/(typename negator<R>::StdNumber)(R)r;}
template <int M, class E, int S, class R> inline
typename CNT<R>::template Result<SymMat<M,E,S> >::Dvd
operator/(const negator<R>& l, const SymMat<M,E,S>& r) {return (typename negator<R>::StdNumber)(R)l/r;}


// Add and subtract are odd as scalar ops. They behave as though the
// scalar stands for a conforming matrix whose diagonal elements are that,
// scalar and then a normal matrix add or subtract is done.

// SCALAR ADD

// m = m+real, real+m 
template <int M, class E, int S> inline
typename SymMat<M,E,S>::template Result<float>::Add
operator+(const SymMat<M,E,S>& l, const float& r)
  { return SymMat<M,E,S>::template Result<float>::AddOp::perform(l,r); }
template <int M, class E, int S> inline
typename SymMat<M,E,S>::template Result<float>::Add
operator+(const float& l, const SymMat<M,E,S>& r) {return r+l;}

template <int M, class E, int S> inline
typename SymMat<M,E,S>::template Result<double>::Add
operator+(const SymMat<M,E,S>& l, const double& r)
  { return SymMat<M,E,S>::template Result<double>::AddOp::perform(l,r); }
template <int M, class E, int S> inline
typename SymMat<M,E,S>::template Result<double>::Add
operator+(const double& l, const SymMat<M,E,S>& r) {return r+l;}

template <int M, class E, int S> inline
typename SymMat<M,E,S>::template Result<long double>::Add
operator+(const SymMat<M,E,S>& l, const long double& r)
  { return SymMat<M,E,S>::template Result<long double>::AddOp::perform(l,r); }
template <int M, class E, int S> inline
typename SymMat<M,E,S>::template Result<long double>::Add
operator+(const long double& l, const SymMat<M,E,S>& r) {return r+l;}

// m = m+int, int+m -- just convert int to m's precision float
template <int M, class E, int S> inline
typename SymMat<M,E,S>::template Result<typename CNT<E>::Precision>::Add
operator+(const SymMat<M,E,S>& l, int r) {return l + (typename CNT<E>::Precision)r;}
template <int M, class E, int S> inline
typename SymMat<M,E,S>::template Result<typename CNT<E>::Precision>::Add
operator+(int l, const SymMat<M,E,S>& r) {return r + (typename CNT<E>::Precision)l;}

// Complex, conjugate, and negator are all easy to templatize.

// m = m+complex, complex+m
template <int M, class E, int S, class R> inline
typename SymMat<M,E,S>::template Result<std::complex<R> >::Add
operator+(const SymMat<M,E,S>& l, const std::complex<R>& r)
  { return SymMat<M,E,S>::template Result<std::complex<R> >::AddOp::perform(l,r); }
template <int M, class E, int S, class R> inline
typename SymMat<M,E,S>::template Result<std::complex<R> >::Add
operator+(const std::complex<R>& l, const SymMat<M,E,S>& r) {return r+l;}

// m = m+conjugate, conjugate+m (convert conjugate->complex)
template <int M, class E, int S, class R> inline
typename SymMat<M,E,S>::template Result<std::complex<R> >::Add
operator+(const SymMat<M,E,S>& l, const conjugate<R>& r) {return l+(std::complex<R>)r;}
template <int M, class E, int S, class R> inline
typename SymMat<M,E,S>::template Result<std::complex<R> >::Add
operator+(const conjugate<R>& l, const SymMat<M,E,S>& r) {return r+(std::complex<R>)l;}

// m = m+negator, negator+m: convert negator to standard number
template <int M, class E, int S, class R> inline
typename SymMat<M,E,S>::template Result<typename negator<R>::StdNumber>::Add
operator+(const SymMat<M,E,S>& l, const negator<R>& r) {return l + (typename negator<R>::StdNumber)(R)r;}
template <int M, class E, int S, class R> inline
typename SymMat<M,E,S>::template Result<typename negator<R>::StdNumber>::Add
operator+(const negator<R>& l, const SymMat<M,E,S>& r) {return r + (typename negator<R>::StdNumber)(R)l;}

// SCALAR SUBTRACT -- careful, not commutative.

// m = m-real, real-m 
template <int M, class E, int S> inline
typename SymMat<M,E,S>::template Result<float>::Sub
operator-(const SymMat<M,E,S>& l, const float& r)
  { return SymMat<M,E,S>::template Result<float>::SubOp::perform(l,r); }
template <int M, class E, int S> inline
typename CNT<float>::template Result<SymMat<M,E,S> >::Sub
operator-(const float& l, const SymMat<M,E,S>& r)
  { return CNT<float>::template Result<SymMat<M,E,S> >::SubOp::perform(l,r); }

template <int M, class E, int S> inline
typename SymMat<M,E,S>::template Result<double>::Sub
operator-(const SymMat<M,E,S>& l, const double& r)
  { return SymMat<M,E,S>::template Result<double>::SubOp::perform(l,r); }
template <int M, class E, int S> inline
typename CNT<double>::template Result<SymMat<M,E,S> >::Sub
operator-(const double& l, const SymMat<M,E,S>& r)
  { return CNT<double>::template Result<SymMat<M,E,S> >::SubOp::perform(l,r); }

template <int M, class E, int S> inline
typename SymMat<M,E,S>::template Result<long double>::Sub
operator-(const SymMat<M,E,S>& l, const long double& r)
  { return SymMat<M,E,S>::template Result<long double>::SubOp::perform(l,r); }
template <int M, class E, int S> inline
typename CNT<long double>::template Result<SymMat<M,E,S> >::Sub
operator-(const long double& l, const SymMat<M,E,S>& r)
  { return CNT<long double>::template Result<SymMat<M,E,S> >::SubOp::perform(l,r); }

// m = m-int, int-m // just convert int to m's precision float
template <int M, class E, int S> inline
typename SymMat<M,E,S>::template Result<typename CNT<E>::Precision>::Sub
operator-(const SymMat<M,E,S>& l, int r) {return l - (typename CNT<E>::Precision)r;}
template <int M, class E, int S> inline
typename CNT<typename CNT<E>::Precision>::template Result<SymMat<M,E,S> >::Sub
operator-(int l, const SymMat<M,E,S>& r) {return (typename CNT<E>::Precision)l - r;}


// Complex, conjugate, and negator are all easy to templatize.

// m = m-complex, complex-m
template <int M, class E, int S, class R> inline
typename SymMat<M,E,S>::template Result<std::complex<R> >::Sub
operator-(const SymMat<M,E,S>& l, const std::complex<R>& r)
  { return SymMat<M,E,S>::template Result<std::complex<R> >::SubOp::perform(l,r); }
template <int M, class E, int S, class R> inline
typename CNT<std::complex<R> >::template Result<SymMat<M,E,S> >::Sub
operator-(const std::complex<R>& l, const SymMat<M,E,S>& r)
  { return CNT<std::complex<R> >::template Result<SymMat<M,E,S> >::SubOp::perform(l,r); }

// m = m-conjugate, conjugate-m (convert conjugate->complex)
template <int M, class E, int S, class R> inline
typename SymMat<M,E,S>::template Result<std::complex<R> >::Sub
operator-(const SymMat<M,E,S>& l, const conjugate<R>& r) {return l-(std::complex<R>)r;}
template <int M, class E, int S, class R> inline
typename CNT<std::complex<R> >::template Result<SymMat<M,E,S> >::Sub
operator-(const conjugate<R>& l, const SymMat<M,E,S>& r) {return (std::complex<R>)l-r;}

// m = m-negator, negator-m: convert negator to standard number
template <int M, class E, int S, class R> inline
typename SymMat<M,E,S>::template Result<typename negator<R>::StdNumber>::Sub
operator-(const SymMat<M,E,S>& l, const negator<R>& r) {return l-(typename negator<R>::StdNumber)(R)r;}
template <int M, class E, int S, class R> inline
typename CNT<R>::template Result<SymMat<M,E,S> >::Sub
operator-(const negator<R>& l, const SymMat<M,E,S>& r) {return (typename negator<R>::StdNumber)(R)l-r;}


// SymMat I/O
template <int M, class E, int RS, class CHAR, class TRAITS> inline
std::basic_ostream<CHAR,TRAITS>&
operator<<(std::basic_ostream<CHAR,TRAITS>& o, const SymMat<M,E,RS>& m) {
    for (int i=0;i<M;++i) {
        o << std::endl << "[";
        for (int j=0; j<=i; ++j)         
            o << (j>0?" ":"") << m(i,j);
        for (int j=i+1; j<M; ++j)
            o << " *";
        o << "]";
    }
    if (M) o << std::endl;
    return o; 
}

template <int M, class E, int RS, class CHAR, class TRAITS> inline
std::basic_istream<CHAR,TRAITS>&
operator>>(std::basic_istream<CHAR,TRAITS>& is, SymMat<M,E,RS>& m) {
    // TODO: not sure how to do input yet
    assert(false);
    return is;
}

} //namespace SimTK


#endif //SimTK_SIMMATRIX_SMALLMATRIX_SYMMAT_H_
