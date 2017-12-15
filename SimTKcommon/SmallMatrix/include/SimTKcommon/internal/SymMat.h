#ifndef SimTK_SIMMATRIX_SMALLMATRIX_SYMMAT_H_
#define SimTK_SIMMATRIX_SMALLMATRIX_SYMMAT_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-13 Stanford University and the Authors.        *
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

/**@file
 This file declares templatized class SymMat for small, fixed-size symmetric
 matrices. **/


#include "SimTKcommon/internal/common.h"

namespace SimTK {

/** @brief This is a small, fixed-size symmetric or Hermitian matrix designed
for no-overhead inline computation.

@ingroup MatVecUtilities

@tparam     M       The dimension of this square matrix.
@tparam     ELT     The element type. Must be a composite numerical type (CNT).
The default is ELT=Real.
@tparam     RS  The spacing from one element to the next in memory, as an
integer number of elements of type ELT. The default is RS=1.

This is logically a square MxM Hermitian matrix, but only the diagonal
and \e lower triangle are stored. Elements above the diagonal are
the Hermitan transpose of their mirrored elements below the diagonal.

The storage is packed by column, with an optional stride RS when going element-
to-element. We use an unconventional storage scheme here which
provides substantial conveniences and some performance advantage over the
conventional (LAPACK) scheme which stores the matrix in packed
columns. Instead, we store the diagonal first, as a Vec<M> and then store
the lower triangle separately in the conventional columnwise form. This
allows easy handling of the diagonal elements separately, which is very
common in a symmetric matrix. The diagonals are "special"; for example
they must be "self-Hermitian" meaning their imaginary parts must be zero.
We could just store real elements for the diagonals, but we don't for
a number of good reasons which I'll leave as an exercise for the reader.
However, we will ensure that their imaginary parts are zero, although
a determined user could screw this up.

Here is an example. Imagine a full 4x4 Hermitian matrix with complex
elements. We'll use a' to mean conjugate(a). Then the elements are
<pre>
    w  a' b' c'
    a  x  d' e'
    b  d  y  f'
    c  e  f  z
</pre>
Note that it must be the case that w=w', x=x', y=y', and z=z' so the
diagonal elements must be real. Here's how we store that as a SymMat:
<pre>
    w x y z a b c d e f
</pre>
These are stored in consecutive memory locations (possibly separated
by meaningless elements as indicated by the row spacing RS). Mere
recasting allows us to view this as a 4x4 matrix, or a 4-element
diagonal and 6-element "lower" vector, or as these vectors
apparently containing the Hermitian types (i.e. one may cast
the "lower" portion into what looks like the "upper" portion).

This is a Composite Numerical Type (CNT).
**/
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

    typedef typename CNT<E>::TSqrt              ESqrt;
    typedef typename CNT<E>::TAbs               EAbs;
    typedef typename CNT<E>::TStandard          EStandard;
    typedef typename CNT<E>::TInvert            EInvert;
    typedef typename CNT<E>::TNormalize         ENormalize;

    typedef typename CNT<E>::Scalar             EScalar;
    typedef typename CNT<E>::ULessScalar        EULessScalar;
    typedef typename CNT<E>::Number             ENumber;
    typedef typename CNT<E>::StdNumber          EStdNumber;
    typedef typename CNT<E>::Precision          EPrecision;
    typedef typename CNT<E>::ScalarNormSq       EScalarNormSq;

public:

    enum {
        NRows               = M,
        NCols               = M,
        NDiagElements       = M,
        NLowerElements      = (M*(M-1))/2,
        NPackedElements     = NDiagElements+NLowerElements,
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
        IsULessScalar       = 0,
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
    typedef SymMat<M,ESqrt,1>           TSqrt;
    typedef SymMat<M,EAbs,1>            TAbs;
    typedef SymMat<M,EStandard,1>       TStandard;
    typedef SymMat<M,EInvert,1>         TInvert;
    typedef SymMat<M,ENormalize,1>      TNormalize;
    typedef SymMat<M,ESqHermT,1>        TSqHermT; // ~Mat*Mat
    typedef SymMat<M,ESqTHerm,1>        TSqTHerm; // Mat*~Mat
    typedef SymMat<M,E,1>               TPacked;  // no extra row spacing for new data
    
    typedef EScalar                     Scalar;
    typedef EULessScalar                ULessScalar;
    typedef ENumber                     Number;
    typedef EStdNumber                  StdNumber;
    typedef EPrecision                  Precision;
    typedef EScalarNormSq               ScalarNormSq;

    static int size() { return (M*(M+1))/2; }
    static int nrow() { return M; }
    static int ncol() { return M; }

    // Scalar norm square is sum( squares of all scalars ). The off-diagonals
    // come up twice.
    ScalarNormSq scalarNormSqr() const { 
        return getDiag().scalarNormSqr() + 2.*getLower().scalarNormSqr();
    }

    // sqrt() is elementwise square root; that is, the return value has the same
    // dimension as this SymMat but with each element replaced by whatever it thinks
    // its square root is.
    TSqrt sqrt() const { 
        return TSqrt(getAsVec().sqrt());
    }

    // abs() is elementwise absolute value; that is, the return value has the same
    // dimension as this SymMat but with each element replaced by whatever it thinks
    // its absolute value is.
    TAbs abs() const { 
        return TAbs(getAsVec().abs());
    }

    TStandard standardize() const {
        return TStandard(getAsVec().standardize());
    }

    EStandard trace() const {return getDiag().sum();}

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

    /// Default construction initializes to NaN when debugging but
    /// is left uninitialized otherwise.
    SymMat(){ 
    #ifndef NDEBUG
        setToNaN();
    #endif
    }

    /// Copy constructor.
    SymMat(const SymMat& src) {
        updAsVec() = src.getAsVec();
    }

    /// Copy assignment; no harm if source and this are the same matrix.
    SymMat& operator=(const SymMat& src) {
        updAsVec() = src.getAsVec();
        return *this;
    }

    /// This is an \e explicit conversion from square Mat of right size, assuming
    /// that the source matrix is symmetric to within a reasonable numerical 
    /// tolerance. In Debug mode we'll test that assumption and throw an exception
    /// if it is wrong. In Release mode you're on your own. All the elements of
    /// the source Mat are used; off-diagonal elements (i,j) are averaged with
    /// their corresponding element (j,i); the imaginary part of the diagonal
    /// is set exactly to zero. If you don't want to spend the flops to average
    /// the off-diagonals, and you're sure the source is symmetric, use either
    /// setFromLower() or setFromUpper() which will just copy the elements.
    /// @see setFromLower(), setFromUpper()
    template <class EE, int CSS, int RSS>
    explicit SymMat(const Mat<M,M,EE,CSS,RSS>& m)
    {   setFromSymmetric(m); }

    /// Create a new SymMat of this type from a square Mat of the right
    /// size, looking only at lower elements and the real part of the
    /// diagonal.
    template <class EE, int CSS, int RSS>
    SymMat& setFromLower(const Mat<M,M,EE,CSS,RSS>& m) {
        this->updDiag() = m.diag().real();
        for (int j=0; j<M; ++j)
            for (int i=j+1; i<M; ++i)
                this->updEltLower(i,j) = m(i,j);
        return *this;
    }

    /// Create a new SymMat of this type from a square Mat of the right
    /// size, looking only at upper elements and the real part of the
    /// diagonal. Note that the SymMat's stored elements are still in
    /// its \e lower triangle; they are just initialized from the Mat's
    /// \e upper triangle. There is no transposing of elements here;
    /// we simply copy the upper elements of the Mat to the corresponding
    /// lower elements of the SymMat.
    template <class EE, int CSS, int RSS>
    SymMat& setFromUpper(const Mat<M,M,EE,CSS,RSS>& m) {
        this->updDiag() = m.diag().real();
        for (int j=0; j<M; ++j)
            for (int i=j+1; i<M; ++i)
                this->updEltLower(i,j) = m(j,i);
        return *this;
    }

    /// Create a new SymMat of this type from a square Mat of the right
    /// size, that is expected to be symmetric (hermitian) to within
    /// a tolerance. All elements are used; we average the upper and
    /// lower elements of the Mat to produce the corresponding element
    /// of the SymMat.
    template <class EE, int CSS, int RSS>
    SymMat& setFromSymmetric(const Mat<M,M,EE,CSS,RSS>& m) {
        SimTK_ERRCHK1(m.isNumericallySymmetric(), "SymMat::setFromSymmetric()",
            "The allegedly symmetric source matrix was not symmetric to within "
            "a tolerance of %g.", m.getDefaultTolerance());
        this->updDiag() = m.diag().real();
        for (int j=0; j<M; ++j)
            for (int i=j+1; i<M; ++i)
                this->updEltLower(i,j) = 
                    (m(i,j) + CNT<EE>::transpose(m(j,i)))/2;
        return *this;
    }

    /// This is an \e implicit conversion from a SymMat of the same length
    /// and element type but with different spacing.
    template <int RSS> SymMat(const SymMat<M,E,RSS>& src) 
      { updAsVec() = src.getAsVec(); }

    /// This is an \e implicit conversion from a SymMat of the same length
    /// and \e negated element type, possibly with different spacing.
    template <int RSS> SymMat(const SymMat<M,ENeg,RSS>& src)
      { updAsVec() = src.getAsVec(); }

    /// Construct a SymMat from a SymMat of the same dimensions, with any
    /// element type and spacing. Works as long as the element types are 
    /// assignment compatible.
    template <class EE, int RSS> explicit SymMat(const SymMat<M,EE,RSS>& src)
      { updAsVec() = src.getAsVec(); }

    // Construction using an element repeats that element on the diagonal
    // but sets the rest of the matrix to zero.
    explicit SymMat(const E& e) {
        updDiag() = CNT<E>::real(e); 
        for (int i=0; i < NLowerElements; ++i) updlowerE(i) = E(0); 
    }

    // Construction using a negated element is just like construction from
    // the element.
    explicit SymMat(const ENeg& e) {
        updDiag() = CNT<ENeg>::real(e); 
        for (int i=0; i < NLowerElements; ++i) updlowerE(i) = E(0); 
    }

    // Given an int, turn it into a suitable floating point number
    // and then feed that to the above single-element constructor.
    explicit SymMat(int i) 
      { new (this) SymMat(E(Precision(i))); }

    /// A bevy of constructors from individual exact-match elements IN ROW ORDER,
    /// giving the LOWER TRIANGLE, like this:
    /// <pre>
    ///     a
    ///     b c
    ///     d e f
    ///     g h i j
    /// </pre>
    /// Note that this will be mapped to our diagonal/lower layout, which in
    /// the above example would be:
    /// <pre>
    ///     [a c f j][b d g e h i]
    /// </pre>

    SymMat(const E& e0,
           const E& e1,const E& e2)
    {   assert(M==2); TDiag& d=updDiag(); TLower& l=updLower();
        d[0]=CNT<E>::real(e0); d[1]=CNT<E>::real(e2); 
        l[0]=e1; }

    SymMat(const E& e0,
           const E& e1,const E& e2,
           const E& e3,const E& e4, const E& e5)
    {   assert(M==3); TDiag& d=updDiag(); TLower& l=updLower();
        d[0]=CNT<E>::real(e0);d[1]=CNT<E>::real(e2);d[2]=CNT<E>::real(e5); 
        l[0]=e1;l[1]=e3;
        l[2]=e4; }

    SymMat(const E& e0,
           const E& e1,const E& e2,
           const E& e3,const E& e4, const E& e5,
           const E& e6,const E& e7, const E& e8, const E& e9)
    {   assert(M==4); TDiag& d=updDiag(); TLower& l=updLower();
        d[0]=CNT<E>::real(e0);d[1]=CNT<E>::real(e2);d[2]=CNT<E>::real(e5);d[3]=CNT<E>::real(e9); 
        l[0]=e1;l[1]=e3;l[2]=e6;
        l[3]=e4;l[4]=e7;
        l[5]=e8; }

    SymMat(const E& e0,
           const E& e1, const E& e2,
           const E& e3, const E& e4,  const E& e5,
           const E& e6, const E& e7,  const E& e8,  const E& e9,
           const E& e10,const E& e11, const E& e12, const E& e13, const E& e14)
    {   assert(M==5); TDiag& d=updDiag(); TLower& l=updLower();
        d[0]=CNT<E>::real(e0);d[1]=CNT<E>::real(e2);d[2]=CNT<E>::real(e5);d[3]=CNT<E>::real(e9);d[4]=CNT<E>::real(e14); 
        l[0]=e1;l[1]=e3;l[2]=e6;l[3]=e10;
        l[4]=e4;l[5]=e7;l[6]=e11;
        l[7]=e8;l[8]=e12;
        l[9]=e13; }

    SymMat(const E& e0,
           const E& e1, const E& e2,
           const E& e3, const E& e4,  const E& e5,
           const E& e6, const E& e7,  const E& e8,  const E& e9,
           const E& e10,const E& e11, const E& e12, const E& e13, const E& e14,
           const E& e15,const E& e16, const E& e17, const E& e18, const E& e19, const E& e20)
    {   assert(M==6); TDiag& d=updDiag(); TLower& l=updLower();
        d[0]=CNT<E>::real(e0);d[1]=CNT<E>::real(e2);d[2]=CNT<E>::real(e5);
        d[3]=CNT<E>::real(e9);d[4]=CNT<E>::real(e14);d[5]=CNT<E>::real(e20); 
        l[0] =e1; l[1] =e3; l[2] =e6;  l[3]=e10; l[4]=e15;
        l[5] =e4; l[6] =e7; l[7] =e11; l[8]=e16;
        l[9] =e8; l[10]=e12;l[11]=e17;
        l[12]=e13;l[13]=e18;
        l[14]=e19; }

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
            updDiag()[i] = CNT<EE>::real(p[rowStart + i]);
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
            updDiag()[i] = CNT<EE>::real(p[rowStart + i]);
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

    // symmetric * symmetric produces a full result
    // m= this * s
    // TODO: this is not a good implementation
    template <class E2, int RS2>
    typename Result<SymMat<M,E2,RS2> >::Mul
    conformingMultiply(const SymMat<M,E2,RS2>& s) const {
        typename Result<SymMat<M,E2,RS2> >::Mul result;
        for (int j=0;j<M;++j)
            for (int i=0;i<M;++i)
                result(i,j) = (*this)[i] * s(j);
        return result;
    }

    // sy= this .* sy
    template <class E2, int RS2> 
    typename EltResult<E2>::Mul
    elementwiseMultiply(const SymMat<M,E2,RS2>& r) const {
        return typename EltResult<E2>::Mul
            (getAsVec().elementwiseMultiply(r.getAsVec()));
    }

    // sy= this ./ sy
    template <class E2, int RS2> 
    typename EltResult<E2>::Dvd
    elementwiseDivide(const SymMat<M,E2,RS2>& r) const {
        return typename EltResult<E2>::Dvd
            (getAsVec().elementwiseDivide(r.getAsVec()));
    }

    // TODO: need the rest of the SymMat operators
    
    // Must be i >= j.
    const E& operator()(int i,int j) const 
      { return i==j ? getDiag()[i] : getEltLower(i,j); }
    E& operator()(int i,int j)
      { return i==j ? updDiag()[i] : updEltLower(i,j); }

    // These are slow for a symmetric matrix, requiring copying and
    // possibly floating point operations for conjugation.
    TRow operator[](int i) const {return row(i);}
    TCol operator()(int j) const {return col(j);}


    // This is the scalar Frobenius norm.
    ScalarNormSq normSqr() const { return scalarNormSqr(); }
    typename CNT<ScalarNormSq>::TSqrt 
        norm() const { return CNT<ScalarNormSq>::sqrt(scalarNormSqr()); }

    /// There is no conventional meaning for normalize() applied to a matrix. We
    /// choose to define it as follows:
    /// If the elements of this SymMat are scalars, the result is what you get by
    /// dividing each element by the Frobenius norm() calculated above. If the elements are
    /// *not* scalars, then the elements are *separately* normalized.
    ///
    /// Normalize returns a matrix of the same dimension but in new, packed storage
    /// and with a return type that does not include negator<> even if the original
    /// SymMat<> does, because we can eliminate the negation here almost for free.
    /// But we can't standardize (change conjugate to complex) for free, so we'll retain
    /// conjugates if there are any.
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

    void setToNaN()  { updAsVec().setToNaN();  }
    void setToZero() { updAsVec().setToZero(); }

    // These assume we are given a pointer to d[0] of a SymMat<M,E,RS> like this one.
    static const SymMat& getAs(const ELT* p)  {return *reinterpret_cast<const SymMat*>(p);}
    static SymMat&       updAs(ELT* p)        {return *reinterpret_cast<SymMat*>(p);}

    // Note packed spacing
    static TPacked getNaN() {
        return TPacked(CNT<typename TPacked::TDiag>::getNaN(),
                       CNT<typename TPacked::TLower>::getNaN());
    }

    /// Return true if any element of this SymMat contains a NaN anywhere.
    bool isNaN() const {return getAsVec().isNaN();}

    /// Return true if any element of this SymMat contains a +Inf
    /// or -Inf somewhere but no element contains a NaN anywhere.
    bool isInf() const {return getAsVec().isInf();}

    /// Return true if no element contains an Infinity or a NaN.
    bool isFinite() const {return getAsVec().isFinite();}

    /// For approximate comparisons, the default tolerance to use for a matrix is
    /// its shortest dimension times its elements' default tolerance.
    static double getDefaultTolerance() {return M*CNT<ELT>::getDefaultTolerance();}

    /// %Test whether this matrix is numerically equal to some other matrix with
    /// the same shape, using a specified tolerance.
    template <class E2, int RS2>
    bool isNumericallyEqual(const SymMat<M,E2,RS2>& m, double tol) const {
        return getAsVec().isNumericallyEqual(m.getAsVec(), tol);
    }

    /// %Test whether this matrix is numerically equal to some other matrix with
    /// the same shape, using a default tolerance which is the looser of the
    /// default tolerances of the two objects being compared.
    template <class E2, int RS2>
    bool isNumericallyEqual(const SymMat<M,E2,RS2>& m) const {
        const double tol = std::max(getDefaultTolerance(),m.getDefaultTolerance());
        return isNumericallyEqual(m, tol);
    }

    /// %Test whether this is numerically a "scalar" matrix, meaning that it is 
    /// a diagonal matrix in which each diagonal element is numerically equal to 
    /// the same scalar, using either a specified tolerance or the matrix's 
    /// default tolerance (which is always the same or looser than the default
    /// tolerance for one of its elements).
    bool isNumericallyEqual
       (const ELT& e,
        double     tol = getDefaultTolerance()) const 
    {
        if (!diag().isNumericallyEqual(e, tol))
            return false;
        return getLower().isNumericallyEqual(ELT(0), tol);
    }

    // Rows and columns have to be copied and Hermitian elements have to
    // be conjugated at a floating point cost. This isn't the best way
    // to work with a symmetric matrix.
    TRow row(int i) const {
        SimTK_INDEXCHECK(i,M,"SymMat::row[i]");
        TRow rowi;
        // Columns left of diagonal are lower.
        for (int j=0; j<i; ++j)
            rowi[j] = getEltLower(i,j);
        rowi[i] = getEltDiag(i);
        for (int j=i+1; j<M; ++j)
            rowi[j] = getEltUpper(i,j); // conversion from EHerm to E may cost flops
        return rowi;
    }

    TCol col(int j) const {
        SimTK_INDEXCHECK(j,M,"SymMat::col(j)");
        TCol colj;
        // Rows above diagonal are upper (with conjugated elements).
        for (int i=0; i<j; ++i)
            colj[i] = getEltUpper(i,j); // conversion from EHerm to E may cost flops
        colj[j] = getEltDiag(j);
        for (int i=j+1; i<M; ++i)
            colj[i] = getEltLower(i,j);
        return colj;
    }

    /// Return a value for \e any element of a symmetric matrix, even those
    /// in the upper triangle which aren't actually stored anywhere. For elements
    /// whose underlying numeric types are complex, this will require computation
    /// in order to return the conjugates. So we always have to copy out the
    /// element, and may also have to conjugate it (one flop per complex number).
    E elt(int i, int j) const {
        SimTK_INDEXCHECK(i,M,"SymMat::elt(i,j)");
        SimTK_INDEXCHECK(j,M,"SymMat::elt(i,j)");
        if      (i>j)  return getEltLower(i,j); // copy element
        else if (i==j) return getEltDiag(i);    // copy element
        else           return getEltUpper(i,j); // conversion from EHerm to E may cost flops 
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

    const E& getEltDiag(int i) const {return getDiag()[i];}
    E&       updEltDiag(int i)       {return updDiag()[i];}

    // must be i > j
    const E& getEltLower(int i, int j) const {return getLower()[lowerIx(i,j)];}
    E&       updEltLower(int i, int j)       {return updLower()[lowerIx(i,j)];}

    // must be i < j
    const EHerm& getEltUpper(int i, int j) const {return getUpper()[lowerIx(j,i)];}
    EHerm&       updEltUpper(int i, int j)       {return updUpper()[lowerIx(j,i)];}
    
    /// Returns a row vector (Row) containing the column sums of this matrix.
    TRow colSum() const {
        TRow temp(~getDiag());
        for (int i = 1; i < M; ++i)
            for (int j = 0; j < i; ++j) {
                const E& value = getEltLower(i, j);;
                temp[j] += value;
                temp[i] += E(reinterpret_cast<const EHerm&>(value));
            }
        return temp;
    }
    /// This is an alternate name for colSum(); behaves like the Matlab
    /// function of the same name.
    TRow sum() const {return colSum();}

    /// Returns a column vector (Vec) containing the row sums of this matrix.
    /// This will be the conjugate transpose of the column sums.
    TCol rowSum() const {
        TCol temp(getDiag());
        for (int i = 1; i < M; ++i)
            for (int j = 0; j < i; ++j) {
                const E& value = getEltLower(i, j);;
                temp[i] += value;
                temp[j] += E(reinterpret_cast<const EHerm&>(value));
            }
        return temp;
    }
private:
    E d[NActualElements];

    // This utility doesn't turn lower or upper into a Vec which could turn
    // out to have zero length if this is a 1x1 matrix.
    const E& getlowerE(int i) const {return d[(M+i)*RS];}
    E& updlowerE(int i) {return d[(M+i)*RS];}

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
operator*(const SymMat<M,E1,S1>& l, const SymMat<M,E2,S2>& r) {
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
