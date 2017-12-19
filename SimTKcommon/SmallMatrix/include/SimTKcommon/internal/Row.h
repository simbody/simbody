#ifndef SimTK_SIMMATRIX_SMALLMATRIX_ROW_H_
#define SimTK_SIMMATRIX_SMALLMATRIX_ROW_H_

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
 * Contributors: Peter Eastman                                                *
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
 * This file declares class Row<NCOLS, ELEMENT_TYPE, STRIDE>.
 */

#include "SimTKcommon/internal/common.h"


namespace SimTK {

// The following functions are used internally by Row.

namespace Impl {

// For those wimpy compilers that don't unroll short, constant-limit loops, 
// Peter Eastman added these recursive template implementations of 
// elementwise add, subtract, and copy. Sherm added multiply and divide.

template <class E1, int S1, class E2, int S2> void
conformingAdd(const Row<1,E1,S1>& r1, const Row<1,E2,S2>& r2, 
              Row<1,typename CNT<E1>::template Result<E2>::Add>& result) {
    result[0] = r1[0] + r2[0];
}
template <int N, class E1, int S1, class E2, int S2> void
conformingAdd(const Row<N,E1,S1>& r1, const Row<N,E2,S2>& r2, 
              Row<N,typename CNT<E1>::template Result<E2>::Add>& result) {
    conformingAdd(reinterpret_cast<const Row<N-1,E1,S1>&>(r1), 
                  reinterpret_cast<const Row<N-1,E2,S2>&>(r2), 
                  reinterpret_cast<Row<N-1,typename CNT<E1>::
                              template Result<E2>::Add>&>(result));
    result[N-1] = r1[N-1] + r2[N-1];
}

template <class E1, int S1, class E2, int S2> void
conformingSubtract(const Row<1,E1,S1>& r1, const Row<1,E2,S2>& r2, 
                   Row<1,typename CNT<E1>::template Result<E2>::Sub>& result) {
    result[0] = r1[0] - r2[0];
}
template <int N, class E1, int S1, class E2, int S2> void
conformingSubtract(const Row<N,E1,S1>& r1, const Row<N,E2,S2>& r2,
                   Row<N,typename CNT<E1>::template Result<E2>::Sub>& result) {
    conformingSubtract(reinterpret_cast<const Row<N-1,E1,S1>&>(r1), 
                       reinterpret_cast<const Row<N-1,E2,S2>&>(r2), 
                       reinterpret_cast<Row<N-1,typename CNT<E1>::
                                   template Result<E2>::Sub>&>(result));
    result[N-1] = r1[N-1] - r2[N-1];
}

template <class E1, int S1, class E2, int S2> void
elementwiseMultiply(const Row<1,E1,S1>& r1, const Row<1,E2,S2>& r2, 
              Row<1,typename CNT<E1>::template Result<E2>::Mul>& result) {
    result[0] = r1[0] * r2[0];
}
template <int N, class E1, int S1, class E2, int S2> void
elementwiseMultiply(const Row<N,E1,S1>& r1, const Row<N,E2,S2>& r2, 
              Row<N,typename CNT<E1>::template Result<E2>::Mul>& result) {
    elementwiseMultiply(reinterpret_cast<const Row<N-1,E1,S1>&>(r1), 
                        reinterpret_cast<const Row<N-1,E2,S2>&>(r2), 
                        reinterpret_cast<Row<N-1,typename CNT<E1>::
                                    template Result<E2>::Mul>&>(result));
    result[N-1] = r1[N-1] * r2[N-1];
}

template <class E1, int S1, class E2, int S2> void
elementwiseDivide(const Row<1,E1,S1>& r1, const Row<1,E2,S2>& r2, 
              Row<1,typename CNT<E1>::template Result<E2>::Dvd>& result) {
    result[0] = r1[0] / r2[0];
}
template <int N, class E1, int S1, class E2, int S2> void
elementwiseDivide(const Row<N,E1,S1>& r1, const Row<N,E2,S2>& r2, 
              Row<N,typename CNT<E1>::template Result<E2>::Dvd>& result) {
    elementwiseDivide(reinterpret_cast<const Row<N-1,E1,S1>&>(r1), 
                        reinterpret_cast<const Row<N-1,E2,S2>&>(r2), 
                        reinterpret_cast<Row<N-1,typename CNT<E1>::
                                    template Result<E2>::Dvd>&>(result));
    result[N-1] = r1[N-1] / r2[N-1];
}

template <class E1, int S1, class E2, int S2> void
copy(Row<1,E1,S1>& r1, const Row<1,E2,S2>& r2) {
    r1[0] = r2[0];
}
template <int N, class E1, int S1, class E2, int S2> void
copy(Row<N,E1,S1>& r1, const Row<N,E2,S2>& r2) {
    copy(reinterpret_cast<Row<N-1,E1,S1>&>(r1), 
         reinterpret_cast<const Row<N-1,E2,S2>&>(r2));
    r1[N-1] = r2[N-1];
}

}

/** @brief This is a fixed-length row vector designed for no-overhead inline
computation.

@ingroup MatVecUtilities

The %Row type is not commonly used in Simbody user programs; the column vector
class Vec is much more common. Typically %Row objects arise either from
transposing a Vec or selecting rows from a Mat.

@tparam     N       The number of columns in the row vector.
@tparam     ELT     The element type. Must be a composite numerical type (CNT).
The default is ELT=Real.
@tparam     STRIDE  The spacing from one element to the next in memory, as an
integer number of elements of type ELT. The default is STRIDE=1.
**/
template <int N, class ELT, int STRIDE> class Row {
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
        NRows               = 1,
        NCols               = N,
        NPackedElements     = N,
        NActualElements     = N * STRIDE,   // includes trailing gap
        NActualScalars      = CNT<E>::NActualScalars * NActualElements,
        RowSpacing          = NActualElements,
        ColSpacing          = STRIDE,
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

    typedef Row<N,E,STRIDE>                 T;
    typedef Row<N,ENeg,STRIDE>              TNeg;
    typedef Row<N,EWithoutNegator,STRIDE>   TWithoutNegator;

    typedef Row<N,EReal,STRIDE*CNT<E>::RealStrideFactor>         
                                            TReal;
    typedef Row<N,EImag,STRIDE*CNT<E>::RealStrideFactor>         
                                            TImag;
    typedef Row<N,EComplex,STRIDE>          TComplex;
    typedef Vec<N,EHerm,STRIDE>             THerm;
    typedef Vec<N,E,STRIDE>                 TPosTrans;
    typedef E                               TElement;
    typedef Row                             TRow;
    typedef E                               TCol;

    // These are the results of calculations, so are returned in new, packed
    // memory. Be sure to refer to element types here which are also packed.
    typedef Vec<N,ESqrt,1>                  TSqrt;      // Note stride
    typedef Row<N,EAbs,1>                   TAbs;       // Note stride
    typedef Row<N,EStandard,1>              TStandard;
    typedef Vec<N,EInvert,1>                TInvert;    // packed
    typedef Row<N,ENormalize,1>             TNormalize;

    typedef SymMat<N,ESqHermT>              TSqHermT;   // result of self outer product
    typedef EScalarNormSq                   TSqTHerm;   // result of self dot product

    // These recurse right down to the underlying scalar type no matter how
    // deep the elements are.
    typedef EScalar                         Scalar;
    typedef EULessScalar                    ULessScalar;
    typedef ENumber                         Number;
    typedef EStdNumber                      StdNumber;
    typedef EPrecision                      Precision;
    typedef EScalarNormSq                   ScalarNormSq;

    static int size() { return N; }
    static int nrow() { return 1; }
    static int ncol() { return N; }


    // Scalar norm square is sum( conjugate squares of all scalars )
    ScalarNormSq scalarNormSqr() const { 
        ScalarNormSq sum(0);
        for(int i=0;i<N;++i) sum += CNT<E>::scalarNormSqr(d[i*STRIDE]);
        return sum;
    }

    // sqrt() is elementwise square root; that is, the return value has the same
    // dimension as this Vec but with each element replaced by whatever it thinks
    // its square root is.
    TSqrt sqrt() const {
        TSqrt rsqrt;
        for(int i=0;i<N;++i) rsqrt[i] = CNT<E>::sqrt(d[i*STRIDE]);
        return rsqrt;
    }

    // abs() is elementwise absolute value; that is, the return value has the same
    // dimension as this Row but with each element replaced by whatever it thinks
    // its absolute value is.
    TAbs abs() const {
        TAbs rabs;
        for(int i=0;i<N;++i) rabs[i] = CNT<E>::abs(d[i*STRIDE]);
        return rabs;
    }

    TStandard standardize() const {
        TStandard rstd;
        for(int i=0;i<N;++i) rstd[i] = CNT<E>::standardize(d[i*STRIDE]);
        return rstd;
    }

    // Sum just adds up all the elements, getting rid of negators and
    // conjugates in the process.
    EStandard sum() const {
        E sum(0);
        for (int i=0;i<N;++i) sum += d[i*STRIDE];
        return CNT<E>::standardize(sum);
    }

    // This gives the resulting rowvector type when (v[i] op P) is applied to each element of v.
    // It is a row of length N, stride 1, and element types which are the regular
    // composite result of E op P. Typically P is a scalar type but it doesn't have to be.
    template <class P> struct EltResult { 
        typedef Row<N, typename CNT<E>::template Result<P>::Mul, 1> Mul;
        typedef Row<N, typename CNT<E>::template Result<P>::Dvd, 1> Dvd;
        typedef Row<N, typename CNT<E>::template Result<P>::Add, 1> Add;
        typedef Row<N, typename CNT<E>::template Result<P>::Sub, 1> Sub;
    };

    // This is the composite result for v op P where P is some kind of appropriately shaped
    // non-scalar type.
    template <class P> struct Result { 
        typedef MulCNTs<1,N,ArgDepth,Row,ColSpacing,RowSpacing,
            CNT<P>::NRows, CNT<P>::NCols, CNT<P>::ArgDepth,
            P, CNT<P>::ColSpacing, CNT<P>::RowSpacing> MulOp;
        typedef typename MulOp::Type Mul;

        typedef MulCNTsNonConforming<1,N,ArgDepth,Row,ColSpacing,RowSpacing,
            CNT<P>::NRows, CNT<P>::NCols, CNT<P>::ArgDepth,
            P, CNT<P>::ColSpacing, CNT<P>::RowSpacing> MulOpNonConforming;
        typedef typename MulOpNonConforming::Type MulNon;


        typedef DvdCNTs<1,N,ArgDepth,Row,ColSpacing,RowSpacing,
            CNT<P>::NRows, CNT<P>::NCols, CNT<P>::ArgDepth,
            P, CNT<P>::ColSpacing, CNT<P>::RowSpacing> DvdOp;
        typedef typename DvdOp::Type Dvd;

        typedef AddCNTs<1,N,ArgDepth,Row,ColSpacing,RowSpacing,
            CNT<P>::NRows, CNT<P>::NCols, CNT<P>::ArgDepth,
            P, CNT<P>::ColSpacing, CNT<P>::RowSpacing> AddOp;
        typedef typename AddOp::Type Add;

        typedef SubCNTs<1,N,ArgDepth,Row,ColSpacing,RowSpacing,
            CNT<P>::NRows, CNT<P>::NCols, CNT<P>::ArgDepth,
            P, CNT<P>::ColSpacing, CNT<P>::RowSpacing> SubOp;
        typedef typename SubOp::Type Sub;
    };

    // Shape-preserving element substitution (always packed)
    template <class P> struct Substitute {
        typedef Row<N,P> Type;
    };

    // Default construction initializes to NaN when debugging but
    // is left uninitialized otherwise.
    Row(){ 
    #ifndef NDEBUG
        setToNaN();
    #endif
    }

    // It's important not to use the default copy constructor or copy
    // assignment because the compiler doesn't understand that we may
    // have noncontiguous storage and will try to copy the whole array.
    Row(const Row& src) {
        Impl::copy(*this, src);
    }
    Row& operator=(const Row& src) {    // no harm if src and 'this' are the same
        Impl::copy(*this, src);
        return *this;
    }

    // We want an implicit conversion from a Row of the same length
    // and element type but with a different stride.
    template <int SS> Row(const Row<N,E,SS>& src) {
        Impl::copy(*this, src);
    }

    // We want an implicit conversion from a Row of the same length
    // and *negated* element type, possibly with a different stride.
    template <int SS> Row(const Row<N,ENeg,SS>& src) {
        Impl::copy(*this, src);
    }

    // Construct a Row from a Row of the same length, with any
    // stride. Works as long as the element types are compatible.
    template <class EE, int SS> explicit Row(const Row<N,EE,SS>& vv) {
        Impl::copy(*this, vv);
    }

    // Construction using an element assigns to each element.
    explicit Row(const E& e)
      { for (int i=0;i<N;++i) d[i*STRIDE]=e; }

    // Construction using a negated element assigns to each element.
    explicit Row(const ENeg& ne)
      { for (int i=0;i<N;++i) d[i*STRIDE]=ne; }

    // Given an int, turn it into a suitable floating point number
    // and then feed that to the above single-element constructor.
    explicit Row(int i) 
      { new (this) Row(E(Precision(i))); }

    // A bevy of constructors for Rows up to length 6.
    Row(const E& e0,const E& e1)
      { assert(N==2);(*this)[0]=e0;(*this)[1]=e1; }
    Row(const E& e0,const E& e1,const E& e2)
      { assert(N==3);(*this)[0]=e0;(*this)[1]=e1;(*this)[2]=e2; }
    Row(const E& e0,const E& e1,const E& e2,const E& e3)
      { assert(N==4);(*this)[0]=e0;(*this)[1]=e1;(*this)[2]=e2;(*this)[3]=e3; }
    Row(const E& e0,const E& e1,const E& e2,const E& e3,const E& e4)
      { assert(N==5);(*this)[0]=e0;(*this)[1]=e1;(*this)[2]=e2;
        (*this)[3]=e3;(*this)[4]=e4; }
    Row(const E& e0,const E& e1,const E& e2,const E& e3,const E& e4,const E& e5)
      { assert(N==6);(*this)[0]=e0;(*this)[1]=e1;(*this)[2]=e2;
        (*this)[3]=e3;(*this)[4]=e4;(*this)[5]=e5; }
    Row(const E& e0,const E& e1,const E& e2,const E& e3,const E& e4,const E& e5,const E& e6)
      { assert(N==7);(*this)[0]=e0;(*this)[1]=e1;(*this)[2]=e2;
        (*this)[3]=e3;(*this)[4]=e4;(*this)[5]=e5;(*this)[6]=e6; }
    Row(const E& e0,const E& e1,const E& e2,const E& e3,const E& e4,const E& e5,const E& e6,const E& e7)
      { assert(N==8);(*this)[0]=e0;(*this)[1]=e1;(*this)[2]=e2;
        (*this)[3]=e3;(*this)[4]=e4;(*this)[5]=e5;(*this)[6]=e6;(*this)[7]=e7; }
    Row(const E& e0,const E& e1,const E& e2,const E& e3,const E& e4,const E& e5,const E& e6,const E& e7,const E& e8)
      { assert(N==9);(*this)[0]=e0;(*this)[1]=e1;(*this)[2]=e2;
        (*this)[3]=e3;(*this)[4]=e4;(*this)[5]=e5;(*this)[6]=e6;(*this)[7]=e7;(*this)[8]=e8; }

    // Construction from a pointer to anything assumes we're pointing
    // at an element list of the right length.
    template <class EE> explicit Row(const EE* p)
      { assert(p); for(int i=0;i<N;++i) d[i*STRIDE]=p[i]; }
    template <class EE> Row& operator=(const EE* p)
      { assert(p); for(int i=0;i<N;++i) d[i*STRIDE]=p[i]; return *this; }

    // Conforming assignment ops.
    template <class EE, int SS> Row& operator=(const Row<N,EE,SS>& vv) {
        Impl::copy(*this, vv);
        return *this;
    }
    template <class EE, int SS> Row& operator+=(const Row<N,EE,SS>& r)
      { for(int i=0;i<N;++i) d[i*STRIDE] += r[i]; return *this; }
    template <class EE, int SS> Row& operator+=(const Row<N,negator<EE>,SS>& r)
      { for(int i=0;i<N;++i) d[i*STRIDE] -= -(r[i]); return *this; }
    template <class EE, int SS> Row& operator-=(const Row<N,EE,SS>& r)
      { for(int i=0;i<N;++i) d[i*STRIDE] -= r[i]; return *this; }
    template <class EE, int SS> Row& operator-=(const Row<N,negator<EE>,SS>& r)
      { for(int i=0;i<N;++i) d[i*STRIDE] += -(r[i]); return *this; }

    // Conforming binary ops with 'this' on left, producing new packed result.
    // Cases: r=r+r, r=r-r, s=r*v r=r*m

    /** Vector addition -- use operator+ instead. **/
    template <class EE, int SS> Row<N,typename CNT<E>::template Result<EE>::Add>
    conformingAdd(const Row<N,EE,SS>& r) const {
        Row<N,typename CNT<E>::template Result<EE>::Add> result;
        Impl::conformingAdd(*this, r, result);
        return result;
    }

    /** Vector subtraction -- use operator- instead. **/
    template <class EE, int SS> Row<N,typename CNT<E>::template Result<EE>::Sub>
    conformingSubtract(const Row<N,EE,SS>& r) const {
        Row<N,typename CNT<E>::template Result<EE>::Sub> result;
        Impl::conformingSubtract(*this, r, result);
        return result;
    }

    /** Same as dot product (s = row*col) -- use operator* or dot() instead. **/
    template <class EE, int SS> typename CNT<E>::template Result<EE>::Mul
    conformingMultiply(const Vec<N,EE,SS>& r) const {
        return (*this)*r;
    }

    /** Row times a conforming matrix, row=row*mat -- use operator* instead. **/
    template <int MatNCol, class EE, int CS, int RS> 
    Row<MatNCol,typename CNT<E>::template Result<EE>::Mul>
    conformingMultiply(const Mat<N,MatNCol,EE,CS,RS>& m) const {
        Row<MatNCol,typename CNT<E>::template Result<EE>::Mul> result;
        for (int j=0;j<N;++j) result[j] = conformingMultiply(m(j));
        return result;
    }

    /** Elementwise multiply (Matlab .* operator). **/
    template <class EE, int SS> Row<N,typename CNT<E>::template Result<EE>::Mul>
    elementwiseMultiply(const Row<N,EE,SS>& r) const {
        Row<N,typename CNT<E>::template Result<EE>::Mul> result;
        Impl::elementwiseMultiply(*this, r, result);
        return result;
    }

    /** Elementwise divide (Matlab ./ operator). **/
    template <class EE, int SS> Row<N,typename CNT<E>::template Result<EE>::Dvd>
    elementwiseDivide(const Row<N,EE,SS>& r) const {
        Row<N,typename CNT<E>::template Result<EE>::Dvd> result;
        Impl::elementwiseDivide(*this, r, result);
        return result;
    }

    const E& operator[](int i) const { assert(0 <= i && i < N); return d[i*STRIDE]; }
    E&       operator[](int i)         { assert(0 <= i && i < N); return d[i*STRIDE]; }
    const E& operator()(int i) const { return (*this)[i]; }
    E&       operator()(int i)         { return (*this)[i]; }

    ScalarNormSq normSqr() const { return scalarNormSqr(); }
    typename CNT<ScalarNormSq>::TSqrt 
        norm() const { return CNT<ScalarNormSq>::sqrt(scalarNormSqr()); }

    // If the elements of this Row are scalars, the result is what you get by
    // dividing each element by the norm() calculated above. If the elements are
    // *not* scalars, then the elements are *separately* normalized. That means
    // you will get a different answer from Row<2,Row3>::normalize() than you
    // would from a Row<6>::normalize() containing the same scalars.
    //
    // Normalize returns a row of the same dimension but in new, packed storage
    // and with a return type that does not include negator<> even if the original
    // Row<> does, because we can eliminate the negation here almost for free.
    // But we can't standardize (change conjugate to complex) for free, so we'll retain
    // conjugates if there are any.
    TNormalize normalize() const {
        if (CNT<E>::IsScalar) {
            return castAwayNegatorIfAny() / (SignInterpretation*norm());
        } else {
            TNormalize elementwiseNormalized;
            for (int j=0; j<N; ++j) 
                elementwiseNormalized[j] = CNT<E>::normalize((*this)[j]);
            return elementwiseNormalized;
        }
    }

    TInvert invert() const {assert(false); return TInvert();} // TODO default inversion

    const Row&   operator+() const { return *this; }
    const TNeg&  operator-() const { return negate(); }
    TNeg&        operator-()       { return updNegate(); }
    const THerm& operator~() const { return transpose(); }
    THerm&       operator~()       { return updTranspose(); }

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


    // These are elementwise binary operators, (this op ee) by default but 
    // (ee op this) if 'FromLeft' appears in the name. The result is a packed 
    // Row<N> but the element type may change. These are mostly used to 
    // implement global operators. We call these "scalar" operators but 
    // actually the "scalar" can be a composite type.

    //TODO: consider converting 'e' to Standard Numbers as precalculation and 
    // changing return type appropriately.
    template <class EE> Row<N, typename CNT<E>::template Result<EE>::Mul>
    scalarMultiply(const EE& e) const {
        Row<N, typename CNT<E>::template Result<EE>::Mul> result;
        for (int j=0; j<N; ++j) result[j] = (*this)[j] * e;
        return result;
    }
    template <class EE> Row<N, typename CNT<EE>::template Result<E>::Mul>
    scalarMultiplyFromLeft(const EE& e) const {
        Row<N, typename CNT<EE>::template Result<E>::Mul> result;
        for (int j=0; j<N; ++j) result[j] = e * (*this)[j];
        return result;
    }

    // TODO: should precalculate and store 1/e, while converting to Standard 
    // Numbers. Note that return type should change appropriately.
    template <class EE> Row<N, typename CNT<E>::template Result<EE>::Dvd>
    scalarDivide(const EE& e) const {
        Row<N, typename CNT<E>::template Result<EE>::Dvd> result;
        for (int j=0; j<N; ++j) result[j] = (*this)[j] / e;
        return result;
    }
    template <class EE> Row<N, typename CNT<EE>::template Result<E>::Dvd>
    scalarDivideFromLeft(const EE& e) const {
        Row<N, typename CNT<EE>::template Result<E>::Dvd> result;
        for (int j=0; j<N; ++j) result[j] = e / (*this)[j];
        return result;
    }

    template <class EE> Row<N, typename CNT<E>::template Result<EE>::Add>
    scalarAdd(const EE& e) const {
        Row<N, typename CNT<E>::template Result<EE>::Add> result;
        for (int j=0; j<N; ++j) result[j] = (*this)[j] + e;
        return result;
    }
    // Add is commutative, so no 'FromLeft'.

    template <class EE> Row<N, typename CNT<E>::template Result<EE>::Sub>
    scalarSubtract(const EE& e) const {
        Row<N, typename CNT<E>::template Result<EE>::Sub> result;
        for (int j=0; j<N; ++j) result[j] = (*this)[j] - e;
        return result;
    }
    template <class EE> Row<N, typename CNT<EE>::template Result<E>::Sub>
    scalarSubtractFromLeft(const EE& e) const {
        Row<N, typename CNT<EE>::template Result<E>::Sub> result;
        for (int j=0; j<N; ++j) result[j] = e - (*this)[j];
        return result;
    }

    // Generic assignments for any element type not listed explicitly, including scalars.
    // These are done repeatedly for each element and only work if the operation can
    // be performed leaving the original element type.
    template <class EE> Row& operator =(const EE& e) {return scalarEq(e);}
    template <class EE> Row& operator+=(const EE& e) {return scalarPlusEq(e);}
    template <class EE> Row& operator-=(const EE& e) {return scalarMinusEq(e);}
    template <class EE> Row& operator*=(const EE& e) {return scalarTimesEq(e);}
    template <class EE> Row& operator/=(const EE& e) {return scalarDivideEq(e);}

    // Generalized scalar assignment & computed assignment methods. These will work
    // for any assignment-compatible element, not just scalars.
    template <class EE> Row& scalarEq(const EE& ee)
      { for(int i=0;i<N;++i) d[i*STRIDE] = ee; return *this; }
    template <class EE> Row& scalarPlusEq(const EE& ee)
      { for(int i=0;i<N;++i) d[i*STRIDE] += ee; return *this; }
    template <class EE> Row& scalarMinusEq(const EE& ee)
      { for(int i=0;i<N;++i) d[i*STRIDE] -= ee; return *this; }
    template <class EE> Row& scalarMinusEqFromLeft(const EE& ee)
      { for(int i=0;i<N;++i) d[i*STRIDE] = ee - d[i*STRIDE]; return *this; }
    template <class EE> Row& scalarTimesEq(const EE& ee)
      { for(int i=0;i<N;++i) d[i*STRIDE] *= ee; return *this; }
    template <class EE> Row& scalarTimesEqFromLeft(const EE& ee)
      { for(int i=0;i<N;++i) d[i*STRIDE] = ee * d[i*STRIDE]; return *this; }
    template <class EE> Row& scalarDivideEq(const EE& ee)
      { for(int i=0;i<N;++i) d[i*STRIDE] /= ee; return *this; }
    template <class EE> Row& scalarDivideEqFromLeft(const EE& ee)
      { for(int i=0;i<N;++i) d[i*STRIDE] = ee / d[i*STRIDE]; return *this; }


    // Specialize for int to avoid warnings and ambiguities.
    Row& scalarEq(int ee)       {return scalarEq(Precision(ee));}
    Row& scalarPlusEq(int ee)   {return scalarPlusEq(Precision(ee));}
    Row& scalarMinusEq(int ee)  {return scalarMinusEq(Precision(ee));}
    Row& scalarTimesEq(int ee)  {return scalarTimesEq(Precision(ee));}
    Row& scalarDivideEq(int ee) {return scalarDivideEq(Precision(ee));}
    Row& scalarMinusEqFromLeft(int ee)  {return scalarMinusEqFromLeft(Precision(ee));}
    Row& scalarTimesEqFromLeft(int ee)  {return scalarTimesEqFromLeft(Precision(ee));}
    Row& scalarDivideEqFromLeft(int ee) {return scalarDivideEqFromLeft(Precision(ee));}

    /** Set every scalar in this %Row to NaN; this is the default initial
    value in Debug builds, but not in Release. **/
    void setToNaN() {
        (*this) = CNT<ELT>::getNaN();
    }

    /** Set every scalar in this %Row to zero. **/
    void setToZero() {
        (*this) = ELT(0);
    }

    /** Extract a const reference to a sub-Row with size known at compile time. 
    This must be called with an explicit template argument for the size, for
    example, getSubRow<3>(j). This is only a recast; no copying or computation
    is performed. The size and index are range checked in Debug builds but
    not in Release builds. **/
    template <int NN>
    const Row<NN,ELT,STRIDE>& getSubRow(int j) const {
        assert(0 <= j && j + NN <= N);
        return Row<NN,ELT,STRIDE>::getAs(&(*this)[j]);
    }
    /** Extract a writable reference to a sub-Row with size known at compile time. 
    This must be called with an explicit template argument for the size, for
    example, updSubRow<3>(j). This is only a recast; no copying or computation
    is performed. The size and index are range checked in Debug builds but
    not in Release builds. **/
    template <int NN>
    Row<NN,ELT,STRIDE>& updSubRow(int j) {
        assert(0 <= j && j + NN <= N);
        return Row<NN,ELT,STRIDE>::updAs(&(*this)[j]);
    }

    /** Extract a subvector of type %Row from a longer one that has the same
    element type and stride, and return a const reference to the selected 
    subsequence. **/
    template <int NN>
    static const Row& getSubRow(const Row<NN,ELT,STRIDE>& r, int j) {
        assert(0 <= j && j + N <= NN);
        return getAs(&r[j]);
    }
    /** Extract a subvector of type %Row from a longer one that has the same
    element type and stride, and return a writable reference to the selected 
    subsequence. **/
    template <int NN>
    static Row& updSubRow(Row<NN,ELT,STRIDE>& r, int j) {
        assert(0 <= j && j + N <= NN);
        return updAs(&r[j]);
    }

    /** Return a row one smaller than this one by dropping the element
    at the indicated position p. The result is a packed copy with the same
    element type as this one. **/
    Row<N-1,ELT,1> drop1(int p) const {
        assert(0 <= p && p < N);
        Row<N-1,ELT,1> out;
        int nxt=0;
        for (int i=0; i<N-1; ++i, ++nxt) {
            if (nxt==p) ++nxt;  // skip the loser
            out[i] = (*this)[nxt];
        }
        return out;
    }

    /** Return a row one larger than this one by adding an element
    to the end. The result is a packed copy with the same element type as
    this one. Works for any assignment compatible element. **/
    template <class EE> Row<N+1,ELT,1> append1(const EE& v) const {
        Row<N+1,ELT,1> out;
        Row<N,ELT,1>::updAs(&out[0]) = (*this);
        out[N] = v;
        return out;
    }


    /** Return a row one larger than this one by inserting an element
    \e before the indicated one. The result is a packed copy with the same 
    element type as this one. Works for any assignment compatible element. The 
    index can be one greater than normally allowed in which case the element
    is appended (but use append1() if you know you're appending). **/
    template <class EE> Row<N+1,ELT,1> insert1(int p, const EE& v) const {
        assert(0 <= p && p <= N);
        if (p==N) return append1(v);
        Row<N+1,ELT,1> out;
        int nxt=0;
        for (int i=0; i<N; ++i, ++nxt) {
            if (i==p) out[nxt++] = v;
            out[nxt] = (*this)[i];
        }
        return out;
    }

    /** Recast an ordinary C++ array E[] to a const %Row<N,E,S>; assumes 
    compatible length, stride, and packing. **/
    static const Row& getAs(const ELT* p)  {return *reinterpret_cast<const Row*>(p);}
    /** Recast a writable ordinary C++ array E[] to a writable %Row<N,E,S>; 
    assumes compatible length, stride, and packing. **/
    static Row&       updAs(ELT* p)        {return *reinterpret_cast<Row*>(p);}

    /** Return a %Row of the same length and element type as this one but
    with all elements set to NaN. The result is packed (stride==1) regardless
    of the stride of this %Row. **/
    static Row<N,ELT,1> getNaN() { return Row<N,ELT,1>(CNT<ELT>::getNaN()); }

    /** Return true if any element of this Row contains a NaN anywhere. **/
    bool isNaN() const {
        for (int j=0; j<N; ++j)
            if (CNT<ELT>::isNaN((*this)[j]))
                return true;
        return false;
    }

    /** Return true if any element of this Row contains a +Infinity
    or -Infinity somewhere but no element contains a NaN anywhere. **/
    bool isInf() const {
        bool seenInf = false;
        for (int j=0; j<N; ++j) {
            const ELT& e = (*this)[j];
            if (!CNT<ELT>::isFinite(e)) {
                if (!CNT<ELT>::isInf(e)) 
                    return false; // something bad was found
                seenInf = true; 
            }
        }
        return seenInf;
    }

    /** Return true if no element of this %Row contains an Infinity or a NaN 
    anywhere. **/
    bool isFinite() const {
        for (int j=0; j<N; ++j)
            if (!CNT<ELT>::isFinite((*this)[j]))
                return false;
        return true;
    }

    /** For approximate comparisons, the default tolerance to use for a vector is
    the same as its elements' default tolerance. **/
    static double getDefaultTolerance() {return CNT<ELT>::getDefaultTolerance();}

    /** %Test whether this row is numerically equal to some other row with
    the same shape, using a specified tolerance. **/
    template <class E2, int CS2>
    bool isNumericallyEqual(const Row<N,E2,CS2>& r, double tol) const {
        for (int j=0; j<N; ++j)
            if (!CNT<ELT>::isNumericallyEqual((*this)(j), r(j), tol))
                return false;
        return true;
    }

    /** %Test whether this row vector is numerically equal to some other row with
    the same shape, using a default tolerance which is the looser of the
    default tolerances of the two objects being compared. **/
    template <class E2, int CS2>
    bool isNumericallyEqual(const Row<N,E2,CS2>& r) const {
        const double tol = std::max(getDefaultTolerance(),r.getDefaultTolerance());
        return isNumericallyEqual(r, tol);
    }

    /** %Test whether every element of this row vector is numerically equal to
    the given element, using either a specified tolerance or the row's 
    default tolerance (which is always the same or looser than the default
    tolerance for one of its elements). **/
    bool isNumericallyEqual
       (const ELT& e,
        double     tol = getDefaultTolerance()) const 
    {
        for (int j=0; j<N; ++j)
            if (!CNT<ELT>::isNumericallyEqual((*this)(j), e, tol))
                return false;
        return true;
    }
private:
    ELT d[NActualElements];    // data
};

/////////////////////////////////////////////
// Global operators involving two rows.    //
//   v+v, v-v, v==v, v!=v                  //
/////////////////////////////////////////////

// v3 = v1 + v2 where all v's have the same length N. 
template <int N, class E1, int S1, class E2, int S2> inline
typename Row<N,E1,S1>::template Result< Row<N,E2,S2> >::Add
operator+(const Row<N,E1,S1>& l, const Row<N,E2,S2>& r) { 
    return Row<N,E1,S1>::template Result< Row<N,E2,S2> >
        ::AddOp::perform(l,r);
}

// v3 = v1 - v2, similar to +
template <int N, class E1, int S1, class E2, int S2> inline
typename Row<N,E1,S1>::template Result< Row<N,E2,S2> >::Sub
operator-(const Row<N,E1,S1>& l, const Row<N,E2,S2>& r) { 
    return Row<N,E1,S1>::template Result< Row<N,E2,S2> >
        ::SubOp::perform(l,r);
}

/// bool = v1[i] == v2[i], for all elements i
template <int N, class E1, int S1, class E2, int S2> inline bool
operator==(const Row<N,E1,S1>& l, const Row<N,E2,S2>& r) { 
    for (int i=0; i < N; ++i) if (l[i] != r[i]) return false;
    return true;
}
/// bool = v1[i] != v2[i], for any element i
template <int N, class E1, int S1, class E2, int S2> inline bool
operator!=(const Row<N,E1,S1>& l, const Row<N,E2,S2>& r) {return !(l==r);} 

/// bool = v1[i] < v2[i], for all elements i
template <int N, class E1, int S1, class E2, int S2> inline bool
operator<(const Row<N,E1,S1>& l, const Row<N,E2,S2>& r) 
{   for (int i=0; i < N; ++i) if (l[i] >= r[i]) return false;
    return true; }
/// bool = v[i] < e, for all elements v[i] and element e
template <int N, class E1, int S1, class E2> inline bool
operator<(const Row<N,E1,S1>& v, const E2& e) 
{   for (int i=0; i < N; ++i) if (v[i] >= e) return false;
    return true; }

/// bool = v1[i] > v2[i], for all elements i
template <int N, class E1, int S1, class E2, int S2> inline bool
operator>(const Row<N,E1,S1>& l, const Row<N,E2,S2>& r) 
{   for (int i=0; i < N; ++i) if (l[i] <= r[i]) return false;
    return true; }
/// bool = v[i] > e, for all elements v[i] and element e
template <int N, class E1, int S1, class E2> inline bool
operator>(const Row<N,E1,S1>& v, const E2& e) 
{   for (int i=0; i < N; ++i) if (v[i] <= e) return false;
    return true; }

/// bool = v1[i] <= v2[i], for all elements i.
/// This is not the same as !(v1>v2).
template <int N, class E1, int S1, class E2, int S2> inline bool
operator<=(const Row<N,E1,S1>& l, const Row<N,E2,S2>& r) 
{   for (int i=0; i < N; ++i) if (l[i] > r[i]) return false;
    return true; }
/// bool = v[i] <= e, for all elements v[i] and element e.
/// This is not the same as !(v1>e).
template <int N, class E1, int S1, class E2> inline bool
operator<=(const Row<N,E1,S1>& v, const E2& e) 
{   for (int i=0; i < N; ++i) if (v[i] > e) return false;
    return true; }

/// bool = v1[i] >= v2[i], for all elements i
/// This is not the same as !(v1<v2).
template <int N, class E1, int S1, class E2, int S2> inline bool
operator>=(const Row<N,E1,S1>& l, const Row<N,E2,S2>& r) 
{   for (int i=0; i < N; ++i) if (l[i] < r[i]) return false;
    return true; }
/// bool = v[i] >= e, for all elements v[i] and element e.
/// This is not the same as !(v1<e).
template <int N, class E1, int S1, class E2> inline bool
operator>=(const Row<N,E1,S1>& v, const E2& e) 
{   for (int i=0; i < N; ++i) if (v[i] < e) return false;
    return true; }

////////////////////////////////////////////////////
// Global operators involving a row and a scalar. //
////////////////////////////////////////////////////

// I haven't been able to figure out a nice way to templatize for the
// built-in reals without introducing a lot of unwanted type matches
// as well. So we'll just grind them out explicitly here.

// SCALAR MULTIPLY

// v = v*real, real*v 
template <int N, class E, int S> inline
typename Row<N,E,S>::template Result<float>::Mul
operator*(const Row<N,E,S>& l, const float& r)
  { return Row<N,E,S>::template Result<float>::MulOp::perform(l,r); }
template <int N, class E, int S> inline
typename Row<N,E,S>::template Result<float>::Mul
operator*(const float& l, const Row<N,E,S>& r) {return r*l;}

template <int N, class E, int S> inline
typename Row<N,E,S>::template Result<double>::Mul
operator*(const Row<N,E,S>& l, const double& r)
  { return Row<N,E,S>::template Result<double>::MulOp::perform(l,r); }
template <int N, class E, int S> inline
typename Row<N,E,S>::template Result<double>::Mul
operator*(const double& l, const Row<N,E,S>& r) {return r*l;}

// v = v*int, int*v -- just convert int to v's precision float
template <int N, class E, int S> inline
typename Row<N,E,S>::template Result<typename CNT<E>::Precision>::Mul
operator*(const Row<N,E,S>& l, int r) {return l * (typename CNT<E>::Precision)r;}
template <int N, class E, int S> inline
typename Row<N,E,S>::template Result<typename CNT<E>::Precision>::Mul
operator*(int l, const Row<N,E,S>& r) {return r * (typename CNT<E>::Precision)l;}

// Complex, conjugate, and negator are all easy to templatize.

// v = v*complex, complex*v
template <int N, class E, int S, class R> inline
typename Row<N,E,S>::template Result<std::complex<R> >::Mul
operator*(const Row<N,E,S>& l, const std::complex<R>& r)
  { return Row<N,E,S>::template Result<std::complex<R> >::MulOp::perform(l,r); }
template <int N, class E, int S, class R> inline
typename Row<N,E,S>::template Result<std::complex<R> >::Mul
operator*(const std::complex<R>& l, const Row<N,E,S>& r) {return r*l;}

// v = v*conjugate, conjugate*v (convert conjugate->complex)
template <int N, class E, int S, class R> inline
typename Row<N,E,S>::template Result<std::complex<R> >::Mul
operator*(const Row<N,E,S>& l, const conjugate<R>& r) {return l*(std::complex<R>)r;}
template <int N, class E, int S, class R> inline
typename Row<N,E,S>::template Result<std::complex<R> >::Mul
operator*(const conjugate<R>& l, const Row<N,E,S>& r) {return r*(std::complex<R>)l;}

// v = v*negator, negator*v: convert negator to standard number
template <int N, class E, int S, class R> inline
typename Row<N,E,S>::template Result<typename negator<R>::StdNumber>::Mul
operator*(const Row<N,E,S>& l, const negator<R>& r) {return l * (typename negator<R>::StdNumber)(R)r;}
template <int N, class E, int S, class R> inline
typename Row<N,E,S>::template Result<typename negator<R>::StdNumber>::Mul
operator*(const negator<R>& l, const Row<N,E,S>& r) {return r * (typename negator<R>::StdNumber)(R)l;}


// SCALAR DIVIDE. This is a scalar operation when the scalar is on the right,
// but when it is on the left it means scalar * pseudoInverse(row), which is 
// a vec.

// v = v/real, real/v 
template <int N, class E, int S> inline
typename Row<N,E,S>::template Result<float>::Dvd
operator/(const Row<N,E,S>& l, const float& r)
  { return Row<N,E,S>::template Result<float>::DvdOp::perform(l,r); }
template <int N, class E, int S> inline
typename CNT<float>::template Result<Row<N,E,S> >::Dvd
operator/(const float& l, const Row<N,E,S>& r)
  { return CNT<float>::template Result<Row<N,E,S> >::DvdOp::perform(l,r); }

template <int N, class E, int S> inline
typename Row<N,E,S>::template Result<double>::Dvd
operator/(const Row<N,E,S>& l, const double& r)
  { return Row<N,E,S>::template Result<double>::DvdOp::perform(l,r); }
template <int N, class E, int S> inline
typename CNT<double>::template Result<Row<N,E,S> >::Dvd
operator/(const double& l, const Row<N,E,S>& r)
  { return CNT<double>::template Result<Row<N,E,S> >::DvdOp::perform(l,r); }

// v = v/int, int/v -- just convert int to v's precision float
template <int N, class E, int S> inline
typename Row<N,E,S>::template Result<typename CNT<E>::Precision>::Dvd
operator/(const Row<N,E,S>& l, int r) {return l / (typename CNT<E>::Precision)r;}
template <int N, class E, int S> inline
typename CNT<typename CNT<E>::Precision>::template Result<Row<N,E,S> >::Dvd
operator/(int l, const Row<N,E,S>& r) {return (typename CNT<E>::Precision)l / r;}


// Complex, conjugate, and negator are all easy to templatize.

// v = v/complex, complex/v
template <int N, class E, int S, class R> inline
typename Row<N,E,S>::template Result<std::complex<R> >::Dvd
operator/(const Row<N,E,S>& l, const std::complex<R>& r)
  { return Row<N,E,S>::template Result<std::complex<R> >::DvdOp::perform(l,r); }
template <int N, class E, int S, class R> inline
typename CNT<std::complex<R> >::template Result<Row<N,E,S> >::Dvd
operator/(const std::complex<R>& l, const Row<N,E,S>& r)
  { return CNT<std::complex<R> >::template Result<Row<N,E,S> >::DvdOp::perform(l,r); }

// v = v/conjugate, conjugate/v (convert conjugate->complex)
template <int N, class E, int S, class R> inline
typename Row<N,E,S>::template Result<std::complex<R> >::Dvd
operator/(const Row<N,E,S>& l, const conjugate<R>& r) {return l/(std::complex<R>)r;}
template <int N, class E, int S, class R> inline
typename CNT<std::complex<R> >::template Result<Row<N,E,S> >::Dvd
operator/(const conjugate<R>& l, const Row<N,E,S>& r) {return (std::complex<R>)l/r;}

// v = v/negator, negator/v: convert negator to number
template <int N, class E, int S, class R> inline
typename Row<N,E,S>::template Result<typename negator<R>::StdNumber>::Dvd
operator/(const Row<N,E,S>& l, const negator<R>& r) {return l/(typename negator<R>::StdNumber)(R)r;}
template <int N, class E, int S, class R> inline
typename CNT<R>::template Result<Row<N,E,S> >::Dvd
operator/(const negator<R>& l, const Row<N,E,S>& r) {return (typename negator<R>::StdNumber)(R)l/r;}


// Add and subtract are odd as scalar ops. They behave as though the
// scalar stands for a vector each of whose elements is that scalar,
// and then a normal vector add or subtract is done.

// SCALAR ADD

// v = v+real, real+v 
template <int N, class E, int S> inline
typename Row<N,E,S>::template Result<float>::Add
operator+(const Row<N,E,S>& l, const float& r)
  { return Row<N,E,S>::template Result<float>::AddOp::perform(l,r); }
template <int N, class E, int S> inline
typename Row<N,E,S>::template Result<float>::Add
operator+(const float& l, const Row<N,E,S>& r) {return r+l;}

template <int N, class E, int S> inline
typename Row<N,E,S>::template Result<double>::Add
operator+(const Row<N,E,S>& l, const double& r)
  { return Row<N,E,S>::template Result<double>::AddOp::perform(l,r); }
template <int N, class E, int S> inline
typename Row<N,E,S>::template Result<double>::Add
operator+(const double& l, const Row<N,E,S>& r) {return r+l;}

// v = v+int, int+v -- just convert int to v's precision float
template <int N, class E, int S> inline
typename Row<N,E,S>::template Result<typename CNT<E>::Precision>::Add
operator+(const Row<N,E,S>& l, int r) {return l + (typename CNT<E>::Precision)r;}
template <int N, class E, int S> inline
typename Row<N,E,S>::template Result<typename CNT<E>::Precision>::Add
operator+(int l, const Row<N,E,S>& r) {return r + (typename CNT<E>::Precision)l;}

// Complex, conjugate, and negator are all easy to templatize.

// v = v+complex, complex+v
template <int N, class E, int S, class R> inline
typename Row<N,E,S>::template Result<std::complex<R> >::Add
operator+(const Row<N,E,S>& l, const std::complex<R>& r)
  { return Row<N,E,S>::template Result<std::complex<R> >::AddOp::perform(l,r); }
template <int N, class E, int S, class R> inline
typename Row<N,E,S>::template Result<std::complex<R> >::Add
operator+(const std::complex<R>& l, const Row<N,E,S>& r) {return r+l;}

// v = v+conjugate, conjugate+v (convert conjugate->complex)
template <int N, class E, int S, class R> inline
typename Row<N,E,S>::template Result<std::complex<R> >::Add
operator+(const Row<N,E,S>& l, const conjugate<R>& r) {return l+(std::complex<R>)r;}
template <int N, class E, int S, class R> inline
typename Row<N,E,S>::template Result<std::complex<R> >::Add
operator+(const conjugate<R>& l, const Row<N,E,S>& r) {return r+(std::complex<R>)l;}

// v = v+negator, negator+v: convert negator to standard number
template <int N, class E, int S, class R> inline
typename Row<N,E,S>::template Result<typename negator<R>::StdNumber>::Add
operator+(const Row<N,E,S>& l, const negator<R>& r) {return l + (typename negator<R>::StdNumber)(R)r;}
template <int N, class E, int S, class R> inline
typename Row<N,E,S>::template Result<typename negator<R>::StdNumber>::Add
operator+(const negator<R>& l, const Row<N,E,S>& r) {return r + (typename negator<R>::StdNumber)(R)l;}

// SCALAR SUBTRACT -- careful, not commutative.

// v = v-real, real-v 
template <int N, class E, int S> inline
typename Row<N,E,S>::template Result<float>::Sub
operator-(const Row<N,E,S>& l, const float& r)
  { return Row<N,E,S>::template Result<float>::SubOp::perform(l,r); }
template <int N, class E, int S> inline
typename CNT<float>::template Result<Row<N,E,S> >::Sub
operator-(const float& l, const Row<N,E,S>& r)
  { return CNT<float>::template Result<Row<N,E,S> >::SubOp::perform(l,r); }

template <int N, class E, int S> inline
typename Row<N,E,S>::template Result<double>::Sub
operator-(const Row<N,E,S>& l, const double& r)
  { return Row<N,E,S>::template Result<double>::SubOp::perform(l,r); }
template <int N, class E, int S> inline
typename CNT<double>::template Result<Row<N,E,S> >::Sub
operator-(const double& l, const Row<N,E,S>& r)
  { return CNT<double>::template Result<Row<N,E,S> >::SubOp::perform(l,r); }

// v = v-int, int-v // just convert int to v's precision float
template <int N, class E, int S> inline
typename Row<N,E,S>::template Result<typename CNT<E>::Precision>::Sub
operator-(const Row<N,E,S>& l, int r) {return l - (typename CNT<E>::Precision)r;}
template <int N, class E, int S> inline
typename CNT<typename CNT<E>::Precision>::template Result<Row<N,E,S> >::Sub
operator-(int l, const Row<N,E,S>& r) {return (typename CNT<E>::Precision)l - r;}


// Complex, conjugate, and negator are all easy to templatize.

// v = v-complex, complex-v
template <int N, class E, int S, class R> inline
typename Row<N,E,S>::template Result<std::complex<R> >::Sub
operator-(const Row<N,E,S>& l, const std::complex<R>& r)
  { return Row<N,E,S>::template Result<std::complex<R> >::SubOp::perform(l,r); }
template <int N, class E, int S, class R> inline
typename CNT<std::complex<R> >::template Result<Row<N,E,S> >::Sub
operator-(const std::complex<R>& l, const Row<N,E,S>& r)
  { return CNT<std::complex<R> >::template Result<Row<N,E,S> >::SubOp::perform(l,r); }

// v = v-conjugate, conjugate-v (convert conjugate->complex)
template <int N, class E, int S, class R> inline
typename Row<N,E,S>::template Result<std::complex<R> >::Sub
operator-(const Row<N,E,S>& l, const conjugate<R>& r) {return l-(std::complex<R>)r;}
template <int N, class E, int S, class R> inline
typename CNT<std::complex<R> >::template Result<Row<N,E,S> >::Sub
operator-(const conjugate<R>& l, const Row<N,E,S>& r) {return (std::complex<R>)l-r;}

// v = v-negator, negator-v: convert negator to standard number
template <int N, class E, int S, class R> inline
typename Row<N,E,S>::template Result<typename negator<R>::StdNumber>::Sub
operator-(const Row<N,E,S>& l, const negator<R>& r) {return l-(typename negator<R>::StdNumber)(R)r;}
template <int N, class E, int S, class R> inline
typename CNT<R>::template Result<Row<N,E,S> >::Sub
operator-(const negator<R>& l, const Row<N,E,S>& r) {return (typename negator<R>::StdNumber)(R)l-r;}


// Row I/O
template <int N, class E, int S, class CHAR, class TRAITS> inline
std::basic_ostream<CHAR,TRAITS>&
operator<<(std::basic_ostream<CHAR,TRAITS>& o, const Row<N,E,S>& v) {
    o << "[" << v[0]; for(int i=1;i<N;++i) o<<','<<v[i]; o<<']'; return o;
}

/** Read a Row from a stream as M elements separated by white space or
by commas, optionally enclosed in () or [] (but no leading "~"). **/
template <int N, class E, int S, class CHAR, class TRAITS> inline
std::basic_istream<CHAR,TRAITS>&
operator>>(std::basic_istream<CHAR,TRAITS>& is, Row<N,E,S>& v) {
    CHAR openBracket, closeBracket;
    is >> openBracket; if (is.fail()) return is;
    if (openBracket==CHAR('('))
        closeBracket = CHAR(')');
    else if (openBracket==CHAR('['))
        closeBracket = CHAR(']');
    else {
        closeBracket = CHAR(0);
        is.unget(); if (is.fail()) return is;
    }

    for (int i=0; i < N; ++i) {
        is >> v[i];
        if (is.fail()) return is;
        if (i != N-1) {
            CHAR c; is >> c; if (is.fail()) return is;
            if (c != ',') is.unget();
            if (is.fail()) return is;
        }
    }

    // Get the closing bracket if there was an opening one. If we don't
    // see the expected character we'll set the fail bit in the istream.
    if (closeBracket != CHAR(0)) {
        CHAR closer; is >> closer; if (is.fail()) return is;
        if (closer != closeBracket) {
            is.unget(); if (is.fail()) return is;
            is.setstate( std::ios::failbit );
        }
    }

    return is;
}

} //namespace SimTK


#endif //SimTK_SIMMATRIX_SMALLMATRIX_ROW_H_
