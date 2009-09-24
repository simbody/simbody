#ifndef SimTK_SIMMATRIX_SMALLMATRIX_VEC_H_
#define SimTK_SIMMATRIX_SMALLMATRIX_VEC_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmatrix(tm)                       *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-9 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
 * Contributors: Peter Eastman                                                *
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
 * Declaration of class Vec<NROWS, ELEMENT_TYPE, STRIDE>.
 */

#include "SimTKcommon/internal/common.h"

namespace SimTK {


// The following functions are used internally by Vec.

namespace Impl {

// For those wimpy compilers that don't unroll short, constant-limit loops, Peter Eastman added these
// recursive template implementations of add and subtract.

template <class E1, int S1, class E2, int S2> void
conformingAdd(const Vec<1,E1,S1>& r1, const Vec<1,E2,S2>& r2, Vec<1,typename CNT<E1>::template Result<E2>::Add>& result) {
    result[0] = r1[0] + r2[0];
}
template <int N, class E1, int S1, class E2, int S2> void
conformingAdd(const Vec<N,E1,S1>& r1, const Vec<N,E2,S2>& r2, Vec<N,typename CNT<E1>::template Result<E2>::Add>& result) {
    conformingAdd(reinterpret_cast<const Vec<N-1,E1,S1>&>(r1), reinterpret_cast<const Vec<N-1,E2,S2>&>(r2), reinterpret_cast<Vec<N-1,typename CNT<E1>::template Result<E2>::Add>&>(result));
    result[N-1] = r1[N-1] + r2[N-1];
}
template <class E1, int S1, class E2, int S2> void
conformingSubtract(const Vec<1,E1,S1>& r1, const Vec<1,E2,S2>& r2, Vec<1,typename CNT<E1>::template Result<E2>::Sub>& result) {
    result[0] = r1[0] - r2[0];
}
template <int N, class E1, int S1, class E2, int S2> void
conformingSubtract(const Vec<N,E1,S1>& r1, const Vec<N,E2,S2>& r2, Vec<N,typename CNT<E1>::template Result<E2>::Sub>& result) {
    conformingSubtract(reinterpret_cast<const Vec<N-1,E1,S1>&>(r1), reinterpret_cast<const Vec<N-1,E2,S2>&>(r2), reinterpret_cast<Vec<N-1,typename CNT<E1>::template Result<E2>::Sub>&>(result));
    result[N-1] = r1[N-1] - r2[N-1];
}
template <class E1, int S1, class E2, int S2> void
copy(Vec<1,E1,S1>& r1, const Vec<1,E2,S2>& r2) {
    r1[0] = r2[0];
}
template <int N, class E1, int S1, class E2, int S2> void
copy(Vec<N,E1,S1>& r1, const Vec<N,E2,S2>& r2) {
    copy(reinterpret_cast<Vec<N-1,E1,S1>&>(r1), reinterpret_cast<const Vec<N-1,E2,S2>&>(r2));
    r1[N-1] = r2[N-1];
}

}

/// Generic Vec
template <int M, class ELT, int STRIDE>
class Vec {
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
        NCols               = 1,
        NPackedElements     = M,
        NActualElements     = M * STRIDE,   // includes trailing gap
        NActualScalars      = CNT<E>::NActualScalars * NActualElements,
        RowSpacing          = STRIDE,
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

    // These are reinterpretations of the current data, so have the
    // same packing (stride).
    typedef Vec<M,E,STRIDE>                 T;
    typedef Vec<M,ENeg,STRIDE>              TNeg;
    typedef Vec<M,EWithoutNegator,STRIDE>   TWithoutNegator;

    typedef Vec<M,EReal,STRIDE*CNT<E>::RealStrideFactor>         
                                            TReal;
    typedef Vec<M,EImag,STRIDE*CNT<E>::RealStrideFactor>         
                                            TImag;
    typedef Vec<M,EComplex,STRIDE>          TComplex;
    typedef Row<M,EHerm,STRIDE>             THerm;
    typedef Row<M,E,STRIDE>                 TPosTrans;
    typedef E                               TElement;
    typedef E                               TRow;
    typedef Vec                             TCol;

    // These are the results of calculations, so are returned in new, packed
    // memory. Be sure to refer to element types here which are also packed.
    typedef Vec<M,ESqrt,1>                  TSqrt;      // Note stride
    typedef Vec<M,EAbs,1>                   TAbs;       // Note stride
    typedef Vec<M,EStandard,1>              TStandard;
    typedef Row<M,EInvert,1>                TInvert;
    typedef Vec<M,ENormalize,1>             TNormalize;

    typedef ESqHermT                        TSqHermT;   // result of self dot product
    typedef SymMat<M,ESqTHerm>              TSqTHerm;   // result of self outer product

    // These recurse right down to the underlying scalar type no matter how
    // deep the elements are.
    typedef EScalar                         Scalar;
    typedef EULessScalar                    ULessScalar;
    typedef ENumber                         Number;
    typedef EStdNumber                      StdNumber;
    typedef EPrecision                      Precision;
    typedef EScalarNormSq                   ScalarNormSq;

    int size()   const  { return M; }
    int nrow()   const  { return M; }
    int ncol()   const  { return 1; }


    // Scalar norm square is sum( conjugate squares of all scalars )
    ScalarNormSq scalarNormSqr() const { 
        ScalarNormSq sum(0);
        for(int i=0;i<M;++i) sum += CNT<E>::scalarNormSqr(d[i*STRIDE]);
        return sum;
    }

    // sqrt() is elementwise square root; that is, the return value has the same
    // dimension as this Vec but with each element replaced by whatever it thinks
    // its square root is.
    TSqrt sqrt() const {
        TSqrt vsqrt;
        for(int i=0;i<M;++i) vsqrt[i] = CNT<E>::sqrt(d[i*STRIDE]);
        return vsqrt;
    }

    // abs() is elementwise absolute value; that is, the return value has the same
    // dimension as this Vec but with each element replaced by whatever it thinks
    // its absolute value is.
    TAbs abs() const {
        TAbs vabs;
        for(int i=0;i<M;++i) vabs[i] = CNT<E>::abs(d[i*STRIDE]);
        return vabs;
    }

    TStandard standardize() const {
        TStandard vstd;
        for(int i=0;i<M;++i) vstd[i] = CNT<E>::standardize(d[i*STRIDE]);
        return vstd;
    }

    // Sum just adds up all the elements, getting rid of negators and
    // conjugates in the process.
    EStandard sum() const {
        E sum(0);
        for (int i=0;i<M;++i) sum += d[i*STRIDE];
        return CNT<E>::standardize(sum);
    }


    // This gives the resulting vector type when (v[i] op P) is applied to each element of v.
    // It is a vector of length M, stride 1, and element types which are the regular
    // composite result of E op P. Typically P is a scalar type but it doesn't have to be.
    template <class P> struct EltResult { 
        typedef Vec<M, typename CNT<E>::template Result<P>::Mul, 1> Mul;
        typedef Vec<M, typename CNT<E>::template Result<P>::Dvd, 1> Dvd;
        typedef Vec<M, typename CNT<E>::template Result<P>::Add, 1> Add;
        typedef Vec<M, typename CNT<E>::template Result<P>::Sub, 1> Sub;
    };

    // This is the composite result for v op P where P is some kind of appropriately shaped
    // non-scalar type.
    template <class P> struct Result { 
        typedef MulCNTs<M,1,ArgDepth,Vec,ColSpacing,RowSpacing,
            CNT<P>::NRows, CNT<P>::NCols, CNT<P>::ArgDepth,
            P, CNT<P>::ColSpacing, CNT<P>::RowSpacing> MulOp;
        typedef typename MulOp::Type Mul;

        typedef MulCNTsNonConforming<M,1,ArgDepth,Vec,ColSpacing,RowSpacing,
            CNT<P>::NRows, CNT<P>::NCols, CNT<P>::ArgDepth,
            P, CNT<P>::ColSpacing, CNT<P>::RowSpacing> MulOpNonConforming;
        typedef typename MulOpNonConforming::Type MulNon;

        typedef DvdCNTs<M,1,ArgDepth,Vec,ColSpacing,RowSpacing,
            CNT<P>::NRows, CNT<P>::NCols, CNT<P>::ArgDepth,
            P, CNT<P>::ColSpacing, CNT<P>::RowSpacing> DvdOp;
        typedef typename DvdOp::Type Dvd;

        typedef AddCNTs<M,1,ArgDepth,Vec,ColSpacing,RowSpacing,
            CNT<P>::NRows, CNT<P>::NCols, CNT<P>::ArgDepth,
            P, CNT<P>::ColSpacing, CNT<P>::RowSpacing> AddOp;
        typedef typename AddOp::Type Add;

        typedef SubCNTs<M,1,ArgDepth,Vec,ColSpacing,RowSpacing,
            CNT<P>::NRows, CNT<P>::NCols, CNT<P>::ArgDepth,
            P, CNT<P>::ColSpacing, CNT<P>::RowSpacing> SubOp;
        typedef typename SubOp::Type Sub;
    };

    // Shape-preserving element substitution (always packed)
    template <class P> struct Substitute {
        typedef Vec<M,P> Type;
    };

    // Default construction initializes to NaN when debugging but
    // is left uninitialized otherwise.
	Vec(){ 
    #ifndef NDEBUG
        setToNaN();
    #endif
    }

    // It's important not to use the default copy constructor or copy
    // assignment because the compiler doesn't understand that we may
    // have noncontiguous storage and will try to copy the whole array.
    Vec(const Vec& src) {
        Impl::copy(*this, src);
    }
    Vec& operator=(const Vec& src) {    // no harm if src and 'this' are the same
        Impl::copy(*this, src);
        return *this;
    }

    // We want an implicit conversion from a Vec of the same length
    // and element type but with a different stride.
    template <int SS> Vec(const Vec<M,E,SS>& src) {
        Impl::copy(*this, src);
    }

    // We want an implicit conversion from a Vec of the same length
    // and *negated* element type (possibly with a different stride).

    template <int SS> Vec(const Vec<M,ENeg,SS>& src) {
        Impl::copy(*this, src);
    }

    // Construct a Vec from a Vec of the same length, with any
    // stride. Works as long as the element types are compatible.
    template <class EE, int SS> explicit Vec(const Vec<M,EE,SS>& vv) {
        Impl::copy(*this, vv);
    }

    // Construction using an element assigns to each element.
    explicit Vec(const ELT& e)
      { for (int i=0;i<M;++i) d[i*STRIDE]=e; }

    // A bevy of constructors for Vecs up to length 6.
    Vec(const E& e0,const E& e1)
      { assert(M==2);(*this)[0]=e0;(*this)[1]=e1; }
    Vec(const E& e0,const E& e1,const E& e2)
      { assert(M==3);(*this)[0]=e0;(*this)[1]=e1;(*this)[2]=e2; }
    Vec(const E& e0,const E& e1,const E& e2,const E& e3)
      { assert(M==4);(*this)[0]=e0;(*this)[1]=e1;(*this)[2]=e2;(*this)[3]=e3; }
    Vec(const E& e0,const E& e1,const E& e2,const E& e3,const E& e4)
      { assert(M==5);(*this)[0]=e0;(*this)[1]=e1;(*this)[2]=e2;
        (*this)[3]=e3;(*this)[4]=e4; }
    Vec(const E& e0,const E& e1,const E& e2,const E& e3,const E& e4,const E& e5)
      { assert(M==6);(*this)[0]=e0;(*this)[1]=e1;(*this)[2]=e2;
        (*this)[3]=e3;(*this)[4]=e4;(*this)[5]=e5; }
    Vec(const E& e0,const E& e1,const E& e2,const E& e3,const E& e4,const E& e5, const E& e6)
      { assert(M==7);(*this)[0]=e0;(*this)[1]=e1;(*this)[2]=e2;
        (*this)[3]=e3;(*this)[4]=e4;(*this)[5]=e5;(*this)[6]=e6; }
    Vec(const E& e0,const E& e1,const E& e2,const E& e3,const E& e4,const E& e5, const E& e6, const E& e7)
      { assert(M==8);(*this)[0]=e0;(*this)[1]=e1;(*this)[2]=e2;
        (*this)[3]=e3;(*this)[4]=e4;(*this)[5]=e5;(*this)[6]=e6;(*this)[7]=e7; }
    Vec(const E& e0,const E& e1,const E& e2,const E& e3,const E& e4,const E& e5, const E& e6, const E& e7, const E& e8)
      { assert(M==9);(*this)[0]=e0;(*this)[1]=e1;(*this)[2]=e2;
        (*this)[3]=e3;(*this)[4]=e4;(*this)[5]=e5;(*this)[6]=e6;(*this)[7]=e7;(*this)[8]=e8; }

    // Construction or assignment from a pointer to anything assumes we're pointing
    // at an element list of the right length.
    template <class EE> explicit Vec(const EE* p)
      { assert(p); for(int i=0;i<M;++i) d[i*STRIDE]=p[i]; }
    template <class EE> Vec& operator=(const EE* p)
      { assert(p); for(int i=0;i<M;++i) d[i*STRIDE]=p[i]; return *this; }

    // Conforming assignment ops.
    template <class EE, int SS> Vec& operator=(const Vec<M,EE,SS>& vv) {
        Impl::copy(*this, vv);
        return *this;
    }
    template <class EE, int SS> Vec& operator+=(const Vec<M,EE,SS>& r)
      { for(int i=0;i<M;++i) d[i*STRIDE] += r[i]; return *this; }
    template <class EE, int SS> Vec& operator+=(const Vec<M,negator<EE>,SS>& r)
      { for(int i=0;i<M;++i) d[i*STRIDE] -= -(r[i]); return *this; }
    template <class EE, int SS> Vec& operator-=(const Vec<M,EE,SS>& r)
      { for(int i=0;i<M;++i) d[i*STRIDE] -= r[i]; return *this; }
    template <class EE, int SS> Vec& operator-=(const Vec<M,negator<EE>,SS>& r)
      { for(int i=0;i<M;++i) d[i*STRIDE] += -(r[i]); return *this; }

    // Conforming binary ops with 'this' on left, producing new packed result.
    // Cases: v=v+v, v=v-v, m=v*r
    template <class EE, int SS> Vec<M,typename CNT<E>::template Result<EE>::Add>
    conformingAdd(const Vec<M,EE,SS>& r) const {
        Vec<M,typename CNT<E>::template Result<EE>::Add> result;
        Impl::conformingAdd(*this, r, result);
        return result;
    }
    template <class EE, int SS> Vec<M,typename CNT<E>::template Result<EE>::Sub>
    conformingSubtract(const Vec<M,EE,SS>& r) const {
        Vec<M,typename CNT<E>::template Result<EE>::Sub> result;
        Impl::conformingSubtract(*this, r, result);
        return result;
    }
    template <class EE, int SS> Mat<M,M,typename CNT<E>::template Result<EE>::Mul>
    conformingMultiply(const Row<M,EE,SS>& r) const {
        Mat<M,M,typename CNT<E>::template Result<EE>::Mul> result;
        for (int j=0;j<M;++j) result(j) = scalarMultiply(r(j));
        return result;
    }

	const E& operator[](int i) const { assert(0 <= i && i < M); return d[i*STRIDE]; }
	E&       operator[](int i)	     { assert(0 <= i && i < M); return d[i*STRIDE]; }
    const E& operator()(int i) const { return (*this)[i]; }
	E&       operator()(int i)	     { return (*this)[i]; }

    ScalarNormSq normSqr() const { return scalarNormSqr(); }
    typename CNT<ScalarNormSq>::TSqrt 
        norm() const { return CNT<ScalarNormSq>::sqrt(scalarNormSqr()); }

    // If the elements of this Vec are scalars, the result is what you get by
    // dividing each element by the norm() calculated above. If the elements are
    // *not* scalars, then the elements are *separately* normalized. That means
    // you will get a different answer from Vec<2,Vec3>::normalize() than you
    // would from a Vec<6>::normalize() containing the same scalars.
    //
    // Normalize returns a vector of the same dimension but in new, packed storage
    // and with a return type that does not include negator<> even if the original
    // Vec<> does, because we can eliminate the negation here almost for free.
    // But we can't standardize (change conjugate to complex) for free, so we'll retain
    // conjugates if there are any.
    TNormalize normalize() const {
        if (CNT<E>::IsScalar) {
            return castAwayNegatorIfAny() / (SignInterpretation*norm());
        } else {
            TNormalize elementwiseNormalized;
            for (int i=0; i<M; ++i) 
                elementwiseNormalized[i] = CNT<E>::normalize((*this)[i]);
            return elementwiseNormalized;
        }
    }

    TInvert invert() const {assert(false); return TInvert();} // TODO default inversion

    const Vec&   operator+() const { return *this; }
    const TNeg&  operator-() const { return negate(); }
    TNeg&        operator-()       { return updNegate(); }
    const THerm& operator~() const { return transpose(); }
    THerm&       operator~()       { return updTranspose(); }

    const TNeg&  negate() const { return *reinterpret_cast<const TNeg*>(this); }
    TNeg&        updNegate()    { return *reinterpret_cast<      TNeg*>(this); }

    const THerm& transpose()    const { return *reinterpret_cast<const THerm*>(this); }
    THerm&       updTranspose()       { return *reinterpret_cast<      THerm*>(this); }

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
    // 'FromLeft' appears in the name. The result is a packed Vec<M> but the element type
    // may change. These are mostly used to implement global operators.
    // We call these "scalar" operators but actually the "scalar" can be a composite type.

    //TODO: consider converting 'e' to Standard Numbers as precalculation and changing
    // return type appropriately.
    template <class EE> Vec<M, typename CNT<E>::template Result<EE>::Mul>
    scalarMultiply(const EE& e) const {
        Vec<M, typename CNT<E>::template Result<EE>::Mul> result;
        for (int i=0; i<M; ++i) result[i] = (*this)[i] * e;
        return result;
    }
    template <class EE> Vec<M, typename CNT<EE>::template Result<E>::Mul>
    scalarMultiplyFromLeft(const EE& e) const {
        Vec<M, typename CNT<EE>::template Result<E>::Mul> result;
        for (int i=0; i<M; ++i) result[i] = e * (*this)[i];
        return result;
    }

    // TODO: should precalculate and store 1/e, while converting to Standard Numbers. Note
    // that return type should change appropriately.
    template <class EE> Vec<M, typename CNT<E>::template Result<EE>::Dvd>
    scalarDivide(const EE& e) const {
        Vec<M, typename CNT<E>::template Result<EE>::Dvd> result;
        for (int i=0; i<M; ++i) result[i] = (*this)[i] / e;
        return result;
    }
    template <class EE> Vec<M, typename CNT<EE>::template Result<E>::Dvd>
    scalarDivideFromLeft(const EE& e) const {
        Vec<M, typename CNT<EE>::template Result<E>::Dvd> result;
        for (int i=0; i<M; ++i) result[i] = e / (*this)[i];
        return result;
    }

    template <class EE> Vec<M, typename CNT<E>::template Result<EE>::Add>
    scalarAdd(const EE& e) const {
        Vec<M, typename CNT<E>::template Result<EE>::Add> result;
        for (int i=0; i<M; ++i) result[i] = (*this)[i] + e;
        return result;
    }
    // Add is commutative, so no 'FromLeft'.

    template <class EE> Vec<M, typename CNT<E>::template Result<EE>::Sub>
    scalarSubtract(const EE& e) const {
        Vec<M, typename CNT<E>::template Result<EE>::Sub> result;
        for (int i=0; i<M; ++i) result[i] = (*this)[i] - e;
        return result;
    }
    template <class EE> Vec<M, typename CNT<EE>::template Result<E>::Sub>
    scalarSubtractFromLeft(const EE& e) const {
        Vec<M, typename CNT<EE>::template Result<E>::Sub> result;
        for (int i=0; i<M; ++i) result[i] = e - (*this)[i];
        return result;
    }

    // Generic assignments for any element type not listed explicitly, including scalars.
    // These are done repeatedly for each element and only work if the operation can
    // be performed leaving the original element type.
    template <class EE> Vec& operator =(const EE& e) {return scalarEq(e);}
    template <class EE> Vec& operator+=(const EE& e) {return scalarPlusEq(e);}
    template <class EE> Vec& operator-=(const EE& e) {return scalarMinusEq(e);}
    template <class EE> Vec& operator*=(const EE& e) {return scalarTimesEq(e);}
    template <class EE> Vec& operator/=(const EE& e) {return scalarDivideEq(e);}

    // Generalized element assignment & computed assignment methods. These will work
    // for any assignment-compatible element, not just scalars.
    template <class EE> Vec& scalarEq(const EE& ee)
      { for(int i=0;i<M;++i) d[i*STRIDE] = ee; return *this; }

    template <class EE> Vec& scalarPlusEq(const EE& ee)
      { for(int i=0;i<M;++i) d[i*STRIDE] += ee; return *this; }

    template <class EE> Vec& scalarMinusEq(const EE& ee)
      { for(int i=0;i<M;++i) d[i*STRIDE] -= ee; return *this; }
    template <class EE> Vec& scalarMinusEqFromLeft(const EE& ee)
      { for(int i=0;i<M;++i) d[i*STRIDE] = ee - d[i*STRIDE]; return *this; }

    template <class EE> Vec& scalarTimesEq(const EE& ee)
      { for(int i=0;i<M;++i) d[i*STRIDE] *= ee; return *this; }
    template <class EE> Vec& scalarTimesEqFromLeft(const EE& ee)
      { for(int i=0;i<M;++i) d[i*STRIDE] = ee * d[i*STRIDE]; return *this; }

    template <class EE> Vec& scalarDivideEq(const EE& ee)
      { for(int i=0;i<M;++i) d[i*STRIDE] /= ee; return *this; }
    template <class EE> Vec& scalarDivideEqFromLeft(const EE& ee)
      { for(int i=0;i<M;++i) d[i*STRIDE] = ee / d[i*STRIDE]; return *this; }

    void setToNaN() {
        (*this) = CNT<ELT>::getNaN();
    }

    void setToZero() {
        (*this) = ELT(0);
    }

    // Extract a sub-Vec with size known at compile time. These have to be
    // called with explicit template arguments, e.g. getSubVec<3>(i).
    template <int MM>
    const Vec<MM,ELT,STRIDE>& getSubVec(int i) const {
        assert(0 <= i && i + MM <= M);
        return Vec<MM,ELT,STRIDE>::getAs(&(*this)[i]);
    }
    template <int MM>
    Vec<MM,ELT,STRIDE>& updSubVec(int i) {
        assert(0 <= i && i + MM <= M);
        return Vec<MM,ELT,STRIDE>::updAs(&(*this)[i]);
    }

    // Return a vector one smaller than this one by dropping the element
    // at the indicated position p. The result is packed but has same
    // element type as this one.
    Vec<M-1,ELT,1> drop1(int p) const {
        assert(0 <= p && p < M);
        Vec<M-1,ELT,1> out;
        int nxt=0;
        for (int i=0; i<M-1; ++i, ++nxt) {
            if (nxt==p) ++nxt;  // skip the loser
            out[i] = (*this)[nxt];
        }
        return out;
    }

    // Return a vector one larger than this one by adding an element
    // to the end. The result is packed but has same element type as
    // this one. Works for any assignment compatible element.
    template <class EE> Vec<M+1,ELT,1> append1(const EE& v) const {
        Vec<M+1,ELT,1> out;
        Vec<M,ELT,1>::updAs(&out[0]) = (*this);
        out[M] = v;
        return out;
    }


    // Return a vector one larger than this one by inserting an element
    // *before* the indicated one. The result is packed but has same element type as
    // this one. Works for any assignment compatible element. The index
    // can be one greater than normally allowed in which case the element
    // is appended.
    template <class EE> Vec<M+1,ELT,1> insert1(int p, const EE& v) const {
        assert(0 <= p && p <= M);
        if (p==M) return append1(v);
        Vec<M+1,ELT,1> out;
        int nxt=0;
        for (int i=0; i<M; ++i, ++nxt) {
            if (i==p) out[nxt++] = v;
            out[nxt] = (*this)[i];
        }
        return out;
    }
            
    // These assume we are given a pointer to d[0] of a Vec<M,E,S> like this one.
    static const Vec& getAs(const ELT* p)  {return *reinterpret_cast<const Vec*>(p);}
    static Vec&       updAs(ELT* p)        {return *reinterpret_cast<Vec*>(p);}

    // Extract a subvector from a longer one. Element type and stride must match.
    template <int MM>
    static const Vec& getSubVec(const Vec<MM,ELT,STRIDE>& v, int i) {
        assert(0 <= i && i + M <= MM);
        return getAs(&v[i]);
    }
    template <int MM>
    static Vec& updSubVec(Vec<MM,ELT,STRIDE>& v, int i) {
        assert(0 <= i && i + M <= MM);
        return updAs(&v[i]);
    }

    static Vec<M,ELT,1> getNaN() { return Vec<M,ELT,1>(CNT<ELT>::getNaN()); }

    /// Return true if any element of this Vec contains a NaN anywhere.
    bool isNaN() const {
        for (int i=0; i<M; ++i)
            if (CNT<ELT>::isNaN((*this)[i]))
                return true;
        return false;
    }

    /// Return true if any element of this Vec contains a +Inf
    /// or -Inf somewhere but no element contains a NaN anywhere.
    bool isInf() const {
        bool seenInf = false;
        for (int i=0; i<M; ++i) {
            const ELT& e = (*this)[i];
            if (!CNT<ELT>::isFinite(e)) {
                if (!CNT<ELT>::isInf(e)) 
                    return false; // something bad was found
                seenInf = true; 
            }
        }
        return seenInf;
    }

    /// Return true if no element contains an Infinity or a NaN.
    bool isFinite() const {
        for (int i=0; i<M; ++i)
            if (!CNT<ELT>::isFinite((*this)[i]))
                return false;
        return true;
    }

    /// For approximate comparisions, the default tolerance to use for a vector is
    /// the same as its elements' default tolerance.
    static double getDefaultTolerance() {return CNT<ELT>::getDefaultTolerance();}

    /// %Test whether this vector is numerically equal to some other vector with
    /// the same shape, using a specified tolerance.
    template <class E2, int RS2>
    bool isNumericallyEqual(const Vec<M,E2,RS2>& v, double tol) const {
        for (int i=0; i<M; ++i)
            if (!CNT<ELT>::isNumericallyEqual((*this)[i], v[i], tol))
                return false;
        return true;
    }

    /// %Test whether this vector is numerically equal to some other vector with
    /// the same shape, using a default tolerance which is the looser of the
    /// default tolerances of the two objects being compared.
    template <class E2, int RS2>
    bool isNumericallyEqual(const Vec<M,E2,RS2>& v) const {
        const double tol = std::max(getDefaultTolerance(),v.getDefaultTolerance());
        return isNumericallyEqual(v, tol);
    }

    /// %Test whether every element of this vector is numerically equal to the given
    /// element, using either a specified tolerance or the vector's 
    /// default tolerance (which is always the same or looser than the default
    /// tolerance for one of its elements).
    bool isNumericallyEqual
       (const ELT& e,
        double     tol = getDefaultTolerance()) const 
    {
        for (int i=0; i<M; ++i)
            if (!CNT<ELT>::isNumericallyEqual((*this)[i], e, tol))
                return false;
        return true;
    }
private:
	ELT d[NActualElements];    // data
};

/////////////////////////////////////////////
// Global operators involving two vectors. //
//   v+v, v-v, v==v, v!=v                  //
/////////////////////////////////////////////

// v3 = v1 + v2 where all v's have the same length M. 
template <int M, class E1, int S1, class E2, int S2> inline
typename Vec<M,E1,S1>::template Result< Vec<M,E2,S2> >::Add
operator+(const Vec<M,E1,S1>& l, const Vec<M,E2,S2>& r) {
    return Vec<M,E1,S1>::template Result< Vec<M,E2,S2> >
        ::AddOp::perform(l,r);
}

// v3 = v1 - v2, similar to +
template <int M, class E1, int S1, class E2, int S2> inline
typename Vec<M,E1,S1>::template Result< Vec<M,E2,S2> >::Sub
operator-(const Vec<M,E1,S1>& l, const Vec<M,E2,S2>& r) { 
    return Vec<M,E1,S1>::template Result< Vec<M,E2,S2> >
        ::SubOp::perform(l,r);
}

/// bool = v1[i] == v2[i], for all elements i
template <int M, class E1, int S1, class E2, int S2> inline bool
operator==(const Vec<M,E1,S1>& l, const Vec<M,E2,S2>& r) 
{   for (int i=0; i < M; ++i) if (l[i] != r[i]) return false;
    return true; }
/// bool = v1[i] != v2[i], for any element i
template <int M, class E1, int S1, class E2, int S2> inline bool
operator!=(const Vec<M,E1,S1>& l, const Vec<M,E2,S2>& r) {return !(l==r);} 

/// bool = v[i] == e, for all elements v[i] and element e
template <int M, class E1, int S1, class E2> inline bool
operator==(const Vec<M,E1,S1>& v, const E2& e) 
{   for (int i=0; i < M; ++i) if (v[i] != e) return false;
    return true; }
/// bool = v[i] != e, for any element v[i] and element e
template <int M, class E1, int S1, class E2> inline bool
operator!=(const Vec<M,E1,S1>& v, const E2& e) {return !(v==e);} 


/// bool = v1[i] < v2[i], for all elements i
template <int M, class E1, int S1, class E2, int S2> inline bool
operator<(const Vec<M,E1,S1>& l, const Vec<M,E2,S2>& r) 
{   for (int i=0; i < M; ++i) if (l[i] >= r[i]) return false;
    return true; }
/// bool = v[i] < e, for all elements v[i] and element e
template <int M, class E1, int S1, class E2> inline bool
operator<(const Vec<M,E1,S1>& v, const E2& e) 
{   for (int i=0; i < M; ++i) if (v[i] >= e) return false;
    return true; }

/// bool = v1[i] > v2[i], for all elements i
template <int M, class E1, int S1, class E2, int S2> inline bool
operator>(const Vec<M,E1,S1>& l, const Vec<M,E2,S2>& r) 
{   for (int i=0; i < M; ++i) if (l[i] <= r[i]) return false;
    return true; }
/// bool = v[i] > e, for all elements v[i] and element e
template <int M, class E1, int S1, class E2> inline bool
operator>(const Vec<M,E1,S1>& v, const E2& e) 
{   for (int i=0; i < M; ++i) if (v[i] <= e) return false;
    return true; }

/// bool = v1[i] <= v2[i], for all elements i.
/// This is not the same as !(v1>v2).
template <int M, class E1, int S1, class E2, int S2> inline bool
operator<=(const Vec<M,E1,S1>& l, const Vec<M,E2,S2>& r) 
{   for (int i=0; i < M; ++i) if (l[i] > r[i]) return false;
    return true; }
/// bool = v[i] <= e, for all elements v[i] and element e.
/// This is not the same as !(v1>e).
template <int M, class E1, int S1, class E2> inline bool
operator<=(const Vec<M,E1,S1>& v, const E2& e) 
{   for (int i=0; i < M; ++i) if (v[i] > e) return false;
    return true; }

/// bool = v1[i] >= v2[i], for all elements i
/// This is not the same as !(v1<v2).
template <int M, class E1, int S1, class E2, int S2> inline bool
operator>=(const Vec<M,E1,S1>& l, const Vec<M,E2,S2>& r) 
{   for (int i=0; i < M; ++i) if (l[i] < r[i]) return false;
    return true; }
/// bool = v[i] >= e, for all elements v[i] and element e.
/// This is not the same as !(v1<e).
template <int M, class E1, int S1, class E2> inline bool
operator>=(const Vec<M,E1,S1>& v, const E2& e) 
{   for (int i=0; i < M; ++i) if (v[i] < e) return false;
    return true; }

///////////////////////////////////////////////////////
// Global operators involving a vector and a scalar. //
///////////////////////////////////////////////////////

// I haven't been able to figure out a nice way to templatize for the
// built-in reals without introducing a lot of unwanted type matches
// as well. So we'll just grind them out explicitly here.

// SCALAR MULTIPLY

// v = v*real, real*v 
template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<float>::Mul
operator*(const Vec<M,E,S>& l, const float& r)
  { return Vec<M,E,S>::template Result<float>::MulOp::perform(l,r); }
template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<float>::Mul
operator*(const float& l, const Vec<M,E,S>& r) {return r*l;}

template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<double>::Mul
operator*(const Vec<M,E,S>& l, const double& r)
  { return Vec<M,E,S>::template Result<double>::MulOp::perform(l,r); }
template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<double>::Mul
operator*(const double& l, const Vec<M,E,S>& r) {return r*l;}

template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<long double>::Mul
operator*(const Vec<M,E,S>& l, const long double& r)
  { return Vec<M,E,S>::template Result<long double>::MulOp::perform(l,r); }
template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<long double>::Mul
operator*(const long double& l, const Vec<M,E,S>& r) {return r*l;}

// v = v*int, int*v -- just convert int to v's precision float
template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<typename CNT<E>::Precision>::Mul
operator*(const Vec<M,E,S>& l, int r) {return l * (typename CNT<E>::Precision)r;}
template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<typename CNT<E>::Precision>::Mul
operator*(int l, const Vec<M,E,S>& r) {return r * (typename CNT<E>::Precision)l;}

// Complex, conjugate, and negator are all easy to templatize.

// v = v*complex, complex*v
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<std::complex<R> >::Mul
operator*(const Vec<M,E,S>& l, const std::complex<R>& r)
  { return Vec<M,E,S>::template Result<std::complex<R> >::MulOp::perform(l,r); }
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<std::complex<R> >::Mul
operator*(const std::complex<R>& l, const Vec<M,E,S>& r) {return r*l;}

// v = v*conjugate, conjugate*v (convert conjugate->complex)
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<std::complex<R> >::Mul
operator*(const Vec<M,E,S>& l, const conjugate<R>& r) {return l*(std::complex<R>)r;}
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<std::complex<R> >::Mul
operator*(const conjugate<R>& l, const Vec<M,E,S>& r) {return r*(std::complex<R>)l;}

// v = v*negator, negator*v: convert negator to standard number
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<typename negator<R>::StdNumber>::Mul
operator*(const Vec<M,E,S>& l, const negator<R>& r) {return l * (typename negator<R>::StdNumber)(R)r;}
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<typename negator<R>::StdNumber>::Mul
operator*(const negator<R>& l, const Vec<M,E,S>& r) {return r * (typename negator<R>::StdNumber)(R)l;}


// SCALAR DIVIDE. This is a scalar operation when the scalar is on the right,
// but when it is on the left it means scalar * pseudoInverse(vec), which is a row.

// v = v/real, real/v 
template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<float>::Dvd
operator/(const Vec<M,E,S>& l, const float& r)
  { return Vec<M,E,S>::template Result<float>::DvdOp::perform(l,r); }
template <int M, class E, int S> inline
typename CNT<float>::template Result<Vec<M,E,S> >::Dvd
operator/(const float& l, const Vec<M,E,S>& r)
  { return CNT<float>::template Result<Vec<M,E,S> >::DvdOp::perform(l,r); }

template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<double>::Dvd
operator/(const Vec<M,E,S>& l, const double& r)
  { return Vec<M,E,S>::template Result<double>::DvdOp::perform(l,r); }
template <int M, class E, int S> inline
typename CNT<double>::template Result<Vec<M,E,S> >::Dvd
operator/(const double& l, const Vec<M,E,S>& r)
  { return CNT<double>::template Result<Vec<M,E,S> >::DvdOp::perform(l,r); }

template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<long double>::Dvd
operator/(const Vec<M,E,S>& l, const long double& r)
  { return Vec<M,E,S>::template Result<long double>::DvdOp::perform(l,r); }
template <int M, class E, int S> inline
typename CNT<long double>::template Result<Vec<M,E,S> >::Dvd
operator/(const long double& l, const Vec<M,E,S>& r)
  { return CNT<long double>::template Result<Vec<M,E,S> >::DvdOp::perform(l,r); }

// v = v/int, int/v -- just convert int to v's precision float
template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<typename CNT<E>::Precision>::Dvd
operator/(const Vec<M,E,S>& l, int r) {return l / (typename CNT<E>::Precision)r;}
template <int M, class E, int S> inline
typename CNT<typename CNT<E>::Precision>::template Result<Vec<M,E,S> >::Dvd
operator/(int l, const Vec<M,E,S>& r) {return (typename CNT<E>::Precision)l / r;}


// Complex, conjugate, and negator are all easy to templatize.

// v = v/complex, complex/v
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<std::complex<R> >::Dvd
operator/(const Vec<M,E,S>& l, const std::complex<R>& r)
  { return Vec<M,E,S>::template Result<std::complex<R> >::DvdOp::perform(l,r); }
template <int M, class E, int S, class R> inline
typename CNT<std::complex<R> >::template Result<Vec<M,E,S> >::Dvd
operator/(const std::complex<R>& l, const Vec<M,E,S>& r)
  { return CNT<std::complex<R> >::template Result<Vec<M,E,S> >::DvdOp::perform(l,r); }

// v = v/conjugate, conjugate/v (convert conjugate->complex)
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<std::complex<R> >::Dvd
operator/(const Vec<M,E,S>& l, const conjugate<R>& r) {return l/(std::complex<R>)r;}
template <int M, class E, int S, class R> inline
typename CNT<std::complex<R> >::template Result<Vec<M,E,S> >::Dvd
operator/(const conjugate<R>& l, const Vec<M,E,S>& r) {return (std::complex<R>)l/r;}

// v = v/negator, negator/v: convert negator to number
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<typename negator<R>::StdNumber>::Dvd
operator/(const Vec<M,E,S>& l, const negator<R>& r) {return l/(typename negator<R>::StdNumber)(R)r;}
template <int M, class E, int S, class R> inline
typename CNT<R>::template Result<Vec<M,E,S> >::Dvd
operator/(const negator<R>& l, const Vec<M,E,S>& r) {return (typename negator<R>::StdNumber)(R)l/r;}


// Add and subtract are odd as scalar ops. They behave as though the
// scalar stands for a vector each of whose elements is that scalar,
// and then a normal vector add or subtract is done.

// SCALAR ADD

// v = v+real, real+v 
template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<float>::Add
operator+(const Vec<M,E,S>& l, const float& r)
  { return Vec<M,E,S>::template Result<float>::AddOp::perform(l,r); }
template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<float>::Add
operator+(const float& l, const Vec<M,E,S>& r) {return r+l;}

template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<double>::Add
operator+(const Vec<M,E,S>& l, const double& r)
  { return Vec<M,E,S>::template Result<double>::AddOp::perform(l,r); }
template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<double>::Add
operator+(const double& l, const Vec<M,E,S>& r) {return r+l;}

template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<long double>::Add
operator+(const Vec<M,E,S>& l, const long double& r)
  { return Vec<M,E,S>::template Result<long double>::AddOp::perform(l,r); }
template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<long double>::Add
operator+(const long double& l, const Vec<M,E,S>& r) {return r+l;}

// v = v+int, int+v -- just convert int to v's precision float
template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<typename CNT<E>::Precision>::Add
operator+(const Vec<M,E,S>& l, int r) {return l + (typename CNT<E>::Precision)r;}
template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<typename CNT<E>::Precision>::Add
operator+(int l, const Vec<M,E,S>& r) {return r + (typename CNT<E>::Precision)l;}

// Complex, conjugate, and negator are all easy to templatize.

// v = v+complex, complex+v
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<std::complex<R> >::Add
operator+(const Vec<M,E,S>& l, const std::complex<R>& r)
  { return Vec<M,E,S>::template Result<std::complex<R> >::AddOp::perform(l,r); }
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<std::complex<R> >::Add
operator+(const std::complex<R>& l, const Vec<M,E,S>& r) {return r+l;}

// v = v+conjugate, conjugate+v (convert conjugate->complex)
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<std::complex<R> >::Add
operator+(const Vec<M,E,S>& l, const conjugate<R>& r) {return l+(std::complex<R>)r;}
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<std::complex<R> >::Add
operator+(const conjugate<R>& l, const Vec<M,E,S>& r) {return r+(std::complex<R>)l;}

// v = v+negator, negator+v: convert negator to standard number
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<typename negator<R>::StdNumber>::Add
operator+(const Vec<M,E,S>& l, const negator<R>& r) {return l + (typename negator<R>::StdNumber)(R)r;}
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<typename negator<R>::StdNumber>::Add
operator+(const negator<R>& l, const Vec<M,E,S>& r) {return r + (typename negator<R>::StdNumber)(R)l;}

// SCALAR SUBTRACT -- careful, not commutative.

// v = v-real, real-v 
template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<float>::Sub
operator-(const Vec<M,E,S>& l, const float& r)
  { return Vec<M,E,S>::template Result<float>::SubOp::perform(l,r); }
template <int M, class E, int S> inline
typename CNT<float>::template Result<Vec<M,E,S> >::Sub
operator-(const float& l, const Vec<M,E,S>& r)
  { return CNT<float>::template Result<Vec<M,E,S> >::SubOp::perform(l,r); }

template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<double>::Sub
operator-(const Vec<M,E,S>& l, const double& r)
  { return Vec<M,E,S>::template Result<double>::SubOp::perform(l,r); }
template <int M, class E, int S> inline
typename CNT<double>::template Result<Vec<M,E,S> >::Sub
operator-(const double& l, const Vec<M,E,S>& r)
  { return CNT<double>::template Result<Vec<M,E,S> >::SubOp::perform(l,r); }

template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<long double>::Sub
operator-(const Vec<M,E,S>& l, const long double& r)
  { return Vec<M,E,S>::template Result<long double>::SubOp::perform(l,r); }
template <int M, class E, int S> inline
typename CNT<long double>::template Result<Vec<M,E,S> >::Sub
operator-(const long double& l, const Vec<M,E,S>& r)
  { return CNT<long double>::template Result<Vec<M,E,S> >::SubOp::perform(l,r); }

// v = v-int, int-v // just convert int to v's precision float
template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<typename CNT<E>::Precision>::Sub
operator-(const Vec<M,E,S>& l, int r) {return l - (typename CNT<E>::Precision)r;}
template <int M, class E, int S> inline
typename CNT<typename CNT<E>::Precision>::template Result<Vec<M,E,S> >::Sub
operator-(int l, const Vec<M,E,S>& r) {return (typename CNT<E>::Precision)l - r;}


// Complex, conjugate, and negator are all easy to templatize.

// v = v-complex, complex-v
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<std::complex<R> >::Sub
operator-(const Vec<M,E,S>& l, const std::complex<R>& r)
  { return Vec<M,E,S>::template Result<std::complex<R> >::SubOp::perform(l,r); }
template <int M, class E, int S, class R> inline
typename CNT<std::complex<R> >::template Result<Vec<M,E,S> >::Sub
operator-(const std::complex<R>& l, const Vec<M,E,S>& r)
  { return CNT<std::complex<R> >::template Result<Vec<M,E,S> >::SubOp::perform(l,r); }

// v = v-conjugate, conjugate-v (convert conjugate->complex)
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<std::complex<R> >::Sub
operator-(const Vec<M,E,S>& l, const conjugate<R>& r) {return l-(std::complex<R>)r;}
template <int M, class E, int S, class R> inline
typename CNT<std::complex<R> >::template Result<Vec<M,E,S> >::Sub
operator-(const conjugate<R>& l, const Vec<M,E,S>& r) {return (std::complex<R>)l-r;}

// v = v-negator, negator-v: convert negator to standard number
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<typename negator<R>::StdNumber>::Sub
operator-(const Vec<M,E,S>& l, const negator<R>& r) {return l-(typename negator<R>::StdNumber)(R)r;}
template <int M, class E, int S, class R> inline
typename CNT<R>::template Result<Vec<M,E,S> >::Sub
operator-(const negator<R>& l, const Vec<M,E,S>& r) {return (typename negator<R>::StdNumber)(R)l-r;}

// Vec I/O
template <int M, class E, int S, class CHAR, class TRAITS> inline
std::basic_ostream<CHAR,TRAITS>&
operator<<(std::basic_ostream<CHAR,TRAITS>& o, const Vec<M,E,S>& v) {
    o << "~[" << v[0]; for(int i=1;i<M;++i) o<<','<<v[i]; o<<']'; return o;
}

template <int M, class E, int S, class CHAR, class TRAITS> inline
std::basic_istream<CHAR,TRAITS>&
operator>>(std::basic_istream<CHAR,TRAITS>& is, Vec<M,E,S>& v) {
    // TODO: not sure how to do Vec input yet
    assert(false);
    return is;
}

} //namespace SimTK


#endif //SimTK_SIMMATRIX_SMALLMATRIX_VEC_H_
