#ifndef SimTK_SIMMATRIX_SMALLMATRIX_MAT_H_
#define SimTK_SIMMATRIX_SMALLMATRIX_MAT_H_

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

/**@file
 * This file declares class Mat<NROWS, NCOLS, ELEMENT_TYPE, COL_SPACING, ROW_SPACING>.
 */

#include "SimTKcommon/internal/common.h"

namespace SimTK {

/// CS is total spacing between columns in memory (default M)
/// RS is total spacing between rows in memory (default 1) 
template <int M, int N, class ELT, int CS, int RS> class Mat {
    typedef ELT                         E;
    typedef typename CNT<E>::TNeg       ENeg;
    typedef typename CNT<E>::TAbs       EAbs;
    typedef typename CNT<E>::TStandard  EStandard;
    typedef typename CNT<E>::TReal      EReal;
    typedef typename CNT<E>::TImag      EImag;
    typedef typename CNT<E>::TComplex   EComplex;
    typedef typename CNT<E>::THerm      EHerm;
    typedef typename CNT<E>::TInvert    EInvert;
    typedef typename CNT<E>::TPosTrans  EPosTrans;
    typedef typename CNT<E>::TSqHermT   ESqHermT;
    typedef typename CNT<E>::TSqTHerm   ESqTHerm;

    typedef typename CNT<E>::Scalar     EScalar;
    typedef typename CNT<E>::Number     ENumber;
    typedef typename CNT<E>::StdNumber  EStdNumber;
    typedef typename CNT<E>::Precision  EPrecision;
    typedef typename CNT<E>::ScalarSq   EScalarSq;

public:

    enum {
        NRows               = M,
        NCols               = N,
        MinDim              = N < M ? N : M,
        MaxDim              = N > M ? N : M,
        RowSpacing          = RS,
        ColSpacing          = CS,
        NPackedElements     = M * N,
        NActualElements     = (N-1)*CS + (M-1)*RS + 1,
        NActualScalars      = CNT<E>::NActualScalars * NActualElements,
        ImagOffset          = NTraits<ENumber>::ImagOffset,
        RealStrideFactor    = 1, // composite types don't change size when
                                 // cast from complex to real or imaginary
        ArgDepth            = ((int)CNT<E>::ArgDepth < (int)MAX_RESOLVED_DEPTH 
                                ? CNT<E>::ArgDepth + 1 
                                : MAX_RESOLVED_DEPTH),
        IsScalar            = 0,
        IsNumber            = 0,
        IsStdNumber         = 0,
        IsPrecision         = 0
    };

    typedef Mat<M,N,E,CS,RS>            T;
    typedef Mat<M,N,ENeg,CS,RS>         TNeg;
    typedef Mat<M,N,EAbs,M,1>           TAbs;       // Note strides are packed
    typedef Mat<M,N,EStandard,M,1>      TStandard;
    typedef Mat<M,N,EReal,CS*CNT<E>::RealStrideFactor,RS*CNT<E>::RealStrideFactor>
                                        TReal;
    typedef Mat<M,N,EImag,CS*CNT<E>::RealStrideFactor,RS*CNT<E>::RealStrideFactor>
                                        TImag;
    typedef Mat<M,N,EComplex,CS,RS>     TComplex;
    typedef Mat<N,M,EHerm,RS,CS>        THerm;
    typedef Mat<N,M,EInvert,N,1>        TInvert;    // like THerm but packed
    typedef Mat<N,M,E,RS,CS>            TPosTrans;
    typedef SymMat<N,ESqHermT>          TSqHermT;   // ~Mat*Mat
    typedef SymMat<M,ESqTHerm>          TSqTHerm;   // Mat*~Mat
    typedef E                           TElement;
    typedef Row<N,E,CS>                 TRow;    
    typedef Vec<M,E,RS>                 TCol;

    typedef EScalar                     Scalar;
    typedef ENumber                     Number;
    typedef EStdNumber                  StdNumber;
    typedef EPrecision                  Precision;
    typedef EScalarSq                   ScalarSq;

    typedef THerm                       TransposeType; // TODO

    typedef Vec<MinDim,E,RS+CS>              TDiag;

    int size() const { return M*N; }
    int nrow() const { return M; }
    int ncol() const { return N; }

    // Scalar norm square is sum( squares of all scalars )
    ScalarSq scalarNormSqr() const { 
        ScalarSq sum(0);
        for(int j=0;j<N;++j) sum += CNT<TCol>::scalarNormSqr((*this)(j));
        return sum;
    }

    // abs() is elementwise absolute value; that is, the return value has the same
    // dimension as this Mat but with each element replaced by whatever it thinks
    // its absolute value is.
    TAbs abs() const { 
        TAbs mabs;
        for(int j=0;j<N;++j) mabs(j) = (*this)(j).abs();
        return mabs;
    }

    TStandard standardize() const {
        TStandard mstd;
        for(int j=0;j<N;++j) mstd(j) = (*this)(j).standardize();
    }

    // This gives the resulting matrix type when (m(i,j) op P) is applied to each element.
    // It is an MxN vector with default column and row spacing, and element types which
    // are the regular composite result of E op P. Typically P is a scalar type but
    // it doesn't have to be.
    template <class P> struct EltResult { 
        typedef Mat<M,N, typename CNT<E>::template Result<P>::Mul, M, 1> Mul;
        typedef Mat<M,N, typename CNT<E>::template Result<P>::Dvd, M, 1> Dvd;
        typedef Mat<M,N, typename CNT<E>::template Result<P>::Add, M, 1> Add;
        typedef Mat<M,N, typename CNT<E>::template Result<P>::Sub, M, 1> Sub;
    };

    // This is the composite result for m op P where P is some kind of appropriately shaped
    // non-scalar type.
    template <class P> struct Result { 
        typedef MulCNTs<M,N,ArgDepth,Mat,ColSpacing,RowSpacing,
            CNT<P>::NRows, CNT<P>::NCols, CNT<P>::ArgDepth,
            P, CNT<P>::ColSpacing, CNT<P>::RowSpacing> MulOp;
        typedef typename MulOp::Type Mul;

        typedef MulCNTsNonConforming<M,N,ArgDepth,Mat,ColSpacing,RowSpacing,
            CNT<P>::NRows, CNT<P>::NCols, CNT<P>::ArgDepth,
            P, CNT<P>::ColSpacing, CNT<P>::RowSpacing> MulOpNonConforming;
        typedef typename MulOpNonConforming::Type MulNon;

        typedef DvdCNTs<M,N,ArgDepth,Mat,ColSpacing,RowSpacing,
            CNT<P>::NRows, CNT<P>::NCols, CNT<P>::ArgDepth,
            P, CNT<P>::ColSpacing, CNT<P>::RowSpacing> DvdOp;
        typedef typename DvdOp::Type Dvd;

        typedef AddCNTs<M,N,ArgDepth,Mat,ColSpacing,RowSpacing,
            CNT<P>::NRows, CNT<P>::NCols, CNT<P>::ArgDepth,
            P, CNT<P>::ColSpacing, CNT<P>::RowSpacing> AddOp;
        typedef typename AddOp::Type Add;

        typedef SubCNTs<M,N,ArgDepth,Mat,ColSpacing,RowSpacing,
            CNT<P>::NRows, CNT<P>::NCols, CNT<P>::ArgDepth,
            P, CNT<P>::ColSpacing, CNT<P>::RowSpacing> SubOp;
        typedef typename SubOp::Type Sub;
    };

    // Shape-preserving element substitution (always packed)
    template <class P> struct Substitute {
        typedef Mat<M,N,P> Type;
    };

    // Default construction initializes to NaN when debugging but
    // is left uninitialized otherwise.
	Mat(){ 
    #ifndef NDEBUG
        setToNaN();
    #endif
    }

    // It's important not to use the default copy constructor or copy
    // assignment because the compiler doesn't understand that we may
    // have noncontiguous storage and will try to copy the whole array.
    Mat(const Mat& src) {
        for (int j=0; j<N; ++j)
            (*this)(j) = src(j);
    }
    Mat& operator=(const Mat& src) {    // no harm if src and 'this' are the same
        for (int j=0; j<N; ++j)
           (*this)(j) = src(j);
        return *this;
    }

    // We want an implicit conversion from a Mat of the same length
    // and element type but with different spacings.
    template <int CSS, int RSS> Mat(const Mat<M,N,E,CSS,RSS>& src) {
        for (int j=0; j<N; ++j)
            (*this)(j) = src(j);
    }

    // We want an implicit conversion from a Mat of the same length
    // and *negated* element type, possibly with different spacings.
    template <int CSS, int RSS> Mat(const Mat<M,N,ENeg,CSS,RSS>& src) {
        for (int j=0; j<N; ++j)
            (*this)(j) = src(j);
    }

    // Construct a Mat from a Mat of the same dimensions, with any
    // spacings. Works as long as the element types are assignment compatible.
    template <class EE, int CSS, int RSS> explicit Mat(const Mat<M,N,EE,CSS,RSS>& mm)
      { for (int j=0;j<N;++j) (*this)(j) = mm(j);}

    // Construction using an element repeats that element on the diagonal
    // but sets the rest of the matrix to zero.
    explicit Mat(const E& e)
      { for (int j=0;j<N;++j) (*this)(j) = E(0); diag()=e; }

    // A bevy of constructors from individual exact-match elements IN ROW ORDER.
    Mat(const E& e0,const E& e1)
      {assert(M*N==2);d[rIx(0)]=e0;d[rIx(1)]=e1;}
    Mat(const E& e0,const E& e1,const E& e2)
      {assert(M*N==3);d[rIx(0)]=e0;d[rIx(1)]=e1;d[rIx(2)]=e2;}
    Mat(const E& e0,const E& e1,const E& e2,const E& e3)
      {assert(M*N==4);d[rIx(0)]=e0;d[rIx(1)]=e1;d[rIx(2)]=e2;d[rIx(3)]=e3;}
    Mat(const E& e0,const E& e1,const E& e2,const E& e3,const E& e4)
      {assert(M*N==5);d[rIx(0)]=e0;d[rIx(1)]=e1;d[rIx(2)]=e2;d[rIx(3)]=e3;d[rIx(4)]=e4;}
    Mat(const E& e0,const E& e1,const E& e2,const E& e3,const E& e4,
        const E& e5)
      {assert(M*N==6);d[rIx(0)]=e0;d[rIx(1)]=e1;d[rIx(2)]=e2;d[rIx(3)]=e3;d[rIx(4)]=e4;
       d[rIx(5)]=e5;}
    Mat(const E& e0,const E& e1,const E& e2,const E& e3,const E& e4,
        const E& e5,const E& e6)
      {assert(M*N==7);d[rIx(0)]=e0;d[rIx(1)]=e1;d[rIx(2)]=e2;d[rIx(3)]=e3;d[rIx(4)]=e4;
       d[rIx(5)]=e5;d[rIx(6)]=e6;}
    Mat(const E& e0,const E& e1,const E& e2,const E& e3,const E& e4,
        const E& e5,const E& e6,const E& e7)
      {assert(M*N==8);d[rIx(0)]=e0;d[rIx(1)]=e1;d[rIx(2)]=e2;d[rIx(3)]=e3;d[rIx(4)]=e4;
       d[rIx(5)]=e5;d[rIx(6)]=e6;d[rIx(7)]=e7;}
    Mat(const E& e0,const E& e1,const E& e2,const E& e3,const E& e4,
        const E& e5,const E& e6,const E& e7,const E& e8)
      {assert(M*N==9);d[rIx(0)]=e0;d[rIx(1)]=e1;d[rIx(2)]=e2;d[rIx(3)]=e3;d[rIx(4)]=e4;
       d[rIx(5)]=e5;d[rIx(6)]=e6;d[rIx(7)]=e7;d[rIx(8)]=e8;}
    Mat(const E& e0,const E& e1,const E& e2,const E& e3,const E& e4,
        const E& e5,const E& e6,const E& e7,const E& e8,const E& e9)
      {assert(M*N==10);d[rIx(0)]=e0;d[rIx(1)]=e1;d[rIx(2)]=e2;d[rIx(3)]=e3;d[rIx(4)]=e4;
       d[rIx(5)]=e5;d[rIx(6)]=e6;d[rIx(7)]=e7;d[rIx(8)]=e8;d[rIx(9)]=e9;}
    Mat(const E& e0,const E& e1,const E& e2,const E& e3,const E& e4,
        const E& e5,const E& e6,const E& e7,const E& e8,const E& e9,
        const E& e10)
      {assert(M*N==11);d[rIx(0)]=e0;d[rIx(1)]=e1;d[rIx(2)]=e2;d[rIx(3)]=e3;d[rIx(4)]=e4;
       d[rIx(5)]=e5;d[rIx(6)]=e6;d[rIx(7)]=e7;d[rIx(8)]=e8;d[rIx(9)]=e9;d[rIx(10)]=e10;}
    Mat(const E& e0,const E& e1,const E& e2,const E& e3,const E& e4,
        const E& e5,const E& e6,const E& e7,const E& e8,const E& e9,
        const E& e10, const E& e11)
      {assert(M*N==12);d[rIx(0)]=e0;d[rIx(1)]=e1;d[rIx(2)]=e2;d[rIx(3)]=e3;d[rIx(4)]=e4;
       d[rIx(5)]=e5;d[rIx(6)]=e6;d[rIx(7)]=e7;d[rIx(8)]=e8;d[rIx(9)]=e9;d[rIx(10)]=e10;
       d[rIx(11)]=e11;}
    Mat(const E& e0,const E& e1,const E& e2,const E& e3,const E& e4,
        const E& e5,const E& e6,const E& e7,const E& e8,const E& e9,
        const E& e10, const E& e11, const E& e12)
      {assert(M*N==13);d[rIx(0)]=e0;d[rIx(1)]=e1;d[rIx(2)]=e2;d[rIx(3)]=e3;d[rIx(4)]=e4;
       d[rIx(5)]=e5;d[rIx(6)]=e6;d[rIx(7)]=e7;d[rIx(8)]=e8;d[rIx(9)]=e9;d[rIx(10)]=e10;
       d[rIx(11)]=e11;d[rIx(12)]=e12;}
    Mat(const E& e0,const E& e1,const E& e2,const E& e3,const E& e4,
        const E& e5,const E& e6,const E& e7,const E& e8,const E& e9,
        const E& e10, const E& e11, const E& e12, const E& e13)
      {assert(M*N==14);d[rIx(0)]=e0;d[rIx(1)]=e1;d[rIx(2)]=e2;d[rIx(3)]=e3;d[rIx(4)]=e4;
       d[rIx(5)]=e5;d[rIx(6)]=e6;d[rIx(7)]=e7;d[rIx(8)]=e8;d[rIx(9)]=e9;d[rIx(10)]=e10;
       d[rIx(11)]=e11;d[rIx(12)]=e12;d[rIx(13)]=e13;}
    Mat(const E& e0,const E& e1,const E& e2,const E& e3,const E& e4,
        const E& e5,const E& e6,const E& e7,const E& e8,const E& e9,
        const E& e10, const E& e11, const E& e12, const E& e13, const E& e14)
      {assert(M*N==15);d[rIx(0)]=e0;d[rIx(1)]=e1;d[rIx(2)]=e2;d[rIx(3)]=e3;d[rIx(4)]=e4;
       d[rIx(5)]=e5;d[rIx(6)]=e6;d[rIx(7)]=e7;d[rIx(8)]=e8;d[rIx(9)]=e9;d[rIx(10)]=e10;
       d[rIx(11)]=e11;d[rIx(12)]=e12;d[rIx(13)]=e13;d[rIx(14)]=e14;}
    Mat(const E& e0,const E& e1,const E& e2,const E& e3,const E& e4,
        const E& e5,const E& e6,const E& e7,const E& e8,const E& e9,
        const E& e10, const E& e11, const E& e12, const E& e13, const E& e14, 
        const E& e15)
      {assert(M*N==16);d[rIx(0)]=e0;d[rIx(1)]=e1;d[rIx(2)]=e2;d[rIx(3)]=e3;d[rIx(4)]=e4;
       d[rIx(5)]=e5;d[rIx(6)]=e6;d[rIx(7)]=e7;d[rIx(8)]=e8;d[rIx(9)]=e9;d[rIx(10)]=e10;
       d[rIx(11)]=e11;d[rIx(12)]=e12;d[rIx(13)]=e13;d[rIx(14)]=e14;d[rIx(15)]=e15;}

    // Construction from 1-6 *exact match* Rows
    explicit Mat(const TRow& r0)
      { assert(M==1); (*this)[0]=r0; }
    Mat(const TRow& r0,const TRow& r1)
      { assert(M==2);(*this)[0]=r0;(*this)[1]=r1; }
    Mat(const TRow& r0,const TRow& r1,const TRow& r2)
      { assert(M==3);(*this)[0]=r0;(*this)[1]=r1;(*this)[2]=r2; }
    Mat(const TRow& r0,const TRow& r1,const TRow& r2,
        const TRow& r3)
      { assert(M==4);(*this)[0]=r0;(*this)[1]=r1;(*this)[2]=r2;(*this)[3]=r3; }
    Mat(const TRow& r0,const TRow& r1,const TRow& r2,
        const TRow& r3,const TRow& r4)
      { assert(M==5);(*this)[0]=r0;(*this)[1]=r1;(*this)[2]=r2;
        (*this)[3]=r3;(*this)[4]=r4; }
    Mat(const TRow& r0,const TRow& r1,const TRow& r2,
        const TRow& r3,const TRow& r4,const TRow& r5)
      { assert(M==6);(*this)[0]=r0;(*this)[1]=r1;(*this)[2]=r2;
        (*this)[3]=r3;(*this)[4]=r4;(*this)[5]=r5; }

    // Construction from 1-6 *compatible* Rows
    template <class EE, int SS> explicit Mat(const Row<N,EE,SS>& r0)
      { assert(M==1); (*this)[0]=r0; }
    template <class EE, int SS> Mat(const Row<N,EE,SS>& r0,const Row<N,EE,SS>& r1)
      { assert(M==2);(*this)[0]=r0;(*this)[1]=r1; }
    template <class EE, int SS> 
    Mat(const Row<N,EE,SS>& r0,const Row<N,EE,SS>& r1,const Row<N,EE,SS>& r2)
      { assert(M==3);(*this)[0]=r0;(*this)[1]=r1;(*this)[2]=r2; }
    template <class EE, int SS> 
    Mat(const Row<N,EE,SS>& r0,const Row<N,EE,SS>& r1,const Row<N,EE,SS>& r2,
        const Row<N,EE,SS>& r3)
      { assert(M==4);(*this)[0]=r0;(*this)[1]=r1;(*this)[2]=r2;(*this)[3]=r3; }
    template <class EE, int SS> 
    Mat(const Row<N,EE,SS>& r0,const Row<N,EE,SS>& r1,const Row<N,EE,SS>& r2,
        const Row<N,EE,SS>& r3,const Row<N,EE,SS>& r4)
      { assert(M==5);(*this)[0]=r0;(*this)[1]=r1;(*this)[2]=r2;
        (*this)[3]=r3;(*this)[4]=r4; }
    template <class EE, int SS> 
    Mat(const Row<N,EE,SS>& r0,const Row<N,EE,SS>& r1,const Row<N,EE,SS>& r2,
        const Row<N,EE,SS>& r3,const Row<N,EE,SS>& r4,const Row<N,EE,SS>& r5)
      { assert(M==6);(*this)[0]=r0;(*this)[1]=r1;(*this)[2]=r2;
        (*this)[3]=r3;(*this)[4]=r4;(*this)[5]=r5; }


    // Construction from 1-6 *exact match* Vecs
    explicit Mat(const TCol& r0)
      { assert(N==1); (*this)(0)=r0; }
    Mat(const TCol& r0,const TCol& r1)
      { assert(N==2);(*this)(0)=r0;(*this)(1)=r1; }
    Mat(const TCol& r0,const TCol& r1,const TCol& r2)
      { assert(N==3);(*this)(0)=r0;(*this)(1)=r1;(*this)(2)=r2; }
    Mat(const TCol& r0,const TCol& r1,const TCol& r2,
        const TCol& r3)
      { assert(N==4);(*this)(0)=r0;(*this)(1)=r1;(*this)(2)=r2;(*this)(3)=r3; }
    Mat(const TCol& r0,const TCol& r1,const TCol& r2,
        const TCol& r3,const TCol& r4)
      { assert(N==5);(*this)(0)=r0;(*this)(1)=r1;(*this)(2)=r2;
        (*this)(3)=r3;(*this)(4)=r4; }
    Mat(const TCol& r0,const TCol& r1,const TCol& r2,
        const TCol& r3,const TCol& r4,const TCol& r5)
      { assert(N==6);(*this)(0)=r0;(*this)(1)=r1;(*this)(2)=r2;
        (*this)(3)=r3;(*this)(4)=r4;(*this)(5)=r5; }

    // Construction from 1-6 *compatible* Vecs
    template <class EE, int SS> explicit Mat(const Vec<M,EE,SS>& r0)
      { assert(N==1); (*this)(0)=r0; }
    template <class EE, int SS> Mat(const Vec<M,EE,SS>& r0,const Vec<M,EE,SS>& r1)
      { assert(N==2);(*this)(0)=r0;(*this)(1)=r1; }
    template <class EE, int SS> 
    Mat(const Vec<M,EE,SS>& r0,const Vec<M,EE,SS>& r1,const Vec<M,EE,SS>& r2)
      { assert(N==3);(*this)(0)=r0;(*this)(1)=r1;(*this)(2)=r2; }
    template <class EE, int SS> 
    Mat(const Vec<M,EE,SS>& r0,const Vec<M,EE,SS>& r1,const Vec<M,EE,SS>& r2,
        const Vec<M,EE,SS>& r3)
      { assert(N==4);(*this)(0)=r0;(*this)(1)=r1;(*this)(2)=r2;(*this)(3)=r3; }
    template <class EE, int SS> 
    Mat(const Vec<M,EE,SS>& r0,const Vec<M,EE,SS>& r1,const Vec<M,EE,SS>& r2,
        const Vec<M,EE,SS>& r3,const Vec<M,EE,SS>& r4)
      { assert(N==5);(*this)(0)=r0;(*this)(1)=r1;(*this)(2)=r2;
        (*this)(3)=r3;(*this)(4)=r4; }
    template <class EE, int SS> 
    Mat(const Vec<M,EE,SS>& r0,const Vec<M,EE,SS>& r1,const Vec<M,EE,SS>& r2,
        const Vec<M,EE,SS>& r3,const Vec<M,EE,SS>& r4,const Vec<M,EE,SS>& r5)
      { assert(N==6);(*this)(0)=r0;(*this)(1)=r1;(*this)(2)=r2;
        (*this)(3)=r3;(*this)(4)=r4;(*this)(5)=r5; }

    // Construction from a pointer to anything assumes we're pointing
    // at a packed element list of the right length, in row order.
    template <class EE> explicit Mat(const EE* p)
      { assert(p); for(int i=0;i<M;++i) (*this)[i]=&p[i*N]; }

    // Assignment works similarly to copy -- if the lengths match,
    // go element-by-element. Otherwise, zero and then copy to each 
    // diagonal element.
    template <class EE, int CSS, int RSS> Mat& operator=(const Mat<M,N,EE,CSS,RSS>& mm) {
        for (int j=0; j<N; ++j) (*this)(j) = mm(j);
        return *this;
    }

    template <class EE> Mat& operator=(const EE* p) {
        assert(p); for(int i=0;i<M;++i) (*this)[i]=&p[i*N];
        return *this;
    }

    // Assignment ops
    template <class EE, int CSS, int RSS> Mat& 
    operator+=(const Mat<M,N,EE,CSS,RSS>& mm) {
        for (int j=0; j<N; ++j) (*this)(j) += mm(j);
        return *this;
    }
    template <class EE, int CSS, int RSS> Mat&
    operator+=(const Mat<M,N,negator<EE>,CSS,RSS>& mm) {
        for (int j=0; j<N; ++j) (*this)(j) -= -(mm(j));
        return *this;
    }

    template <class EE, int CSS, int RSS> Mat&
    operator-=(const Mat<M,N,EE,CSS,RSS>& mm) {
        for (int j=0; j<N; ++j) (*this)(j) -= mm(j);
        return *this;
    }
    template <class EE, int CSS, int RSS> Mat&
    operator-=(const Mat<M,N,negator<EE>,CSS,RSS>& mm) {
        for (int j=0; j<N; ++j) (*this)(j) += -(mm(j));
        return *this;
    }

    // In place matrix multiply can only be done when the RHS matrix is square of dimension
    // N (i.e., number of columns), and the elements are also *= compatible.
    template <class EE, int CSS, int RSS> Mat&
    operator*=(const Mat<N,N,EE,CSS,RSS>& mm) {
        const Mat t(*this);
        for (int j=0; j<N; ++j)
            for (int i=0; i<M; ++i)
                (*this)(i,j) = t[i] * mm(j);
        return *this;
    }

    // Conforming binary ops with 'this' on left, producing new packed result.
    // Cases: m=m+-m, m=m+-sy, m=m*m, m=m*sy, v=m*v

    // m= this + m
    template <class E2, int CS2, int RS2> 
    typename Result<Mat<M,N,E2,CS2,RS2> >::Add
    conformingAdd(const Mat<M,N,E2,CS2,RS2>& r) const {
        typename Result<Mat<M,N,E2,CS2,RS2> >::Add result;
        for (int j=0;j<N;++j) result(j) = (*this)(j) + r(j);
        return result;
    }
    // m= this - m
    template <class E2, int CS2, int RS2> 
    typename Result<Mat<M,N,E2,CS2,RS2> >::Sub
    conformingSubtract(const Mat<M,N,E2,CS2,RS2>& r) const {
        typename Result<Mat<M,N,E2,CS2,RS2> >::Sub result;
        for (int j=0;j<N;++j) result(j) = (*this)(j) - r(j);
        return result;
    }
    // m= m - this
    template <class E2, int CS2, int RS2> 
    typename Mat<M,N,E2,CS2,RS2>::template Result<Mat>::Sub
    conformingSubtractFromLeft(const Mat<M,N,E2,CS2,RS2>& l) const {
        return l.conformingSubtract(*this);
    }

    // We always punt to the SymMat since it knows better what to do.
    // m = this + sym
    template <class E2, int RS2> 
    typename Result<SymMat<M,E2,RS2> >::Add
    conformingAdd(const SymMat<M,E2,RS2>& sy) const {
        assert(M==N);
        return sy.conformingAdd(*this);
    }
    // m= this - sym
    template <class E2, int RS2> 
    typename Result<SymMat<M,E2,RS2> >::Sub
    conformingSubtract(const SymMat<M,E2,RS2>& sy) const {
        assert(M==N);
        return sy.conformingSubtractFromLeft(*this);
    }
    // m= sym - this
    template <class E2, int RS2> 
    typename SymMat<M,E2,RS2>::template Result<Mat>::Sub
    conformingSubtractFromLeft(const SymMat<M,E2,RS2>& sy) const {
        assert(M==N);
        return sy.conformingSubtract(*this);
    }

    // m= this * m
    template <int N2, class E2, int CS2, int RS2>
    typename Result<Mat<N,N2,E2,CS2,RS2> >::Mul
    conformingMultiply(const Mat<N,N2,E2,CS2,RS2>& m) const {
        typename Result<Mat<N,N2,E2,CS2,RS2> >::Mul result;
        for (int j=0;j<N2;++j)
            for (int i=0;i<M;++i)
                result(i,j) = (*this)[i].conformingMultiply(m(j)); // i.e., dot()
        return result;
    }
    // m= m * this
    template <int M2, class E2, int CS2, int RS2>
    typename Mat<M2,M,E2,CS2,RS2>::template Result<Mat>::Mul
    conformingMultiplyFromLeft(const Mat<M2,M,E2,CS2,RS2>& m) const {
        return m.conformingMultiply(*this);
    }

    // m= this / m = this * m.invert()
    template <int M2, class E2, int CS2, int RS2> 
    typename Result<Mat<M2,N,E2,CS2,RS2> >::Dvd
    conformingDivide(const Mat<M2,N,E2,CS2,RS2>& m) const {
        return conformingMultiply(m.invert());
    }
    // m= m / this = m * this.invert()
    template <int M2, class E2, int CS2, int RS2> 
    typename Mat<M2,N,E2,CS2,RS2>::template Result<Mat>::Dvd
    conformingDivideFromLeft(const Mat<M2,N,E2,CS2,RS2>& m) const {
        return m.conformingMultiply((*this).invert());
    }
    
    const TRow& operator[](int i) const { return row(i); }
    TRow&       operator[](int i)       { return row(i); }    
    const TCol& operator()(int j) const { return col(j); }
    TCol&       operator()(int j)       { return col(j); }
    
    const E& operator()(int i,int j) const { return d[i*RS+j*CS]; }
    E&       operator()(int i,int j)       { return d[i*RS+j*CS]; }

    // This is the scalar Frobenius norm.
    ScalarSq normSqr() const { return scalarNormSqr(); }
    ScalarSq norm()    const { return std::sqrt(scalarNormSqr()); }

    // Default inversion. Assume full rank if square, otherwise return
    // pseudoinverse. (Mostly TODO)
    TInvert invert() const;

    const Mat&   operator+() const { return *this; }
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

    const TRow& row(int i) const 
      { assert(0<=i&&i<M); return *reinterpret_cast<const TRow*>(&d[i*RS]); }
    TRow&       row(int i)       
      { assert(0<=i&&i<M); return *reinterpret_cast<      TRow*>(&d[i*RS]); }

    const TCol& col(int j) const 
      { assert(0<=j&&j<N); return *reinterpret_cast<const TCol*>(&d[j*CS]); }
    TCol&       col(int j)       
      { assert(0<=j&&j<N); return *reinterpret_cast<      TCol*>(&d[j*CS]); }    


    const TDiag& diag() const { return *reinterpret_cast<const TDiag*>(d); }
    TDiag&       diag()       { return *reinterpret_cast<TDiag*>(d); }

    StdNumber trace() const {return diag().sum();}

    // These are elementwise binary operators, (this op ee) by default but (ee op this) if
    // 'FromLeft' appears in the name. The result is a packed Mat<M,N> but the element type
    // may change. These are mostly used to implement global operators.
    // We call these "scalar" operators but actually the "scalar" can be a composite type.

    //TODO: consider converting 'e' to Standard Numbers as precalculation and changing
    // return type appropriately.
    template <class EE> Mat<M,N, typename CNT<E>::template Result<EE>::Mul>
    scalarMultiply(const EE& e) const {
        Mat<M,N, typename CNT<E>::template Result<EE>::Mul> result;
        for (int j=0; j<N; ++j) result(j) = (*this)(j).scalarMultiply(e);
        return result;
    }
    template <class EE> Mat<M,N, typename CNT<EE>::template Result<E>::Mul>
    scalarMultiplyFromLeft(const EE& e) const {
        Mat<M,N, typename CNT<EE>::template Result<E>::Mul> result;
        for (int j=0; j<N; ++j) result(j) = (*this)(j).scalarMultiplyFromLeft(e);
        return result;
    }

    // TODO: should precalculate and store 1/e, while converting to Standard Numbers. Note
    // that return type should change appropriately.
    template <class EE> Mat<M,N, typename CNT<E>::template Result<EE>::Dvd>
    scalarDivide(const EE& e) const {
        Mat<M,N, typename CNT<E>::template Result<EE>::Dvd> result;
        for (int j=0; j<N; ++j) result(j) = (*this)(j).scalarDivide(e);
        return result;
    }
    template <class EE> Mat<M,N, typename CNT<EE>::template Result<E>::Dvd>
    scalarDivideFromLeft(const EE& e) const {
        Mat<M,N, typename CNT<EE>::template Result<E>::Dvd> result;
        for (int j=0; j<N; ++j) result(j) = (*this)(j).scalarDivideFromLeft(e);
        return result;
    }

    template <class EE> Mat<M,N, typename CNT<E>::template Result<EE>::Add>
    scalarAdd(const EE& e) const {
        Mat<M,N, typename CNT<E>::template Result<EE>::Add> result;
        for (int j=0; j<N; ++j) result(j) = (*this)(j).scalarAdd(e);
        return result;
    }
    // Add is commutative, so no 'FromLeft'.

    template <class EE> Mat<M,N, typename CNT<E>::template Result<EE>::Sub>
    scalarSubtract(const EE& e) const {
        Mat<M,N, typename CNT<E>::template Result<EE>::Sub> result;
        for (int j=0; j<N; ++j) result(j) = (*this)(j).scalarSubtract(e);
        return result;
    }
    template <class EE> Mat<M,N, typename CNT<EE>::template Result<E>::Sub>
    scalarSubtractFromLeft(const EE& e) const {
        Mat<M,N, typename CNT<EE>::template Result<E>::Sub> result;
        for (int j=0; j<N; ++j) result(j) = (*this)(j).scalarSubtractFromLeft(e);
        return result;
    }

    // Generic assignments for any element type not listed explicitly, including scalars.
    // These are done repeatedly for each element and only work if the operation can
    // be performed leaving the original element type.
    template <class EE> Mat& operator =(const EE& e) {return scalarEq(e);}
    template <class EE> Mat& operator+=(const EE& e) {return scalarPlusEq(e);}
    template <class EE> Mat& operator-=(const EE& e) {return scalarMinusEq(e);}
    template <class EE> Mat& operator*=(const EE& e) {return scalarTimesEq(e);}
    template <class EE> Mat& operator/=(const EE& e) {return scalarDivideEq(e);}

    // Generalized scalar assignment & computed assignment methods. These will work
    // for any assignment-compatible element, not just scalars.
    template <class EE> Mat& scalarEq(const EE& ee)
      { for(int j=0; j<N; ++j) (*this)(j).scalarEq(EE(0)); diag().scalarEq(ee); return *this; }

    template <class EE> Mat& scalarPlusEq(const EE& ee)
      { for(int j=0; j<N; ++j) diag().scalarPlusEq(ee); return *this; }

    template <class EE> Mat& scalarMinusEq(const EE& ee)
      { for(int j=0; j<N; ++j) diag().scalarMinusEq(ee); return *this; }
    template <class EE> Mat& scalarMinusEqFromLeft(const EE& ee)
      { for(int j=0; j<N; ++j) diag().scalarMinusEqFromLeft(ee); return *this; }

    template <class EE> Mat& scalarTimesEq(const EE& ee)
      { for(int j=0; j<N; ++j) (*this)(j).scalarTimesEq(ee); return *this; }
    template <class EE> Mat& scalarTimesEqFromLeft(const EE& ee)
      { for(int j=0; j<N; ++j) (*this)(j).scalarTimesEqFromLeft(ee); return *this; } 

    template <class EE> Mat& scalarDivideEq(const EE& ee)
      { for(int j=0; j<N; ++j) (*this)(j).scalarDivideEq(ee); return *this; }
    template <class EE> Mat& scalarDivideEqFromLeft(const EE& ee)
      { for(int j=0; j<N; ++j) (*this)(j).scalarDivideEqFromLeft(ee); return *this; } 

    void setToNaN() {
        for (int j=0; j<N; ++j)
            (*this)(j).setToNaN();
    }

    // Extract a sub-Mat with size known at compile time. These have to be
    // called with explicit template arguments, e.g. getSubMat<3,4>(i,j).

    template <int MM, int NN> struct SubMat {
        typedef Mat<MM,NN,ELT,CS,RS> Type;
    };

    template <int MM, int NN>
    const typename SubMat<MM,NN>::Type& getSubMat(int i, int j) const {
        assert(0 <= i && i + MM <= M);
        assert(0 <= j && j + NN <= N);
        return SubMat<MM,NN>::Type::getAs(&(*this)(i,j));
    }
    template <int MM, int NN>
    typename SubMat<MM,NN>::Type& updSubMat(int i, int j) {
        assert(0 <= i && i + MM <= M);
        assert(0 <= j && j + NN <= N);
        return SubMat<MM,NN>::Type::updAs(&(*this)(i,j));
    }

    // These assume we are given a pointer to d[0] of a Mat<M,N,E,CS,RS> like this one.
    static const Mat& getAs(const ELT* p)  {return *reinterpret_cast<const Mat*>(p);}
    static Mat&       updAs(ELT* p)        {return *reinterpret_cast<Mat*>(p);}

    // Note packed spacing
    static Mat<M,N,ELT,M,1> getNaN() { 
        Mat<M,N,ELT,M,1> m;
        m.setToNaN();
        return m;
    }

private:
    E d[NActualElements];

    // This permits running through d as though it were stored
    // in row order with packed elements. Pass in the k'th value
    // of the row ordering and get back the index into d where
    // that element is stored.
    int rIx(int k) const {
        const int row = k / N;
        const int col = k % N; // that's modulus, not cross product!
        return row*RS + col*CS;
    }
};

//////////////////////////////////////////////
// Global operators involving two matrices. //
//   m+m, m-m, m*m, m==m, m!=m              //
//////////////////////////////////////////////

template <int M, int N, class EL, int CSL, int RSL, class ER, int CSR, int RSR> inline 
typename Mat<M,N,EL,CSL,RSL>::template Result<Mat<M,N,ER,CSR,RSR> >::Add
operator+(const Mat<M,N,EL,CSL,RSL>& l, const Mat<M,N,ER,CSR,RSR>& r) { 
    return Mat<M,N,EL,CSL,RSL>::template Result<Mat<M,N,ER,CSR,RSR> >
        ::AddOp::perform(l,r);
}

template <int M, int N, class EL, int CSL, int RSL, class ER, int CSR, int RSR> inline
typename Mat<M,N,EL,CSL,RSL>::template Result<Mat<M,N,ER,CSR,RSR> >::Sub
operator-(const Mat<M,N,EL,CSL,RSL>& l, const Mat<M,N,ER,CSR,RSR>& r) { 
    return Mat<M,N,EL,CSL,RSL>::template Result<Mat<M,N,ER,CSR,RSR> >
        ::SubOp::perform(l,r);
}

// Matrix multiply of an MxN by NxP to produce a packed MxP.
template <int M, int N, class EL, int CSL, int RSL, int P, class ER, int CSR, int RSR> inline
typename Mat<M,N,EL,CSL,RSL>::template Result<Mat<N,P,ER,CSR,RSR> >::Mul
operator*(const Mat<M,N,EL,CSL,RSL>& l, const Mat<N,P,ER,CSR,RSR>& r) { 
    return Mat<M,N,EL,CSL,RSL>::template Result<Mat<N,P,ER,CSR,RSR> >
        ::MulOp::perform(l,r);
}

// Non-conforming matrix multiply of an MxN by MMxNN; will be a scalar multiply if one
// has scalar elements and the other has composite elements.
template <int M, int N, class EL, int CSL, int RSL, int MM, int NN, class ER, int CSR, int RSR> inline
typename Mat<M,N,EL,CSL,RSL>::template Result<Mat<MM,NN,ER,CSR,RSR> >::MulNon
operator*(const Mat<M,N,EL,CSL,RSL>& l, const Mat<MM,NN,ER,CSR,RSR>& r) { 
    return Mat<M,N,EL,CSL,RSL>::template Result<Mat<MM,NN,ER,CSR,RSR> >
                ::MulOpNonConforming::perform(l,r);
}

template <int M, int N, class EL, int CSL, int RSL, class ER, int CSR, int RSR> inline
bool operator==(const Mat<M,N,EL,CSL,RSL>& l, const Mat<M,N,ER,CSR,RSR>& r) { 
    for (int j=0; j<N; ++j)
        if (l(j) != r(j)) return false;
    return true;
}
template <int M, int N, class EL, int CSL, int RSL, class ER, int CSR, int RSR> inline
bool operator!=(const Mat<M,N,EL,CSL,RSL>& l, const Mat<M,N,ER,CSR,RSR>& r) { 
    return !(l==r);
}


///////////////////////////////////////////////////////
// Global operators involving a matrix and a scalar. //
///////////////////////////////////////////////////////

// SCALAR MULTIPLY

// m = m*real, real*m 
template <int M, int N, class E, int CS, int RS> inline
typename Mat<M,N,E,CS,RS>::template Result<float>::Mul
operator*(const Mat<M,N,E,CS,RS>& l, const float& r)
  { return Mat<M,N,E,CS,RS>::template Result<float>::MulOp::perform(l,r); }
template <int M, int N, class E, int CS, int RS> inline
typename Mat<M,N,E,CS,RS>::template Result<float>::Mul
operator*(const float& l, const Mat<M,N,E,CS,RS>& r) {return r*l;}

template <int M, int N, class E, int CS, int RS> inline
typename Mat<M,N,E,CS,RS>::template Result<double>::Mul
operator*(const Mat<M,N,E,CS,RS>& l, const double& r)
  { return Mat<M,N,E,CS,RS>::template Result<double>::MulOp::perform(l,r); }
template <int M, int N, class E, int CS, int RS> inline
typename Mat<M,N,E,CS,RS>::template Result<double>::Mul
operator*(const double& l, const Mat<M,N,E,CS,RS>& r) {return r*l;}

template <int M, int N, class E, int CS, int RS> inline
typename Mat<M,N,E,CS,RS>::template Result<long double>::Mul
operator*(const Mat<M,N,E,CS,RS>& l, const long double& r)
  { return Mat<M,N,E,CS,RS>::template Result<long double>::MulOp::perform(l,r); }
template <int M, int N, class E, int CS, int RS> inline
typename Mat<M,N,E,CS,RS>::template Result<long double>::Mul
operator*(const long double& l, const Mat<M,N,E,CS,RS>& r) {return r*l;}

// m = m*int, int*m -- just convert int to m's precision float
template <int M, int N, class E, int CS, int RS> inline
typename Mat<M,N,E,CS,RS>::template Result<typename CNT<E>::Precision>::Mul
operator*(const Mat<M,N,E,CS,RS>& l, int r) {return l * (typename CNT<E>::Precision)r;}
template <int M, int N, class E, int CS, int RS> inline
typename Mat<M,N,E,CS,RS>::template Result<typename CNT<E>::Precision>::Mul
operator*(int l, const Mat<M,N,E,CS,RS>& r) {return r * (typename CNT<E>::Precision)l;}

// Complex, conjugate, and negator are all easy to templatize.

// m = m*complex, complex*m
template <int M, int N, class E, int CS, int RS, class R> inline
typename Mat<M,N,E,CS,RS>::template Result<std::complex<R> >::Mul
operator*(const Mat<M,N,E,CS,RS>& l, const std::complex<R>& r)
  { return Mat<M,N,E,CS,RS>::template Result<std::complex<R> >::MulOp::perform(l,r); }
template <int M, int N, class E, int CS, int RS, class R> inline
typename Mat<M,N,E,CS,RS>::template Result<std::complex<R> >::Mul
operator*(const std::complex<R>& l, const Mat<M,N,E,CS,RS>& r) {return r*l;}

// m = m*conjugate, conjugate*m (convert conjugate->complex)
template <int M, int N, class E, int CS, int RS, class R> inline
typename Mat<M,N,E,CS,RS>::template Result<std::complex<R> >::Mul
operator*(const Mat<M,N,E,CS,RS>& l, const conjugate<R>& r) {return l*(std::complex<R>)r;}
template <int M, int N, class E, int CS, int RS, class R> inline
typename Mat<M,N,E,CS,RS>::template Result<std::complex<R> >::Mul
operator*(const conjugate<R>& l, const Mat<M,N,E,CS,RS>& r) {return r*(std::complex<R>)l;}

// m = m*negator, negator*m: convert negator to standard number
template <int M, int N, class E, int CS, int RS, class R> inline
typename Mat<M,N,E,CS,RS>::template Result<typename negator<R>::StdNumber>::Mul
operator*(const Mat<M,N,E,CS,RS>& l, const negator<R>& r) {return l * (typename negator<R>::StdNumber)(R)r;}
template <int M, int N, class E, int CS, int RS, class R> inline
typename Mat<M,N,E,CS,RS>::template Result<typename negator<R>::StdNumber>::Mul
operator*(const negator<R>& l, const Mat<M,N,E,CS,RS>& r) {return r * (typename negator<R>::StdNumber)(R)l;}


// SCALAR DIVIDE. This is a scalar operation when the scalar is on the right,
// but when it is on the left it means scalar * pseudoInverse(mat), 
// which is a matrix whose type is like the matrix's Hermitian transpose.

// m = m/real, real/m 
template <int M, int N, class E, int CS, int RS> inline
typename Mat<M,N,E,CS,RS>::template Result<float>::Dvd
operator/(const Mat<M,N,E,CS,RS>& l, const float& r)
  { return Mat<M,N,E,CS,RS>::template Result<float>::DvdOp::perform(l,r); }
template <int M, int N, class E, int CS, int RS> inline
typename CNT<float>::template Result<Mat<M,N,E,CS,RS> >::Dvd
operator/(const float& l, const Mat<M,N,E,CS,RS>& r)
  { return CNT<float>::template Result<Mat<M,N,E,CS,RS> >::DvdOp::perform(l,r); }

template <int M, int N, class E, int CS, int RS> inline
typename Mat<M,N,E,CS,RS>::template Result<double>::Dvd
operator/(const Mat<M,N,E,CS,RS>& l, const double& r)
  { return Mat<M,N,E,CS,RS>::template Result<double>::DvdOp::perform(l,r); }
template <int M, int N, class E, int CS, int RS> inline
typename CNT<double>::template Result<Mat<M,N,E,CS,RS> >::Dvd
operator/(const double& l, const Mat<M,N,E,CS,RS>& r)
  { return CNT<double>::template Result<Mat<M,N,E,CS,RS> >::DvdOp::perform(l,r); }

template <int M, int N, class E, int CS, int RS> inline
typename Mat<M,N,E,CS,RS>::template Result<long double>::Dvd
operator/(const Mat<M,N,E,CS,RS>& l, const long double& r)
  { return Mat<M,N,E,CS,RS>::template Result<long double>::DvdOp::perform(l,r); }
template <int M, int N, class E, int CS, int RS> inline
typename CNT<long double>::template Result<Mat<M,N,E,CS,RS> >::Dvd
operator/(const long double& l, const Mat<M,N,E,CS,RS>& r)
  { return CNT<long double>::template Result<Mat<M,N,E,CS,RS> >::DvdOp::perform(l,r); }

// m = m/int, int/m -- just convert int to m's precision float
template <int M, int N, class E, int CS, int RS> inline
typename Mat<M,N,E,CS,RS>::template Result<typename CNT<E>::Precision>::Dvd
operator/(const Mat<M,N,E,CS,RS>& l, int r) {return l / (typename CNT<E>::Precision)r;}
template <int M, int N, class E, int CS, int RS> inline
typename CNT<typename CNT<E>::Precision>::template Result<Mat<M,N,E,CS,RS> >::Dvd
operator/(int l, const Mat<M,N,E,CS,RS>& r) {return (typename CNT<E>::Precision)l / r;}


// Complex, conjugate, and negator are all easy to templatize.

// m = m/complex, complex/m
template <int M, int N, class E, int CS, int RS, class R> inline
typename Mat<M,N,E,CS,RS>::template Result<std::complex<R> >::Dvd
operator/(const Mat<M,N,E,CS,RS>& l, const std::complex<R>& r)
  { return Mat<M,N,E,CS,RS>::template Result<std::complex<R> >::DvdOp::perform(l,r); }
template <int M, int N, class E, int CS, int RS, class R> inline
typename CNT<std::complex<R> >::template Result<Mat<M,N,E,CS,RS> >::Dvd
operator/(const std::complex<R>& l, const Mat<M,N,E,CS,RS>& r)
  { return CNT<std::complex<R> >::template Result<Mat<M,N,E,CS,RS> >::DvdOp::perform(l,r); }

// m = m/conjugate, conjugate/m (convert conjugate->complex)
template <int M, int N, class E, int CS, int RS, class R> inline
typename Mat<M,N,E,CS,RS>::template Result<std::complex<R> >::Dvd
operator/(const Mat<M,N,E,CS,RS>& l, const conjugate<R>& r) {return l/(std::complex<R>)r;}
template <int M, int N, class E, int CS, int RS, class R> inline
typename CNT<std::complex<R> >::template Result<Mat<M,N,E,CS,RS> >::Dvd
operator/(const conjugate<R>& l, const Mat<M,N,E,CS,RS>& r) {return (std::complex<R>)l/r;}

// m = m/negator, negator/m: convert negator to a standard number
template <int M, int N, class E, int CS, int RS, class R> inline
typename Mat<M,N,E,CS,RS>::template Result<typename negator<R>::StdNumber>::Dvd
operator/(const Mat<M,N,E,CS,RS>& l, const negator<R>& r) {return l/(typename negator<R>::StdNumber)(R)r;}
template <int M, int N, class E, int CS, int RS, class R> inline
typename CNT<R>::template Result<Mat<M,N,E,CS,RS> >::Dvd
operator/(const negator<R>& l, const Mat<M,N,E,CS,RS>& r) {return (typename negator<R>::StdNumber)(R)l/r;}


// Add and subtract are odd as scalar ops. They behave as though the
// scalar stands for a conforming matrix whose diagonal elements are that,
// scalar and then a normal matrix add or subtract is done.

// SCALAR ADD

// m = m+real, real+m 
template <int M, int N, class E, int CS, int RS> inline
typename Mat<M,N,E,CS,RS>::template Result<float>::Add
operator+(const Mat<M,N,E,CS,RS>& l, const float& r)
  { return Mat<M,N,E,CS,RS>::template Result<float>::AddOp::perform(l,r); }
template <int M, int N, class E, int CS, int RS> inline
typename Mat<M,N,E,CS,RS>::template Result<float>::Add
operator+(const float& l, const Mat<M,N,E,CS,RS>& r) {return r+l;}

template <int M, int N, class E, int CS, int RS> inline
typename Mat<M,N,E,CS,RS>::template Result<double>::Add
operator+(const Mat<M,N,E,CS,RS>& l, const double& r)
  { return Mat<M,N,E,CS,RS>::template Result<double>::AddOp::perform(l,r); }
template <int M, int N, class E, int CS, int RS> inline
typename Mat<M,N,E,CS,RS>::template Result<double>::Add
operator+(const double& l, const Mat<M,N,E,CS,RS>& r) {return r+l;}

template <int M, int N, class E, int CS, int RS> inline
typename Mat<M,N,E,CS,RS>::template Result<long double>::Add
operator+(const Mat<M,N,E,CS,RS>& l, const long double& r)
  { return Mat<M,N,E,CS,RS>::template Result<long double>::AddOp::perform(l,r); }
template <int M, int N, class E, int CS, int RS> inline
typename Mat<M,N,E,CS,RS>::template Result<long double>::Add
operator+(const long double& l, const Mat<M,N,E,CS,RS>& r) {return r+l;}

// m = m+int, int+m -- just convert int to m's precision float
template <int M, int N, class E, int CS, int RS> inline
typename Mat<M,N,E,CS,RS>::template Result<typename CNT<E>::Precision>::Add
operator+(const Mat<M,N,E,CS,RS>& l, int r) {return l + (typename CNT<E>::Precision)r;}
template <int M, int N, class E, int CS, int RS> inline
typename Mat<M,N,E,CS,RS>::template Result<typename CNT<E>::Precision>::Add
operator+(int l, const Mat<M,N,E,CS,RS>& r) {return r + (typename CNT<E>::Precision)l;}

// Complex, conjugate, and negator are all easy to templatize.

// m = m+complex, complex+m
template <int M, int N, class E, int CS, int RS, class R> inline
typename Mat<M,N,E,CS,RS>::template Result<std::complex<R> >::Add
operator+(const Mat<M,N,E,CS,RS>& l, const std::complex<R>& r)
  { return Mat<M,N,E,CS,RS>::template Result<std::complex<R> >::AddOp::perform(l,r); }
template <int M, int N, class E, int CS, int RS, class R> inline
typename Mat<M,N,E,CS,RS>::template Result<std::complex<R> >::Add
operator+(const std::complex<R>& l, const Mat<M,N,E,CS,RS>& r) {return r+l;}

// m = m+conjugate, conjugate+m (convert conjugate->complex)
template <int M, int N, class E, int CS, int RS, class R> inline
typename Mat<M,N,E,CS,RS>::template Result<std::complex<R> >::Add
operator+(const Mat<M,N,E,CS,RS>& l, const conjugate<R>& r) {return l+(std::complex<R>)r;}
template <int M, int N, class E, int CS, int RS, class R> inline
typename Mat<M,N,E,CS,RS>::template Result<std::complex<R> >::Add
operator+(const conjugate<R>& l, const Mat<M,N,E,CS,RS>& r) {return r+(std::complex<R>)l;}

// m = m+negator, negator+m: convert negator to standard number
template <int M, int N, class E, int CS, int RS, class R> inline
typename Mat<M,N,E,CS,RS>::template Result<typename negator<R>::StdNumber>::Add
operator+(const Mat<M,N,E,CS,RS>& l, const negator<R>& r) {return l + (typename negator<R>::StdNumber)(R)r;}
template <int M, int N, class E, int CS, int RS, class R> inline
typename Mat<M,N,E,CS,RS>::template Result<typename negator<R>::StdNumber>::Add
operator+(const negator<R>& l, const Mat<M,N,E,CS,RS>& r) {return r + (typename negator<R>::StdNumber)(R)l;}

// SCALAR SUBTRACT -- careful, not commutative.

// m = m-real, real-m 
template <int M, int N, class E, int CS, int RS> inline
typename Mat<M,N,E,CS,RS>::template Result<float>::Sub
operator-(const Mat<M,N,E,CS,RS>& l, const float& r)
  { return Mat<M,N,E,CS,RS>::template Result<float>::SubOp::perform(l,r); }
template <int M, int N, class E, int CS, int RS> inline
typename CNT<float>::template Result<Mat<M,N,E,CS,RS> >::Sub
operator-(const float& l, const Mat<M,N,E,CS,RS>& r)
  { return CNT<float>::template Result<Mat<M,N,E,CS,RS> >::SubOp::perform(l,r); }

template <int M, int N, class E, int CS, int RS> inline
typename Mat<M,N,E,CS,RS>::template Result<double>::Sub
operator-(const Mat<M,N,E,CS,RS>& l, const double& r)
  { return Mat<M,N,E,CS,RS>::template Result<double>::SubOp::perform(l,r); }
template <int M, int N, class E, int CS, int RS> inline
typename CNT<double>::template Result<Mat<M,N,E,CS,RS> >::Sub
operator-(const double& l, const Mat<M,N,E,CS,RS>& r)
  { return CNT<double>::template Result<Mat<M,N,E,CS,RS> >::SubOp::perform(l,r); }

template <int M, int N, class E, int CS, int RS> inline
typename Mat<M,N,E,CS,RS>::template Result<long double>::Sub
operator-(const Mat<M,N,E,CS,RS>& l, const long double& r)
  { return Mat<M,N,E,CS,RS>::template Result<long double>::SubOp::perform(l,r); }
template <int M, int N, class E, int CS, int RS> inline
typename CNT<long double>::template Result<Mat<M,N,E,CS,RS> >::Sub
operator-(const long double& l, const Mat<M,N,E,CS,RS>& r)
  { return CNT<long double>::template Result<Mat<M,N,E,CS,RS> >::SubOp::perform(l,r); }

// m = m-int, int-m // just convert int to m's precision float
template <int M, int N, class E, int CS, int RS> inline
typename Mat<M,N,E,CS,RS>::template Result<typename CNT<E>::Precision>::Sub
operator-(const Mat<M,N,E,CS,RS>& l, int r) {return l - (typename CNT<E>::Precision)r;}
template <int M, int N, class E, int CS, int RS> inline
typename CNT<typename CNT<E>::Precision>::template Result<Mat<M,N,E,CS,RS> >::Sub
operator-(int l, const Mat<M,N,E,CS,RS>& r) {return (typename CNT<E>::Precision)l - r;}


// Complex, conjugate, and negator are all easy to templatize.

// m = m-complex, complex-m
template <int M, int N, class E, int CS, int RS, class R> inline
typename Mat<M,N,E,CS,RS>::template Result<std::complex<R> >::Sub
operator-(const Mat<M,N,E,CS,RS>& l, const std::complex<R>& r)
  { return Mat<M,N,E,CS,RS>::template Result<std::complex<R> >::SubOp::perform(l,r); }
template <int M, int N, class E, int CS, int RS, class R> inline
typename CNT<std::complex<R> >::template Result<Mat<M,N,E,CS,RS> >::Sub
operator-(const std::complex<R>& l, const Mat<M,N,E,CS,RS>& r)
  { return CNT<std::complex<R> >::template Result<Mat<M,N,E,CS,RS> >::SubOp::perform(l,r); }

// m = m-conjugate, conjugate-m (convert conjugate->complex)
template <int M, int N, class E, int CS, int RS, class R> inline
typename Mat<M,N,E,CS,RS>::template Result<std::complex<R> >::Sub
operator-(const Mat<M,N,E,CS,RS>& l, const conjugate<R>& r) {return l-(std::complex<R>)r;}
template <int M, int N, class E, int CS, int RS, class R> inline
typename CNT<std::complex<R> >::template Result<Mat<M,N,E,CS,RS> >::Sub
operator-(const conjugate<R>& l, const Mat<M,N,E,CS,RS>& r) {return (std::complex<R>)l-r;}

// m = m-negator, negator-m: convert negator to standard number
template <int M, int N, class E, int CS, int RS, class R> inline
typename Mat<M,N,E,CS,RS>::template Result<typename negator<R>::StdNumber>::Sub
operator-(const Mat<M,N,E,CS,RS>& l, const negator<R>& r) {return l-(typename negator<R>::StdNumber)(R)r;}
template <int M, int N, class E, int CS, int RS, class R> inline
typename CNT<R>::template Result<Mat<M,N,E,CS,RS> >::Sub
operator-(const negator<R>& l, const Mat<M,N,E,CS,RS>& r) {return (typename negator<R>::StdNumber)(R)l-r;}


// Mat I/O
template <int M, int N, class E, int CS, int RS, class CHAR, class TRAITS> inline
std::basic_ostream<CHAR,TRAITS>&
operator<<(std::basic_ostream<CHAR,TRAITS>& o, const Mat<M,N,E,CS,RS>& m) {
    for (int i=0;i<M;++i) {
        o << std::endl << "[";
        for (int j=0;j<N;++j)         
            o << (j>0?",":"") << m(i,j);
        o << "]";
    }
    if (M) o << std::endl;
    return o; 
}

template <int M, int N, class E, int CS, int RS, class CHAR, class TRAITS> inline
std::basic_istream<CHAR,TRAITS>&
operator>>(std::basic_istream<CHAR,TRAITS>& is, Mat<M,N,E,CS,RS>& m) {
    // TODO: not sure how to do Vec input yet
    assert(false);
    return is;
}

} //namespace SimTK


#endif //SimTK_SIMMATRIX_SMALLMATRIX_MAT_H_
