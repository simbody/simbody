#ifndef SimTK_SIMMATRIX_SMALLMATRIX_MAT_H_
#define SimTK_SIMMATRIX_SMALLMATRIX_MAT_H_

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

/**@file
 * This file declares class Mat<NROWS, NCOLS, ELEMENT_TYPE, COL_SPACING, ROW_SPACING>.
 */

#include "SimTKcommon/internal/common.h"

namespace SimTK {

/** @brief This class represents a small matrix whose size is known at compile 
time, containing elements of any Composite Numerical Type (CNT) and engineered 
to have no runtime overhead whatsoever. 

@ingroup MatVecUtilities

Memory layout defaults to packed, column-ordered storage but can be specified to
have any regular row and column spacing. A %Mat object is itself a Composite 
Numerical Type and can thus be the element type for other matrix and vector 
types. Some common use cases are provided below.

@tparam M   The number of rows in this matrix (no default).
@tparam N   The number of columns in this matrix (no default).
@tparam ELT The element type; default is Real.
@tparam CS  Column spacing in memory as a multiple of element size (default M).
@tparam RS  %Row spacing in memory as a multiple of element size (default 1).

<b>Construction</b>

A 3x3 identity matrix can be constructed in the following ways:
\code
Mat<3,3,Real>(1,0,0, 0,1,0, 0,0,1);          //row-major
Mat33(1,0,0,0,1,0,0,0,1);
Mat33(Vec3(1,0,0),Vec3(0,1,0),Vec3(0,0,1));  //column-by-column
Mat33(1);
\endcode
Note that the default element type is Real, and that Mat33 is a typedef for
%Mat<3,3>; analogous typedefs exist for matrices of up to 9x9 elements.

<b>Manipulation</b>

Standard arithmetic operators can be used, as well as methods like %trace() and
%transpose(). Here are some usage examples, each of which prints a 2x2 identity
matrix:
\code
Mat23 myMat(2,1,0, 3,0,1);
std::cout << myMat.getSubMat<2,2>(0,1) << std::endl;
std::cout << myMat.dropCol(0) << std::endl;
std::cout << Mat12(1,0).appendRow(~Vec2(0,1)) << std::endl;
std::cout << Mat21(0,1).insertCol(0,Vec2(1,0)) << std::endl;
\endcode

<b>Conversion</b>

It may be necessary to convert between a %Mat and a Matrix (to interface with
%FactorQTZ, for instance). In the example below, we print a Mat33 created from
a 3x3 Matrix:
\code
Matrix myMatrix(3,3);
for (int i=0; i<9; ++i) myMatrix[i/3][i%3] = i+1;
Mat33 myMat33 = Mat33(&myMatrix[0][0]).transpose();
std::cout << myMat33 << std::endl;
\endcode
Converting from a Mat33 to a Matrix is straightforward:
\code
Mat33 myMat33(1,2,3, 4,5,6, 7,8,9);
std::cout << Matrix(myMat33) << std::endl;
\endcode

@see Matrix_ for handling of large or variable-sized matrices.
@see SymMat, Vec, Row
**/
template <int M, int N, class ELT, int CS, int RS> class Mat {
public:
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

    /** Every Composite Numerical Type (CNT) must define these values. **/
    #ifndef SWIG
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
        IsULessScalar       = 0,
        IsNumber            = 0,
        IsStdNumber         = 0,
        IsPrecision         = 0,
        SignInterpretation  = CNT<E>::SignInterpretation
    };
    #endif

    typedef Mat<M,N,E,CS,RS>                T;
    typedef Mat<M,N,ENeg,CS,RS>             TNeg;
    typedef Mat<M,N,EWithoutNegator,CS,RS>  TWithoutNegator;

    typedef Mat<M,N,EReal,CS*CNT<E>::RealStrideFactor,RS*CNT<E>::RealStrideFactor>
                                            TReal;
    typedef Mat<M,N,EImag,CS*CNT<E>::RealStrideFactor,RS*CNT<E>::RealStrideFactor>
                                            TImag;
    typedef Mat<M,N,EComplex,CS,RS>         TComplex;
    typedef Mat<N,M,EHerm,RS,CS>            THerm;
    typedef Mat<N,M,E,RS,CS>                TPosTrans;
    typedef E                               TElement;
    typedef Row<N,E,CS>                     TRow;    
    typedef Vec<M,E,RS>                     TCol;
    typedef Vec<MinDim,E,RS+CS>             TDiag;

    // These are the results of calculations, so are returned in new, packed
    // memory. Be sure to refer to element types here which are also packed.
    typedef Mat<M,N,ESqrt,M,1>              TSqrt;      // Note strides are packed
    typedef Mat<M,N,EAbs,M,1>               TAbs;       // Note strides are packed
    typedef Mat<M,N,EStandard,M,1>          TStandard;
    typedef Mat<N,M,EInvert,N,1>            TInvert;    // like THerm but packed
    typedef Mat<M,N,ENormalize,M,1>         TNormalize;

    typedef SymMat<N,ESqHermT>              TSqHermT;   // ~Mat*Mat
    typedef SymMat<M,ESqTHerm>              TSqTHerm;   // Mat*~Mat

    // Here the elements are copied unchanged but the result matrix
    // is an ordinary packed, column order matrix.
    typedef Mat<M,N,E,M,1>                  TPacked;
    typedef Mat<M-1,N,E,M-1,1>              TDropRow;
    typedef Mat<M,N-1,E,M,1>                TDropCol;
    typedef Mat<M-1,N-1,E,M-1,1>            TDropRowCol;
    typedef Mat<M+1,N,E,M+1,1>              TAppendRow;
    typedef Mat<M,N+1,E,M,1>                TAppendCol;
    typedef Mat<M+1,N+1,E,M+1,1>            TAppendRowCol;

    typedef EScalar                         Scalar;
    typedef EULessScalar                    ULessScalar;
    typedef ENumber                         Number;
    typedef EStdNumber                      StdNumber;
    typedef EPrecision                      Precision;
    typedef EScalarNormSq                   ScalarNormSq;

    typedef THerm                           TransposeType; // TODO

    /** Return the total number of elements M*N contained in this Mat. **/
    static int size() { return M*N; }
    /** Return the number of rows in this Mat, echoing the value supplied
    for the template parameter \a M. **/
    static int nrow() { return M; }
    /** Return the number of columns in this Mat, echoing the value supplied
    for the template parameter \a N. **/
    static int ncol() { return N; }

    /** Scalar norm square is the sum of squares of all the scalars that 
    comprise the value of this Mat. For Mat objects with composite element
    types, this is defined recursively as the sum of the scalar norm squares 
    of all the elements, where the scalar norm square of a scalar is just that
    scalar squared. **/
    ScalarNormSq scalarNormSqr() const { 
        ScalarNormSq sum(0);
        for(int j=0;j<N;++j) sum += CNT<TCol>::scalarNormSqr((*this)(j));
        return sum;
    }

    /** Elementwise square root; that is, the return value has the same
    dimensions as this Mat but with each element replaced by whatever it thinks
    its square root is. **/
    TSqrt sqrt() const { 
        TSqrt msqrt;
        for(int j=0;j<N;++j) msqrt(j) = (*this)(j).sqrt();
        return msqrt;
    }

    /** Elementwise absolute value; that is, the return value has the same
    dimensions as this Mat but with each element replaced by whatever it thinks
    its absolute value is. **/
    TAbs abs() const { 
        TAbs mabs;
        for(int j=0;j<N;++j) mabs(j) = (*this)(j).abs();
        return mabs;
    }

    TStandard standardize() const {
        TStandard mstd;
        for(int j=0;j<N;++j) mstd(j) = (*this)(j).standardize();
        return mstd;
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

    /** Default construction initializes to NaN when debugging but is left 
    uninitialized otherwise to ensure that there is no overhead. **/
    Mat(){ 
    #ifndef NDEBUG
        setToNaN();
    #endif
    }

    // It's important not to use the default copy constructor or copy
    // assignment because the compiler doesn't understand that we may
    // have noncontiguous storage and will try to copy the whole array.

    /** Copy constructor copies only the elements that are present and does
    not touch any unused memory space between them if they are not packed. **/
    Mat(const Mat& src) {
        for (int j=0; j<N; ++j)
            (*this)(j) = src(j);
    }
    /** Copy assignment copies only the elements that are present and does
    not touch any unused memory space between them if they are not packed. 
    Works correctly even if source and destination are the same object. **/
    Mat& operator=(const Mat& src) {    
        for (int j=0; j<N; ++j)
           (*this)(j) = src(j); // no harm if src and 'this' are the same
        return *this;
    }

    /** Explicit construction of a Mat from a SymMat (symmetric/Hermitian 
    matrix). Note that a SymMat is a Hermitian matrix when the elements are 
    complex, so in that case the resulting Mat's upper triangle values are 
    complex conjugates of the lower triangle ones. **/
    explicit Mat(const SymMat<M, ELT>& src) {
        updDiag() = src.diag();
        for (int j = 0; j < M; ++j)
            for (int i = j+1; i < M; ++i) {
                (*this)(i, j) = src.getEltLower(i, j);
                (*this)(j, i) = src.getEltUpper(j, i);
            }
    }

    /** This provides an \e implicit conversion from a Mat of the same 
    dimensions and element type but with different element spacing. 
    @tparam CSS Column spacing of the source Mat.
    @tparam RSS %Row spacing of the source Mat. **/
    template <int CSS, int RSS> 
    Mat(const Mat<M,N,E,CSS,RSS>& src) {
        for (int j=0; j<N; ++j)
            (*this)(j) = src(j);
    }

    /** This provides an \e implicit conversion from a Mat of the same 
    dimensions and \e negated element type, possibly with different element 
    spacing.
    @tparam CSS Column spacing of the source Mat.
    @tparam RSS %Row spacing of the source Mat. **/
    template <int CSS, int RSS> 
    Mat(const Mat<M,N,ENeg,CSS,RSS>& src) {
        for (int j=0; j<N; ++j)
            (*this)(j) = src(j);
    }

    /** Explicit construction of a Mat from a source Mat of the same 
    dimensions and an assignment-compatible element type, with any element 
    spacing allowed.
    @tparam EE  The element type of the source Mat; must be assignment 
                compatible with element type E of this Mat.
    @tparam CSS Column spacing of the source Mat.
    @tparam RSS %Row spacing of the source Mat. **/
    template <class EE, int CSS, int RSS> 
    explicit Mat(const Mat<M,N,EE,CSS,RSS>& mm)
      { for (int j=0;j<N;++j) (*this)(j) = mm(j);}

    /** Explicit construction from a single element \a e of this Mat's element
    type E sets all the main diagonal elements to \a e but sets the rest of 
    the elements to zero. **/
    explicit Mat(const E& e)
      { for (int j=0;j<N;++j) (*this)(j) = E(0); diag()=e; }

    /** Explicit construction from a single element \a e whose type is
    negator<E> (abbreviated ENeg here) where E is this Mat's element
    type sets all the main diagonal elements to \a e but sets the rest of 
    the elements to zero. **/
    explicit Mat(const ENeg& e)
      { for (int j=0;j<N;++j) (*this)(j) = E(0); diag()=e; }

    /** Explicit construction from an int value means we convert the int into
    an object of this Mat's element type E, and then apply the single-element
    constructor above which sets the Mat to zero except for its main diagonal
    elements which will all be set to the given value. To convert an int to
    an element, we first turn it into the appropriate-precision floating point 
    number, and then call E's constructor that takes a single scalar. **/
    explicit Mat(int i) 
      { new (this) Mat(E(Precision(i))); }

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

    // m= this .* m
    template <class E2, int CS2, int RS2> 
    typename EltResult<E2>::Mul
    elementwiseMultiply(const Mat<M,N,E2,CS2,RS2>& r) const {
        typename EltResult<E2>::Mul result;
        for (int j=0;j<N;++j) 
            result(j) = (*this)(j).elementwiseMultiply(r(j));
        return result;
    }

    // m= this ./ m
    template <class E2, int CS2, int RS2> 
    typename EltResult<E2>::Dvd
    elementwiseDivide(const Mat<M,N,E2,CS2,RS2>& r) const {
        typename EltResult<E2>::Dvd result;
        for (int j=0;j<N;++j) 
            result(j) = (*this)(j).elementwiseDivide(r(j));
        return result;
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
    
    const E& operator()(int i,int j) const { return elt(i,j); }
    E&       operator()(int i,int j)       { return elt(i,j); }

    // This is the scalar Frobenius norm.
    ScalarNormSq normSqr() const { return scalarNormSqr(); }
    typename CNT<ScalarNormSq>::TSqrt 
        norm() const { return CNT<ScalarNormSq>::sqrt(scalarNormSqr()); }

    // There is no conventional meaning for normalize() applied to a matrix. We
    // choose to define it as follows:
    // If the elements of this Mat are scalars, the result is what you get by
    // dividing each element by the Frobenius norm() calculated above. If the elements are
    // *not* scalars, then the elements are *separately* normalized. That means
    // you will get a different answer from Mat<2,2,Mat33>::normalize() than you
    // would from a Mat<6,6>::normalize() containing the same scalars.
    //
    // Normalize returns a matrix of the same dimension but in new, packed storage
    // and with a return type that does not include negator<> even if the original
    // Mat<> does, because we can eliminate the negation here almost for free.
    // But we can't standardize (change conjugate to complex) for free, so we'll retain
    // conjugates if there are any.
    TNormalize normalize() const {
        if (CNT<E>::IsScalar) {
            return castAwayNegatorIfAny() / (SignInterpretation*norm());
        } else {
            TNormalize elementwiseNormalized;
            // punt to the column Vec to deal with the elements
            for (int j=0; j<N; ++j) 
                elementwiseNormalized(j) = (*this)(j).normalize();
            return elementwiseNormalized;
        }
    }

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

    // If the underlying scalars are complex or conjugate, we can return a
    // reference to the real part just by recasting the data to a matrix of
    // the same dimensions but containing real elements, with the scalar
    // spacing doubled.
    const TReal& real() const { return *reinterpret_cast<const TReal*>(this); }
    TReal&       real()       { return *reinterpret_cast<      TReal*>(this); }

    // Getting the imaginary part is almost the same as real, but we have
    // to shift the starting address of the returned object by 1 real-size
    // ("precision") scalar so that the first element is the imaginary part
    // of the original first element.
    // TODO: should blow up or return a reference to a zero matrix if called
    // on a real object.
    // Had to contort these routines to get them through VC++ 7.net
    const TImag& imag()    const { 
        const int offs = ImagOffset;
        const Precision* p = reinterpret_cast<const Precision*>(this);
        return *reinterpret_cast<const TImag*>(p+offs);
    }
    TImag& imag() { 
        const int offs = ImagOffset;
        Precision* p = reinterpret_cast<Precision*>(this);
        return *reinterpret_cast<TImag*>(p+offs);
    }

    const TWithoutNegator& castAwayNegatorIfAny() const {return *reinterpret_cast<const TWithoutNegator*>(this);}
    TWithoutNegator&       updCastAwayNegatorIfAny()    {return *reinterpret_cast<TWithoutNegator*>(this);}

    const TRow& row(int i) const { 
        SimTK_INDEXCHECK(i,M, "Mat::row[i]");
        return *reinterpret_cast<const TRow*>(&d[i*RS]); 
    }
    TRow& row(int i) { 
        SimTK_INDEXCHECK(i,M, "Mat::row[i]");
        return *reinterpret_cast<TRow*>(&d[i*RS]); 
    }

    const TCol& col(int j) const { 
        SimTK_INDEXCHECK(j,N, "Mat::col(j)");
        return *reinterpret_cast<const TCol*>(&d[j*CS]); 
    }
    TCol& col(int j) { 
        SimTK_INDEXCHECK(j,N, "Mat::col(j)");
        return *reinterpret_cast<TCol*>(&d[j*CS]); 
    }    
    
    const E& elt(int i, int j) const {
        SimTK_INDEXCHECK(i,M, "Mat::elt(i,j)");
        SimTK_INDEXCHECK(j,N, "Mat::elt(i,j)");
        return d[i*RS+j*CS]; 
    }
    E& elt(int i, int j) { 
        SimTK_INDEXCHECK(i,M, "Mat::elt(i,j)");
        SimTK_INDEXCHECK(j,N, "Mat::elt(i,j)");
        return d[i*RS+j*CS]; 
    }

    /// Select main diagonal (of largest leading square if rectangular) and
    /// return it as a read-only view (as a Vec) of the diagonal elements 
    /// of this Mat.
    const TDiag& diag() const { return *reinterpret_cast<const TDiag*>(d); }
    /// Select main diagonal (of largest leading square if rectangular) and
    /// return it as a writable view (as a Vec) of the diagonal elements 
    /// of this Mat.
    TDiag&       updDiag()    { return *reinterpret_cast<TDiag*>(d); }
    /// This non-const version of diag() is an alternate name for updDiag()
    /// available for historical reasons.
    TDiag&       diag()       { return *reinterpret_cast<TDiag*>(d); }

    EStandard trace() const {return diag().sum();}

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

    // Additive operators for scalars operate only on the diagonal.
    template <class EE> Mat<M,N, typename CNT<E>::template Result<EE>::Add>
    scalarAdd(const EE& e) const {
        Mat<M,N, typename CNT<E>::template Result<EE>::Add> result(*this);
        result.diag() += e;
        return result;
    }
    // Add is commutative, so no 'FromLeft'.

    template <class EE> Mat<M,N, typename CNT<E>::template Result<EE>::Sub>
    scalarSubtract(const EE& e) const {
        Mat<M,N, typename CNT<E>::template Result<EE>::Sub> result(*this);
        result.diag() -= e;
        return result;
    }
    // Should probably do something clever with negation here (s - m)
    template <class EE> Mat<M,N, typename CNT<EE>::template Result<E>::Sub>
    scalarSubtractFromLeft(const EE& e) const {
        Mat<M,N, typename CNT<EE>::template Result<E>::Sub> result(-(*this));
        result.diag() += e; // yes, add
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
      { for(int j=0; j<N; ++j) (*this)(j).scalarEq(EE(0)); 
        diag().scalarEq(ee); 
        return *this; }

    template <class EE> Mat& scalarPlusEq(const EE& ee)
      { diag().scalarPlusEq(ee); return *this; }

    template <class EE> Mat& scalarMinusEq(const EE& ee)
      { diag().scalarMinusEq(ee); return *this; }
    // m = s - m; negate m, then add s
    template <class EE> Mat& scalarMinusEqFromLeft(const EE& ee)
      { scalarTimesEq(E(-1)); diag().scalarAdd(ee); return *this; }

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

    void setToZero() {
        for (int j=0; j<N; ++j)
            (*this)(j).setToZero();
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
    template <int MM, int NN>
    void setSubMat(int i, int j, const typename SubMat<MM,NN>::Type& value) {
        assert(0 <= i && i + MM <= M);
        assert(0 <= j && j + NN <= N);
        SubMat<MM,NN>::Type::updAs(&(*this)(i,j)) = value;
    }

    /// Return a matrix one row smaller than this one by dropping row
    /// i. The result is packed but has same element type as this one.
    TDropRow dropRow(int i) const {
        assert(0 <= i && i < M);
        TDropRow out;
        for (int r=0, nxt=0; r<M-1; ++r, ++nxt) {
            if (nxt==i) ++nxt;  // skip the loser
            out[r] = (*this)[nxt];
        }
        return out;
    }

    /// Return a matrix one column smaller than this one by dropping column
    /// j. The result is packed but has same element type as this one.
    TDropCol dropCol(int j) const {
        assert(0 <= j && j < N);
        TDropCol out;
        for (int c=0, nxt=0; c<N-1; ++c, ++nxt) {
            if (nxt==j) ++nxt;  // skip the loser
            out(c) = (*this)(nxt);
        }
        return out;
    }

    /// Return a matrix one row and one column smaller than this one by 
    /// dropping row i and column j. The result is packed but has same 
    /// element type as this one.
    TDropRowCol dropRowCol(int i, int j) const {
        assert(0 <= i && i < M);
        assert(0 <= j && j < N);
        TDropRowCol out;
        for (int c=0, nxtc=0; c<N-1; ++c, ++nxtc) { 
            if (nxtc==j) ++nxtc;
            for (int r=0, nxtr=0; r<M-1; ++r, ++nxtr) {
                if (nxtr==i) ++nxtr;
                out(r,c) = (*this)(nxtr,nxtc);
            }
        }
        return out;
    }

    /// Return a matrix one row larger than this one by adding a row
    /// to the end. The result is packed but has same element type as
    /// this one. Works for any assignment compatible row.
    template <class EE, int SS> 
    TAppendRow appendRow(const Row<N,EE,SS>& row) const {
        TAppendRow out;
        out.template updSubMat<M,N>(0,0) = (*this);
        out[M] = row;
        return out;
    }

    /// Return a matrix one column larger than this one by adding a column
    /// to the end. The result is packed but has same element type as
    /// this one. Works for any assignment compatible column.
    template <class EE, int SS> 
    TAppendCol appendCol(const Vec<M,EE,SS>& col) const {
        TAppendCol out;
        out.template updSubMat<M,N>(0,0) = (*this);
        out(N) = col;
        return out;
    }

    /// Return a matrix one row and one column larger than this one by 
    /// adding a row to the bottom and a column to the right. The final
    /// element of the row is ignored; that value is taken from the 
    /// final element of the column instead. The result is packed
    /// but has same element type as this one. Works for any assignment 
    /// compatible row and column.
    template <class ER, int SR, class EC, int SC> 
    TAppendRowCol appendRowCol(const Row<N+1,ER,SR>& row,
                               const Vec<M+1,EC,SC>& col) const 
    {
        TAppendRowCol out;
        out.template updSubMat<M,N>(0,0) = (*this);
        out[M].template updSubRow<N>(0) = 
            row.template getSubRow<N>(0); // ignore last element
        out(N) = col;
        return out;
    }

    /// Return a matrix one row larger than this one by inserting a row
    /// *before* row i. The result is packed but has same element type as
    /// this one. Works for any assignment compatible row. The index
    /// can be one greater than normally allowed in which case the row
    /// is appended.
    template <class EE, int SS> 
    TAppendRow insertRow(int i, const Row<N,EE,SS>& row) const {
        assert(0 <= i && i <= M);
        if (i==M) return appendRow(row);
        TAppendRow out;
        for (int r=0, nxt=0; r<M; ++r, ++nxt) {
            if (nxt==i) out[nxt++] = row;
            out[nxt] = (*this)[r];
        }
        return out;
    }

    /// Return a matrix one column larger than this one by inserting a column
    /// *before* column j. The result is packed but has same element type as
    /// this one. Works for any assignment compatible column. The index
    /// can be one greater than normally allowed in which case the column
    /// is appended.
    template <class EE, int SS> 
    TAppendCol insertCol(int j, const Vec<M,EE,SS>& col) const {
        assert(0 <= j && j <= N);
        if (j==N) return appendCol(col);
        TAppendCol out;
        for (int c=0, nxt=0; c<N; ++c, ++nxt) {
            if (nxt==j) out(nxt++) = col;
            out(nxt) = (*this)(c);
        }
        return out;
    }

    /// Return a matrix one row and one column larger than this one by 
    /// inserting a row *before* row i and a column *before* column j. 
    /// The intersecting element of the row is ignored; that element is
    /// taken from the column. The result is packed but has same element 
    /// type as this one. Works for any assignment compatible row and 
    /// column. The indices can be one greater than normally allowed 
    /// in which case the row or column is appended.
    template <class ER, int SR, class EC, int SC>
    TAppendRowCol insertRowCol(int i, int j, const Row<N+1,ER,SR>& row,
                                             const Vec<M+1,EC,SC>& col) const {
        assert(0 <= i && i <= M);
        assert(0 <= j && j <= N);
        TAppendRowCol out;
        for (int c=0, nxtc=0; c<N; ++c, ++nxtc) { 
            if (nxtc==j) ++nxtc;   // leave room
            for (int r=0, nxtr=0; r<M; ++r, ++nxtr) {
                if (nxtr==i) ++nxtr;
                out(nxtr,nxtc) = (*this)(r,c);
            }
        }
        out[i] = row;
        out(j) = col; // overwrites row's j'th element
        return out;
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

    /// Return true if any element of this Mat contains a NaN anywhere.
    bool isNaN() const {
        for (int j=0; j<N; ++j)
            if (this->col(j).isNaN())
                return true;
        return false;
    }

    /// Return true if any element of this Mat contains a +Inf
    /// or -Inf somewhere but no element contains a NaN anywhere.
    bool isInf() const {
        bool seenInf = false;
        for (int j=0; j<N; ++j) {
            if (!this->col(j).isFinite()) {
                if (!this->col(j).isInf()) 
                    return false; // something bad was found
                seenInf = true; 
            }
        }
        return seenInf;
    }

    /// Return true if no element contains an Infinity or a NaN.
    bool isFinite() const {
        for (int j=0; j<N; ++j)
            if (!this->col(j).isFinite())
                return false;
        return true;
    }

    /// For approximate comparisons, the default tolerance to use for a matrix is
    /// its shortest dimension times its elements' default tolerance.
    static double getDefaultTolerance() {return MinDim*CNT<ELT>::getDefaultTolerance();}

    /// %Test whether this matrix is numerically equal to some other matrix with
    /// the same shape, using a specified tolerance.
    template <class E2, int CS2, int RS2>
    bool isNumericallyEqual(const Mat<M,N,E2,CS2,RS2>& m, double tol) const {
        for (int j=0; j < N; ++j)
            if (!(*this)(j).isNumericallyEqual(m(j), tol))
                return false;
        return true;
    }

    /// %Test whether this matrix is numerically equal to some other matrix with
    /// the same shape, using a default tolerance which is the looser of the
    /// default tolerances of the two objects being compared.
    template <class E2, int CS2, int RS2>
    bool isNumericallyEqual(const Mat<M,N,E2,CS2,RS2>& m) const {
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
        for (int i=0; i<M; ++i)
            for (int j=0; j<N; ++j) {
                if (i==j) {
                    if (!CNT<ELT>::isNumericallyEqual((*this)(i,i), e, tol))
                        return false;
                } else {
                    // off-diagonals must be zero
                    if (!CNT<ELT>::isNumericallyEqual((*this)(i,j), ELT(0), tol))
                        return false;
                }
            }
        return true;
    }

    /// A Matrix is symmetric (actually Hermitian) if it is square and each 
    /// element (i,j) is the Hermitian transpose of element (j,i). Here we
    /// are testing for numerical symmetry, meaning that the symmetry condition
    /// is satisified to within a tolerance (supplied or default). This is 
    /// a relatively expensive test since all elements must be examined but
    /// can be very useful in Debug mode to check assumptions.
    /// @see isExactlySymmetric() for a rarely-used exact equality test
    bool isNumericallySymmetric(double tol = getDefaultTolerance()) const {
        if (M != N) return false; // handled at compile time
        for (int j=0; j<M; ++j)
            for (int i=j; i<M; ++i)
                if (!CNT<ELT>::isNumericallyEqual(elt(j,i), CNT<ELT>::transpose(elt(i,j)), tol))
                    return false;
        return true;
    }

    /// A Matrix is symmetric (actually Hermitian) if it is square and each 
    /// element (i,j) is the Hermitian (conjugate) transpose of element (j,i). This
    /// method tests for exact (bitwise) equality and is too stringent for most 
    /// purposes; don't use it unless you know that the corresponding elements
    /// should be bitwise conjugates, typically because you put them there directly.
    /// @see isNumericallySymmetric() for a more useful method
    bool isExactlySymmetric() const {
        if (M != N) return false; // handled at compile time
        for (int j=0; j<M; ++j)
            for (int i=j; i<M; ++i)
                if (elt(j,i) != CNT<ELT>::transpose(elt(i,j)))
                    return false;
        return true;
    }
    
    /// Returns a row vector (Row) containing the column sums of this matrix.
    TRow colSum() const {
        TRow temp;
        for (int j = 0; j < N; ++j)
            temp[j] = col(j).sum();
        return temp;
    }
    /// This is an alternate name for colSum(); behaves like the Matlab
    /// function of the same name.
    TRow sum() const {return colSum();}

    /// Returns a column vector (Vec) containing the row sums of this matrix.
    TCol rowSum() const {
        TCol temp;
        for (int i = 0; i < M; ++i)
            temp[i] = row(i).sum();
        return temp;
    }

    // Functions to be used for Scripting in MATLAB and languages that do not support operator overloading
    /** toString() returns a string representation of the Mat. Please refer to operator<< for details. **/
    std::string toString() const {
        std::stringstream stream;
        stream <<  (*this) ;
        return stream.str(); 
    }
    /** Variant of indexing operator that's scripting friendly to get entry (i, j) **/
    const ELT& get(int i,int j) const { return elt(i,j); }
    /** Variant of indexing operator that's scripting friendly to set entry (i, j) **/
    void       set(int i,int j, const ELT& value)       { elt(i,j)=value; }

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
// TODO: for now it is just going to call mat.invert() which will fail on
// singular matrices.

// m = m/real, real/m 
template <int M, int N, class E, int CS, int RS> inline
typename Mat<M,N,E,CS,RS>::template Result<float>::Dvd
operator/(const Mat<M,N,E,CS,RS>& l, const float& r)
{   return Mat<M,N,E,CS,RS>::template Result<float>::DvdOp::perform(l,r); }

template <int M, int N, class E, int CS, int RS> inline
typename CNT<float>::template Result<Mat<M,N,E,CS,RS> >::Dvd
operator/(const float& l, const Mat<M,N,E,CS,RS>& r)
{   return l * r.invert(); }

template <int M, int N, class E, int CS, int RS> inline
typename Mat<M,N,E,CS,RS>::template Result<double>::Dvd
operator/(const Mat<M,N,E,CS,RS>& l, const double& r)
{   return Mat<M,N,E,CS,RS>::template Result<double>::DvdOp::perform(l,r); }

template <int M, int N, class E, int CS, int RS> inline
typename CNT<double>::template Result<Mat<M,N,E,CS,RS> >::Dvd
operator/(const double& l, const Mat<M,N,E,CS,RS>& r)
{   return l * r.invert(); }

// m = m/int, int/m -- just convert int to m's precision float
template <int M, int N, class E, int CS, int RS> inline
typename Mat<M,N,E,CS,RS>::template Result<typename CNT<E>::Precision>::Dvd
operator/(const Mat<M,N,E,CS,RS>& l, int r) 
{   return l / (typename CNT<E>::Precision)r; }

template <int M, int N, class E, int CS, int RS> inline
typename CNT<typename CNT<E>::Precision>::template Result<Mat<M,N,E,CS,RS> >::Dvd
operator/(int l, const Mat<M,N,E,CS,RS>& r) 
{   return (typename CNT<E>::Precision)l / r; }


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
